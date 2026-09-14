# Design: a repeatable mutation-testing harness

**Status:** proposal, not implemented. **Verified against:** `main` @ `f9971d0` (2026-08-27);
the Cython-extension prerequisite (below) additionally re-verified against `main` @ `ad89735`.

## Problem

Seven tests in this repo passed regardless of whether the code they claimed to protect
worked, and all seven were found by hand-run mutation testing, not by reading:

- `ba60b80`: `TestSeFilterGating`'s only test passed TSV input, which a neutralizer
  always clears `single_end`/`se_max_fragment_length` on before the SE-length gate at
  `fragments_h5.py:782` is reached. The gate had zero coverage.
- `214b214`: `test_gs_url_takes_remote_branch` asserted `_is_remote_url("gs://...") == True`
  and never drove the code path that uses it.
- `98a0c32`/current `test_contig_with_zero_mapped_reads_skipped`: asserted a contig was
  absent from output, true whether or not the skip logic ran (its own docstring now says so).
- `d6c1d16`: an all-`255` GC sentinel assertion was vacuous because `255` decodes to `NaN` on read.
- `5262f3a` finding 1: a `pytest.raises(KeyError)` passed via an unrelated downstream lookup
  raising the same exception type, independent of the gate under test.
- `bb6f0d9`: `_build_code_revision`'s write site had no test; replacing the write with `pass`
  left all 17 tests in `test_build_provenance.py` green.
- `f9971d0`: `test_revision_dist_editable_branch` asserted a property of the *machine*
  (that its install happened to be editable), not of the code.

Mutation testing (mutate one line, confirm the target test goes red, restore) is the one
technique with a track record here. It has also failed, by hand, twice:

- An auditor ran a mutated copy under `cp -r src/ /tmp && PYTHONPATH=/tmp/src ...`. This
  displaces `__file__`; `_resolve_build_code_revision`'s tracked-file gate
  (`git -C dirname(__file__) ls-files --error-unmatch ...`) then hits "not a git repository"
  and falls through to a different branch than production ever takes for anyone with a real
  checkout. Confirmed live in this session: `git ls-files` on a plain `cp` copy exits 128;
  the intended branch is silently unreachable. The test went red for a reason unrelated to
  the mutation — a false positive, caught only because someone happened to check why.
- A mutation targeted `fragments_h5.py:1023`; per `fragment_selection_and_build_provenance.md:815-819`,
  that address had already drifted to `:1021` after one earlier edit in the same review round.
- Two design docs in this repo hit 854 and 1023 lines; a round-2 reviewer found a Critical
  defect "hundreds of lines above the prescription that contradicted it" and could only
  confirm it "by patching a throwaway copy and running it" — reading did not surface it.

This doc designs a small, committed tool so the check stops depending on memory and luck.

## Decision: isolate via `git worktree`, never mutate the shared tree

Every prior attempt mutated `src/` in place (or a manual `cp`) and relied on remembering to
restore it. `try/finally` does not run under `SIGKILL`; a mutation left in the tree that
another session reads is worse than no harness, and this tree is shared by more than one
active session right now (`.claude/worktrees/debug-fragments-h5-parallel-gc` and
`tsv-fragment-input` both exist alongside this one).

The harness instead runs `git worktree add <ephemeral-path> HEAD` and mutates only the file
inside that linked worktree. This one decision resolves three of the six problems at once:

1. **Restoration under kill.** The shared tree is never written to, so there is nothing to
   restore in it. A `SIGKILL` mid-run leaves an abandoned worktree directory and, at worst,
   a stale entry in `git worktree list` — inert bookkeeping, cleaned by `git worktree prune`,
   not a corrupted tracked file. This is a strictly stronger guarantee than "restore in a
   finally block," which a kill signal bypasses.
2. **Concurrent sessions.** Each invocation gets its own path (`.mutation-tmp/<uuid8>/`,
   repo-root-local). Two simultaneous runs never touch the same files; `git worktree add`
   with distinct paths does not conflict. **Adding `.mutation-tmp/` to `.gitignore` is an
   implementation prerequisite, not cosmetic — it is not there today.** A killed run
   leaves an untracked `.mutation-tmp/` directory behind, and `make require-clean-tree`
   (which gates `make tag`, `make conda-build`, and `make docker-build`) fails any
   untracked file via `git ls-files --others --exclude-standard` — verified against this
   repo's `Makefile`. Without the gitignore entry, one killed mutation run silently blocks
   the next release until someone works out why.
3. **The exact false positive above stops being possible.** A linked worktree is a real git
   checkout (`git -C <worktree> ls-files` succeeds), unlike a `cp`. Verified live in a
   throwaway repo this session (see Worked example): the tracked-file gate passes in a
   worktree and fails in a plain copy, for the identical file.

## Prerequisite: the compiled Cython extension

`git worktree add` checks out tracked files only; `sequence*.so` is gitignored, so a fresh
worktree has no compiled extension, and the import chain is universal —
`fragments_h5/__init__.py` → `fragments_h5.py` → `fragment.py`'s `from fragments_h5.sequence
import one_hot_encode_sequences` — so **every** test file fails at collection. Verified
against a real worktree of this repo: `ERROR tests/test_build_provenance.py` /
`Interrupted: 1 error during collection`, `ModuleNotFoundError: No module named
'fragments_h5.sequence'`; copying `sequence*.so` from the main tree's `src/fragments_h5/`
into the worktree's `src/fragments_h5/` before the baseline run fixes it (18 tests then
collect). The harness must do this copy as a setup step.

That reuse is only conditionally valid: a `.so` built from a different `sequence.pyx` is
silently wrong — the harness would then pin mutations against compiled code that disagrees
with the source, the exact class of error it exists to catch. The cheap check is mtimes:
before creating the worktree, compare `sequence.pyx`'s mtime against its `.so`'s mtime in
the *main* tree (older `.so` than `.pyx` means not rebuilt since the last edit). If `.pyx`
is newer, hard-error and refuse to proceed — same fail-loud posture as the anchor guard —
telling the caller to rebuild first. Warning and continuing would reintroduce the exact
silent-wrongness this check exists to prevent.

## Mutation spec: content anchor, not line number, with a hard zero/multi guard

A mutation is `(file, anchor, replacement, test_nodeid)`. `anchor` is an exact source
substring (a full line, or a few lines) copied from the file — the same unit already used
by hand in `ba60b80`'s commit message (`` `single_end and` ``). The harness applies it as
`content.replace(anchor, replacement, 1)` after checking `content.count(anchor)`:

- **0 matches → hard error, non-zero exit, no test is run.** Silently doing nothing and
  then reporting "the test still passed against the mutant" would manufacture exactly the
  false confidence this tool exists to detect. This is the one failure mode that must never
  degrade into a result.
- **>1 matches → hard error.** An ambiguous anchor could silently mutate the wrong
  occurrence; require the caller to add enough surrounding context to be unique.
- Anchors rot too (a line can be deleted or reworded), but they fail loudly at the point of
  use instead of silently targeting the wrong line, which is what a line number does.

## What is asserted — the crux

"The test failed" is not "the test failed for the right reason." The false positive above
came from exactly this conflation. The harness must distinguish three outcomes, not two:

1. **Baseline, in the worktree, before mutating:** run only `test_nodeid` with
   `--junitxml`. Require exactly one `<testcase>` with no `<failure>` and no `<error>`
   child (pytest's own outcome categories — no custom parsing needed). If baseline does not
   cleanly pass, abort with "baseline failed in isolation" and never touch the mutation —
   this is what would have caught the editable-install false failure in `f9971d0` early,
   and what distinguishes "isolation broke something" from "the mutation had an effect."
2. **Mutant:** apply the mutation, run the identical `test_nodeid` the identical way.
   Require exactly one `<testcase>` with a `<failure>` child. An `<error>` (collection or
   fixture failure) is reported separately and is **not** treated as caught — a syntactically
   broken mutation would trivially "fail" every test in the file for a reason that has
   nothing to do with the test's logic.
3. Report both raw outcomes. "1 passed → 1 failed" is the only claim the harness makes; it
   does not attempt traceback provenance (confirming the failure line is inside the target
   test body, as opposed to a fixture) — pytest's failure/error split already rules out the
   collection-error class of false positive, and going further adds real complexity for a
   marginal gain given how narrow this tool's job is. Out of scope, noted here on purpose.
   One variant this misses concretely: a mutation inside a shared test helper (e.g.
   `_make_simple_bam()`) can turn the target test red for a reason unrelated to its own
   assertions — a genuine failure, not a collection error, so the failed/error split does
   not catch it either. Mutations should target production code, not shared test
   infrastructure, for this reason.

Restricting both runs to a single `test_nodeid` (not the file, not the suite) is what makes
this check about *that* assertion protecting *that* line, rather than "something somewhere
changed" — the ambiguity that caused several of the seven original gaps.

## Form, and CI

A standalone script, `scripts/mutation_check.py`, invoked directly or via a Makefile target
(`make mutation-check FILE=... ANCHOR=... REPLACEMENT=... TEST=...`), matching this repo's
existing `scripts/*.sh` + `Makefile` convention. Not a pytest plugin — there is nothing to
hook into pytest's collection for; it shells out to a fresh `pytest` process per run by
design, so a crash in the mutant's test process can't take the harness down with it.

**No CI exists in this repo** — no `.github/workflows`, no `.gitlab-ci.yml`, no Jenkinsfile
(checked directly). This tool is therefore audit tooling, run on demand by whoever is
reviewing test coverage, the same way the manual recipe in `ba60b80`/`bb6f0d9`/`5262f3a` was
already being run — not a merge gate. Wiring it into CI is out of scope until CI exists.

## Manifest: a home for the seven known pairs

Four positional arguments per invocation have nowhere to live between runs. Re-checking
the seven pairs in the Problem section by hand-typing `mutation-check` seven times, with
multi-line anchors, will not happen — and an unused harness is worse than the documented
recipe it replaces. v1 therefore adds a manifest, `tests/mutation_findings.yaml`: a list
of `{file, anchor, replacement, test_nodeid, note}` entries, one per finding above, plus
a `make mutation-check-all` target that reads it and invokes `mutation_check.py` once per
entry, printing a caught/not-caught line per pair and exiting non-zero if any pair is no
longer caught. `make mutation-check FILE=... ANCHOR=... REPLACEMENT=... TEST=...` remains
the primitive for adding an eighth pair by hand; the manifest is what turns re-checking
the first seven into one command instead of seven.

## Scope vs. `mutmut` / `cosmic-ray`

Neither is installed (`ModuleNotFoundError` for both, checked directly), and prior work in
this repo already noted their absence (`worker_args_refactor.md:236`). Both are built for
*discovering* unknown gaps: they generate the full space of AST mutations for a module and
run the *entire* suite against each one. That is the wrong shape for this repo's actual,
demonstrated need — pinning seven *already-known* (mutation, test) pairs found by audit, and
re-checking them cheaply as the code changes. Running either tool's default sweep against
this suite (92 tests, 331s wall time measured this session, includes multiprocessing stress
tests) per generated mutant would be slow and would not let a reviewer say "this specific
assertion protects this specific line" — the actual crux above. Adopting a general framework
for a one-test-one-line check would also reintroduce the sprawl this repo's own design docs
have already been penalized for. Not adopted. If the goal ever shifts from *pinning known
gaps* to *discovering unknown ones*, mutmut is worth a second look then, as a separate tool.

## Worked example (executed this session, not transcribed)

In a throwaway repo (`/tmp/mh_demo`, unrelated to this working tree, no writing git command
touched `fragments_h5`):

```
module.py:  if flag and x is not None and x > 0: return x * 2   else return x
test:       assert gated_value(False, 5) == 5
```

1. `git -C /tmp/mh_demo_copy ls-files --error-unmatch module.py` on a plain `cp` (no
   `.git`) → **exit 128**. The same command against a `git worktree add` checkout →
   **exit 0**. This is the displacement bug and its fix, both reproduced live.
2. Anchor guard: `mutate.py` on a 0-match anchor → `SystemExit: ANCHOR NOT FOUND`. On a
   2-match anchor (`return 1` appearing in two functions) → `SystemExit: ANCHOR AMBIGUOUS`.
   Both refuse to proceed rather than mutating the wrong (or no) line.
3. Baseline, worktree, `--junitxml`: `1 passed`, `failures="0" errors="0"`.
4. Applied mutation (anchor `if flag and x is not None and x > 0:` →
   `if x is not None and x > 0:`), re-ran the identical nodeid in the same worktree:
   `1 failed`, `failures="1"`, `AssertionError: assert 10 == 5`.
5. `git -C /home/nathanboley/src/fragments_h5 status --short` throughout: empty. The shared
   tree was never touched.
6. `git worktree remove --force` cleaned up the demo worktree (in the throwaway repo).

## Reference sketch (not implemented)

```python
# scripts/mutation_check.py — sketch, ~60 lines at full size
def run(file, anchor, replacement, test_nodeid, repo_root):
    require_so_not_stale(repo_root)                             # .pyx mtime > .so mtime -> hard error
    wt = repo_root / ".mutation-tmp" / uuid.uuid4().hex[:8]     # gitignored, repo-local
    subprocess.run(["git", "worktree", "add", "-q", "--detach", str(wt), "HEAD"], check=True)
    try:
        copy_so_files(repo_root / "src/fragments_h5", wt / "src/fragments_h5")
        baseline = pytest_junit(wt, test_nodeid)                # PYTHONPATH=wt/src
        require(baseline.passed == 1 and baseline.failed == 0 and baseline.errors == 0,
                "baseline did not cleanly pass in isolation")
        apply_anchor_mutation(wt / file, anchor, replacement)    # 0/>1 match -> hard error
        mutant = pytest_junit(wt, test_nodeid)
        require(mutant.failed == 1 and mutant.errors == 0,
                "mutation NOT caught (test still green, or errored instead of failing)")
        print(f"CAUGHT: {test_nodeid} pins {file}:{anchor!r}")
    finally:
        subprocess.run(["git", "worktree", "remove", "--force", str(wt)])
```

## Self-assessment

- The worktree mechanism (tracked-file gate + anchor guard + red/green via junit outcome
  categories) was executed end-to-end this session, in an isolated throwaway repo, and every
  claim above about its behavior is a transcript of a real run, not a prediction.
- The Cython-extension prerequisite was subsequently verified against a real `git worktree
  add` of this repo (cleaned up with `git worktree remove` afterward): a fresh worktree
  reproduces `ERROR tests/test_build_provenance.py` / `Interrupted: 1 error during
  collection` exactly as prescribed, and copying `sequence*.so` in fixes it (18 tests
  collect). The anchor-guard/junit red-green mechanism itself is still only verified in the
  throwaway demo repo, not against a real anchor in `fragments_h5`'s own source — that
  remains the first thing an implementer should do.
- The "no traceback provenance check" scope cut (end of the Assertion section) is a real
  limitation, not a false one: a mutation that changes *which* line raises inside the same
  test function, while still landing on an assertion failure, would be reported as caught
  without checking it's the intended assertion. Judged acceptable because it matches the
  granularity of every fix in the Problem section above, none of which needed that.
- Length: 248 lines (`wc -l`) after addressing design review, against the well-under-250
  target — closer to the target than the original 197, but the additions were the minimum
  needed to close two required conditions plus two smaller items, not new scope.

## Addendum: script vs. parametrised pytest test — measured (2026-08-27, `main` @ `c47ae57`)

**Recommendation: hybrid, not a pick.** Keep `scripts/mutation_check.py` (`run_pair()`) as
the primitive an auditor calls by hand when adding an eighth pair — authoring a new
anchor/replacement needs iteration a test suite is a bad venue for. Add a thin
`tests/test_mutations.py` that imports the same `run_pair()` and parametrizes over
`tests/mutation_findings.yaml`, running under ordinary `pytest`, not opt-in. A pure script
reintroduces exactly the review's "nothing makes anyone re-run it" gap; a pure test with no
authoring primitive duplicates work the script already does. One implementation, two callers.

**Cost, measured.** Baseline this session: 92 passed, 258s (biomarker_env pin, no PATH
prefix) — the design's "331s" and this number both land in "a few minutes" but are not the
same measurement; suite timing is noisy across sessions, not a fixed constant. Prototyped
`run_pair()` end-to-end (worktree add → copy `.so` → baseline junit run → mutate → mutant
junit run → `worktree remove`) against two real pairs: `ba60b80`'s SE gate
(`test_gate_blocks_filter_when_not_single_end`, BAM build run twice) — **28.7s**; `bb6f0d9`'s
provenance write (`test_revision_written_into_real_build`) — **8.6s**. Worktree-add +
`.so`-copy is <0.3s combined in both; the cost is almost entirely the two pytest
startups/collections. Extrapolating the ~18.6s/pair average to all seven: **~130s**, growing
the suite from ~258s to **~390s (~1.5×)** if the pinned test runs unconditionally. Real, not
prohibitive.

**Nested pytest is sound, with one real gotcha.** Ran the prototype sequentially and twice
**concurrently**; `subprocess.run(capture_output=True)` isolates stdout cleanly, exit codes
are the expected 0/1. Surprise: pytest 9.x's junit XML root is `<testsuites><testsuite>...`,
not the bare `<testsuite>` the design's "no custom parsing needed" line implies — both forms
must handle the wrapper. `-p no:cacheprovider` is hygiene only. Outer `PYTEST_CURRENT_TEST`
leaks into the inner subprocess's env unless stripped; harmless here, worth stripping on
principle.

**Isolation and concurrency both hold from inside pytest.** The resolver returns
`git:v2.12.1-22-gc47ae57` (not `dist-editable:...`) whether driven directly or via the
prototype's `subprocess.run` call — `PYTHONPATH` must point at the worktree's `src` since
`fragments_h5` is installed editable against the *main* tree in `biomarker_env`; getting it
wrong silently re-tests the main tree, not the mutant, with no error. Two simultaneous
meta-test runs (uuid'd `.mutation-tmp/<id>` paths) both passed, no leftover worktrees after,
~8% wall-time slowdown from CPU contention only.

**On `pytest.mark.slow` as the answer to speed:** with no CI (verified: no
`.github/workflows`, no `.gitlab-ci.yml`, no Jenkinsfile), marking it slow/excluded-by-default
reintroduces the review's own failure mode one step removed — "remember not to exclude the
marker" instead of "remember to run the script." Recommend including it in the default
`pytest` run at ~1.5× suite time; revisit only if that becomes a real complaint, as an
explicit later decision, not a default now.
