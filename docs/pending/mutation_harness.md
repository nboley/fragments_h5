# Design: a repeatable mutation-testing harness

**Status:** proposal, not implemented. **Verified against:** `main` @ `f9971d0` (2026-08-27).

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
   repo-root-local, to be gitignored). Two simultaneous runs never touch the same files;
   `git worktree add` with distinct paths does not conflict.
3. **The exact false positive above stops being possible.** A linked worktree is a real git
   checkout (`git -C <worktree> ls-files` succeeds), unlike a `cp`. Verified live in a
   throwaway repo this session (see Worked example): the tracked-file gate passes in a
   worktree and fails in a plain copy, for the identical file.

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
    wt = repo_root / ".mutation-tmp" / uuid.uuid4().hex[:8]     # gitignored, repo-local
    subprocess.run(["git", "worktree", "add", "-q", "--detach", str(wt), "HEAD"], check=True)
    try:
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
- **Not executed inside `fragments_h5` itself**, deliberately — `git worktree add` was kept
  out of this shared tree per this task's constraint on writing git commands here. The
  `.mutation-tmp/` path and the real `SubBuildArgs`/SE-gate-shaped mutation are therefore
  reasoned from the validated mechanism, not independently run against this repo's own
  fixtures. That is the honest gap in this design's verification, and the first thing an
  implementer should do is rerun the worked example against a real anchor in this repo.
- The "no traceback provenance check" scope cut (end of the Assertion section) is a real
  limitation, not a false one: a mutation that changes *which* line raises inside the same
  test function, while still landing on an assertion failure, would be reported as caught
  without checking it's the intended assertion. Judged acceptable because it matches the
  granularity of every fix in the Problem section above, none of which needed that.
- Length: 197 lines (`wc -l`), against the well-under-250 target.
