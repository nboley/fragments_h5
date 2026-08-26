# Build version provenance: labelling which oracle `_build_code_revision` trusted

Status: **implemented** (branch `build-revision-provenance`, commits `0f787af`..`bb6f0d9`).
Scope: `_build_code_revision` resolver, `_build_version` removal, and clean-tree gate.
`_build_argv` is correct and out of scope.
Not related to `docs/pending/build_provenance_metadata.md`, which is REJECTED — do not extend that design.

## Problem

2.12.0 writes `_build_version` from `importlib.metadata.version("fragments-h5")`, which reads
**installed dist-info, not the working tree**. Under an editable install — the standard dev setup — those
diverge silently. Reproduced on this machine today:

```
$ .../biomarker_env/bin/python -c "import importlib.metadata as m; print(m.version('fragments-h5'))"  # 2.11.0
$ grep -m1 -i "^version" pyproject.toml                                                # version = "2.12.1"
```

`__editable__.fragments_h5-2.11.0.pth` points at the live tree (`fragments_h5.__file__` →
`.../src/fragments_h5/__init__.py`), so imports execute 2.12.1 code; the `fragments_h5-2.11.0.dist-info/` snapshot
is frozen at install time (METADATA mtime 3 days older than pyproject.toml) and never rewritten. A file built from
the current clean tree recorded `_build_version = '2.11.0'` — two releases stale on a *clean* tree: not dirty-tree
drift, just `pip install -e .` not re-run.

Same label-disagrees-with-artifact defect that motivated provenance: `v2.5.0` declares 2.3.0, `v2.6.0` declares
2.5.0, `v2.10.1` declares 2.10.0, and image `ghcr.io/nboley/fragments-h5:2.10.1` ships `--min-mapq`, a flag absent
from `git show v2.10.1:src/fragments_h5/main.py` (added later in `f45b6f6`). A field that can be wrong invites trust.

## Current state

- Write: `fragments_h5.py`, the `f.attrs["_build_version"] = importlib.metadata.version(...)` block —
  unconditional, guarded only against `PackageNotFoundError` (absence = unknown).
- Read: `FragmentsH5.__init__`, the `.attrs.get("_build_argv")` / `.attrs.get("_build_version")` pair → `None` when absent.

**Live demonstration, measured 2026-08-24 in the development environment:** `pyproject.toml` declares
`2.12.1`, `git describe` reports `v2.12.1-6-g020fb8d`, and `importlib.metadata.version("fragments-h5")`
returns **`2.11.0`** — the editable install's dist-info is frozen two releases back. An h5 built on this
machine right now stamps itself `_build_version = "2.11.0"`. This is not a hypothesis about editable
installs; it is the current state of the machine these releases were cut from.
- Version is declared once, statically, at `pyproject.toml:7`. No setuptools-scm, no `__version__` in `src/`, no
  `FORMAT_VERSION` constant, and nothing anywhere branches on `_build_version`.
- Measured git availability: `git describe --tags --always --dirty` succeeds inside the repo and returns a
  value of the form `v2.12.1-2-g704a630`. **Every `v2.12.1-…` string in this document is an illustrative
  snapshot of the format, not a current value** — it was accurate when written and was already stale a day
  later (the tree reached `v2.12.1-6-g020fb8d`). Read them as shape, never as fact; a document arguing that
  recorded versions go stale should not ask you to trust its own. Outside the repo, e.g. `/tmp`, the command
  exits 128 without raising (measured). `.dockerignore` excludes `.git`, so git is structurally unavailable
  in the container *at runtime*.

**Two failure modes, both fixable — at different times.** Editable install: git works at runtime, so
`_build_code_revision` is correct there directly. Non-editable install (container, `pip install .`): git is absent
at *runtime* but fully available at *image-build* time, on the host, inside the repo (`Makefile:94-97` runs
`docker build` from the checkout), so the revision can be **baked in** rather than resolved at runtime.

## Recommendation (as originally approved)

> Keep `_build_version` exactly as is.

**Superseded (2026-08-25, commit `2d71995`).** The user decided after an EM critical review to
stop writing `_build_version` entirely. A locally built h5 carried
`_build_version='2.11.0'` (stale dist-info) right next to a correct
`_build_code_revision='git:v2.12.1-11-g0f787af'` — two provenance fields disagreeing, where the
one with the more authoritative-sounding name was the less trustworthy one. `_build_version` is
now read-only: `FragmentsH5.build_version` still returns the value for files written by 2.12.0
and 2.12.1, and returns `None` for newly built files. The reasoning below for *adding*
`_build_code_revision` remains valid; only the recommendation to *also keep writing*
`_build_version` is withdrawn.

Add one new flat attribute, `_build_code_revision` (accessor
`build_code_revision`) — *revision*, not *version*, so it cannot be confused with the adjacent, less
trustworthy `_build_version`: a single **self-labeling** string naming both the value and its oracle.

Resolution order, computed from the package's own `__file__` (not cwd); every git call is
`git -C <dir of fragments_h5.py> ...`:

1. Git, **gated** on the running file being tracked in whatever repo git finds (`ls-files --error-unmatch
   <basename(__file__)>` must exit 0) — ungated, a venv inside the checkout labels a frozen `pip install .`
   snapshot with the *live* tree's describe, and site-packages sitting inside an unrelated repo picks up that
   repo. Then `describe --tags --always --dirty`; "succeeds" means `returncode == 0 and stdout.strip()` (a
   nonzero exit, e.g. rc 128 outside a repo — measured — does not raise and must not be treated as success) →
   `git:v2.12.1-2-g704a630`, plus an `-untracked` suffix when
   `git -C <dir> ls-files --others --exclude-standard ':/'` prints anything.
2. Else, a revision **baked in at image-build time** → `baked:v2.12.1-2-g704a630`. `docker-build` passes
   `--build-arg BUILD_CODE_REVISION="$$(git describe --tags --always --dirty)"` — a **Makefile recipe line**, so
   `$$(` is mandatory: a literal `$(` is a *make* expansion and silently bakes an empty string (measured: `$(` →
   `--build-arg BUILD_CODE_REVISION=`; `$$(` → `...=v2.12.1-2-g704a630`). The Dockerfile declares
   `ARG BUILD_CODE_REVISION` **inside the builder stage** (after `FROM ... AS builder` at `:2`, else invisible)
   and writes `src/fragments_h5/_build_revision.py` (`REVISION = "<value>"`) in a new `RUN` before `:16` —
   `:16-17` is a single `RUN` (`:16` ends in `&& \`), so nothing can be inserted between them. The resolver reads
   it as `try: from ._build_revision import REVISION` / `except ImportError: REVISION = None`, leaving dev trees,
   which have no such module, unaffected. Packaging needs no change: `pyproject.toml:36-37`
   (`[tool.setuptools.packages.find]`, `where = ["src"]`) makes `build_py` sweep every `.py` in the package, so
   the unlisted generated module ships in the wheel, and `MANIFEST.in:5` covers the sdist. Neither setuptools-scm
   nor un-ignoring `.git` is needed.
3. Else, install is not marked editable → `dist:2.12.1` from `importlib.metadata`. Reflects the tree's declared
   version at install time; cannot detect that that tree was dirty or ahead of its tag.
4. Else (editable per PEP 610, but no git) → `dist-editable:2.11.0`: same number as today, flagged as possibly
   stale. Probe: `json.loads(importlib.metadata.distribution("fragments-h5").read_text("direct_url.json"))`
   `["dir_info"]["editable"]` — verified → `True` here. When the file is absent `read_text` returns `None` and
   `json.loads(None)` raises `TypeError` (verified on `pytest`), so wrap it and read *any* exception as "not
   editable", falling to step 3.
5. All fail → omit the attribute; reader sees `None`.

`--dirty` detects modifications to **tracked** files only: in a scratch repo an untracked file still describes as
a clean `v1.0.0`, and this repo describes as `v2.12.1-2-g704a630` while two untracked files sit in it — so a new
untracked module under `src/fragments_h5/` would otherwise get a clean-looking `git:` label on code in no commit,
hence `-untracked`. That suffix must come from `ls-files --others`, **not** `status --porcelain`: status is also
non-empty for *modified tracked* files, so it stamps `-untracked` on a tree with zero untracked files (measured:
one modified `README.md`, nothing untracked → `git:v1.0.0-1-ge6a1fb3-dirty-untracked`), a claim both false and,
`--dirty` having already fired, redundant. `':/'` is mandatory because `ls-files` is cwd-scoped where `status` is
not (measured from `src/pkg` with `docs/NOTES.md` untracked: bare → empty; `':/'` → `../../docs/NOTES.md`).
`dist:`/`dist-editable:` cannot detect a dirty or uncommitted tree at all; the prefix says which guarantee
applies, making a dirty build *labeled*, not reproducible.

**Where it is determined: in the library, not at the CLI boundary.** The `sys.argv` precedent (`main.py:177`;
`fragment_selection_and_build_provenance.md` §5: a library caller would record `['pytest', '-v', ...]`) does
**not** transfer — argv is a property of the *invocation*, known only to the CLI, while code revision is a property
of the *installed package*, derivable from `__file__`; at the CLI the ~30 library/test callers would get nothing.

**Compatibility.** No format version bump — consistent with the provenance decisions in
`fragment_selection_and_build_provenance.md`. 2.12.0/2.12.1 files carry `_build_version` but
lack `_build_code_revision`; post-`2d71995` files carry `_build_code_revision` but lack
`_build_version`. Both are read via `.attrs.get()`, so absence is unremarkable in either
direction.

**Container.** Runtime git is absent, but `make docker-build` runs on the host inside the checkout, so the
`baked:` step closes the gap without un-ignoring `.git`.

> **Original text (superseded 2026-08-25, commit `3eb4c59`).** This section said
> "`Makefile:94-97` has no cleanliness check at all" and prescribed a whole-tree gate
> on `docker-build` only. The gate was implemented first as a per-target check on
> `docker-build` (commit `f054e6a`), then refactored into a shared `require-clean-tree`
> prerequisite guarding `conda-build`, `docker-build`, and `tag` (commit `3eb4c59`).
> The broader scope was a user decision: leaving `conda-build` unguarded would have
> reproduced the exact asymmetry — one artifact type gated, another not — that caused
> the original v2.10.1 incident. `tag` additionally keeps its narrower
> `check-pyproject-clean` prerequisite ordered first, so the tailored diagnostic
> ("`pyproject.toml` has uncommitted changes; refusing to tag") still fires when
> relevant.

**Priority** (lettered, to stay distinct from the numbered resolution order above): **A.** the whole-tree
cleanliness gate on all artifact-producing targets — highest value, lowest cost, and it directly addresses the
motivating defect;
**B.** baking `git describe` in via `--build-arg`; **C.** the resolver itself. All three shipped.

## Rejected alternatives

Replacing `_build_version` with a git SHA (breaks non-editable; changes a released attribute's meaning);
setuptools-scm (contradicts `AGENT_CONTEXT.md:960` and the just-merged single-version-source work); parsing
`pyproject.toml` at runtime (not shipped into site-packages — `Dockerfile:12` copies then installs);
`__version__` in `__init__.py` (a second declaration site, still blind to a dirty tree); a structured/nested
provenance object (over-engineered per the accepted design; stays flat).

**A new attribute holding only a git SHA, absent when git is unavailable** — the obvious simpler design, and
the one to beat. Rejected because it is absent in exactly the case that matters most: `.dockerignore:2`
excludes `.git`, and the runtime image has no git binary, so every container-built h5 — the artifacts whose
provenance is least traceable and which caused the motivating incident — would carry nothing at all. The
prefix scheme exists so that case yields `baked:` rather than silence. A field that is empty precisely when
you most need it is not a simpler solution to the same problem; it is a solution to a different, easier one.

## Implementation

Single phase; no ordering dependencies.

1. `src/fragments_h5/fragments_h5.py` — add `import subprocess` (not currently imported; the stdlib imports are
   at `:140-146`) and a module-level `_resolve_build_code_revision() -> str | None` implementing the five-step
   order above. Step 1 is now three `subprocess.run(["git", "-C", os.path.dirname(__file__), ...],
   capture_output=True, timeout=5, text=True)` calls — `ls-files --error-unmatch`, `describe`, then
   `ls-files --others` — successful only when `returncode == 0` (and non-empty `stdout` for `describe`;
   `--others` legitimately returns empty). **Wrap step 1 in its own `try/except Exception`**, not merely the whole
   resolver: a missing `git` binary raises `FileNotFoundError` from `subprocess.run`, so a single outer catch
   unwinds past steps 2-4 and returns `None`, discarding `baked:` — the very container case step 2 exists for.
   Measured with `git` off `PATH`: step-1-local catch → `baked:v2.12.1-2-g704a630` when a baked module is present,
   `dist-editable:2.11.0` when not; single outer catch → `None` in **both**. The outer catch remains as a backstop
   but must never be the path a missing binary takes. The three calls together measure ~0.09s (7-run median
   0.0890s, range 0.0858-0.0947s; individually 0.0084 / 0.0681 / 0.0144s), once per built file — negligible.
2. `fragments_h5.py`, the `try/except PackageNotFoundError` block that writes
   `f.attrs["_build_version"] = importlib.metadata.version("fragments-h5")` — **removed** (commit
   `2d71995`; see "Recommendation" above). The original proposal said "leave it untouched"; that was
   superseded. A conditional write of `_build_code_revision` when the resolver returns non-`None`
   was added in its place, mirroring the `f.attrs["_build_argv"]` write immediately above it.
3. `fragments_h5.py`, in `FragmentsH5.__init__`, the line
   `self.build_version = self._f.attrs.get("_build_version")` — append
   `self.build_code_revision = self._f.attrs.get("_build_code_revision")` immediately after it.

   > **Locate both sites by grepping for the quoted attribute names, not by line number.** These
   > references were accurate when written and had drifted by 2 and 50 lines within a day, purely from
   > unrelated commits on a sibling branch — a reviewer duly flagged them as stale while they were still
   > correct on `main`. Line numbers in this repo rot faster than the prose around them; a review of an
   > earlier design here raised the identical finding. Symbol names do not rot.
4. `docs/architecture/fragment_selection_and_build_provenance.md` §5 — add a short subsection recording that
   revision determination is deliberately library-side, unlike argv, with the reason above.
5. `Makefile` / `Dockerfile` — in scope here: a `--build-arg` plus a generated module is bounded, not
   a release-process redesign. This step also carries priority **A**, which no other step covers.

   **As implemented (commits `f054e6a`, `3eb4c59`):** a shared `require-clean-tree` prerequisite
   gates on `git diff --quiet HEAD` (tracked, staged and unstaged) **and** an empty
   `git ls-files --others --exclude-standard` (untracked), aborting otherwise. It guards
   `conda-build`, `docker-build`, and `tag` — broader than the original proposal, which
   specified only `docker-build`. See the "Container" section above for the rationale.

## Testing

Extend `tests/test_build_provenance.py` (currently 6 tests):

- Prefix contract: build a file; assert `_build_version` is **absent** (no longer written per commit `2d71995`;
  the original proposal said "unchanged") and `_build_code_revision` is either
  **absent** (step 5 permits omission — sdist install, no git, no dist-info) or a `str` starting with one of
  `git:` / `baked:` / `dist:` / `dist-editable:`. As an unconditional key lookup it errors on that legal case.
- `git:` branch: assert the prefix, skipping on **the gate itself** (`ls-files --error-unmatch <basename>` exiting
  nonzero), *not* on `.git` presence — they diverge exactly where the gate matters. A venv/site-packages copy
  inside the checkout has a `.git` above it, so a `.git`-keyed skip does not fire, yet the gate rightly rejects
  and the resolver returns `dist-editable:`, failing `startswith("git:")` spuriously (measured in a scratch repo:
  frozen copy → gate exit 1 → `dist-editable:2.11.0`; live tracked tree → `git:v1.0.0-1-g5a6c241-dirty-untracked`).
  A skip is needed at all *not* for CI (this repo has none — no `.github/`, travis, circle, gitlab, jenkins or
  buildkite config) but because a non-editable install or an sdist has no git-tracked package dir.
- Gate: point the resolver at a package copy inside a scratch repo but *not* added to it; assert the result does
  **not** start with `git:`. Without this the gate ships unverified.
- `baked:` branch: force git to fail, inject `_build_revision` with `REVISION = "vX"`, assert `== "baked:vX"` —
  the newest, least-reviewed step, and one no dev machine reaches by default.
- `dist:` branch: force git to fail **and** force non-editable (monkeypatch the `direct_url.json` probe); assert
  `startswith("dist:")` exactly. Failing git alone lands on `dist-editable:` because this install *is* editable
  (verified: real tree → `git:v2.12.1-2-g704a630-untracked`; `-C /tmp` and a missing git binary both →
  `dist-editable:2.11.0`), so a loose `startswith("dist")` passes while testing the wrong branch.
- `dist-editable:` branch: force git to fail only; assert `startswith("dist-editable:")` — the branch this machine
  actually lands on, and the one a loose `dist` assertion conflates with step 3. Every branch test must force its
  own branch rather than read the ambient environment, else `git:` is exercised only on a dev editable checkout
  and the rest nowhere, shipping three of four prefixes of a trustworthiness feature untested.
- Total failure: force git to fail and `importlib.metadata.version` to raise `PackageNotFoundError`; assert the
  attribute is omitted and `build_code_revision is None`.
- Seam: `subprocess.run` sits inline in a module-level function, so the seam is
  `monkeypatch.setattr(fragments_h5.fragments_h5.subprocess, "run", ...)`; it intercepts every subprocess call in
  the test, not just git's, which is acceptable. The resolver must run in the parent process, and does: the write
  site `:1188-1191` is in the post-worker merge block, so the patch is not lost across a fork.
- Old-file degradation: alongside `test_file_without_provenance_opens:94`, assert a file with
  `_build_code_revision` deleted opens and exposes `None` (empirically confirmed for the existing two attrs).
- `test_package_not_found_omits_build_version` — **removed** (commit `2d71995`). It existed solely
  to cover the `try/except PackageNotFoundError` write block, which was deleted in the same commit.
  The broader property it tested (the resolver must not crash the build on `PackageNotFoundError`)
  remains covered by `test_revision_total_failure_omits_attribute`.

Baseline, measured on branch `build-revision-provenance` at `bb6f0d9` (all implementation and
this design doc's corrections applied):
→ **`92 passed, 2 warnings, 0 errors`**. The 6 CLI-script errors from the original baseline
(console script not on PATH) were resolved by the implementation work. "No regression" means
these counts hold; `tests/test_docker_build.py` has zero `def test_` functions — a standalone
manual script, not a pytest suite — so it contributes nothing regardless.

### Coverage gap found post-review (2026-08-25)

A coverage gap was found *after* implementation review passed: replacing the
`f.attrs["_build_code_revision"] = _code_rev` write with `pass` — leaving the
`if _code_rev is not None:` guard intact — left all 17 provenance tests green. Every test either
drove the resolver directly (the `test_revision_*_branch` tests) or asserted the `None` case
(backward-compat tests for files lacking the attribute), and none pinned that a real
`build_fragments_h5()` call actually writes the attribute to the output file.

Closed in commit `bb6f0d9` by `test_revision_written_into_real_build`, which monkeypatches
`_resolve_build_code_revision` to a fixed sentinel, runs a real build, and asserts the sentinel
appears in the output file's attrs. The monkeypatch avoids a rotting `git:...` string (this
project has been bitten by hardcoded references three times).

This is the sixth such gap found on this project where a write-site mutation left the suite
green. Every one was caught by mutation testing (manual or automated), not by reading the tests.

## Risks

- `--dirty`/`-untracked` reflect the whole repo, so an edited README marks a build dirty; false positives are the safe direction, accepted. A missing `git` binary or slow/NFS repo degrades via the timeout and broad catch.
- `dist-editable:` depends on PEP 610 `direct_url.json`, which pip writes but not every installer guarantees;
  absent it, an editable install without git degrades to a plain `dist:` — the status quo.
- The prefix vocabulary is a new micro-convention; nothing branches on it today and nothing should.
