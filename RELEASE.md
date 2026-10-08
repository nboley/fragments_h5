# Release Guide for fragments-h5

This guide explains how to build and push Docker images and packages after making changes.

## Prerequisites

1. **Docker** installed and running
2. **GitHub CLI** (`gh`) installed and authenticated:
   ```bash
   gh auth login
   ```

Conda packaging was **retired 2026-10-07**. The `conda-build` / `conda` / `conda-login`
Makefile targets and this guide's conda sections are gone. Docker + git tag are the only
release artifacts.

## Current Version

The version is automatically read from `pyproject.toml` (currently **2.16.0**).

This line was stale at **2.12.1** through the v2.13.0–v2.14.0 releases. It is a hand-edited
duplicate of a value that `pyproject.toml` already owns, so it rots silently. Treat
`pyproject.toml` as the only source of truth and verify after tagging with
`git show v<VERSION>:pyproject.toml | grep '^version'`.

## Changelog

### v2.16.0 (2026-10-08)

**Changed:**
- `build_fragments_h5` no longer reads the reference for chunks that hold no fragments. The
  three fragment generators (`bam_to_fragments`, `single_end_bam_to_fragments`,
  `tsv_to_fragments`) now build the GC cumsum on the first fragment instead of before
  iterating. **Output is byte-identical:** empty chunks already wrote nothing and were
  dropped before the merge. Verified on the RD-56804 simulator BED (7/7 datasets) and on
  `fragments.bam` chr21 (8/8 datasets, 1,626,404 fragments).
- **One observable behaviour change, and the reason this is a minor and not a patch bump:**
  a FASTA *fetch* error now surfaces only in a chunk that has a fragment. Opening the FASTA
  is still eager, so a missing or unreadable file still raises. This is a relaxation —
  nothing that previously succeeded now fails.

**Performance** (warm page cache; cold could not be measured without root):
- Sparse input (322 fragments, 1 occupied chunk of chr1's 25): **7.1 s → ~1 s**.
  `get_g_or_c_cumsum` calls **25 → 1**; it was 89% of the build.
- Dense input (chr21): no change expected — every chunk is occupied. Wall-clock runs were
  too noisy to resolve the per-fragment `is None` check (estimated <0.1%); outputs identical.
- A fragment reaching into the next chunk still makes that chunk load its FASTA (~0.3 s);
  the caller drops it afterwards.

**Testing:** +6 tests in `tests/test_empty_chunk_skip.py`. They fail on an eager revert, on a
`not gc_offset` sentinel (a missing FASTA contig returns `(None, 0)` and must count as
loaded), and on removing GC. Region queries across skipped chunks are checked too.

### v2.15.0 (2026-10-07)

**Changed:**
- `FragmentsH5.has_methyl`, `has_strand`, `has_gc` and `has_fragment_end_clipped` are now
  `functools.cached_property` instead of `property`. Each answer is computed at most once
  per open handle. The `any(...)` bodies are byte-identical, so **what every property
  returns for every file is unchanged**. Production diff is one import plus four decorator
  lines.
- **One observable behaviour change, and the reason this is a minor and not a patch
  bump:** on a *closed* handle, an answer that was already computed now returns instead of
  raising. Pre-change all four raised `ValueError: Invalid group (or file) id` (measured
  directly). An answer that was never computed still raises. This is a strict relaxation —
  nothing that previously succeeded now fails — and matches `contig_lengths`,
  `max_fragment_length` and `n_fragments`, which were already readable after `close()`.

**Performance** (measured on a 497 MB, 195-contig production h5; 300 per-region `chr1`
fetches; h5py call counts quoted because they do not depend on host load):
- 300 `fetch_array` calls: **5.49–5.98 s → 1.20–1.33 s**
- `has_methyl`: **300 calls / 6.369 s cumulative → 1 call / 0.231 s**
- `Group.__contains__`: **59,700 calls → 1,096**; `Group.__getitem__`: **61,500 → 2,896**
- `read_direct` (the actual data read) unchanged at 900 calls, 0.607 → 0.580 s, and is now
  the dominant cost, as it should be
- `fetch_array` consults `has_methyl` and `has_fragment_end_clipped` on every call, because
  their `return_*` arguments default to `None` meaning "ask the file" — that is why a caller
  doing per-region fetches previously paid a full per-contig scan per region.

**Rejected on measurement — do not reintroduce:** an eager single-pass scan in `__init__`.
It looks tidier and is slower. A combined pass must resolve the worst-case flag, which
destroys the per-flag short-circuit that `has_gc`/`has_strand`/`has_fragment_end_clipped`
rely on (true at contig #1), because `has_methyl` is typically false and walks every contig.
Measured: open cost **20 ms → 256 ms**, and `open + has_gc only` **41 ms → 256 ms**, on a
195-contig file. Lazy per-flag caching strictly dominates it in every scenario measured.

**Testing:** suite **191 → 212 passed** / 3 skipped. The implementation never changed after
the first commit; the tests were wrong three times. Five mutants survived successive
versions of the suite and were each found by *execution*, not review — `has_methyl` hardcoded
false; `has_strand` hardcoded true; dropping only the legacy `len(shape) == 1` strand guard;
`any()` → `all()` on all four; and `any(...)` → `any([...])`, which returns the identical
value for every file while visiting every contig and so silently removes the short-circuit.
Root cause each time was fixture coverage, not assertion logic: no fixture disagreed with
the mutant. Four fixtures were added to force disagreement (`methyl_h5_path` via the `YM`
tag, `no_strand_h5_path`, `two_bit_strand_h5_path` and `non_uniform_h5_path`, the last two
forged with raw h5py because the builder cannot emit those layouts), plus per-flag
`*_differs_between_fixtures` guards and a first-scan cost assertion.

**Known and unchanged:** the `any()` semantics remain a live hazard — a file carrying a
dataset on some contigs but not others reports `True`, and a consumer that then iterates
every contig gets a `KeyError` (see `AGENT_CONTEXT.md` §7.3). Caching did not cause or
worsen this; the answer is identical, merely computed once. The behaviour is now pinned by
`test_has_properties_use_any_not_all_semantics`, so changing it to `all()` requires a
deliberate test change rather than happening by accident.

### v2.11.0 (unreleased)

**Added:**
- `--se-max-fragment-length` CLI flag: maximum fragment length filter for single-end mode.
  Required with `--single-end` for BAM input. Range: 1–65535.
- `--min-mapq` CLI flag: minimum mapping quality filter (default: 0 = keep all).
- TSV/BED safety: `--single-end`, `--se-max-fragment-length`, and `--min-mapq` are each
  warned about and neutralized for TSV/BED input (BAM-only flags).
- CLI validation: range checks, mutual requirement of `--single-end` and
  `--se-max-fragment-length` (BAM only).
- First CLI-level tests; mutation-verified coverage for the SE filter gate.
- Build provenance: `_build_argv` and `_build_code_revision` h5 attributes record
  CLI arguments (JSON) and a self-labeling code revision string. Exposed as
  `FragmentsH5.build_argv`, `.build_code_revision`, and `.build_version` (legacy
  read-only). `_build_version` is no longer written to new files (see note below);
  `build_argv` is recorded only for CLI builds — library callers get
  `_build_code_revision` only.
- `numpy>=1.24` dependency floor in `pyproject.toml`, ensuring out-of-range
  uint16 assignment always raises (closes environment-dependent failure mode).

**Changed:**
- Secondary alignments (`is_secondary`) are now excluded in both paired-end
  and single-end filters. Unconditional, no flag. Measured impact on current
  data: zero — 0 secondary alignments in ~61k sampled reads, because
  `bwa-mem2`/`bowtie2` are not given `-a`/`-k`. Re-check trigger: if an
  aligner config ever gains those flags, this becomes material.
- Single-end over-length spans (>`65535`) now raise `ValueError` with contig,
  position, read name, and CIGAR when `se_max_fragment_length` is unset.
  Previously raised an opaque `OverflowError` from inside a multiprocessing
  worker. When `se_max_fragment_length` is set, over-long spans are still
  silently skipped (unchanged behavior).
- `num_mapped` (now `num_mapped_alignments`) no longer halves the BAM index
  alignment count with `// 2`. A single-end contig with exactly one mapped
  read is no longer silently dropped from the output.
- Remote URL detection now uses a generic scheme regex instead of a prefix
  list, covering `gs://`, `ftp://`, etc. in addition to `s3://` and `http(s)://`.

**Fixed:**
- `--read-methyl` help text: corrected "YN tag" to "YM tag" (code always read "YM").
- S3 input: `os.path.abspath` was mangling `s3://b/k.bam` into
  `/cwd/s3:/b/k.bam`. Remote URLs are now left untouched; local paths are
  still absolutized for worker CWD independence.

**Correction (added 2026-08-24):** The merge commit for this release (`aa753c7`)
stated that `build_se_fragment_h5s.nf`'s container "has neither" flag and that
`errorStrategy 'ignore'` made the resulting argparse failure silent, so "those
samples simply produced no h5." Both claims are false, per direct measurement:
the `ghcr.io/nboley/fragments-h5:2.10.1` *image* does have both flags (it was
built from a tree ahead of the `v2.10.1` git tag, which lacks them); and all 48
expected h5 files for the affected project exist in S3, built successfully.
Separately, `build_se_fragment_h5s.config`'s `standard` profile sets
`errorStrategy = 'terminate'` (loud); only the `remote` profile retries twice
then falls back to `'ignore'`. The recurring lesson: a git tag is not evidence
of what a container contains — verify a container by running it, not by
reading the tag it was built from. What the release did correctly: the CLI
flags genuinely were unreachable from `main.py` before this work, and exposing
them was the right fix.

**Tag `v2.10.1` deleted (2026-08-24), local and origin.** It pointed at commit
`dbed0ae` (a merge dated 2026-06-08), which declares `version = "2.10.0"` — no
`v2.10.0` tag ever existed — and whose source cannot build an h5 at all:
`total_bases = sum(a[3] - a[2] ...)` computed `chunk_start - output_contig`,
raising `TypeError: unsupported operand type(s) for -: 'int' and 'str'`. That
accessor drifted when `output_contig` was inserted at tuple index 2; the pack and
unpack sites were both updated correctly and this third reader, 390 lines away,
was not. Verified by execution at `num_processes` 1, 2 and 4.

The SHA is recorded here deliberately. The 48 h5 files referenced above were
built by the `ghcr.io/nboley/fragments-h5:2.10.1` **image**, which works and is
unaffected by the tag's removal; with the tag gone, this note is the only
remaining git-side anchor for that artifact. The tag was deleted because a label
that points at unbuildable source is worse than no label — but the information it
implied is preserved here rather than destroyed.

**Worker-args refactor: `SubBuildArgs` replaces the positional tuple (2026-08-25, branch
`worker-args-refactor`, merged `9430e40`).** `build_sub_fragments_h5` took a single positional
17-element tuple; it now takes a module-scope `@dataclass(frozen=True, slots=True)`,
`SubBuildArgs`, constructed with keyword arguments. Motivation: the tuple shipped a total
failure in the (deleted, see above) `v2.10.1` tag — inserting `output_contig` at index 2 was
correctly reflected at the pack and unpack sites, but not at a *third*, derived reader,
`total_bases = sum(a[3] - a[2] for a in args)`, ~370 lines away, which then computed
`chunk_start - output_contig` (`int - str`), raising `TypeError` on every build at
`num_processes` 1, 2, and 4. Restoring positional access now raises
`TypeError: 'SubBuildArgs' object is not subscriptable` — the defect class is structurally
unreachable, not merely absent. Shipped alongside: `--contig-name-map` test coverage (zero to
seven tests, including the multiprocessing path — it is the flag that makes `output_contig`
differ from `bam_contig`, the exact field whose insertion caused the defect); and the
`target_h5_path` CLI fixture switched from a bare `build-fragments-h5` under `shell=True` to
`sys.executable -m fragments_h5.main`, fixing six tests that had been silently erroring (exit
127) whenever pytest was launched by absolute interpreter path. Known and accepted, not
defects: keyword construction removes ordering errors but not wrong-value binding between six
adjacent booleans; 8 of the 17 fields are per-build invariants resent with every chunk (an
invariant/config split is deferred); `max_tlen=1000` in `single_end_bam_to_fragments` is dead
in the body but must not be removed (a shared call passes it unconditionally). No version bump
— internal-only change, no external callers. See `docs/architecture/worker_args_refactor.md`.
This changelog entry completes a documentation gate that was missed when `9430e40` merged; it
lands here, on `build-revision-provenance`, because that branch already edits this file and the
gap was found while doing so.

**`_build_version` no longer written (2026-08-25).** Decided by the user
after an EM critical review of the 2.12.0/2.12.1 build-provenance work above.
`_build_version` (from installed dist-info) and `_build_code_revision` (from
`git describe`, primarily) could disagree — measured on this machine, a
locally built h5 carried a stale `_build_version` next to a correct
`_build_code_revision`. `_build_code_revision` is now the sole authoritative
field for code identity; `_build_version` is retained read-only
(`FragmentsH5.build_version`) for backward compatibility with files written
by 2.12.0 and 2.12.1, but is no longer written to new files. No format
version bump. See
`docs/architecture/fragment_selection_and_build_provenance.md`'s 2026-08-25
addendum for the full rationale.

**`require-clean-tree` guards all artifact-producing targets (2026-08-25).**
`conda-build`, `docker-build`, and `tag` now all depend on a shared
`require-clean-tree` Makefile prerequisite. It refuses to proceed when the
working tree has tracked changes (staged or unstaged) or untracked files.
`tag` additionally keeps its `check-pyproject-clean` prerequisite (ordered
first for a tailored diagnostic). This closes the asymmetry where some
artifact types were gated and others were not — the root cause of the
v2.10.1 image/tag mismatch.

### v2.7.0 (2026-02-11)

**Fixed:**
- **Multiprocessing hang with small BAMs**: Switched from `forkserver` to `fork` start method
  - Previous implementation caused race conditions when workers completed before forkserver initialization
  - Fork is safe here because output HDF5 opened after all workers complete
  - Added stress tests with 8 workers on single-contig BAMs to prevent regression

**Added:**
- Multiprocessing stress tests with timeouts (`test_multiprocessing_with_small_bam`, `test_multiprocessing_stress_test`)
- pytest-timeout dependency for catching hangs in tests
- Default 300-second timeout for all tests

**Changed:**
- Multiprocessing start method: `forkserver` → `fork`
- Updated documentation to explain fork safety

## Complete Release Workflow

**Tag first, then build the image. The order is not cosmetic.**

```bash
# 1. Commit the version bump in pyproject.toml (nothing below will run otherwise)
# 2. Create and push the git tag
make tag

# 3. Build and push the Docker image
make docker-push

# 4. Verify by RUNNING the image -- see Verification below
```

Or in one step, which orders this correctly for you:
```bash
make all   # = login tag docker clean
```

### Why tag before build

`docker-build` bakes `BUILD_CODE_REVISION` from
`git describe --tags --always --dirty`. If you build *before* tagging, the newest tag
reachable is the *previous* release, so an image labelled `2.15.0` self-reports
`v2.14.0-7-gcb6be7a`. The artifact then disagrees with its own label — which is exactly
how `v2.10.1` went wrong (see the correction note in the changelog above).

Tag first and `git describe` returns `v2.15.0`, so the image reports precisely the release
it is. Verified on the v2.15.0 release: baked revision `v2.15.0`.

Earlier revisions of this guide instructed docker-then-tag. That was wrong, and the
Makefile's own `all` target always disagreed with it.

## Building and Pushing the Docker Image

Images are pushed to `ghcr.io/$(GITHUB_USER)/fragments-h5:$(VERSION)` and `:latest`, where
`VERSION` comes from `pyproject.toml`.

```bash
make docker-build   # build locally only
make docker-push    # build, authenticate to GHCR via `gh auth token`, tag, push
make docker         # both of the above
```

There is **no `make push` target**. Earlier versions of this guide told you to run it;
it never existed.

### Custom Configuration

```bash
GITHUB_USER=your-org make docker-push   # different GitHub org/user
VERSION=2.6.0 make docker-push          # override version (default: pyproject.toml)
```

## Creating the Git Tag

```bash
make tag
```

Creates and pushes `v$(VERSION)`. It refuses to run if the tag already exists, if
`pyproject.toml` is uncommitted (`check-pyproject-clean`), or if the tree has any tracked
change or untracked file (`require-clean-tree`). Those guards are what stop a tag from
pointing at a commit that declares a different version — the `v2.10.1` failure.

An untracked file from an unrelated session is enough to block this. Commit it, stash it,
or gitignore it; do not work around the gate.

## Verification

After pushing, verify the Docker image:
```bash
VERSION=$(grep 'version = ' pyproject.toml | head -1 | sed 's/.*"\(.*\)".*/\1/')

# pull FRESH, so you test what is in the registry and not your local build cache
docker rmi ghcr.io/nboley/fragments-h5:$VERSION 2>/dev/null
docker pull ghcr.io/nboley/fragments-h5:$VERSION
docker run --rm ghcr.io/nboley/fragments-h5:$VERSION build-fragments-h5 --help

# and confirm the image contains the code it claims to:
docker run --rm ghcr.io/nboley/fragments-h5:$VERSION python -c "
import importlib.metadata as md
import fragments_h5._build_revision as br
print('dist version  :', md.version('fragments_h5'))
print('baked revision:', br.BUILD_CODE_REVISION)"
```

This is not optional flourish: an image can be built from a tree ahead of its
git tag (see the v2.11.0 correction note above). Treat the tag as a label, not
evidence — confirm what a container contains by running it.

## Troubleshooting

- **Docker push fails**: Ensure `gh auth login` is completed
- **Version mismatch**: `pyproject.toml` is the single source of truth. Both
  artifact-producing targets (`docker-build`, `tag`) depend on `require-clean-tree`, which
  refuses to proceed when the tree has tracked changes or untracked files. `tag`
  additionally has `check-pyproject-clean` with a tailored diagnostic. Verify after
  tagging: `git show v<VERSION>:pyproject.toml | grep '^version'`.
- **Version appears in more than one place**: `pyproject.toml` owns it, but
  `AGENT_CONTEXT.md` carries two hand-copied `**Version:**` lines and this guide's
  "Current Version" line is a third copy. They rot — that line sat at `2.12.1` from before
  v2.13.0 until v2.15.0. Update all of them, or trust only `pyproject.toml`.
- **`require-clean-tree` blocks on a file you did not create**: on a shared host another
  session may leave an untracked file in the tree. Commit it, stash it, or gitignore it.
  Do not bypass the gate — it is the only thing preventing an artifact from disagreeing
  with the commit it claims to be built from.
