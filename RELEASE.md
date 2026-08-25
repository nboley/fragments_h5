# Release Guide for fragments-h5

This guide explains how to build and push Docker images and packages after making changes.

## Prerequisites

1. **Docker** installed and running
2. **GitHub CLI** (`gh`) installed and authenticated:
   ```bash
   gh auth login
   ```
3. **Conda** installed (for conda package builds)

## Current Version

The version is automatically read from `pyproject.toml` (currently **2.12.1**).

## Changelog

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

## Building and Pushing Docker Image

The Docker image will be pushed to `ghcr.io/nboley/fragments-h5:2.11.0` and `ghcr.io/nboley/fragments-h5:latest`.

### Quick Command
```bash
cd /home/nathanboley/src/fragments_h5
make push
```

This will:
1. Build the Docker image locally
2. Authenticate with GHCR using GitHub CLI
3. Tag the image with version and `latest`
4. Push both tags to GHCR

### Step-by-Step

1. **Build Docker image:**
   ```bash
   make docker
   ```
   This creates: `fragments-h5:2.7.0` and `fragments-h5:latest`

2. **Push to GHCR:**
   ```bash
   make push
   ```
   This authenticates, tags, and pushes to `ghcr.io/nboley/fragments-h5:2.11.0`

### Custom Configuration

Override defaults with environment variables:
```bash
GITHUB_USER=your-org make push  # Use different GitHub org/user
VERSION=2.6.0 make push         # Override version (defaults to pyproject.toml)
```

## Building Conda Package

### Build Only
```bash
make conda
```
Output: `conda-build-output/`

### Build and Upload to JFrog (if configured)
```bash
./scripts/build_conda_package.sh --upload
```

Requires environment variables:
- `JFROG_URL` - Artifactory URL
- `JFROG_REPO` - Repository name
- `JFROG_ACCESS_TOKEN` or `JFROG_USER`/`JFROG_PASSWORD`

## Creating Git Tag

After building and pushing, create a git tag:
```bash
make tag
```

This creates and pushes tag `v2.5.0` to the repository.

## Complete Release Workflow

To do everything at once:
```bash
# 1. Build conda package
make conda

# 2. Build and push Docker image
make push

# 3. Create and push git tag
make tag
```

Or use the `all` target (builds conda + docker + push, but not tag):
```bash
make all
make tag  # Still need to tag separately
```

## Verification

After pushing, verify the Docker image:
```bash
docker pull ghcr.io/nboley/fragments-h5:2.11.0
docker run --rm ghcr.io/nboley/fragments-h5:2.11.0 build-fragments-h5 --help
```

This is not optional flourish: an image can be built from a tree ahead of its
git tag (see the v2.11.0 correction note above). Treat the tag as a label, not
evidence — confirm what a container contains by running it.

## Troubleshooting

- **Docker push fails**: Ensure `gh auth login` is completed
- **Version mismatch**: `pyproject.toml` is the single source of truth — the conda recipe receives the version at build time via `--variant pkg_version=$(VERSION)`, and all artifact-producing targets (`conda-build`, `docker-build`, `tag`) depend on `require-clean-tree`, which refuses to proceed when the working tree has tracked or untracked changes. `tag` additionally has a `check-pyproject-clean` prerequisite with a tailored diagnostic. Verify after tagging: `git show v<VERSION>:pyproject.toml | grep '^version'`.
- **Conda build fails**: Ensure conda-forge and bioconda channels are available
