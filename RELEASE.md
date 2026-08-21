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

The version is automatically read from `pyproject.toml` (currently **2.11.0**).

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
- Build provenance: `_build_argv` and `_build_version` h5 attributes record CLI
  arguments (JSON) and the installed package version. Exposed as
  `FragmentsH5.build_argv` and `.build_version`; both return `None` for files
  built before this feature. `build_argv` is recorded only for CLI builds —
  library callers get `_build_version` only.
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

## Troubleshooting

- **Docker push fails**: Ensure `gh auth login` is completed
- **Version mismatch**: `pyproject.toml` is the single source of truth — the conda recipe receives the version at build time via `--variant pkg_version=$(VERSION)`, and the `tag` target refuses to tag when `pyproject.toml` is uncommitted. However, `docker-build` has no such check, so a Docker image *can* be built from a dirty tree with a stale version. Verify after tagging: `git show v<VERSION>:pyproject.toml | grep '^version'`.
- **Conda build fails**: Ensure conda-forge and bioconda channels are available
