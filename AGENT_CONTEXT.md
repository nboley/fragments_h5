# fragments-h5 Agent Context Document

**Last Updated:** 2026-08-21
**Version:** 2.11.0
**Project Location:** `/home/nathanboley/src/fragments_h5`  
**Repository:** https://github.com/nboley/fragments_h5

---

## 1. Project Overview

### 1.1 Purpose and Domain

**fragments-h5** (PyPI: `fragments-h5`) is a Python library for storing and querying DNA sequencing fragments in HDF5 format. It was developed for the analysis of cell-free DNA (cfDNA) fragments using fragmentomics approaches.

**Key Features:**
- Fast and memory-efficient fragment storage (~10× faster retrieval than BED files)
- Interval query support with bisection algorithm
- Multiprocessing build pipeline (using fork for HDF5 safety)
- tqdm progress bar during build (weighted by contig length)
- S3/GS/HTTP streaming support for BAM and FASTA files (any htslib-supported scheme)
- Optional metadata: GC content, strand, mapping quality, methylation, fragment clipping status
- Build provenance: `_build_argv` (CLI invocations) and `_build_version` recorded in the h5

**Target Use Cases:**
- Cell-free DNA analysis and fragmentomics research
- Fragment-level genomic analysis requiring fast interval queries
- Large-scale genomic pipelines with S3 integration

**Related Paper:** https://www.nature.com/articles/s41416-021-01635-z

### 1.2 Current Status

- **Version:** 2.11.0
- **License:** GPL-3.0-or-later
- **Python Support:** 3.10+
- **Build System:** pip (setuptools + Cython), conda (rattler-build), Docker
- **Deployment:** GHCR (Docker), JFrog Artifactory (conda)

---

## 2. Architecture and Design

### 2.1 Core Data Model

**Storage Strategy:**
- **Starts + Lengths:** Fragment starts stored as `int32`, lengths as `uint16` (not starts + stops)
- **Stops Computed:** `stops = starts + lengths` (calculated on-the-fly)
- **Sorted Starts:** Enables bisection algorithm for fast interval overlap queries
- **Index Structure:** Block-based index (INDEX_BLOCK_SIZE = 5000) for efficient disk access

**HDF5 Layout:**
```
/
├── attrs:
│   ├── index_block_size          [int]
│   ├── max_fragment_length       [int, always 65535]
│   ├── _bam_header               [str, "" for TSV input]
│   ├── _source_format            [str, "BAM" or "TSV"]
│   ├── _contig_lengths_str       [str, Python dict literal]
│   ├── _build_argv               [str, JSON array] (optional, CLI builds only)
│   └── _build_version            [str] (optional, absent if package not installed)
├── data/
│   └── {contig}/
│       ├── starts          [int32, sorted]
│       ├── lengths         [uint16]
│       ├── mapq            [uint8, 2 columns] (optional)
│       ├── gc              [uint8] (optional, if FASTA provided)
│       ├── strand          [char] (optional)
│       ├── methyl_*        [uint8] (optional, if --read-methyl)
│       └── fragment_end_clipped [uint8] (optional, 0/1/255)
├── index/
│   └── {contig}            [uint32, genomic_pos -> array_index mapping]
└── fragment_length_counts  [uint16 histogram, 0 to 65535]
```

**Compression:** HDF5 datasets use gzip level 4

### 2.2 Interval Overlap Algorithm

**Bisection Search:**
```
Goal: Find all fragments overlapping [region_start, region_stop)
Criteria: f_stop > region_start AND f_start < region_stop

With max_frag_len known:
  (q1) Lower bound: searchsorted(starts, region_start - max_frag_len - 1, 'left')
  (q2) Upper bound: searchsorted(starts, region_stop - 1, 'right')

Optimization: Use block index to reduce searchsorted range
  index_lb = INDEX[search_pos // INDEX_BLOCK_SIZE]
  index_ub = INDEX[1 + search_pos // INDEX_BLOCK_SIZE]
  searchsorted(starts[index_lb:index_ub], search_pos, side)
```

**Performance Notes:**
- Majority of time spent in HDF5 array access (not algorithm)
- Uses `h5py.read_direct()` for efficient reads
- Potential 2× speedup with memory-mapped arrays (not implemented due to maintenance cost)

### 2.3 Multiprocessing Architecture

**Build Pipeline (Chunk-Based Parallelization, since v2.8.0):**
1. **Fork:** Uses `fork` start method for minimal overhead (safe because output HDF5 opened after workers complete)
2. **Per-Worker Temp Dirs:** Each worker gets isolated temp directory for S3 index caching
3. **Chunk-Based Workers:** Each contig is split into fixed-size genomic chunks (GENOMIC_CHUNK_SIZE = 10M bases). Each chunk is a work unit processed independently. Fragments are assigned to the chunk containing their start position. Workers fetch only the relevant BAM/FASTA region (with MAX_FRAGMENT_LENGTH buffer for GC).
4. **Main Process Merge:** Chunks are grouped by contig, sorted by chunk_start, and arrays are concatenated to produce the final per-contig datasets.
5. **Signal Handling:** Workers ignore SIGINT/SIGTERM; main process handles cleanup
6. **Skip Chunking:** `--skip-chunking` flag reverts to whole-contig processing (each contig = one work unit)
7. **Single Process:** When num_processes=1, work is executed directly in the main process (no fork overhead)

**Why Fork is Safe Here:**
- Output HDF5 file opened AFTER all workers complete (no open file handles at fork)
- Workers open independent BAM/FASTA file handles
- No shared locks or state during fork
- Each worker writes to isolated temp files

**Previous Implementation Note:**
- Originally used `forkserver` for extra safety
- Caused race conditions and hangs with small BAMs (workers completed before server initialization)
- Switched to `fork` in v2.7.0

**Worker Workflow:**
```python
def build_sub_fragments_h5(args):
    # args = (bam_path, contig, chunk_start, chunk_stop, contig_length, ...)
    # Worker fetches BAM region [chunk_start, chunk_stop) and filters fragments by start position
    # FASTA fetched for region [chunk_start, chunk_stop + MAX_FRAG_LENGTH) for GC calculation
    # Returns (contig, chunk_start, chunk_stop, temp_h5_path)
```

---

## 3. Codebase Structure

### 3.1 Directory Layout

```
fragments_h5/
├── src/fragments_h5/           # Main package
│   ├── __init__.py             # Exports FragmentsH5
│   ├── main.py                 # CLI entry point
│   ├── fragments_h5.py         # Core: FragmentsH5 class, build logic
│   ├── fragment.py             # Fragment dataclass, BAM parsing
│   ├── sequence.pyx            # Cython: GC calculation, one-hot encoding
│   └── _logging.py             # Logging configuration
├── tests/                      # pytest test suite
│   ├── data/                   # Test BAMs, FASTA files
│   ├── test_fragments_h5.py    # Core library tests
│   ├── test_fragment.py        # Fragment dataclass tests
│   ├── test_s3_bam.py          # S3 streaming tests (conditional)
│   ├── test_create_duplicate_sam.py
│   ├── test_docker_build.py    # Manual Docker test
│   └── specialized/            # Standalone comparison/integration tests
│       └── compare_chunked_vs_unchunked.py  # Chunked vs unchunked + Docker comparison
├── scripts/
│   ├── build_conda_package.sh
│   └── publish_conda_package.sh
├── conda-recipe/               # rattler-build conda package
│   ├── recipe.yaml
│   ├── variant_config.yaml
│   └── conda_build_config.yaml
├── docs/
│   ├── architecture/
│   │   └── fragment_selection_and_build_provenance.md
│   └── plan_chunk_based_parallelization.md
├── pyproject.toml              # pip package metadata
├── setup.py                    # Cython extension build
├── Makefile                    # Build/release automation
├── Dockerfile
├── environment.yml
├── pytest.ini
├── README.md
└── RELEASE.md
```

### 3.2 Core Modules

#### `fragments_h5.py` (Core Logic)

**FragmentsH5 Class:**
```python
class FragmentsH5:
    def __init__(self, fname, mode='r', cache_pointers=False):
        # cache_pointers currently disabled (no-op)
        
    def fetch_array(self, contig, region_start=None, region_stop=None,
                   max_frag_len=None, return_mapqs=False, return_gc=False,
                   return_strand=False, return_methyl=False,
                   return_fragment_end_clipped=None,
                   filter_to_midpoint_frags=False):
        """Returns (starts, stops, supp_data) tuples"""
        
    def fetch(self, ...):
        """Iterator over Fragment objects (slow, for inspection)"""
        
    def fetch_counts(self, contig, region_start, region_stop):
        """Count fragments in region"""
        
    # Properties: filename, name, has_methyl, has_strand,
    #            has_fragment_end_clipped, max_fragment_length,
    #            fragment_length_counts, n_fragments,
    #            build_argv, build_version, source_format
```

**Build Functions:**
```python
def build_fragments_h5(input_fname, ofname, fasta_filename=None,
                      allowed_contigs=None, set_mapq_255_to_none=False,
                      read_strand=True, read_methyl=False, single_end=False,
                      se_max_fragment_length=None, num_processes=None,
                      include_duplicates=False, store_fragment_end_clipped=True,
                      skip_chunking=False, contig_name_map=None, min_mapq=None,
                      *, build_argv=None):
    """Main entry point for building fragment H5 from BAM or TSV/BED"""

def build_sub_fragments_h5(args):
    """Worker function — processes one chunk (contig, chunk_start, chunk_stop)"""
```

#### `fragment.py` (Data Classes and BAM Parsing)

**Fragment Dataclass:**
```python
@dataclass(frozen=True, slots=True)
class Fragment:
    chrom: str
    start: int
    stop: int
    mapq1: int | None  # -1 for None
    mapq2: int | None  # -1 for None
    gc: float | None = None
    strand: str | None = None
    cell_barcode: str | None = None
    methyl_counts: MethylCounts | None = None
    fragment_end_clipped: int | None = None  # 0/1/255
    
    @property
    def length(self) -> int:
        return self.stop - self.start
    
    @property
    def tlen(self) -> int:
        """Template length (signed)"""
    
    @property
    def midpoint(self) -> int:
        return (self.start + self.stop) // 2
```

**BAM Parsing:**
```python
def bam_to_fragments(bam, contig, start, stop, fasta_file, max_tlen,
                     min_mapq, include_duplicates, ...):
    """Paired-end BAM → Fragment iterator
    Filters: is_qcfail, is_secondary, is_supplementary, is_duplicate (unless included),
             is_unmapped, mate_is_unmapped, tlen==0, abs(tlen)>max_tlen, mapq<min_mapq"""

def single_end_bam_to_fragments(bam, ..., se_max_fragment_length=None):
    """Single-end BAM → Fragment iterator
    Filters: is_qcfail, is_secondary, is_supplementary, is_duplicate (unless included),
             is_unmapped, mapq<min_mapq
    Over-length: if se_max_fragment_length is None and span > 65535, raises ValueError;
                 if set, the skip happens at the chunk level in build_sub_fragments_h5"""
```

#### `sequence.pyx` (Cython Performance-Critical Code)

```python
def one_hot_encode_sequences(sequences: list[str]) -> np.ndarray:
    """One-hot encode DNA sequences for GC calculation
    Returns: [N, 4, L] array (A, C, G, T)"""
    
# Used by fragment.py to calculate GC content efficiently
```

#### `main.py` (CLI Entry Point)

**CLI Arguments:**
- **Required:** `input_bam`, `output_frags_h5`
- **Optional:**
  - `--fasta`: Reference FASTA (local or S3, required for GC content)
  - `--contigs`: Restrict to specific contigs
  - `--set-mapq-255-to-none`: Map MAPQ 255 to None
  - `--exclude-strand`: Don't store strand info
  - `--read-methyl`: Parse methylation from YM tag
  - `--single-end`: Single-end sequencing mode
  - `--se-max-fragment-length`: Max fragment length for single-end mode (required with `--single-end` for BAM input, range 1–65535)
  - `--min-mapq`: Minimum mapping quality filter (default: 0 = keep all, must be >= 0)
  - `--contig-name-map`: TSV file mapping input contig names to output names
  - `--include-duplicates`: Include duplicate-marked reads (default: exclude)
  - `--no-store-fragment-end-clipped`: Don't store clipping info
  - `--num-processes`: Number of workers (default: 1, use 'all' for all cores)
  - `--skip-chunking`: Disable chunk-based parallelization, process each contig as a whole
  - Logging: `--quiet`, `--verbose`, `--debug`, `--log-filename`, etc.

**Validation:**
- Checks if output file already exists (fails early with error)
- Range validation: `--se-max-fragment-length` must be 1–65535; `--min-mapq` must be >= 0
- For BAM input: `--se-max-fragment-length` required with `--single-end`, rejected without it
- For TSV/BED input: `--single-end`, `--se-max-fragment-length`, and `--min-mapq` are warned about and neutralized (these flags are BAM-only)
- Checks BAM index exists (runs `samtools index` if local, fails if remote)
- Validates FASTA accessibility for S3 URLs (requires `.fai`, `.gzi` for compressed)

---

## 4. Dependencies and Build System

### 4.1 Dependencies

**Runtime:**
- `numpy >= 1.24` - Array operations (floor ensures out-of-range integer assignment raises rather than wraps)
- `h5py` - HDF5 file I/O
- `pysam` - BAM/FASTA reading (requires htslib with S3 support for S3 URLs)
- `tqdm` - Progress bars

**Build-time:**
- `setuptools >= 61.0`
- `Cython >= 3.0`
- `numpy` (for Cython extension compilation)

**Development:**
- `pytest`, `pytest-timeout` - Testing framework
- `rattler-build` - Conda package builds

### 4.2 Build Systems

#### Pip (PyPI)

```bash
# Install from PyPI
pip install fragments-h5

# Build from source
python -m build
pip install dist/fragments_h5-*.whl
```

**Build Process:**
1. `setup.py` builds Cython extension `fragments_h5.sequence` from `sequence.pyx`
2. `setuptools` packages wheel with compiled extension

#### Conda (JFrog Artifactory)

```bash
# Build conda package
make conda-build
# or
rattler-build build \
    --recipe conda-recipe/recipe.yaml \
    --output-dir conda-build-output \
    --channel conda-forge \
    --channel bioconda \
    --variant-config conda-recipe/variant_config.yaml

# Publish to JFrog
make conda-publish
# or
./scripts/publish_conda_package.sh
```

**Conda Recipe:** `conda-recipe/recipe.yaml`
- Platform: linux-64
- Python: 3.13 (variant_config.yaml), 3.10-3.13 (conda_build_config.yaml)
- Build: `python setup.py build_ext --inplace && pip install .`

#### Docker (GHCR)

```bash
# Build Docker image
make docker-build

# Push to GitHub Container Registry
make docker-push

# Both
make docker
```

**Image:** `ghcr.io/nboley/fragments-h5:2.11.0` and `:latest`

**Dockerfile Highlights:**
- Base: Python 3.10+
- Installs pysam with S3 support
- Compiles Cython extensions
- Entry point: `build-fragments-h5` CLI

### 4.3 Makefile Targets

| Target | Description |
|--------|-------------|
| `conda-build` | Build conda package with rattler-build |
| `conda-publish` | Upload to JFrog Artifactory (requires credentials) |
| `conda` | Build + publish conda |
| `docker-build` | Build Docker image |
| `docker-push` | Push to GHCR (requires gh auth) |
| `docker` | Build + push Docker |
| `tag` | Create and push git tag v$(VERSION) |
| `all` | conda + docker + tag + clean |
| `clean` | Remove build artifacts |

---

## 5. Testing Infrastructure

### 5.1 Test Framework

**Framework:** pytest  
**Config:** `pytest.ini` (log_level=INFO)  
**Location:** `tests/`

### 5.2 Test Files

| File | Focus |
|------|-------|
| `test_fragments_h5.py` | Core library: build, fetch, fetch_array, counts, cache, regions, pickle, include_duplicates, fragment_end_clipped, output-exists check, multiprocessing stress tests |
| `test_fragment.py` | Fragment dataclass: equality, length, midpoint, mapq, strand normalization |
| `test_s3_bam.py` | S3 BAM/FASTA streaming (skipped when S3 not available) |
| `test_create_duplicate_sam.py` | Generate test data for duplicate tests |
| `test_docker_build.py` | Manual Docker vs local build comparison (not run via pytest) |

### 5.3 Test Data

**Location:** `tests/data/`

| File | Description |
|------|-------------|
| `small.chr6.bam` | Basic test BAM |
| `scATAC_breast_v1_chr6_*.bam` | scATAC test data |
| `test_duplicates.bam` | Duplicate-marked reads |
| `GRCh38.p12.genome.chr6_99110000_99130000.fa.gz` | Reference FASTA segment |
| `.fai`, `.gzi` | FASTA index files |

### 5.4 Key Test Scenarios

**Core Functionality:**
- Build from BAM (with/without FASTA)
- Fetch fragments by region
- Fetch with filters (max_frag_len, midpoint filtering)
- Return supplementary data (mapqs, gc, strand, methyl, fragment_end_clipped)
- Count fragments in region

**Edge Cases:**
- Empty regions
- Single-fragment regions
- Overlapping vs midpoint filtering
- MAPQ 255 handling
- Duplicate read handling (include vs exclude)

**S3 Integration:**
- S3 BAM streaming
- S3 FASTA streaming (gzipped)
- Index file requirements

---

## 6. S3 Integration and Remote Files

### 6.1 S3 Support

**Supported URLs (any scheme htslib recognizes):**
- `s3://bucket/path/to/file.bam`
- `gs://bucket/path/to/file.bam`
- `http(s)://url/to/file.bam`
- `ftp://url/to/file.bam`

The remote URL predicate uses a generic scheme regex (`^[a-z][a-z0-9+.\-]*://`), so any scheme htslib supports now or gains later is handled without code changes.

**Requirements:**
- pysam/htslib built with S3 support (libcurl)
- AWS credentials configured (for non-public buckets)
- Index files must exist at same S3 path:
  - BAM: `.bai` file (e.g., `s3://bucket/sample.bam.bai`)
  - FASTA: `.fai` file (e.g., `s3://bucket/ref.fa.gz.fai`)
  - Compressed FASTA: `.gzi` file (e.g., `s3://bucket/ref.fa.gz.gzi`)

### 6.2 Multiprocessing with S3

**Challenge:** pysam downloads S3 indexes to local temp files, causing conflicts when multiple processes access the same S3 URL.

**Solution:** Per-worker temporary working directories (`docs/plan_chunk_based_parallelization.md`)

```python
@contextmanager
def _temporary_working_directory():
    """Create isolated temp dir for worker to avoid S3 index conflicts"""
    orig_cwd = os.getcwd()
    temp_dir = tempfile.mkdtemp(prefix='fragments_h5_worker_')
    try:
        os.chdir(temp_dir)
        yield temp_dir
    finally:
        os.chdir(orig_cwd)
        shutil.rmtree(temp_dir)
```

**Performance Note:** For genome-wide BAMs, the tool fetches complete chromosome sequences from FASTA, which is efficient (minimizes S3 requests).

---

## 7. Known Issues and Limitations

### 7.1 Disabled Features

**cache_pointers Parameter:**
- **Status:** Currently disabled (no-op)
- **Original Purpose:** Load entire index into memory for faster queries
- **Reason for Disabling:** Memory overhead concerns
- **TODO:** Re-implement with more memory-efficient approach

### 7.2 Current Limitations

**Fragment Length:**
- Maximum: 65535 (uint16 limit)
- In single-end mode without `--se-max-fragment-length`: over-long spans raise `ValueError` naming the contig, position, read, and CIGAR
- In single-end mode with `--se-max-fragment-length`: over-long spans are silently skipped
- In paired-end mode: `max_tlen` filter silently skips reads with `abs(tlen) > MAX_FRAG_LENGTH`

**Mapping Quality:**
- Stored as uint8 (0-255)
- MAPQ 255 can be mapped to None via `--set-mapq-255-to-none`
- Unknown MAPQ stored as -1 in fetch results

**Methylation:**
- Only supports YM tag format
- Limited to uint8 counts (0-255 CpGs per fragment)

**S3 Performance:**
- S3 index files downloaded per worker (temp dir overhead)
- No streaming for index data (entire contig index loaded)

### 7.3 Resolved Issues

**Multiprocessing Test Hang (Fixed in v2.7.0):**
- **Previous Issue:** Tests hung when num_processes > 1 with small BAMs
- **Root Cause:** forkserver initialization race with fast job completion
- **Resolution:** Switched to fork method, added timeout tests
- **Validation:** Stress tests with 8 workers on 1-contig BAM

### 7.4 Build System Quirks

**rattler-build Cleanup:**
- rattler-build may exit with non-zero status during cleanup even when build succeeds
- Makefile checks for presence of `.conda` files to determine success
- Known rattler-build issue, does not affect package quality

**Cython Compilation:**
- `sequence.c` generated from `sequence.pyx` during build
- Not checked into git (excluded in .gitignore)
- Requires Cython 3+ at build time

---

## 8. Usage Examples

### 8.1 CLI Usage

**Basic Build:**
```bash
build-fragments-h5 input.bam output.fragments.h5
```

**With GC Content (requires FASTA):**
```bash
build-fragments-h5 input.bam output.fragments.h5 \
    --fasta reference.fa.gz
```

**S3 Inputs:**
```bash
build-fragments-h5 \
    s3://bucket/sample.bam \
    ./output.fragments.h5 \
    --fasta s3://bucket/reference.fa.gz
```

**Multiprocessing:**
```bash
build-fragments-h5 input.bam output.fragments.h5 \
    --fasta reference.fa.gz \
    --num-processes all
```

**Include Duplicates:**
```bash
build-fragments-h5 input.bam output.fragments.h5 \
    --include-duplicates
```

### 8.2 Python API Usage

**Load and Query:**
```python
from fragments_h5 import FragmentsH5

# Open file
fh5 = FragmentsH5("sample.fragments.h5")

# Fetch all fragments in region (fast)
starts, stops, supp_data = fh5.fetch_array('chr6', 1000000, 2000000)

# With supplementary data
starts, stops, supp = fh5.fetch_array(
    'chr6', 1000000, 2000000,
    max_frag_len=500,
    return_mapqs=True,
    return_gc=True,
    return_strand=True,
    filter_to_midpoint_frags=True
)

# Access supplementary data
mapqs = supp['mapq']  # [N, 2] array of (mapq1, mapq2)
gc = supp['gc']       # [N] array of GC fractions
strand = supp['strand']  # [N] array of '+'/'-'

# Iterate over Fragment objects (slow, for inspection)
for frag in fh5.fetch('chr6', 1000000, 2000000):
    print(f"{frag.chrom}:{frag.start}-{frag.stop} "
          f"mapq={frag.mapq1},{frag.mapq2} gc={frag.gc}")

# Count fragments
n = fh5.fetch_counts('chr6', 1000000, 2000000)

# Properties
print(fh5.filename)
print(fh5.max_fragment_length)
print(fh5.n_fragments)
print(fh5.has_gc)
print(fh5.has_strand)
print(fh5.has_methyl)

# Build provenance (None for files built before this feature)
print(fh5.build_argv)      # list of CLI args, or None
print(fh5.build_version)   # package version string, or None
```

**Build Programmatically:**
```python
from fragments_h5.fragments_h5 import build_fragments_h5

build_fragments_h5(
    'input.bam',
    'output.fragments.h5',
    fasta_filename='reference.fa.gz',
    allowed_contigs=['chr1', 'chr2'],
    set_mapq_255_to_none=False,
    read_strand=True,
    read_methyl=False,
    single_end=True,
    se_max_fragment_length=120,
    num_processes=8,
    include_duplicates=False,
    store_fragment_end_clipped=True,
    min_mapq=30,
    # build_argv is keyword-only; omit for library callers (records _build_version only)
    # CLI passes build_argv=sys.argv automatically
)
```

---

## 9. Integration with Fragmentomics Pipeline

### 9.1 Nextflow Integration

**Related Files:**
- `omni/pipeline/fragmentomics/workflow.nf` - Defines `build_fragments_h5` process
- `omni/pipeline/fragmentomics/specialized_workflows/rebuild_all_frag_h5s.nf` - Batch rebuild workflow
- `omni/pipeline/fragmentomics/tests/functional/test_build_fragments_h5.nf` - Integration tests

**Nextflow Process Example:**
```groovy
process build_fragments_h5 {
    container "ghcr.io/nboley/fragments-h5:2.11.0"
    
    input:
    tuple val(sample_id), path(bam), path(bai), 
          path(fasta), path(fasta_gzi), path(fasta_fai)
    
    output:
    tuple val(sample_id), path("${sample_id}.fragments.h5")
    
    script:
    """
    build-fragments-h5 \\
        ${bam} \\
        ${sample_id}.fragments.h5 \\
        --fasta ${fasta} \\
        --num-processes ${task.cpus}
    """
}
```

### 9.2 Known Integration Issues

**BAM Naming Convention:**
- Original prompt specified `human.mapped.bam`
- Actual implementation uses `mapped.sorted.bam`
- See: `omni/pipeline/fragmentomics/docs/plans/REBUILD_FRAG_H5S_FROM_BAMS_PROMPT.md`

**Index Type:**
- CSI indexes (`.csi`) used instead of BAI (`.bai`)
- Path: `bam_path + '.csi'`

**Output Path:**
- Hardcoded in `rebuild_all_frag_h5s.nf`
- Cannot be overridden by `-profile remote`

---

## 10. Release Process

### 10.1 Version Management

**Version Location:** `pyproject.toml`
```toml
[project]
version = "2.7.1"
```

**Version source:** `pyproject.toml` only. The conda recipe carries no version
literal — it receives one at build time via `--variant pkg_version=$(VERSION)`,
and the Makefile derives `VERSION` from `pyproject.toml`. The git tag is the one
remaining value that must be set to match.

### 10.2 Release Workflow

**Standard Release:**
```bash
# 1. Update version in pyproject.toml
vim pyproject.toml

# 2. Verify credentials
make login

# 3. Build and publish conda package
make conda

# 4. Build and push Docker image
make docker

# 5. Create and push git tag
make tag

# 6. Clean up
make clean
```

**Or use all-in-one:**
```bash
make all
```

### 10.3 Publishing Targets

**Conda:**
- **Target:** JFrog Artifactory (`karius.jfrog.io/artifactory/karius-conda`)
- **Credentials:** Environment variables or pip.conf
  - `ARTIFACTORY_HOST`, `ARTIFACTORY_USER`, `ARTIFACTORY_TOKEN`
- **Script:** `scripts/publish_conda_package.sh`

**Docker:**
- **Target:** GitHub Container Registry (`ghcr.io/nboley/fragments-h5`)
- **Auth:** GitHub CLI (`gh auth login`)
- **Tags:** `2.6.0` and `latest`

**Git:**
- **Tag Format:** `v2.7.1`
- **Remote:** `origin`

---

## 11. Performance Characteristics

### 11.1 Profiling Results (from README.md)

**Test Setup:**
- Benchmark: chr10 from two BAMs (capture and WGS)
- Cache: Data cached in OS (disk I/O not bottleneck)
- Iterations: 2500 queries over multiple regions

**Best Configuration:**
```
Data Model: int32 starts + uint16 lengths + uint8 mapqs/gc
INDEX_BLOCK_SIZE: 5000
CHUNK_SIZE: 50000
Optimization: read_direct + contiguous arrays
```

**Results:**

| Scenario | WGS Time | Capture Time | WGS Size | Capture Size |
|----------|----------|--------------|----------|--------------|
| max_frag_len=511 | 1.45s | 1.71s | 68M | 9.8M |
| max_frag_len=511 + return mapqs + gc | 3.13s | 5.18s | 68M | 9.8M |

**Key Insights:**
- ~10× faster than iterating BED files
- Majority of time in HDF5 array access (not algorithm)
- Potential 2× speedup with memory-mapped arrays (not implemented)
- Smaller max_frag_len improves query speed (~40% faster for 511 vs 65535)

### 11.2 File Size Comparison

**WGS BAM (chr10):**
- Original BAM: ~1-5 GB (estimated)
- Fragment H5: 68M (with mapqs + GC)

**Capture BAM (chr10):**
- Original BAM: ~100-500 MB (estimated)
- Fragment H5: 9.8M (with mapqs + GC)

**Space Efficiency:**
- ~10-20× smaller than BAM
- Stores only fragment-level information (no full reads)

---

## 12. Debugging and Development

### 12.1 Common Debug Workflows

**Run Tests:**
```bash
cd fragments_h5
pytest tests/
pytest tests/test_fragments_h5.py::test_fetch_array -v
```

**Test S3 Support:**
```bash
pytest tests/test_s3_bam.py -v -s
# Skipped if S3 not configured
```

**Manual Testing:**
```bash
# Build test fragment H5
build-fragments-h5 \
    tests/data/small.chr6.bam \
    test_output.fragments.h5 \
    --fasta tests/data/GRCh38.p12.genome.chr6_99110000_99130000.fa.gz \
    --verbose

# Inspect with h5ls
h5ls -r test_output.fragments.h5

# Test in Python
python -c "
from fragments_h5 import FragmentsH5
fh5 = FragmentsH5('test_output.fragments.h5')
print(fh5.fetch_array('chr6'))
"
```

**Debug Multiprocessing:**
```bash
# Single process (easier to debug)
build-fragments-h5 input.bam output.h5 --num-processes 1 --debug

# Multi-process with logging
build-fragments-h5 input.bam output.h5 --num-processes all \
    --debug --log-filename build.log
```

### 12.2 Common Issues

**Issue: "BAM has no index"**
- **Cause:** BAM file not indexed
- **Fix (local):** Tool auto-runs `samtools index`
- **Fix (S3):** Ensure `.bai` exists at same S3 path

**Issue: "Failed to open remote FASTA"**
- **Cause:** Missing FASTA index files
- **Fix:** Ensure `.fai` and `.gzi` (for compressed) exist at same S3 path

**Issue: Multiprocess S3 conflicts**
- **Cause:** Workers sharing S3 index temp files
- **Status:** Fixed via per-worker temp dirs
- **Check:** `docs/plan_chunk_based_parallelization.md`

**Issue: ImportError for `fragments_h5.sequence`**
- **Cause:** Cython extension not compiled
- **Fix:** `python setup.py build_ext --inplace`

**Issue: MAPQ values appear as -1**
- **Cause:** MAPQ 255 stored as None (represented as -1 in fetch results)
- **Fix:** Use `--set-mapq-255-to-none` to explicitly map 255 → None

### 12.3 Logging Configuration

**CLI Verbosity:**
```bash
--quiet    # ERROR and above
--verbose  # INFO and above (default)
--debug    # DEBUG and above
```

**File Logging:**
```bash
--log-filename build.log
--log-file-verbosity-level DEBUG
--log-file-format "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
```

**Python API:**
```python
import fragments_h5._logging as logging
import argparse

parser = argparse.ArgumentParser(parents=[logging.build_log_parser()])
args = parser.parse_args()
logging.configure_root_logger_from_args(args)
```

---

## 13. Development Environment

### 13.1 Setup

**Create Conda Environment:**
```bash
cd fragments_h5
conda env create -f environment.yml
conda activate fragments_h5
```

**Or with pip:**
```bash
python -m venv venv
source venv/bin/activate
pip install -e .[dev]
```

### 13.2 Development Workflow

**Edit Cython Code:**
```bash
vim src/fragments_h5/sequence.pyx
python setup.py build_ext --inplace
pytest tests/
```

**Edit Python Code:**
```bash
vim src/fragments_h5/fragments_h5.py
pytest tests/test_fragments_h5.py -v
```

**Update Tests:**
```bash
vim tests/test_fragments_h5.py
pytest tests/test_fragments_h5.py::test_new_feature -v
```

### 13.3 Pre-Release Checklist

- [ ] Update version in `pyproject.toml` (the only place it is declared)
- [ ] Update `RELEASE.md` if needed
- [ ] Run full test suite: `pytest tests/`
- [ ] Test conda build: `make conda-build`
- [ ] Test Docker build: `make docker-build`
- [ ] Test CLI: `build-fragments-h5 tests/data/small.chr6.bam test.h5`
- [ ] Update `AGENT_CONTEXT.md` (this file) if architecture changed

---

## 14. Related Projects and Context

### 14.1 Fragmentomics Pipeline

**Location:** `/home/nathanboley/src/omni/pipeline/fragmentomics/`

**Relationship:**
- Nextflow pipeline uses `fragments-h5` as a process
- Outputs used by downstream fragmentomics analysis
- Integration tests in `omni/pipeline/fragmentomics/tests/functional/`

**Key Documents:**
- `omni/pipeline/fragmentomics/agent_context/RebuildFragH5sAndRetry.2026-01-29.projectcontext.md`
- `omni/pipeline/fragmentomics/docs/plans/REBUILD_FRAG_H5S_FROM_BAMS_PROMPT.md`

### 14.2 Test Installation Copy

**Location:** `/home/nathanboley/src/test_fragmentomics_tools_install/fragments_h5/`

**Purpose:** Snapshot of fragments_h5 for integration testing with fragmentomics pipeline

---

## 15. Quick Reference

### 15.1 File Locations

| What | Path |
|------|------|
| Main code | `src/fragments_h5/` |
| Tests | `tests/` |
| CLI entry point | `src/fragments_h5/main.py` |
| Core logic | `src/fragments_h5/fragments_h5.py` |
| Fragment class | `src/fragments_h5/fragment.py` |
| Cython extension | `src/fragments_h5/sequence.pyx` |
| Conda recipe | `conda-recipe/recipe.yaml` |
| Docker | `Dockerfile` |
| Build automation | `Makefile` |
| Package metadata | `pyproject.toml` |
| Documentation | `README.md`, `RELEASE.md` |

### 15.2 Key Commands

```bash
# Build conda package
make conda-build

# Publish conda package
make conda-publish

# Build Docker image
make docker-build

# Push Docker image
make docker-push

# Run tests
pytest tests/

# Build from source
python setup.py build_ext --inplace
pip install -e .

# Use CLI
build-fragments-h5 input.bam output.h5 --fasta ref.fa.gz --num-processes all
```

### 15.3 Key Constants

```python
INDEX_BLOCK_SIZE = 5000       # Genomic positions per index block
CHUNK_SIZE = 1000000          # Array resize chunk size
MAX_FRAGMENT_LENGTH = 65535   # uint16 limit
GENOMIC_CHUNK_SIZE = 10000000 # 10M bases per parallelization chunk
```

### 15.4 Key Classes and Functions

**Public API:**
- `FragmentsH5` class
- `Fragment` dataclass
- `build_fragments_h5()` function

**Internal:**
- `bam_to_fragments()` - BAM → Fragment iterator
- `single_end_bam_to_fragments()` - Single-end BAM → Fragment iterator
- `one_hot_encode_sequences()` - Cython GC calculation

---

## 16. Next Steps and TODOs

### 16.1 Known TODOs (from code comments)

**cache_pointers:**
- Re-implement with memory-efficient approach
- Original version disabled due to memory overhead

**num_processes validation:**
- Add user-friendly error message for invalid values
- Current: converts to int without validation

**Test Coverage Gaps:**
- MethylCounts and YM tag parsing (fragment.py)
- Single-end BAM processing end-to-end (basic SE filter gating is now tested)

### 16.2 Potential Improvements

**Performance:**
- Investigate memory-mapped arrays for 2× speedup
- Optimize S3 index caching across workers
- Reduce HDF5 read overhead

**Features:**
- Support CSI indexes (in addition to BAI)
- Support alternative methylation tag formats

**Usability:**
- Validate FASTA/BAM compatibility (contig names)
- Add dry-run mode to estimate output size

**Testing:**
- Expand S3 test coverage
- Add performance regression tests
- Test with real-world large BAMs

---

## 17. Contact and Support

**Maintainer:** Nathan Boley (npboley@gmail.com)  
**Repository:** https://github.com/nboley/fragments_h5  
**Issues:** https://github.com/nboley/fragments_h5/issues

---

### Unreleased (fragment-selection-and-provenance branch)
- **Secondary alignment exclusion:** `is_secondary` now excluded in both paired-end (`bam_to_align`) and single-end (`single_end_bam_to_fragments`) filters. Unconditional, no flag. Measured impact on current data: zero (0 secondary alignments in ~61k sampled reads across fixtures and two production BAMs, because `bwa-mem2`/`bowtie2` are not given `-a`/`-k`). If an aligner config ever gains `-a`/`-k`, this becomes material.
- **SE over-length span raises `ValueError`:** When `se_max_fragment_length` is unset, a single-end read whose reference span exceeds 65535 now raises `ValueError` with contig, position, read name, and CIGAR — instead of an opaque `OverflowError` from inside a multiprocessing worker. When `se_max_fragment_length` is set, over-long spans are still silently skipped.
- **`num_mapped` fix:** `num_mapped_alignments` (formerly `num_mapped`) no longer halves the alignment count with `// 2`. A single-end contig with exactly one mapped read is no longer silently dropped.
- **S3 input fix:** `os.path.abspath` was mangling `s3://b/k.bam` into `/cwd/s3:/b/k.bam`. Remote URLs are now detected by a generic scheme regex and left untouched. Also covers `gs://`, `https://`, `ftp://`, etc.
- **Build provenance:** New h5 attributes `_build_argv` (JSON array, CLI builds only) and `_build_version` (package version string). Exposed as `FragmentsH5.build_argv` and `FragmentsH5.build_version`. Both return `None` on files that predate this feature. `build_argv` is keyword-only in `build_fragments_h5()`; library callers get `_build_version` only.
- **`numpy>=1.24` floor:** Added to `pyproject.toml` dependencies. Ensures out-of-range uint16 assignment always raises, closing the last environment-dependent failure mode.

### v2.11.0 Changelog
- **`--se-max-fragment-length` CLI flag:** Maximum fragment length filter for single-end mode. Required with `--single-end` for BAM input. Range: 1–65535 (uint16 `lengths_arr` limit). Fragments with alignment span exceeding this value are excluded.
- **`--min-mapq` CLI flag:** Minimum mapping quality filter. Fragments below this threshold are excluded at build time. Default: 0 (no filtering). Range: >= 0.
- **TSV/BED input safety:** `--single-end`, `--se-max-fragment-length`, and `--min-mapq` are each individually warned about and neutralized when input is TSV/BED, since these flags are BAM-only concepts.
- **CLI validation:** Range checks for both new flags; `--se-max-fragment-length` requires `--single-end` and vice versa (BAM only).
- **Test coverage:** First CLI-level tests in this repo; SE filter gate at `fragments_h5.py:782` now has mutation-verified test coverage.
- **Pre-existing fix:** `--read-methyl` help text corrected from "YN tag" to "YM tag" (code always read "YM").

**Correction (added 2026-08-24):** The merge commit for this release (`aa753c7`)
claimed the `build_se_fragment_h5s.nf` container "has neither" flag and that
`errorStrategy 'ignore'` made the resulting argparse failure silent, so
"those samples simply produced no h5." Both claims are false, per direct
measurement: `ghcr.io/nboley/fragments-h5:2.10.1` (the image) has both flags —
it was built from a tree ahead of the `v2.10.1` git tag, which lacks them —
and all 48 expected h5 files for the affected project exist in S3, built
successfully. Also, `build_se_fragment_h5s.config`'s `standard` profile sets
`errorStrategy = 'terminate'` (loud); only the `remote` profile retries twice
then falls back to `'ignore'`. Recurring lesson for this repo: a git tag is
not evidence of what a container contains — verify by running the image, not
by reading the tag it was built from. What was correct: the CLI flags were
genuinely unreachable from `main.py` before this release, and exposing them
was the right fix.

### v2.8.0 Changelog
- **Chunk-based parallelization:** Contigs split into 10M-base chunks for balanced multiprocessing load. Eliminates bottleneck where large contigs (e.g., chr1 ~249M) serialize work.
- **Float64 GC cumsum:** GC content cumulative sum uses float64 instead of float32, fixing precision loss on large chromosomes.
- **`--skip-chunking` flag:** CLI option to disable chunking and revert to whole-contig processing.
- **Single-process optimization:** When num_processes=1, work runs in-process without forking.
- **Bug fix:** `contig_lengths` computation with `--contigs` filter was pairing contig names with wrong lengths.

**Document Version:** 1.3
**Last Updated:** 2026-08-21
**Generated for:** Debugging and development assistance
