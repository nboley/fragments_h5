# fragments-h5 Agent Context Document

**Last Updated:** 2026-03-08
**Version:** 2.7.1
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
- S3 streaming support for BAM and FASTA files
- Optional metadata: GC content, strand, mapping quality, methylation, fragment clipping status

**Target Use Cases:**
- Cell-free DNA analysis and fragmentomics research
- Fragment-level genomic analysis requiring fast interval queries
- Large-scale genomic pipelines with S3 integration

**Related Paper:** https://www.nature.com/articles/s41416-021-01635-z

### 1.2 Current Status

- **Version:** 2.7.1
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

**Build Pipeline:**
1. **Fork:** Uses `fork` start method for minimal overhead (safe because output HDF5 opened after workers complete)
2. **Per-Worker Temp Dirs:** Each worker gets isolated temp directory for S3 index caching
3. **Per-Contig Workers:** Workers process individual contigs, write to temp HDF5 files
4. **Main Process Merge:** Main process opens output HDF5 after workers finish, copies contig data
5. **Signal Handling:** Workers ignore SIGINT/SIGTERM; main process handles cleanup

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
def worker_func(contig, bam_path, fasta_path, temp_output_path):
    # Worker gets isolated temp dir for S3 indexes
    with _temporary_working_directory():
        fragments = list(bam_to_fragments(bam_path, contig, ...))
        _write_contig_h5(temp_output_path, contig, fragments)
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
│   └── test_docker_build.py    # Manual Docker test
├── scripts/
│   ├── build_conda_package.sh
│   └── publish_conda_package.sh
├── conda-recipe/               # rattler-build conda package
│   ├── recipe.yaml
│   ├── variant_config.yaml
│   └── conda_build_config.yaml
├── docs/
│   └── plan_s3_multiprocess_worker_cwd.md
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
    #            fragment_length_counts, n_fragments
```

**Build Functions:**
```python
def build_fragments_h5(bam_fname, output_fname, fasta_filename=None,
                      allowed_contigs=None, set_mapq_255_to_none=False,
                      read_strand=True, read_methyl=False, single_end=False,
                      num_processes=1, include_duplicates=False,
                      store_fragment_end_clipped=True):
    """Main entry point for building fragment H5 from BAM"""
    
def build_sub_fragments_h5(args):
    """Worker function for multiprocessing (per-contig builds)"""
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
def bam_to_fragments(bam, contig=None, region_start=None, region_stop=None,
                     fasta=None, set_mapq_255_to_none=False, read_strand=True,
                     read_methyl=False, include_duplicates=False,
                     store_fragment_end_clipped=True):
    """Paired-end BAM → Fragment iterator"""
    
def single_end_bam_to_fragments(bam, ...):
    """Single-end BAM → Fragment iterator"""
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
  - `--include-duplicates`: Include duplicate-marked reads (default: exclude)
  - `--no-store-fragment-end-clipped`: Don't store clipping info
  - `--num-processes`: Number of workers (default: 1, use 'all' for all cores)
  - Logging: `--quiet`, `--verbose`, `--debug`, `--log-filename`, etc.

**Validation:**
- Checks if output file already exists (fails early with error)
- Checks BAM index exists (runs `samtools index` if local, fails if remote)
- Validates FASTA accessibility for S3 URLs (requires `.fai`, `.gzi` for compressed)

---

## 4. Dependencies and Build System

### 4.1 Dependencies

**Runtime:**
- `numpy` - Array operations
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

**Image:** `ghcr.io/nboley/fragments-h5:2.7.1` and `:latest`

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

**Supported URLs:**
- `s3://bucket/path/to/file.bam`
- `http(s)://url/to/file.bam`

**Requirements:**
- pysam/htslib built with S3 support (libcurl)
- AWS credentials configured (for non-public buckets)
- Index files must exist at same S3 path:
  - BAM: `.bai` file (e.g., `s3://bucket/sample.bam.bai`)
  - FASTA: `.fai` file (e.g., `s3://bucket/ref.fa.gz.fai`)
  - Compressed FASTA: `.gzi` file (e.g., `s3://bucket/ref.fa.gz.gzi`)

### 6.2 Multiprocessing with S3

**Challenge:** pysam downloads S3 indexes to local temp files, causing conflicts when multiple processes access the same S3 URL.

**Solution:** Per-worker temporary working directories (`docs/plan_s3_multiprocess_worker_cwd.md`)

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
- Longer fragments silently truncated or cause overflow

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
    single_end=False,
    num_processes=8,
    include_duplicates=False,
    store_fragment_end_clipped=True
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
    container "ghcr.io/nboley/fragments-h5:2.7.1"
    
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

**Sync Requirements:**
- `pyproject.toml` version
- `conda-recipe/recipe.yaml` version
- Git tag (`v2.7.1`)

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
- **Check:** `docs/plan_s3_multiprocess_worker_cwd.md`

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

- [ ] Update version in `pyproject.toml`
- [ ] Update version in `conda-recipe/recipe.yaml`
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
INDEX_BLOCK_SIZE = 5000      # Genomic positions per index block
CHUNK_SIZE = 1000000         # Array resize chunk size
MAX_FRAGMENT_LENGTH = 65535  # uint16 limit
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
- Single-end BAM processing (single_end_bam_to_fragments)

### 16.2 Potential Improvements

**Performance:**
- Investigate memory-mapped arrays for 2× speedup
- Optimize S3 index caching across workers
- Reduce HDF5 read overhead

**Features:**
- Support CSI indexes (in addition to BAI)
- Support alternative methylation tag formats
- Add fragment filtering during build (e.g., min/max length)

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

**Document Version:** 1.1
**Last Updated:** 2026-03-08
**Generated for:** Debugging and development assistance
