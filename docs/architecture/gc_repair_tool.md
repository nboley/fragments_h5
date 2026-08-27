# GC Repair Tool and Known Fragment H5 Defects

The console script `repair-fragments-h5-gc` (`src/fragments_h5/repair.py`) repairs two
independent defects in fragment H5 files. It exists and has a full test suite but has
**not been run against production data**. The full design document with correctness
arguments and review history is at `docs/pending/gc_repair_tool.md`.

## The two defects

### 1. GC float32 cumsum saturation

**Introduced:** `1c550f0` (2025-11-19). **Fixed:** `95c76f5` (2026-03-09).
**Affected tagged releases:** v2.2.1 through v2.6.0 (6 tags). Fixed from v2.7.2 onward.

`get_g_or_c_cumsum` (`fragment.py:404`) builds a per-base cumulative count of G/C over
a contig. Before the fix the cumsum was accumulated in the one-hot encoder's native
**float32** (`sequence.pyx:12`, `DTYPE = np.float32`):

```python
g_or_c_cumsum = seq[:, (1,2)].sum(axis=1).cumsum()          # pre-fix: float32
```

The fix added `.astype(numpy.float64)` before `.cumsum()`.

float32 represents integers exactly only up to 2\*\*24 = 16,777,216. Once the running
G/C count reaches 2\*\*24, adding 1.0 rounds back to 2\*\*24 and the accumulator sticks
permanently. There is also an intermediate band: N bases contribute +0.5 to the
accumulator, and float32 half-integers are exact only up to 2\*\*23 = 8,388,608. So
there are three regions per contig:

| Region (G/C cumsum) | Approx. offset (hg38) | Behaviour |
|---|---|---|
| < 2\*\*23 | first ~20 Mb | Fully exact |
| 2\*\*23 .. 2\*\*24 | ~20–40 Mb | C/G exact; N's +0.5 is an exact tie, discarded inside N runs |
| >= 2\*\*24 | past ~40 Mb | Fully saturated; `gc` stored as 0 |

The offsets are genome-average estimates at accumulator density 0.41470 (not the GC
fraction 0.4102 — the accumulator also gains +0.5 per N base). Per-contig positions vary
with local GC and N content.

**How to detect from a file.** Download the first ~256 KB of data for one long contig
(e.g. chr1). If `gc` bytes past ~40 Mb are all 0, the file is affected. A `gc` value of
0 is a valid encoding (`0 / 254.0 = 0.0` GC fraction) but biologically implausible at
scale. See `gc` encoding below.

### 2. Phantom trailing fragment (padding row)

**Introduced:** `caddb89` (2025-11-19). **Fixed:** `778f4d1` (2025-12-17).
**No tagged release carries this defect** — the earliest tag (v2.2.1, 2025-12-18)
postdates the fix by one day. Affected files came from unreleased `main`.

`mk_dataset` (`fragments_h5.py`, the `mk_dataset` closure inside `build_sub_fragments_h5`)
truncated each per-contig array with `data[: ff + 1]` where `ff` was already the fragment
count, so every dataset was written one element too long. The extra slot carries the zero
fill from the preallocated buffer.

This affects **every** per-contig dataset: `starts`, `lengths`, `mapq`, `gc`, `strand`,
methylation arrays, and `fragment_end_clipped`. Note `mapq` is 2-D (per-mate); the
phantom row is along axis 0.

Three consequences:
1. **`starts` is not sorted** — the phantom has `start = 0` appended after the last real
   fragment's start. This is a latent hazard (formally undefined `searchsorted` input),
   not a live incorrect result: the overlap mask at `fragments_h5.py:538-540` unconditionally
   filters the phantom, and the index entries yield correct empty results.
2. **The index counts the phantom** in its sentinel value. Benign today for the same reason.
3. **`fragment_length_counts` is inflated** — the phantom adds one count at bin 0 per
   contig group. `FragmentsH5.n_fragments` (which is `fragment_length_counts.sum()`) is
   therefore over-reported by ~one per contig group (~189 per file). **This is a live wrong
   value.**

## The 218 affected production files

218 fragment H5 files under `s3://fragmentomics.kariusdx.com/nboley/ibd_v2/build_frag_h5s/`
carry **both** defects. They were built 2025-11-24 from unreleased `main`, between the
introduction of both bugs and the fix of either. All 218 are fingerprinted to the REF-P12
reference (`GRCh38.p12.genome.fa.gz`).

**Backups exist:** all 218 are copied to
`s3://fragmentomics.kariusdx.com/nboley/gc_repair_backup_2026-08-25/` and verified by
CRC64NVME. The bucket is unversioned; the backup is the only rollback.

**Nothing currently reads `gc`.** `return_gc=False` is the default and no code in the
consuming project passes `True`.

## The repair tool

### What it does

Per target file, in order:

1. Download the S3 object (or accept a local file via `--local-file`)
2. Preflight: reference safety checks (§5 of the design doc)
3. Detect and truncate the padding row from every contig; rebuild the index and
   `fragment_length_counts`
4. Recompute all `gc` values from the reference FASTA
5. Validate: per-region diff rules, 255↔non-255 transition checks, pre-saturation oracle
6. If `--dry-run` (the default): stop and report. Nothing has been written.
7. Write the truncated datasets, rebuilt index, new `gc`, and `_repair_history` attr
8. Re-verify structural invariants
9. Upload and verify (CRC64NVME, single-part `put_object`)
10. Append one ledger record

### Key properties

- **`--dry-run` is the default.** There is no way to write by accident.
- **`--apply` requires** `--target-list`, `--expect-fasta-sha256`, `--ledger`, and
  `--max-files` (default 1).
- **Idempotent in datasets.** Running the tool on its own output changes no dataset byte
  and truncates no row. Only `_repair_history` grows, so the file sha256 changes.
- **No production identifiers in the package.** No bucket names, FASTA URIs, or file lists
  are hardcoded in `src/`.

### CLI

```
repair-fragments-h5-gc \
  --fasta /path/to/GRCh38.p12.genome.fa.gz \
  --local-file /path/to/file.h5 \
  --dry-run
```

For S3 targets:
```
repair-fragments-h5-gc \
  --fasta /path/to/ref.fa.gz \
  --target-list targets.txt \
  --bucket fragmentomics.kariusdx.com \
  --apply \
  --expect-fasta-sha256 <hex> \
  --ledger repair_ledger.jsonl \
  --max-files 1
```

### Key files

| File | Role |
|---|---|
| `src/fragments_h5/repair.py` | Tool implementation (CLI + repair pipeline) |
| `tests/test_gc_repair.py` | Test suite (54 tests incl. rounding, region classification, e2e) |
| `tests/test_gc_cumsum_overflow.py` | Regression test for the float32 saturation bug |
| `docs/pending/gc_repair_tool.md` | Full design document with correctness arguments |

## Reference safety

Two non-equivalent hg38 FASTAs were in use across the GC-bug window:
- **REF-P12:** `GRCh38.p12.genome.fa.gz` — used before 2025-12-16 and for the rebuild workflow
- **REF-ASSETS:** `pipeline-assets/hg38/sequence.fa.bgz` — used after 2025-12-16

They have identical contig names and lengths on all 25 primaries but differ by ~11.9 Mb
of hard-masking across 18 contigs. Supplying the wrong reference silently writes **127**
(the N value: `round(0.5, 5) * 254 = 127`) over real GC across the differing blocks.

**`_contig_lengths_str` cannot distinguish them** — on these files it is derived from the
BAM header, not the FASTA.

The tool enforces five reference-safety layers (§5 of the design doc), all failing closed:
1. `--fasta` is required (no default, no env fallback)
2. FASTA sha256 is pinned via `--expect-fasta-sha256` in `--apply` mode
3. Per-contig geometry checks and all-255 assertion for scaffolds absent from the FASTA
4. `get_g_or_c_cumsum` returning `None` is treated as an internal error, not "no GC"
5. Alphabet gate: refuses any FASTA with sequence bytes outside `{A,C,G,T,N}`

## GC encoding reference

`gc` is a `uint8` dataset. The encoding (`fragments_h5.py:884`):

| Stored byte | Meaning |
|---|---|
| 0–254 | GC fraction = `byte / 254.0` |
| 255 | Missing (no FASTA provided, contig absent from FASTA, or zero-length fragment) |

N bases encode to `[0.25, 0.25, 0.25, 0.25]` in the one-hot encoder (`sequence.pyx:40`),
contributing 0.5 to C+G. An all-N fragment therefore stores `round(0.5, 5) * 254 = 127`,
**not 0**. A stored value of 0 in the post-saturation region is the bug signature, not a
legitimate encoding of N-content.

## `has_gc` defect (live, unfixed)

`FragmentsH5.has_gc` (`fragments_h5.py:363`) uses `any()`:

```python
@property
def has_gc(self):
    return any("gc" in self.data[contig] for contig in self.data.keys())
```

A file with `gc` on some contigs but not others reports `has_gc == True`. A consumer
iterating all contigs and accessing `gc` on each will hit a `KeyError` on the contig that
lacks it. The repair tool handles this correctly (`repair.py:96-111`, `check_gc_presence`)
by checking every contig individually and aborting on partial presence. But the library
property itself is unchanged.

Note: `has_strand` (`fragments_h5.py:356`), `has_methyl` (`:346`), and
`has_fragment_end_clipped` (`:370`) all use `any()` with the same semantics, so the
partial-presence hazard applies to all of them — it is not specific to `gc`.

## Provenance

The 218 files carry **no build provenance** — `_build_version` / `_build_argv` landed in
v2.12.0 and these files predate it. The repair tool writes `_repair_history` (a JSON array
appended to on each repair, §9.1 of the design doc) but does **not** backfill
`_build_version` or `_build_argv`, because doing so would fabricate provenance for a build
whose parameters are unknown.

`_repair_history` will be dropped by downstream paths that rebuild H5s from a hardcoded
attribute allowlist (documented in
`docs/architecture/fragment_selection_and_build_provenance.md:720-739`).

## Test coverage

The test suite at `tests/test_gc_repair.py` covers:
- Padding row detection and truncation (including 2-D `mapq`)
- Index rebuild with key-set and sortedness invariants
- `fragment_length_counts` exact arithmetic identity
- GC recomputation with a naive reference encoder (no shared code with `repair.py`)
- Rounding agreement: exhaustive over 160,800 `(num, den)` classes with `den <= 400`,
  plus random fragments and engineered half-ulp cases
- Region classification (`< T23`, `[T23, T24)`, `>= T24`) in both directions
- Float32 accumulator simulation
- End-to-end saturating-contig integration test (~17.2 Mb, crosses both thresholds)
- Idempotency (second repair pass produces byte-identical datasets)
- `--expect-fasta-sha256` enforcement
- Partial-`gc` abort
- Files without `gc` dataset

**Honest residual:** no byte in the test suite comes from a genuinely damaged production
file. The saturating fixture reconstructs the pre-fix builder's behaviour rather than
observing it.
