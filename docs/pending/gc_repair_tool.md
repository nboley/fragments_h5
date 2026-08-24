# GC Repair Tool

> **Status: PROPOSAL, awaiting approval. Nothing here has been implemented and nothing has been
> run against production data.** Target release: **v2.13.0**. Execution against the 218 production
> files requires a *separate* explicit approval beyond approving this design — see §10.
>
> Version bookkeeping: `main` is tagged v2.12.1 but `pyproject.toml:7` in this worktree reads
> `2.12.0`. The implementing PR must reconcile that and set it to `2.13.0`; the tool reads its own
> version with `importlib.metadata.version("fragments-h5")` and never hardcodes it.

---

## 1. Goal

Make the `gc` dataset in 218 corrupted fragment H5 files equal to what a correct build would have
written, and be able to demonstrate that it does.

Secondary goal, stated by the user and binding: the tool must be reusable for future
repair/backfill work on fragment H5s, *without* being over-engineered. §11 resolves that tension
concretely.

## 2. Problem statement

### 2.1 Mechanism

`get_g_or_c_cumsum` (`src/fragments_h5/fragment.py:404`) builds a per-base cumulative count of G/C
over a contig. Before commit `95c76f5` (v2.7.2, 2026-03-09) the cumsum was accumulated in the
one-hot encoder's native **float32** (`sequence.pyx:12`, `DTYPE = np.float32`). The pre-fix line is
`fragment.py:423` at `95c76f5^`:

```python
g_or_c_cumsum = seq[:, (1,2)].sum(axis=1).cumsum()          # pre-fix: cumsum in float32
```

and the fix inserted one call (`fragment.py:442` at HEAD):

```python
g_or_c_cumsum = seq[:, (1,2)].sum(axis=1).astype(numpy.float64).cumsum()
```

float32 represents integers exactly only up to 2**24 = 16,777,216. Once the running G/C count
reaches 2**24, adding 1.0 rounds back to 2**24 (round-half-to-even; the next representable float32
is 2**24 + 2). The accumulator **sticks permanently**.

Consequently, for every contig:

- fragments entirely **before** the saturation point: correct;
- fragments **spanning** the saturation point: `cumsum[stop] - cumsum[start]` is nonzero but
  **undercounted** — a real number, silently wrong;
- fragments entirely **after** it: `cumsum[stop] - cumsum[start] == 0.0` → `gc` stored as `0`.

The accumulator gains `+1.0` per C/G **and `+0.5` per N**, so the rate at which it climbs is not the
GC fraction but the *accumulator density* `(C + G + 0.5·N)/total` = **0.41470** (§4.2). At that
density 2**24 is reached at **~40.45 Mb**, which matches the observed 34–43 Mb per-contig
breakpoints. This and every other offset in this document is a **genome-average estimate**; the
actual per-contig position varies with local GC and N content, which is why the tool classifies
fragments by cumsum value rather than by coordinate (§5b).

**There is an intermediate band, and it is populated.** float32's ulp is 1.0 throughout
`[2**23, 2**24)`, so the representable set there is *exactly the integers*. `C`/`G` contribute
`+1.0`, which is still exact; `N` contributes `+0.5`, which lands on an **exact tie** and is
resolved round-half-to-even. REF-P12 contains **161,331,703 `N` bases** (§4.2), so above
2**23 = 8,388,608 — **~20.23 Mb** at density 0.41470 — N contributions become unreliable while C/G
stay exact. The behaviour in that band is not a slow random walk; it is sharply directional:

- inside a **contiguous `N` run** the accumulator parks on an even integer and every subsequent
  `+0.5` is **discarded outright** (even + 0.5 ties down to even). A megabase N block contributes
  essentially **nothing** instead of 500,000 — a systematic, near-total loss;
- an **isolated `N`** alternates by parity: at even `x` the `+0.5` is lost (−0.5 of error), at odd
  `x` it ties **up** to `x + 1` (+0.5 of error).

Three regions, not two:

| region (G/C cumsum) | approx. offset | float32 behaviour |
|---|---|---|
| `< 2**23` | first ~20.23 Mb | fully exact; original build was correct |
| `2**23 .. 2**24` | ~20.23–40.45 Mb | C/G exact; `N`'s `+0.5` is an exact tie → discarded inside N runs, ±0.5 by parity for isolated N |
| `>= 2**24` | past ~40.45 Mb | fully saturated; `+1.0` is a no-op |

This distinction is load-bearing twice: §7.2 uses 2**23 as its oracle bound, and §5b/§7.3 must
*permit* middle-band corrections on N-containing spans rather than treating them as hard errors.
Those corrections are **few but not small** — see §5b.

**This is why "repair only the values that look broken" is wrong.** The spanning fragments are the
minority but they are the ones a `gc == 0` filter misses. The tool recomputes *everything*
(§3.1).

### 2.2 What is and is not affected

Only `gc`. `starts`, `lengths`, `mapq`, `strand`, methylation, `fragment_end_clipped`,
`index/<contig>`, and all file-level attrs are derived from the BAM and never touch
`get_g_or_c_cumsum`. The GC value is computed in `fragment.py:492-496` / `:603-606` / `:736-743`
and consumed only by the quantizer at `fragments_h5.py:841-846`. Nothing else reads it during a
build. Confidence: **high** — `95c76f5` changed only the accumulator dtype and added region
support.

### 2.3 Blast radius

218 `.h5` objects under `s3://fragmentomics.kariusdx.com/nboley/ibd_v2/build_frag_h5s/`
(prefix total measured: 437 objects / 91.13 GB — 218 `.h5`, 218 zero-byte folder markers, 1 `.gz`),
all from a single build run on 2025-11-24, all fingerprinted to REF-P12 (218/218). They are
consumed by the CD (281) and UC (315) IBD cohorts; 123 CD + 95 UC of the 218 are consumed, and all
218 fall on the *train* side of the split.

**No clean substitute exists.** Every one of the 218 has exactly one alternate copy elsewhere and
every one of those predates the GC feature entirely (no `gc` dataset). Recomputation is the only
option.

### 2.4 Disposition (fixed by the user)

Back up all 218 to a backup prefix, **verify the backups**, then **overwrite the originals in
place**. Chosen because `fragmentomics.kariusdx.com` has **no object versioning** — an overwrite is
terminal. Established empirically: `head_object` on one of the 218 returns no `VersionId`, while the
same call against `karius-biomarker-data-assets` (versioning Enabled) does.
`GetBucketVersioning` and `ListObjectVersions` are both `AccessDenied` for our credentials, so the
`head_object` behaviour is the only evidence — but it is the behaviour that actually matters.

---

## 3. Design

### 3.1 What the tool does, end to end

New console script `repair-fragments-h5-gc` (a second `[project.scripts]` entry alongside
`build-fragments-h5 = "fragments_h5.main:main"`; the repo has no `scripts/` dir and no subcommand
structure, so a flat second entry point is the least-surprising shape).

Per target file, strictly in this order. **No step before 6 touches S3 in any way that writes**, so
`--dry-run` (the default) provably performs zero writes:

1. **Download** the S3 object to local scratch (`/tmp`, 549 GB total / 414 GB free; largest
   observed file ~740 MB). Compute `sha256_local` of the downloaded bytes and record the source
   `ETag` returned by the GET.
2. **Preflight / reference safety** (§5) on the local copy. Abort the file on any mismatch. **While
   the file is still pristine, record a per-dataset raw-chunk `sha256` for every dataset that is
   not a `gc` dataset, and a `sha256` of the serialized attr set excluding `_repair_history`.**
   These are the *only* operand step 8 has for its byte-identity check (the backup is in S3 and the
   pre-mutation local bytes will not survive step 7), so they must be taken here.
3. **Recompute** all `gc` values for every contig present in both the H5 and the FASTA, into
   in-memory `uint8` arrays (§3.5 specifies exactly how). Do **not** write yet.
4. **Pre-saturation oracle gate** (§7.2) and **diff and gate**: compare recomputed vs stored, per
   contig, per region (§2.1). Report counts. The `< T23` byte-identity requirement of §5b **is**
   the §7.2 oracle — one check, stated once in §5b's region table, not two.
5. If `--dry-run` (the default), **stop here** and emit a report — which includes the read-only
   `head_object` probe of the backup key (§6.5, §10.2), so that the "backup exists with different
   content" halt condition is discovered in dry run rather than first seen under `--apply`. Nothing
   local has been modified and nothing in S3 has been written.
6. **Back up** via server-side `copy_object` to the backup prefix and **verify the copy** (§6).
   Backup keys are **write-once** (§6.5). *Nothing is written to the original key until this step
   succeeds for that file.*
7. **Write** the new `gc` arrays into the *local* copy, opened `h5py.File(path, "r+")`, **and append
   the `_repair_history` element (§9.1) in the same open handle**, before the file is closed. The
   ordering is fixed here because `sha256_repaired` is taken at step 8 and must cover the history
   attr.
8. **Re-verify** the local file: reopen read-only, confirm every `gc` dataset equals the recomputed
   array bit-for-bit, confirm every non-`gc` dataset — and every attr **except `_repair_history`** —
   still hashes to the value recorded in step 2, confirm `_repair_history` parses as JSON and gained
   exactly one element, confirm the file opens through `FragmentsH5`. Then compute
   `sha256_repaired` over the closed file.
9. **Append a `status: "started"` ledger record** (§8.3), naming the key, the verified
   `backup_uri`, `sha256_original`, and `sha256_repaired`. This record exists *before* the upload
   so that a crash mid-upload is distinguishable from a crash before backup.
10. **Upload** the local file over the original key with a single-part `put_object` (§6.4).
11. **Verify the upload** (§6.4) and append the `status: "ok"` ledger record.

**Recompute all, not just the broken ones.** Three reasons: (a) spanning fragments are corrupt
without being zero (§2.1); (b) a full recompute makes the clean-file case a *provable `gc` no-op*
(§5b), which is the strongest safety property available here; (c) it removes an entire class of
"which values did we decide to touch" bugs. Cost is bounded and acceptable (§3.4).

Confidence: **high** on the shape; **medium** on step 2's per-dataset raw-chunk hashing being cheap
to express in h5py without decompressing everything — see §12.8.

### 3.2 In-place `r+` vs full-copy rewrite

**Decision: in-place `r+` on the local copy.**

The decisive argument is *fidelity*, not size. Requirement: "must not alter anything currently
correct." An in-place dataset assignment provably cannot touch another dataset's bytes, chunk
layout, filter pipeline, or attrs. A full-copy rewrite re-creates every dataset and every attr,
and any divergence in chunk shape, gzip level, attr dtype (`|S` vs variable-length str,
`_contig_lengths_str` is a `str(dict)` read back with `eval()` at `fragments_h5.py:305`) is a
silent regression we would have to test for exhaustively.

The size objection is empirically weak. Measured on a synthetic chunked+gzip `uint8` dataset
(~40% real beta-distributed GC / ~60% zeroed, rewritten with fully-real data):

| variant | bytes |
|---|---|
| corrupted original | 1,985,045 |
| repaired in place (`r+`) | 4,947,843 |
| built correctly from scratch | 4,935,336 |

Delta vs a from-scratch build: **~0.25%**. HDF5's free-space manager reuses the freed chunk
allocation when the replacement chunk lands in the same size class; five repeated same-shape
rewrites of random data left file size exactly flat at 2,011,192 bytes.

**But note the two different baselines, because only one of them is small.** Relative to a
*from-scratch correct build* the growth is ~0.25%. Relative to the *corrupted original* it is
1,985,045 → 4,947,843 = **~2.5×** on the `gc` dataset, and that is inherent, not a defect: a
zero-filled `gc` dataset gzips to almost nothing while real GC bytes do not. Repaired files will
therefore be **materially larger than the objects they replace**, by roughly the compressed size of
one real `gc` dataset per file. That growth must be carried in three budgets: the upload leg (§3.4),
`/tmp` headroom, and the S3 footprint (§2.3, §6.1). §7.3's "within 5% of expectation" is measured
against the from-scratch build — the right baseline for detecting HDF5 bloat — but it deliberately
says nothing about this absolute growth, which is expected. Report **both** numbers.

**Evidence level: weak-to-moderate.** This was one write pass on a synthetic single-dataset file,
because *no real fragment H5 exists on this machine* (checked `/tmp`, `$HOME`, `tests/`). Real
files have ~25 contig groups × 6–10 datasets each plus an index, and a more fragmented free-space
table. The physics should carry, but **the first real file repaired must have its size compared
against a from-scratch rebuild expectation**, and if in-place growth exceeds **5%**, fall back to
running `h5repack` on the output (note: `h5repack` availability is unverified here) or revisit.
That check is a gate in the graduated rollout (§7.4), not an afterthought.

Precedent for in-place mutation exists (`build_fragments_h5` reopens its own output `"r+"` at
`fragments_h5.py:1248` to add `fragment_length_counts`; `tests/test_build_provenance.py:109` opens
`h5py.File(..., "r+")`), but **overwriting an existing dataset's contents is a new pattern in this
repo**. Treat it as such: it gets its own test.

### 3.3 Where the work happens

Local disk, one file at a time per worker: download → repair → verify → upload. Not streaming, not
in-place-over-S3.

This has a load-bearing safety consequence: **an interrupted write can never produce a corrupted
object in S3.** The `r+` write happens on a throwaway local copy; S3 sees only a complete,
verified `put_object`. The "neither original nor repaired" failure mode is confined
to `/tmp`, where the recovery is `rm` and retry. This alone justifies the download-repair-upload
shape over anything cleverer.

**The upload is a single-part `put_object`, never `upload_file`.** This is not a style preference;
see §6.4. `boto3`'s `upload_file` switches to multipart above an 8 MB threshold, and every one of
the 218 is far above it, which would make the returned `ChecksumSHA256` a composite value that can
never be compared against a full-object hash. Single-part `put_object` is available because every
file is far below the 5 GB single-part limit. A consequence worth stating: because nothing here
uses multipart, **there is no orphaned-multipart-part liability** and no lifecycle rule is needed.

`awscli` is not installed; `boto3` 1.43.74 is. Credentials are the personal IAM user
`arn:aws:iam::573640641260:user/nathan.boley`. Server-side `copy_object` is confirmed working with
these credentials (tested against a zero-byte folder marker, verified with `head_object`, deleted,
deletion confirmed; no `.h5` was touched). `copy_object` has a 5 GB single-call limit — all files
here are far below it, but the tool must **explicitly refuse** any object > 5 GB rather than
silently truncating or erroring late.

### 3.4 Parallelism and memory

`get_g_or_c_cumsum` materializes the whole contig as a Python `str`, then an `(L, 4)` float32
one-hot array, before reducing. For chr1 (248,956,422 bp) that is ~4 GB transient for the one-hot
plus ~2 GB for the float64 cumsum: **~6 GB peak per concurrent contig**. There is **no caching** —
every call re-fetches and re-encodes. Naively, 218 files × 25 contigs = 5,450 full-genome encodes.

**Proposal: a one-shot per-contig cumsum cache.** Preflight computes each contig's float64 cumsum
once and writes it to `--cumsum-cache DIR/<sha256-of-fasta>/<contig>.npy`. Workers open them with
`numpy.load(..., mmap_mode="r")`. Cost: ~3.25e9 bases × 8 B ≈ **26 GB on disk** (fine on the
549 GB `/tmp` overlay, 414 GB free), shared across all workers through the page cache, per-worker
resident memory ≈ 0. Access is near-sequential because fragments are stored sorted by start within
a contig. The cache directory is keyed by the FASTA's sha256, so a different FASTA can never hit a
stale cache.

This is ~30 lines and it deletes the entire memory-pressure question, so it clears the
"not over-engineered" bar. If it is rejected, the fallback is: no cache, `--num-processes 4`
(4 × 6 GB = 24 GB peak against the **96 GB available** on this 16-CPU host, of 124 GB total), which
also works but re-encodes the genome 218 times.

Default `--num-processes 4`; one *file* per worker (not one contig per worker) so that the
download/backup/repair/upload/ledger sequence for a file is owned end-to-end by a single process
and resumability stays simple.

**Runtime — an estimate, not a measurement.** Basis: a full gzip-decompress + per-byte histogram
pass over this FASTA took **~5 min** in pure Python and **~2 min** numpy-vectorized (measured, §4.2).
Cache build is therefore ~5 min one-shot (decompress dominates; the cumsum is one vectorized pass
per contig). Per file: ~740 MB down + **more than 740 MB up** — the repaired file is larger than the
original by roughly the compressed size of a fully-populated `gc` dataset (§3.2), so budget the
upload leg at ~1.2–1.5× the download leg until stage 2 measures it — plus a memory-mapped scan of
the fragment arrays. Call it **2–5 min/file** dominated by network. At `--num-processes 4` and
218 files that is roughly **2–5 hours wall-clock** for the full apply, plus ~1–2 hours for the
§7.4 stage-3 dry-run gate across all 218 (download-bound, no upload). **Every number in this
paragraph is an estimate**; the stage-2 rehearsal produces the first real datum and the runbook
should record it.

### 3.5 Recompute: the vectorized path, and the `round(x, 5)` hazard

Step 3 must be vectorized. A Python per-fragment loop over O(10^8) fragments × 218 files is not
viable. But the naive vectorization is **not** guaranteed to reproduce Stage A, and the difference
is observable in the stored byte:

> `numpy.round(a, 5)` and Python's `round(x, 5)` are different functions. NumPy multiplies by 1e5,
> `rint`s, and divides back — a sequence of three inexact float64 operations. CPython's `round(x, 5)`
> uses correctly-rounded decimal conversion. They disagree on some inputs by 1 ulp, and Stage B
> then multiplies by 254, so a 1-ulp disagreement can straddle a `k + 0.5` boundary and change the
> stored `uint8` by one.

Stage B is safe to vectorize: `int(round(y))` (CPython, on a float) and `numpy.rint(y)` both round
half-to-even on the binary value and agree elementwise. Stage A's `round(x, 5)` does not have that
property, so it is the one that must be pinned.

**Decision: reproduce CPython `round(x, 5)` semantics exactly, and prove it by test.** Concretely:

- Compute `q = (cumsum[stops] - cumsum[starts]) / lengths` as a float64 array (numerator exact per
  §4.2), with `lengths == 0` masked out (below).
- Apply a decimal-correct round-to-5-places, `q5`. The implementation is one of:
  (a) an elementwise loop calling CPython's `round` — the semantics are correct by construction;
  benchmark it, and if it fits the §3.4 runtime budget, **use it and stop**;
  (b) a vectorized approximation, permitted **only** if (a) misses the budget, and **only** with
  the agreement test below green.
  Pick (b) on measured need, never on aesthetics.
- `gc_u8 = numpy.rint(q5 * 254).astype(numpy.uint8)`.

**Required test, blocking:** elementwise equality between the chosen path and a plain Python-loop
reference implementing `int(round(round(float(num)/float(den), 5) * 254))`, over **>= 10^7 real
fragments** drawn from an actual fragment H5, *plus* a synthetic set engineered so that `q5 * 254`
lands within 1 ulp of `k + 0.5` for many `k`. "Spot-checked on a few thousand" is not acceptable —
this is the last bit of every stored byte.

**Zero-length fragments (`length == 0`).** Stage A emits `None` when `g_or_c_cumsum is None` **or**
`frag_stop == frag_start` (`fragment.py:492-493`, `:603-604`, `:736-737`), and Stage B maps
`None`/NaN → 255 (`fragments_h5.py:844-846`). `lengths` is `uint16`, so zero is representable and
must be **checked, not assumed absent**. The rule:

> `length == 0` → force the recomputed value to **255**, exclude the fragment from the division
> (no divide-by-zero, no NaN entering the quantizer), and exclude it from the diff histogram
> entirely.

Without this rule an unguarded vectorized divide produces a NaN that quantizes unpredictably, and
§5b's "`255 → x` is a hard error" would abort the file. Note this is the *only* legitimate source
of 255 on a contig that **is** present in the FASTA; every other 255 means the contig is absent,
which §5 layer 3 handles separately.

---

## 4. Correctness argument

**Claim to evaluate:** recomputation is bit-exact, not approximate, and independent of chunking.

### 4.1 The argument, restated precisely

Stage A (`fragment.py:492-496`, unchanged since before the bug):

```python
gc = round(float(cumsum[frag_stop - gc_offset] - cumsum[frag_start - gc_offset])
           / float(frag_stop - frag_start), 5)
```

Stage B (`fragments_h5.py:841-846`, also unchanged) — quoted with the range guard immediately above
it at `:842-843`, which is part of the contract the recompute must satisfy:

```python
if not ((fragment.gc is None or numpy.isnan(fragment.gc)) or 0 <= fragment.gc <= 1):
    raise ValueError(f"{fragment.gc} is invalid for {fragment}")
gc_arr[ff] = 255 if (gc is None or numpy.isnan(gc)) else int(round(gc * 254))
```

The guard is why a recomputed `q5` outside `[0, 1]` is an internal error rather than something to
clamp: the original build would have raised, so any such value means the recompute diverged.

Both stages are deterministic IEEE-754 / CPython operations with no threading, no `-ffast-math`
path, and no dependence on array size, so **as written they are not where the risk lives** — the
risk in Stage A is in *reproducing* it, since a vectorized `round(x, 5)` is not the same function
(§3.5). Taking faithful reproduction as given, the remaining question for this section is whether
`cumsum[stop] - cumsum[start]` is exact and chunk-invariant.

The exactness of the sum depends entirely on which per-base C+G values can actually occur, and that
depends on **which bytes appear in the FASTA** — a property of the reference, not of the library.
The encoder table (`sequence.pyx:24-41`) admits values that are *not* multiples of 0.5:

```
'B': [0, 1/3, 1/3, 1/3]   -> C+G = 2/3
'V': [1/3, 1/3, 1/3, 0]   -> C+G = 2/3
'H': [1/3, 1/3, 0, 1/3]   -> C+G = 1/3
'D': [1/3, 0, 1/3, 1/3]   -> C+G = 1/3
```

The rest of the table: A/C/G/T give 0 or 1; K/M/R/Y/S/W give 0, 0.5 or 1; N and X give
`[.25,.25,.25,.25]` → C+G = 0.5 exactly (confirming the "all-N fragment stores 127, not 0" fact);
any unlisted byte maps to `[0,0,0,0]` (`sequence.pyx:45-51`), and the `'a'` pad prepended at
`fragment.py:441` maps to `[1,0,0,0]`, contributing 0. Lowercase is aliased to uppercase at
`sequence.pyx:42-43`, so soft-masking is irrelevant to the value.

So the argument cannot be made in the abstract. It has to be made against a **measured** alphabet.

### 4.2 The argument, made against the measured alphabet of REF-P12

**Measurement.** A full byte histogram over
`s3://fragmentomics.kariusdx.com/nboley/resources/GRCh38.p12.genome.fa.gz` (938,986,662 bytes
gzipped), counting sequence-line bytes separately from header-line bytes:

```
=== SEQUENCE bytes ===
  'A'       910,083,024
  'C'       632,693,920
  'G'       635,330,310
  'N'       161,331,703
  'T'       912,769,936

B / V / H / D (the 1/3 codes):            0
K / M / R / Y / S / W (half-integer):     0
lowercase:                                0
```

Total 3,252,208,893 bases. **GC fraction excluding N = 0.4102.** N is 4.96% of the file.

**The accumulator density is not the GC fraction.** `get_g_or_c_cumsum` adds `+1.0` per C/G *and*
`+0.5` per N, so the quantity that determines where the float32 thresholds are crossed is

```
(C + G + 0.5·N) / total
  = (632,693,920 + 635,330,310 + 80,665,851.5) / 3,252,208,893
  = 1,348,690,081.5 / 3,252,208,893
  = 0.41470
```

Every threshold offset in this document (§2.1, §7.2, §7.2.1) is derived from **0.41470**, not from
0.4102, and every one of them is a **genome-average estimate** — per-contig positions vary with
local GC and N content. Nothing downstream depends on the exact offset, because the tool classifies
by cumsum value, never by coordinate.

The only
B/H/V/D bytes anywhere in the file are inside FASTA *header description text*
(`{'B':27,'D':4,'H':857,'V':33}`), which the encoder never sees. REF-P12 is entirely uppercase:
all of its masking is hard-masking via `N`, and there is no soft-masked sequence at all.

**Consequence 1 — the reachable value set is exactly `{0.0, 0.5, 1.0}`.** `A`/`T` → 0, `N` → 0.5
(from `0.25 + 0.25`), `C`/`G` → 1.0. Every reachable per-base C+G is an exact multiple of 0.5.

**Consequence 2 — the per-base float32 add is exact.** `seq[:, (1,2)].sum(axis=1)` at
`fragment.py:442` runs the pairwise C+G addition in **float32**, *before* the
`.astype(numpy.float64)`; only the `cumsum` is float64. That pairwise add is `0+0`, `0+1.0`, or
`0.25+0.25` — all three exactly representable in float32, all three exactly computed. Obligation
discharged, not glossed.

**Consequence 3 — the float64 cumsum is exact with enormous margin.** A float64 sum of multiples of
0.5 is exact while twice the running total stays below 2**53, i.e. below **2**52 ≈ 4.5e15**. The
largest per-contig G/C cumsum here is ~0.41 × 249 Mb ≈ **1.0e8**. The margin is ~4×10^7. There is
no margin question to argue about.

Given exactness, every partial sum is represented exactly, so `cumsum[stop] - cumsum[start]` is
identical whether the cumsum was built over the whole contig or over any sub-region containing
`[start, stop)`. **Chunk-invariance and process-count-invariance follow.** Stage A and Stage B then
carry it through deterministically (subject to §3.5's `round(x, 5)` care).

**Consequence 4 — this also fixes the float32 bounds used elsewhere in the doc.** For multiples of
0.5, float32 is exact only to 2**23 = 8,388,608, which at accumulator density 0.41470 is
**~20.23 Mb** — this is the §7.2 oracle bound, now confirmed rather than assumed. Full saturation
(where `+1.0` stops changing the value) is at 2**24 = 16,777,216 → **~40.45 Mb**, consistent with
the observed 34–43 Mb breakpoints. And because there are 161.3M `N` bases, the intermediate band of
§2.1 is genuinely populated.

**Consequence 5 — the cumsum error is piecewise constant, and this is what licenses §5b's
middle-band rule.** In `[2**23, 2**24)` the `+1.0` steps at C/G positions are exact, so the
difference `err(i) = float32_cumsum(i) − true_cumsum(i)` **cannot change at a C/G or A/T position**;
it changes only where an `N` is consumed. Therefore for any span `[start, stop)` containing no `N`,
`err(stop) − err(start) = 0`, and the recomputed `gc` for that fragment must equal the stored `gc`
exactly. That is a proof, not a heuristic, and it is precisely §5b's predicate: in the middle band,
**a change is permitted only if the span contains an `N`.**

**This is conditional on the reference, and the tool enforces the condition.** `sequence.pyx` still
maps B/V/H/D to thirds; a different FASTA could reach them. The argument above is valid only on the
alphabet `{A,C,G,T,N}`, so the tool **refuses to run** on any FASTA whose sequence lines contain
anything else (§5 layer 5). That gate is load-bearing, not ceremonial: see §7.2 for why B/V/H/D
would have collapsed the oracle to nothing.

**Grade: the conclusion holds and the proof now rests on a measurement rather than on an
assumption about the reference.**

### 4.3 Residual assumptions, and how to discharge each

| # | Assumption | Status | How to discharge |
|---|---|---|---|
| R1(a) | The 2025-11-24 build used **REF-P12, not REF-ASSETS** | **Decidable, and the tool decides it.** Both candidates are in hand, so this is not an assumption — it is a computation. | **§7.2.1's two-reference byte-diff**: recompute the written region under both references and count differing bytes. Zero ⇒ reference choice provably does not affect any byte written to that file, and R1(a) is **closed** for it. Nonzero ⇒ exact blast radius reported before any write, human decides at §7.4 stage 3. |
| R1(b) | The REF-P12 object today is **byte-identical** to the REF-P12 object on 2025-11-24 | **Not established, and not establishable.** Fingerprinting is a 5-contig missingness test, not byte equality. The bucket has no versioning and `ListObjectVersions` is `AccessDenied`, so there is no earlier object to diff against. | **The pre-saturation prefix check (§7.2)** is the only evidence: strong corroboration on the real target files, over the region that needs no repair. It is **corroboration, not proof**, and this row stays open permanently (§12.1). The backstop is that the failure is recoverable from backup (§8.5). |
| R2 | `one_hot_encode_sequences` emits only float32 values, and the reachable per-base C+G set is `{0, 0.5, 1.0}` | **DISCHARGED for REF-P12 by measurement** (§4.2): sequence bytes are exactly `{A,C,G,T,N}`, zero B/V/H/D, zero K/M/R/Y/S/W, zero lowercase. Table at `sequence.pyx:24-41` verified by reading. | Made a **mandatory preflight gate** (§5 layer 5): the tool rescans the supplied FASTA's sequence-line character inventory and refuses to run if anything outside `{A,C,G,T,N}` (case-insensitive) appears. R2 becomes a per-run fact, not an assumption. |
| R3 | `frag_stop == start + length` | **DISCHARGED from source.** `Fragment.length` → `tlen` → `self.stop - self.start` (`fragment.py:244-250`); the writer stores `lengths_arr[ff] = fragment.length` (`fragments_h5.py:817`); the reader reconstructs `stops = starts + lengths` (`fragments_h5.py:537`). | Nothing further needed; §7.1/§7.2 would also catch an off-by-one loudly. |
| R4 | Stage A/B code paths are unchanged between the 2025-11-24 build and HEAD | **High confidence.** `95c76f5` touched only dtype + region support; the three Stage A sites and the Stage B quantizer are textually unchanged. | Confirmed again by §7.1 passing on a 2026-05-22-built file. |
| R5 | `numpy.cumsum` on float64 is a plain left-to-right sequential scan (no pairwise/SIMD reassociation) | Reassociation would not matter **because the sum is exact** — exact addition is associative. R5 is therefore not load-bearing. | Nothing needed. |
| R6 | The 218 files were built without `--contig-name-map` | Established. | n/a |
| R7 | The three Stage A sites compute the same thing | **DISCHARGED from source.** `fragment.py:496`, `:606`, `:739-743` are arithmetically identical: same expression, same `g_or_c_cumsum is None or frag_stop == frag_start → None` guard, same `gc_offset` subtraction. The only difference is line wrapping. | Nothing further needed. It does not matter which of the three built a given file. |

The honest summary: **the arithmetic is now settled by measurement** — a measured alphabet, an exact
float32 pairwise add, and a float64 margin of ~4×10^7. **Reference identity is half settled**: the
two-reference question R1(a) is decided per file by §7.2.1's byte-diff, while R1(b), arbitrary drift
of the REF-P12 object, has no possible oracle and remains the real residual risk — §7.2 constrains
it only over the prefix that was never corrupted.

---

## 5. Reference safety

Using REF-ASSETS where REF-P12 was used (or vice-versa) would **silently write 127 (the N value)
over real GC** across ~11.9 Mb of differently-masked sequence. The two references have identical
contig names and lengths on all 25 primaries and agree on sequence for only 7 of 25, so nothing
about the H5 header or contig table catches this. `_contig_lengths_str` is useless for the purpose
**on these files**: they were built from BAM input, where the attr is derived from the BAM header
rather than the FASTA. (The scoping matters — for TSV/BED input the contig lengths *are* taken from
the FASTA, at `fragments_h5.py:1024-1037`. That path just is not the one that produced the 218.)

Five layers, all mandatory, all failing **closed**:

1. **`--fasta` is required and has no default.** No environment fallback, no "guess from the
   manifest", no hardcoded S3 URI anywhere in the package.
2. **FASTA identity is recorded and pinned.** Preflight computes `sha256` of the FASTA (and of its
   `.fai`), prints it, and writes it into the ledger and into the repaired file's provenance (§9).
   `--expect-fasta-sha256 <hex>` is **required in `--apply` mode** and aborts the whole run on
   mismatch, so the runbook can pin the exact bytes.
3. **Per-file structural fingerprint, before any write.** For every contig group under `data/` in
   the H5:
   - contig **present** in the FASTA → **check its geometry, then** recompute its `gc` in full:
     - `fasta.get_reference_length(c) == h5.contig_lengths[c]` — equal lengths, or **abort the
       file**;
     - `max(starts + lengths) <= fasta.get_reference_length(c)` — every stored fragment lies inside
       the contig, or **abort the file**.

       This is not defensive padding. Stage A indexes `g_or_c_cumsum[frag_stop - gc_offset]` into
       an array of length `len(seq) + 1` (`fragment.py:441-442`); a FASTA whose contig is shorter
       than the H5's stored fragments gives an `IndexError` at best and, with negative-index
       wraparound, silently wrong GC at worst. Both checks are O(1) and O(n) respectively over
       arrays the tool already reads.
   - contig **absent** from the FASTA → **assert its stored `gc` is entirely 255**; if so, leave it
     completely untouched; if not, **abort the file**.

   This rule does the rest of the work. REF-P12 names primaries `chr*` but scaffolds by accession
   while the H5s use UCSC scaffold names, so `chrUn_*` / `chr*_random` legitimately store 255
   throughout — correct, expected behaviour for these files, which the rule reproduces exactly by
   not touching them. And if someone supplies REF-ASSETS (UCSC-named, 195 contigs), those contigs
   *are* present, so they would be recomputed from 255 to real values — a `255 → x` transition,
   which §5b makes a **hard error that aborts the file**. There is no override flag; see §10.3 and
   §11 for why this tool has no bypasses at all.
4. **`get_g_or_c_cumsum` returns `(None, 0)` silently when `chrom not in fasta_file`
   (`fragment.py:404`+). The tool must never call it in a way where that silence is possible.**
   Membership is resolved in preflight (layer 3) and the repair path treats a `None` return as an
   internal error, not as "no GC available". Writing 255 over a contig we merely failed to *find*
   is the single easiest way to destroy data here.
5. **Alphabet gate — the tool refuses to run on a FASTA it cannot reason about.** Preflight scans
   the supplied FASTA's **sequence lines** (header lines excluded — REF-P12's headers contain
   `B`/`D`/`H`/`V` in description text that the encoder never sees) and builds a character
   histogram. If **any** byte outside `{A,C,G,T,N}`, case-insensitively, appears, the tool
   **aborts the entire run** with the offending characters and their counts.

   This is not a formality. §4.2's exactness argument — and therefore §7.2's oracle — is valid
   *only* on that alphabet. Had B/V/H/D been present, the usable pre-saturation prefix would have
   been **zero**, not merely shorter: `float32(1/3)` is not on any half-integer grid, so three
   consecutive `H` bases already push the accumulator off the exactly-representable set. K/M/R/Y/S/W
   would be tolerable in principle (they are multiples of 0.5) but they are measured absent, and
   admitting them would mean re-deriving the bounds; the gate keeps one bound for one alphabet.

   The result is cached as `<fasta-sha256>.alphabet.json` alongside the cumsum cache (§3.4), so the
   scan is one pass per FASTA per machine, not one per file. The measured REF-P12 histogram is
   committed under `docs/pending/gc_repair/` as the expected value.

Additionally: the 5-contig scaffold-missingness classifier from `/tmp/fasta_fingerprint.py`
(`chrUn_KI270442v1`, `chr16_KI270728v1_random`, …) is subsumed by layer 3, which is strictly more
general — it checks *every* contig instead of five. The prototype's value was as a cheap S3-side
probe; the repair tool has the file locally and can do better.

Confidence: **high** that these layers prevent the wrong-FASTA accident. They do not by themselves
prove the *right* FASTA (that is R1 / §7.2 / §7.2.1).

### 5b. The no-op guard as a safety property

A correct repair of a correct file leaves every **dataset** byte-identical. This is not a
nice-to-have; it is the tool's central invariant, and it generalises:

> **Invariant.** For any file, `repair(file)` must produce output byte-identical to the input on
> every dataset the tool did not intend to change, and on every attr except `_repair_history`, and
> the set of `gc` bytes that change must be a subset of what the failure mode predicts.

**Scope note, load-bearing:** the invariant is about *datasets and non-history attrs*, not about the
whole file. `--apply` always appends a `_repair_history` element (§9.1), so the **file sha256 always
changes**, even on a clean file. `gc` is idempotent; the file hash is not. §8.3 depends on this
distinction.

Operationally:

- `--dry-run` is the **default**. It performs the full recompute and reports a per-contig diff
  histogram (`unchanged`, `0 → nonzero`, `nonzero → nonzero` (spanning fragments),
  `255 → x` / `x → 255` (should be **zero** — any occurrence is a hard error)) without opening
  the file for writing and without touching S3.
- On a CLEAN file the expected report is *literally all zeros in every change bucket*. A single
  changed `gc` byte on a known-clean file falsifies the method, immediately, with no reference
  oracle needed.
- On a corrupted file, `255 → x` or `x → 255` transitions must be exactly zero (with the single
  exception of `length == 0` fragments, which are excluded from the histogram entirely — §3.5).

**Where changes are allowed: one rule per region (§2.1).** The doc previously said "anything before
the saturation point is a hard error" with saturation at 2**24, which would have fired on every
file: legitimate `N`-band corrections live *below* 2**24. Define the cumsum thresholds
`T23 = 2**23` and `T24 = 2**24` and classify each fragment by its **`frag_stop` cumsum value**.

**Classify by the simulated float32 accumulator, not by the true float64 cumsum.** This matters and
is easy to get wrong. The 2025-11-24 build crossed `T24` where *its* float32 accumulator reached
2**24 — and per §2.1 that accumulator systematically *under*-counts, discarding the `+0.5` inside
every N run. Across a megabase-scale N block the two positions can diverge by hundreds of thousands
of counts, i.e. by megabases of coordinate. Classifying by the float64 cumsum would put fragments
that the original build already computed under saturation into the restrictive `[T23, T24)` band and
**abort otherwise-repairable files with spurious hard errors**. So the tool simulates the pre-fix
accumulator explicitly — one extra vectorized pass over the same cached per-contig data, replaying
`sum(axis=1).cumsum()` in float32 — and classifies on it. (Widening the band by a measured drift
margin would also work, but the simulation is exact, cheap, and needs no margin to be argued.)

| region (simulated float32 `cumsum(frag_stop)`) | rule | on violation |
|---|---|---|
| `< T23` | **byte-identical required** | hard error — abort the file |
| `T23 <= … < T24` | change permitted **only** if the fragment's span `[start, stop)` contains at least one `N` base | hard error — abort the file |
| `>= T24` | change freely | — |

The `< T23` row **is** the §7.2 pre-saturation oracle — the same check under two names in earlier
drafts. It is stated once, here, and §7.2 explains why the bound is 2**23 rather than 2**24.

The middle-band predicate is cheap: the tool already has the FASTA and `N` is the only non-ACGT
code (§4.2), so a per-contig boolean N-mask plus its prefix-sum answers "does this span contain an
N" in O(1) per fragment. It is also *sound* rather than merely plausible — §4.2 Consequence 5 proves
the float32 cumsum error is piecewise constant, changing only at `N` positions, so an N-free span
has zero error difference and must be byte-identical.

**Do not expect the middle-band changes to be tiny where they occur.** The *count* is small, because
reads essentially never align inside hard-masked blocks. But the *magnitude* is not: where a
fragment straddles an N-block edge the error is up to `0.5 · (N bases in span) / length`, which for
a fragment substantially overlapping an N run approaches 0.5 in GC units — **tens of uint8 units**,
not one. A middle-band diff of 40 is expected behaviour, not a red flag; the flag is a middle-band
diff on an N-free span, of any size.

Because `--apply` is just `--dry-run` plus a write, every production run is preceded by its own
dry run on the same bytes.

---

## 6. Backup and verification

### 6.1 Destination

```
s3://fragmentomics.kariusdx.com/nboley/backups/gc_repair_2026-08/<original key>
```

i.e. the original key path is preserved verbatim under a single backup root, so restore is a
mechanical prefix swap. Backups are deletable once the repair is accepted; the root carries a date
so a future repair gets its own root and can never collide.

**The date in that string is part of a literal, not a computation.** The root is written out
verbatim in the runbook and passed as `--backup-prefix`; the tool **must never derive it from the
current date**, or from anything else evaluated at runtime. A root computed from `date` resolves to
`gc_repair_2026-08` in August and `gc_repair_2026-09` in September, which would make a resumed run
address a *different, empty* backup namespace — and §6.5's write-once protection, which is keyed on
"does this backup key already exist", would then silently pass and let the run copy repaired bytes
into a fresh "backup". §6.6 turns this into an enforced invariant rather than a convention.

**Footprint.** The 218 originals are 91.13 GB. Backups are a full second copy (91.13 GB), and the
repaired originals are larger than what they replace by roughly one compressed real `gc` dataset
each (§3.2). Budget **~190–200 GB** of additional S3 for the duration, released when backups are
deleted after §7.4 stage 7.

### 6.2 Copy mechanism

Server-side `copy_object` — no bytes traverse the client, and it is confirmed working with our
credentials (§3.3). Refuse any object > 5 GB rather than attempt multipart copy; none of the 218
approach that, and building multipart-copy support for a case that does not occur is exactly the
over-engineering §11 refuses.

### 6.3 What counts as "backup verified"

**ETag is not usable as a checksum here.** Objects uploaded via multipart have composite ETags of
the form `<md5-of-md5s>-<partcount>`, which do **not** equal the object MD5; and a single-part
server-side copy of a multipart source produces a *different* (true-MD5) ETag on the destination.
Comparing source ETag to destination ETag would produce spurious failures, and comparing either to
a locally computed MD5 would produce spurious failures in the other direction.

**Use S3 additional checksums instead.** The proposed procedure, per file:

1. Compute `sha256_local` over the downloaded original bytes (we already have them on disk — free),
   and keep the source `ETag` returned by the GET.
2. `copy_object(..., ChecksumAlgorithm="SHA256", CopySourceIfMatch=<source ETag>)`, which makes S3
   compute a SHA-256 on the destination server-side. **`CopySourceIfMatch` is mandatory**: the copy
   is a server-side copy of the *live* object, while `sha256_local` was computed on bytes we
   downloaded earlier. Without the precondition nothing pins those two together, and an object that
   changed in between would be backed up as content we never hashed. With it, a mismatch fails the
   copy with `PreconditionFailed` and the file aborts.
3. `head_object(..., ChecksumMode="ENABLED")` on the destination; compare `ChecksumSHA256`
   (base64) against `sha256_local`, and compare `ContentLength` against the local file size.
4. Backup is **verified** iff both match.

**Must-verify (§12):** that a single-part `copy_object` with `ChecksumAlgorithm="SHA256"` yields a
*full-object* SHA-256 rather than a composite `...-N` value, and that our IAM policy permits the
checksum parameters. Test this on **one** file before the run.

**If the checksum path does not work, the fallback is to re-download and hash *every* backup.** Not
a sample. The requirement is verify-then-overwrite, per file, and a 10% audit would leave ~193 of
218 backups verified by byte count alone before their originals were terminally overwritten — byte
count does not distinguish "correct backup" from "739 MB of the wrong thing". 91 GB of re-download
is affordable (§2.3); it is merely pointless *if* the cheap path works. There is no sampled mode.

### 6.4 Ordering guarantee, and how the upload is verified

Per file, in one worker, no cross-file interleaving of these steps:

```
download → sha256_local + source ETag → preflight → recompute → diff/oracle gates
        → [--dry-run stops here; nothing written anywhere]
        → copy_object(backup, CopySourceIfMatch) → verify backup
        → [barrier: abort file unless verified]
        → write local r+ → verify local → ledger "started"
        → put_object over original → verify upload → ledger "ok"
```

**No original is ever overwritten before its own backup has been verified.** Not "before backups
have started", not "before the backup phase completed" — before *its own* backup verified. The
barrier is per-file and unconditional. Note that the backup now happens *after* the recompute
rather than before it; this is what makes `--dry-run` a genuine no-write mode (§3.1), and it does
not weaken the barrier, which is still immediately upstream of the only write to the original.

**Upload verification cannot reuse the copy path's assumptions.** This is the write over the
irreplaceable original, so it is stated exactly:

- The upload is a **single-part `put_object`** with `ChecksumAlgorithm="SHA256"` and the
  precomputed `ChecksumSHA256` of the repaired local file. Every file is far below the 5 GB
  single-part limit (§3.3).
- `boto3`'s `upload_file` is **forbidden here.** It goes multipart above an 8 MB threshold, and a
  subsequent `head_object` then returns a *composite* `ChecksumSHA256` of the form `<base64>-N`,
  which never equals the full-object SHA-256 — the same trap §6.3 identifies for ETags, one layer
  up. If multipart ever becomes unavoidable, the only acceptable form is multipart with
  `ChecksumType=FULL_OBJECT`; that is not needed at these sizes and is not being built.
- Verification is `head_object(..., ChecksumMode="ENABLED")`, asserting that `ChecksumSHA256`
  contains **no `-N` suffix** and equals `sha256_repaired`, and that `ContentLength` matches.

A `status: "ok"` ledger entry is written only after that assertion passes.

### 6.5 Backup keys are write-once

**The failure mode this closes.** Suppose the upload (§3.1 step 10) succeeds and the process dies
before the ledger append. The key now holds *repaired* bytes and has no `"ok"` record. A resumed
run would treat it as eligible, download it, and — under the original design — unconditionally
`copy_object` it over the backup key. The pristine original would be overwritten by the repaired
bytes and **gone permanently**, because the bucket has no versioning (§2.4). The one file in
flight at crash time is the only dangerous one, and it is exactly the one a "completed files are
fine" argument overlooks.

**The rule:** before any `copy_object` to a backup key, `head_object(..., ChecksumMode="ENABLED")`
that key. **`ChecksumMode="ENABLED"` is not optional** — without it S3 simply omits `ChecksumSHA256`
from the response, and every existing backup would fall into the "no checksum" row below.

| backup key state | action |
|---|---|
| does not exist | perform the copy, then verify (§6.3) |
| exists, `ChecksumSHA256` == `sha256_local` | **skip the copy**; the backup is already correct and verified. Proceed. |
| exists, `ChecksumSHA256` present and differs | **abort this file and stop the run.** Require human review. Do not overwrite, do not proceed. |
| exists, **no `ChecksumSHA256` in the response** | **re-download the backup object and hash it locally**, then apply the matching or differing row above. Never treat an absent checksum as an absent backup. |

The last row is not hypothetical: it is exactly the state produced by §6.3's own fallback, which
creates backups by re-download-and-hash without an S3-side checksum. Without this row a strict
reading sends the fallback path straight to "abort → human review" on every resume, making the
fallback **unresumable**. The re-download costs one object read and is only paid on resume.

The differing-checksum row is the important one: a pre-existing backup with different content means
either a crash-window resume (the local download is the repaired object, and the backup is the true
original — which we must keep) or a key collision (which means the run is misconfigured). Both
demand a human. The §3.1 step-9 `"started"` record makes the first case diagnosable rather than
merely detectable: its presence tells a resumed run "this key was mid-flight; the backup is
authoritative, the current object may already be repaired."

Backup keys are never deleted, never overwritten, and never written to by any path other than this
one.

**This protection is conditional on addressing the same key.** Everything above assumes the resumed
run resolves the *same* backup key as the original run. If it does not, the `head_object` finds
nothing, row 1 fires, and the tool copies the repaired object it just downloaded into a pristine-
looking "backup" — destroying the original silently and terminally. §6.6 closes that.

### 6.6 Run identity: source key, backup prefix and ledger are pinned as a triple

Write-once is **not** a property of a backup key in isolation. It is a property of the triple

> `(source key, backup prefix, ledger)`

held constant across every resume of a run. One mistyped or re-dated `--backup-prefix`, or one fresh
`--ledger` path, bypasses both of §8.2's mechanisms at once. §6.1 pins the prefix as a literal; the
following make it enforceable rather than merely advised:

1. **The ledger has a header record naming its backup prefix.** The first line of a ledger file is
   `{"record": "header", "backup_prefix": "...", "fasta_sha256": "...", "target_list_sha256": "...",
   "tool_version": "...", "created_utc": "..."}`. On resume, if the header's `backup_prefix` differs
   from the `--backup-prefix` passed on the command line, **refuse to start the run.** Not per-file
   abort — the whole run, before the first `head_object`.
2. **Cross-check every prior record, not just the header.** Before any write, scan the ledger for
   records whose `key` appears in the current target list. If any such record's `backup_uri` has a
   root that is not the current `--backup-prefix`, **refuse to start the run.** This catches a
   hand-edited or concatenated ledger that the header alone would not.
3. **The converse case: a populated backup prefix with an empty ledger.** This means a prior run
   very likely existed and its ledger was lost — the single most dangerous state, because every key
   then looks eligible while its backup already holds the pristine original. If a `LIST` of the
   backup prefix returns **any** object and the ledger has **no** records, **refuse to start** and
   emit the count and a sample of the keys found. Proceeding requires the operator to pass
   `--ledger` pointing at the recovered ledger, or to acknowledge explicitly in the runbook and
   re-derive a ledger by listing the backup prefix. There is no flag for this; it is a human step,
   for the same reason §10.3 has no `--yes`.
4. **Refusal is not per-file.** All three checks run once, at startup, before any file is processed.
   A run that would address the wrong backup namespace must not repair even one file first.

Confidence: **high** that these close the resolve-a-different-key hole; **medium** on the exact
shape of the header record, which is a bikeshed and should be settled in the PR.

---

## 7. Validation plan

### 7.1 Clean-file replay (the user-specified test)

Run the repair, in `--dry-run`, against files already known CLEAN and confirm the recomputed `gc`
comes out **byte-identical**. This deliberately avoids needing BAM access: if recomputation
reproduces existing correct values exactly, the method is proven.

**Primary target** (verified by byte-range probe + scaffold fingerprint):

```
s3://fragmentomics.kariusdx.com/nboley/pipeline/output_rebuild_frag_h5s/NC-6724/LB-103925-Lib1/frag_h5s/LB-103925-Lib1.fragments.h5
```
738,378,949 B, LastModified 2026-05-22T18:58:27Z. Probe verdict **CLEAN** (test_n = 10,835,
test_zero_frac = 0.0, control_nonzero_frac = 0.9999); fingerprint **ACCESSION_NAMED → REF-P12**.
This is the only file found that is simultaneously statistically-solid CLEAN, fingerprinted
REF-P12, and full-genome scale. It satisfies both halves of the requirement at once.

**Secondary targets** (same prefix, all fingerprinted REF-P12):
- `.../NC-6724/LB-103927-Lib1/frag_h5s/...` (273,208,801 B, 2026-05-22T17:56:21Z)
- `.../DC4-16772/AC-125674/frag_h5s/AC-125674.fragments.h5` (462,828,323 B, 2026-08-04T22:24:11Z)
- `.../NC-6724/LB-103918-Lib1/frag_h5s/...` (69,303,053 B)
- `.../NC-6724/LB-103913-Lib1/frag_h5s/...` (11,936,498 B — small panel file, probe returned
  INSUFFICIENT_DATA but test_zero_frac = 0.0, directionally clean)

**Pass criteria:**
1. Every recomputed `gc` byte equals the stored byte, on every contig, in every file. Zero
   tolerance — not "99.99%", not "within 1 unit".
2. Every contig absent from REF-P12 (`chrUn_*`, `chr*_random`) is left untouched and was verified
   all-255 by the §5 layer-3 rule. **Scaffold-255 preservation is a pass criterion, not a
   footnote** — a run that "fixes" those has used the wrong reference.
3. The change histogram (§5b) is all zeros.

Note the prefix `s3://fragmentomics.kariusdx.com/nboley/pipeline/output_rebuild_frag_h5s/` holds
489 `.h5` files spanning 2026-01-29 .. 2026-08-04 — **471 pre-fix, 18 post-fix** — so it is also a
convenient source of additional pre-fix (corrupted) test subjects that are *not* in the production
218, if we want a rehearsal target we are willing to break. Recommended: rehearse the full
apply path end-to-end on one of these, in a scratch bucket prefix, before touching the 218.

### 7.2 Pre-saturation prefix agreement (the only available check on R1(b) — and see its limits)

The corrupted files are corrupt only *past* the saturation point. **Everything before it is
correct**, and it was produced by the exact FASTA bytes used on 2025-11-24. That gives us a
per-file reference oracle on the actual targets:

> For each contig, let `S` be the first index where the correct float64 cumsum reaches
> **2**23 = 8,388,608**. Every fragment with `frag_stop <= S` had a bit-exact float32 computation in
> the original build. Recompute those fragments' `gc` and require **byte-identity** with the stored
> values.

This is the same check as the `< T23` row of §5b's region table; it is *implemented* once, as that
row. §7.2 exists to explain the bound and to name the gate.

Why 2**23 and not 2**24: float32 represents integers exactly to 2**24, but *half*-integers only to
2**23, and `N` contributes 0.5 (§4.2). 2**23 is therefore the correct bound for the measured
`{A,C,G,T,N}` alphabet — and the §5 layer-5 gate guarantees that is the alphabet. At accumulator
density 0.41470, 2**23 corresponds to roughly the first **20.23 Mb** of each contig — about
**500 Mb of sequence per file**, across all 25 primaries.

Had B/V/H/D been present, this bound would have been **zero**, not merely smaller: `float32(1/3)`
sits on no half-integer grid, so a handful of such bases destroys exactness immediately and the
oracle would cover nothing at all. That is why §5 layer 5 is a gate and not a log line.

This check:
- **strongly corroborates reference identity (R1(b))** against the 2025-11-24 build, on the real
  target files, which the §7.1 clean replay cannot do (its subjects were built 2026-05-22, six
  months on the wrong side);
- simultaneously validates the coordinate convention, the `'a'` pad offset, Stage A rounding and
  Stage B quantization, per file;
- costs nothing extra — the cumsum is already computed.

**It runs on all 218 as a mandatory preflight gate**, and a single mismatched byte in the
pre-saturation region **aborts that file**.

#### 7.2.1 What this check does *not* prove

It proves agreement over the prefix — roughly the first 20.23 Mb of GC-cumsum per contig — which is
precisely **the region that needs no repair**. The bytes the tool actually writes are derived from
sequence past ~40.45 Mb: the untested complement. The two candidate references differ by ~11.9 Mb of
hard-masking (§5), and **that masking could lie entirely past 20.23 Mb**, in which case all 218 files
pass the oracle and the repair still writes 127-over-real-GC across the differing blocks. This is
not hypothetical enough to ignore; the oracle is corroboration, not proof of reference identity.

R1 in §4.3 is actually **two questions**, and only one of them is answerable here:

- **(a) the two-reference confusion** — did the 2025-11-24 build use REF-P12 or REF-ASSETS? Both
  candidates are in hand, so this is decidable, and the check below decides it *outright* rather
  than narrowing it.
- **(b) arbitrary drift of the REF-P12 object itself** — is the object at
  `nboley/resources/GRCh38.p12.genome.fa.gz` today byte-identical to what it was on 2025-11-24? The
  bucket has no versioning, so there is no candidate to diff against and **no check in this document
  closes (b)**. The pre-saturation prefix oracle is the only evidence, and it is corroboration, not
  proof. (b) remains open and §12.1 says so.

#### The check that closes (a): recompute under both references and diff

The ~11.9 Mb / 18-contig differential between REF-P12 and REF-ASSETS is a **measured, positionally
known** quantity, and both FASTAs are available. So do not test a proxy — compute the answer:

> For every contig the tool intends to write, recompute the `gc` bytes of the written region under
> **REF-P12** and under **REF-ASSETS**, and report **the number of differing bytes**.

- If the count is **0** for a file, then reference choice **provably does not affect any byte this
  tool writes to that file**, and question (a) is closed for that file. Not narrowed — closed.
- If the count is nonzero, the report names the exact blast radius, per contig, *before anything is
  written*, and the operator decides with a number rather than a hunch.

Cost is one extra float32/float64 cumsum pass over the same cached per-contig data (§3.4) plus one
extra FASTA in the cache — affordable at the scale of §3.4's estimates, and it produces a count, not
a heuristic. This is the **primary form** of the check and is a mandatory preflight report in
`--dry-run`.

**Supporting per-contig statistics**, reported alongside the byte-diff because they localize a
nonzero result:

  (a) differing bases in the pre-saturation prefix that have **nonzero fragment coverage** — i.e.
      how much of the oracle's region actually discriminates between the two references;
  (b) differing bases inside the region the tool actually writes;
  (c) **the number of fragments in the written region that overlap a differing position** — the
      denominator of the exposure.

If (c) == 0 for a file, the byte-diff above is necessarily 0 and the two results corroborate.

**What was here before, and why it is gone.** An earlier draft required that each of the 18
disagreeing contigs have "at least one hard-masked `N` block in its pre-saturation prefix, with the
recomputed `gc` agreeing on fragments overlapping it." That check is unsound in two ways and has
been deleted rather than caveated: it tested `N`-ness in *one* FASTA rather than *difference between
the two*, so a shared telomeric N block at 0–10 kb satisfies it while every differing block sits at
60 Mb; and it imposed no minimum fragment count, so — since reads essentially never align inside
hard-masked blocks — the expected number of overlapping fragments is ~0 and zero comparisons
"agree" vacuously. It could pass on literally no evidence.

Contigs with a nonzero byte-diff are listed by name and are a **blocking** item for the §7.4 stage-3
gate: a human decides whether to proceed, with the exposure quantified.

Confidence: **high** that the prefix check is sound over its region; **high** that the two-reference
byte-diff decides question (a) exactly; **medium** that everything passes cleanly on all 218 first
try; and **explicitly not** a claim that the repaired region is verified against the 2025-11-24
build — question (b) has no oracle, and if the prefix check passes we have learned less than "the
FASTA is right".

### 7.3 Post-repair checks on a corrupted file

After writing, for each repaired file. **These are not all the same kind of check**, and the
previous version of this list did not distinguish them — "the GC histogram looks unimodal" and
"zero `255 → x` transitions" are not both grounds to stop a rollout. Split accordingly.

**Blocking — any failure stops the rollout and triggers restore-from-backup:**
- reopen through `FragmentsH5` succeeds (exercises the reader, `eval()` of `_contig_lengths_str` at
  `fragments_h5.py:305`, the `"gc" not in ...` guard at `:562-564`, and the `/254.0` inversion plus
  255→NaN remap at `:570-574`);
- `gc` dataset shapes and dtypes unchanged (`uint8`, same length as `starts`);
- **zero** `255 → x` or `x → 255` transitions, excluding `length == 0` fragments (§3.5);
- every change obeys the per-region rule of §5b: none below `2**23`; those in
  `[2**23, 2**24)` only on fragments whose span contains an `N`; unrestricted at or above `2**24`;
- every non-`gc` dataset, and every attr **except `_repair_history`**, byte-identical to the
  pre-mutation hashes recorded at §3.1 step 2;
- `_repair_history` parses as JSON and has gained **exactly one** element;
- the §7.2 pre-saturation oracle still passes on the repaired file (it must, since the tool does
  not write there — a failure means it wrote somewhere it should not have);
- re-running the tool on the repaired file changes **no `gc` byte** (§5b) — the cheapest possible
  idempotency proof. Note the precise claim: **the `gc` datasets are idempotent, the file sha256 is
  not**, because a re-run under `--apply` appends a second `_repair_history` element. Run this check
  in `--dry-run`, which appends nothing.

**Advisory — recorded in the report and reviewed by a human, but do not by themselves halt:**
- `0 → nonzero` fraction consistent with the predicted post-saturation fragment count;
- resulting genome-wide GC distribution is unimodal near the measured **0.4102** (§4.2) with no
  residual spike at 0 — a distribution sanity check, not a correctness proof;
- file size within 5% of a **from-scratch correct build** (§3.2) — *not* within 5% of the corrupted
  original, which it will exceed substantially and by design. Report both the from-scratch delta and
  the absolute growth over the replaced object. Advisory because bloat is a storage-cost problem,
  not a data-integrity problem; but exceeding the 5% threshold triggers the `h5repack` decision in
  §3.2 before the next stage, so it gates *progression* even though it does not gate *this file*.

### 7.4 Graduated rollout

| Stage | Target | Gate to proceed |
|---|---|---|
| 0 | Unit + regression tests, incl. a new "overwrite a dataset in `r+`" test, the §3.5 rounding-agreement test over >=10^7 fragments, a §6.5 write-once-backup test, a §6.6 refuse-to-start test (mismatched ledger header, differing `backup_uri` root, populated-prefix-with-empty-ledger), and a §5b float32-accumulator-simulation test asserting the simulated band boundary lands where the pre-fix code saturates; the existing overflow regression test at `8820299` (`tests/test_gc_cumsum_overflow.py`, ~676 MB RSS) must still pass | green CI |
| 1 | `--dry-run` clean replay on the 5 §7.1 targets | §7.1 criteria, all zeros |
| 2 | Full `--apply` rehearsal on **one pre-fix file from `output_rebuild_frag_h5s/`**, written to a **scratch prefix**, not over the original | §7.3 blocking checks all pass; file-size advisory reviewed; backup/restore round-trip exercised; first real runtime datum recorded (§3.4) |
| 3 | `--dry-run` §7.2 pre-saturation gate across **all 218**, incl. the §7.2.1 two-reference byte-diff report | 218/218 pass the prefix oracle; every file with a nonzero REF-P12/REF-ASSETS byte-diff named, quantified, and accepted by a human |
| 4 | `--apply` on **1** of the 218 | §7.3 all pass; manual review of the diff report |
| 5 | `--apply` on **10** | §7.3 all pass |
| 6 | `--apply` on the remaining 207 | — |
| 7 | Post-run: re-probe all 218 with the audit detector; spot-restore 3 files from backup and confirm they match the recorded `sha256_local` | before backups are deleted |

Stage 3 is deliberately placed before *any* production write: it is cheap (recompute only, no
writes, and the oracle needs only ~20.23 Mb of each contig's fragments compared) and it is the check
most likely to surface an unmodelled difference.

A note on the existing regression test: `OVERFLOW_THRESHOLD = 2**24`
(`tests/test_gc_cumsum_overflow.py:37`) is the right constant *for that test* and needs no change —
its contig is all-`G` with a short `A` patch, so it contains no half-integer bases and the
`[2**23, 2**24)` band of §2.1 simply does not arise. Its docstring should say so, so that a future
reader does not "fix" it to 2**23 or conclude that 2**24 is the universal bound.

---

## 8. Failure and rollback

### 8.1 Partial failure within one file

Cannot reach S3. The `r+` write is on a local scratch copy; the original object is untouched until
a complete, verified file is uploaded. An interrupted write leaves a corrupt file in `/tmp`;
recovery is `rm` + retry. This is the whole reason for §3.3's shape.

### 8.2 Partial failure across the 218

Each file is independent. A run that dies after 137 files leaves 137 repaired with `ok` records and
81 untouched with no records.

**The completed files are not the interesting case; the one in flight is.** A crash can land
between the upload (§3.1 step 10) and the `"ok"` append (step 11), leaving a key that holds
repaired bytes and advertises nothing. Two mechanisms cover it, and neither is optional:

- the `status: "started"` record written at step 9, *before* the upload, which names the verified
  `backup_uri` and both hashes — so a resumed run can tell "mid-flight" from "never started";
- **write-once backup keys** (§6.5), which mean that even if the `"started"` record is itself lost,
  a resumed run that addresses the **same backup key** cannot overwrite a pristine backup with the
  repaired object it just downloaded.

**That second mechanism is conditional, and the condition has to be enforced.** Write-once protects
a *key*, and the key is `<backup-prefix>/<source key>`. A resume that supplies a different
`--backup-prefix` — a re-dated root, a typo, a copy-pasted line from a different runbook — resolves
a key that does not exist, §6.5 row 1 fires, and the repaired bytes are copied into a fresh backup
namespace. The pristine original is then gone, silently and terminally. Likewise a fresh `--ledger`
path erases the `"started"` record's protection.

So the guarantee, stated exactly: **write-once holds for the triple (source key, backup prefix,
ledger) held constant across resumes.** §6.6 makes that constancy an enforced startup precondition —
ledger header pins the prefix, prior records are cross-checked, and a populated backup prefix with
an empty ledger refuses to start. Without §6.6, one flag value bypasses both mechanisms above at
once.

Without all three, a resume would silently destroy the only surviving copy of the original. The
earlier claim that "nothing is left half-done" was reasoning about completed files and ignoring the
in-flight one, which is the only dangerous one.

### 8.3 Ledger, idempotency and resumability

Append-only JSONL at a `--ledger PATH` on **NFS home** (persistent), not `/tmp`. The **first** line
is the §6.6 header record binding this ledger to one backup prefix, one FASTA and one target list;
a run whose `--backup-prefix` disagrees with the header refuses to start. Thereafter **two** records
per file: `"started"` immediately after the backup is verified and the local file is repaired but
*before* the upload, and `"ok"` after upload verification.

```json
{"record": "header",
 "backup_prefix": "s3://fragmentomics.kariusdx.com/nboley/backups/gc_repair_2026-08/",
 "fasta_sha256": "...", "target_list_sha256": "...",
 "tool_version": "2.13.0", "created_utc": "..."}
```

```json
{"key": "nboley/ibd_v2/build_frag_h5s/.../X.fragments.h5",
 "status": "started",
 "sha256_original": "...", "sha256_repaired": "...",
 "backup_uri": "s3://.../gc_repair_2026-08/...", "backup_verified": true,
 "fasta_sha256": "...", "tool_version": "2.13.0",
 "started_utc": "..."}
```

```json
{"key": "nboley/ibd_v2/build_frag_h5s/.../X.fragments.h5",
 "status": "ok",
 "sha256_original": "...", "sha256_repaired": "...",
 "backup_uri": "s3://.../gc_repair_2026-08/...",
 "fasta_sha256": "...", "tool_version": "2.13.0",
 "contigs_repaired": 25, "contigs_skipped_absent_from_fasta": 170,
 "bytes_changed": 12345678, "repair_history_len": 1,
 "started_utc": "...", "finished_utc": "..."}
```

`tool_version` is read from `importlib.metadata.version("fragments-h5")` at runtime, never
hardcoded — a hardcoded string is provenance that can silently disagree with the code that wrote it.

**How a resumed run classifies each key:**

| ledger state | resume behaviour |
|---|---|
| one or more `"ok"` records | re-check the live object's SHA-256 against the `sha256_repaired` of the **most recent** `"ok"` record for that key. Match → skip. Mismatch → **stop the whole run**, human review. |
| `"started"` but no `"ok"` | the crash window. The backup is verified and authoritative; the live object may be original or repaired. **Stop and require human review** — do not auto-resume this key. |
| `"failed"`, or no record | eligible — but §6.5's write-once backup check still runs first and will halt if a backup already exists with different content. |

**Why "most recent" and not "the" `"ok"` record.** An `--apply` re-run over an already-repaired key
appends a second `_repair_history` element, which changes the file sha256 even though it changes no
`gc` byte (§5b, §7.3). If the resume logic compared against the *first* `"ok"` record's
`sha256_repaired`, a legitimate second repair would permanently mismatch and escalate to "stop the
whole run" — the tool bricking its own run with a hash it caused itself. Comparing against the
latest record makes the ledger track the file's actual current state. Two guards keep this honest:
`--target-list` rejects duplicate keys at load time (§10.2), so no two workers can race the same
key within a run; and each `"ok"` record carries `repair_history_len`, so a resumed run can tell
"this key was repaired twice deliberately" from "this key was written by something else."

Verification uses `head_object(..., ChecksumMode="ENABLED")`. **There is no ETag fallback.** The
ETag is unsound here for exactly the multipart reason §6.3 gives, and `Range: bytes=0-` is a
full-object GET, not a free metadata call. If the checksum path is unavailable, re-download and
hash — same rule as §6.3, and for the same reason.

This is deliberately stricter than "skip anything in the ledger": the ledger asserts a fact about
S3, and we verify the fact rather than trusting it.

Belt-and-braces: the repair is naturally idempotent **in the datasets** — running it twice on a
repaired file changes no `gc` byte (§7.3), because the recompute is total rather than conditional on
the current values. It is **not** idempotent in the file hash: a second `--apply` appends a second
`_repair_history` element and therefore a different `sha256_repaired`. That is by design (the
history is a list precisely so a second repair is recorded, §9.1) and is why the `"ok"` row above
compares against the latest record.

### 8.4 Restore from backup

**There is no `--restore` flag.** Restore is a documented runbook procedure, not a feature.

The reasoning is the doc's own: a restore path that depends on the tool that caused the problem is
not a restore path, so the manual `boto3` procedure has to be written down and rehearsed
*regardless*. Once it is written down, it is ten lines — `copy_object` from `backup_uri` back over
`key` for each ledger record, then `head_object` and compare against `sha256_original` — and
shipping a second implementation of those ten lines inside the tool adds a code path that is
exercised only in an emergency. Ship the runbook; the runbook is the tested artifact.

Because backup keys mirror original keys under one pinned root (§6.1, §6.6), the restore is a
mechanical prefix swap, and the ledger — whose header names that root — supplies the root, the key
list and the expected hashes. §7.4 stage 2
rehearses it end-to-end and stage 7 spot-checks it on 3 files.

Backups are retained until §7.4 stage 7 completes and the cohort owners have re-run at least one
downstream analysis. Deletion is a separate, explicitly-approved action.

### 8.5 One thing we cannot roll back

If the FASTA is wrong in a way that §7.2 does not catch, we overwrite 218 files with wrong-but-
plausible GC and restore from backup — which puts back the *corrupted* originals. That is
recoverable. The unrecoverable case is losing the backups, which is why §6.4's per-file barrier is
unconditional, why §6.5 makes backup keys write-once, and why §6.6 pins the backup prefix and ledger
so that write-once cannot be bypassed by a resume addressing a different key.

Note that §7.2 is *weaker* than this section previously assumed: it constrains the reference only
over the unrepaired prefix (§7.2.1), so "the FASTA is wrong in a way §7.2 does not catch" is a live
possibility. But the specific worry that the ~11.9 Mb of differential masking sits past the oracle's
reach is now **decidable**, not merely mitigated: §7.2.1's two-reference byte-diff recomputes the
written region under both candidates and reports how many bytes actually differ. Where that count is
zero, the REF-P12/REF-ASSETS question is closed for that file. What remains unrecoverable-by-check
is arbitrary drift of the REF-P12 object itself, for which the backstop is simply that this case
*is* recoverable from backup.

---

## 9. Provenance

The 218 files carry **no build provenance at all** — `_build_version` / `_build_argv` landed in
v2.12.0 (`a3be97b`, 2026-08-21) and these files predate it. The current authority,
`docs/architecture/fragment_selection_and_build_provenance.md`, records those as flat root attrs
and **defines no repair-provenance field**. (Note: `docs/pending/build_provenance_metadata.md` is
marked REJECTED/superseded and should not be used as a reference.)

### 9.1 Proposed attr

A single flat root attr, JSON-encoded, consistent with the existing `_build_argv`-is-JSON
convention:

```
_repair_history : str   # JSON array of objects, appended to on each repair
```

Each element:

```json
{"tool": "repair-fragments-h5-gc",
 "version": "2.13.0",
 "argv": ["repair-fragments-h5-gc", "--fasta", "...", "--apply", ...],
 "timestamp_utc": "2026-08-25T17:04:11Z",
 "reason": "gc-float32-cumsum-saturation; see docs/pending/gc_repair_tool.md; fixed by 95c76f5",
 "datasets": ["data/*/gc"],
 "fasta_uri": "s3://fragmentomics.kariusdx.com/nboley/resources/GRCh38.p12.genome.fa.gz",
 "fasta_sha256": "...",
 "source_sha256": "...",
 "backup_uri": "s3://fragmentomics.kariusdx.com/nboley/backups/gc_repair_2026-08/..."}
```

`version` is `importlib.metadata.version("fragments-h5")`, read at runtime.

A *list* rather than scalar attrs, because "this file has been repaired more than once" is a real
future state and rewriting scalars would erase history. This is the main piece of forward-looking
generality being built in, and it costs ~5 lines.

**When it is written: §3.1 step 7**, in the same `r+` handle as the `gc` datasets, before the file
is closed. The ordering is not incidental — `sha256_repaired` is taken at step 8 over the closed
file and recorded at step 9, so it must already cover the history attr. Two consequences follow and
are stated wherever they bite: `_repair_history` is **excluded by name** from every "byte-identical
to the original" assertion (§3.1 step 8, §5b, §7.3), and the file sha256 is **not** idempotent
across repeated `--apply` runs even though the `gc` datasets are (§8.3).

### 9.1.1 Known limitation: this attr does not survive derivation

`_repair_history` will be **dropped by real downstream paths**, and we know this in advance rather
than discovering it later. `docs/architecture/fragment_selection_and_build_provenance.md:720-739`
documents that two consumers rebuild fragment H5s from a hardcoded 5-attribute allowlist and
therefore "will silently drop `_build_argv` and `_build_version`" — one of them,
`liftover_porcine_fragments.py:193-197`, is a real analysis path, not a test fixture. A sixth attr
is dropped by exactly the same mechanism.

Two consequences, neither quietly absorbed here:

1. **Derived H5s already produced from the 218 are not repaired by this tool and will not advertise
   that they need to be.** They carry corrupted GC and no `_repair_history`, and they are
   indistinguishable from a clean derivation. Finding them is a separate task with a separate
   inventory; the tool's target list covers the 218 originals only. Tracked as **§12.9**.
2. Any derivation performed *after* the repair loses the record that the repair happened. Widening
   the allowlist is the fix and it is **out of scope for this PR and not tracked in §12 at all** —
   it belongs as a filed issue against
   `docs/architecture/fragment_selection_and_build_provenance.md`, which owns the allowlist.

(§12.10, the downstream-impact bucket, is about already-published CD/UC *results*, which is a third
and separate thing.)

We are not widening the attr surface to work around this. The honest position is that the
provenance is correct where it is written and lossy where someone else copies it.

### 9.2 Do not backfill `_build_version` / `_build_argv`

**No.** We do not know the original build's argv or version. Writing plausible-looking values would
be **fabricating provenance** — worse than absent provenance, because absent provenance is honestly
absent whereas fabricated provenance is indistinguishable from real. Downstream code that reasons
"`_build_version >= 2.7.2` therefore GC is trustworthy" would then be reasoning from a value we
invented.

The correct signal is the opposite one: `_repair_history` present ⇒ this file's `gc` was produced
by v2.13.0's recompute path, and `_build_version` absent ⇒ we genuinely do not know how the rest
of it was built. If a downstream consumer needs "is the GC trustworthy", it should read
`_repair_history`, and we should say so in the architecture doc.

**Do not add a `_gc_source` scalar.** An earlier draft proposed one as a cheap alternative for
consumers that do not want to parse JSON. It is two sources of truth for one fact, and the failure
mode is the pair disagreeing after a second repair — precisely the situation `_repair_history` was
made a list to handle. One field, one truth. Consumers parse the JSON.

Confidence: **high** on the do-not-fabricate call; **medium** on the exact attr names, which should
be reviewed against `docs/architecture/fragment_selection_and_build_provenance.md` and updated
there in the same PR.

---

## 10. Execution-approval gate

Policy is not a mechanism. The mechanism:

1. **`--dry-run` is the default.** Omitting all flags produces a report and touches nothing. There
   is no way to write by accident. **`--dry-run` exercises the S3 read paths it will depend on**:
   it performs the read-only `head_object` probe of each backup key (§6.5) and the §6.6 startup
   cross-checks, and reports what it found. Without this, "backup key exists with different
   content" — the condition that halts a production run — would first be discovered mid-`--apply`,
   which is the worst possible time to learn it.
2. **`--apply` is required to write**, and `--apply` additionally requires **all** of:
   - `--target-list FILE` — an explicit newline-delimited list of full S3 keys. **Prefix expansion
     / globbing is disabled in `--apply` mode.** You cannot say "everything under this prefix"
     while writing; you must have materialised the list, reviewed it, and checked it in.
     **Duplicate keys are rejected at load time**, in both modes: two workers on one key produce two
     `_repair_history` timestamps, and the second upload invalidates the first ledger record's
     `sha256_repaired`, which §8.3's `"ok"` row escalates to "stop the whole run". The list is a
     set; enforce it where it is read.
   - `--expect-fasta-sha256 <hex>` — pins the reference bytes (§5).
   - `--backup-prefix <s3-uri>` — no default, so a run with no backup destination cannot start. It
     is a **literal from the runbook, never computed** (§6.1), and it is cross-checked against the
     ledger header on every resume (§6.6). Getting this value wrong is the one remaining way to
     lose data, so it is validated before the first file rather than trusted.
   - `--ledger PATH` — the same ledger across all resumes of a run; §6.6 makes a mismatched or
     empty-but-should-not-be ledger a refusal to start.
   - `--max-files N` — hard cap, **default 1**. Scaling from 1 → 10 → 218 requires deliberately
     raising this on the command line, which makes the rollout stages of §7.4 mechanical rather
     than aspirational.
3. **Interactive confirmation, with no bypass.** Before the first write, print the plan (file count,
   total bytes, FASTA sha256, backup prefix, and the sha256 of the target list itself) and require
   the operator to type the literal count of files, e.g. `218`. **There is no `--yes` flag and no
   `--i-know-this-is-production` escape hatch.** An earlier draft had both, gated so that `--yes`
   was refused on the production bucket unless the second flag was also passed; that is two flags
   and a conditional to express "always prompt", and every bypass eventually gets pasted into a
   shell history and reused. One unconditional prompt is stronger *and* simpler. The rehearsal in a
   scratch bucket (§7.4 stage 2) types the number like everyone else — it is one file.
4. **No production identifiers in the package.** No bucket names, prefixes, FASTA URIs, or file
   lists are hardcoded in `src/`. They live in the target list and the runbook. This is what makes
   the tool "general" in the only sense that matters here — it has no opinion about which files
   exist.
5. **Human gate.** The target list file and the runbook are committed to the repo and approved in
   the PR. The user's approval to *run* is recorded there, separately from approval of this design.

---

## 11. Scope: what we deliberately do not build

The user's instruction — *"a general repair/update script going forward. Don't over-engineer but
also do not treat as a throw-away"* — resolves as: **general in interface, specific in behaviour.**

**Generality we are building in** (all cheap, all load-bearing for the current task anyway):
- A flat, config-free CLI with no hardcoded targets (§10.4), so pointing it at a different file set
  is a command-line change.
- The download → backup → verify → mutate → verify → upload → ledger pipeline as the reusable
  skeleton. The *only* thing specific to this incident is one function, "recompute `gc` for one
  contig", plus its preflight rule.
- `_repair_history` as an append-only list (§9.1), so a second repair does not erase the first.
- The ledger + resume machinery (§8), including write-once backup keys (§6.5) and the pinned
  `(source key, backup prefix, ledger)` triple (§6.6), which any future in-place S3 mutation will
  need.
- The §5b no-op invariant expressed as `--dry-run`-by-default, which is a property of the skeleton,
  not of GC.

**Refused, explicitly:**
- **A `fragments_h5.gc_audit` module in `src/`.** An earlier draft proposed promoting the
  `/tmp/gcpilot/probe.py` statistics into the package. Detection is a **finished job** — it has
  already run, and its output is a committed TSV. Shipping and maintaining a module for a job that
  will not run again is over-engineering by the plain meaning of the word. If a future incident
  needs detection, it will need a *different* detector.
- **A `--restore` flag** (§8.4) — the runbook is the restore path.
- **A `--yes` / `--i-know-this-is-production` bypass** (§10.3) — one unconditional prompt.
- **An override for the `255 → x` hard error** (§5, §5b). An earlier draft let `--apply` proceed
  past a scaffold-wide `255 → x` diff "with an explicit override flag". That transition means the
  wrong reference was supplied; there is no legitimate use for continuing. Deleted, not gated.
- **A `_gc_source` scalar attr** (§9.2) — one source of truth.
- A plugin/registry system for "repair kinds". There is exactly one repair. When there is a second,
  the honest refactor is to add a second function, not a framework.
- A generic HDF5 schema-migration engine.
- Support for repairing datasets other than `gc`. The recompute function is isolated so a future
  one can be added; nothing is parameterised in advance.
- Multipart `copy_object` / >5 GB support (§6.2) — refused with a clear error instead.
- Re-deriving anything from BAMs. The whole point of the clean-replay validation is that we never
  need BAM access.
- Bucket scanning / automatic corruption discovery. Detection already happened; its output is a
  TSV. The tool consumes an explicit list.
- AWS Batch / distributed execution. 91 GB and 16 CPUs is a laptop-scale job.
- Any emulation of object versioning, retention policy automation, or lifecycle rules.
- A progress UI beyond structured log lines.

**Prototype disposition** — everything in `/tmp` is not durable and `/tmp` will be cleared. The
rule applied: **commit evidence, not scratch code.**
- `/tmp/ibd_manifest_gc_audit.tsv` (763 per-file verdicts) — **irreplaceable evidence**, cost real
  S3 reads to produce. **Commit it** under `docs/pending/gc_repair/ibd_manifest_gc_audit.tsv`.
- The measured REF-P12 sequence-byte histogram (§4.2) — likewise irreplaceable-ish (a ~5 min
  decompress pass over 939 MB) and load-bearing for the entire §4.2 argument and the §5 layer-5
  gate. **Commit it** under `docs/pending/gc_repair/ref_p12_sequence_histogram.json`, keyed by the
  FASTA's sha256.
- `/tmp/gcpilot/{probe.py,s3h5.py,fallback.py}` and `/tmp/fasta_fingerprint.py` — **do not commit.**
  Scratch stays scratch. They produced the two artifacts above and have no further use: the repair
  tool operates on local files and needs no byte-range S3 reads, and the fingerprint script is
  superseded by the strictly more general §5 layer-3 contig-membership rule. Committing them under
  `prototypes/` with a "these are unmaintained" README just creates code that a future reader has
  to evaluate and discard. The evidence they produced is what has value, and it is committed.

---

## 12. Open questions and unverified claims

Ordered roughly by how much they could hurt. Entries that a previous draft listed here and that have
since been **discharged from source or by measurement** are not restated as "resolved risks" — they
are gone, and the evidence is in the §4.3 table. Specifically: R2 (B/V/H/D presence), R3
(`stop = start + length`), and the three-Stage-A-sites question are all closed. Padding this
register with non-risks makes the real entries harder to weight.

1. **R1(b): arbitrary drift of the REF-P12 object itself is unprovable, and §7.2 does not close it.**
   The bucket has no versioning and `ListObjectVersions` is `AccessDenied`, so the object's history
   is unrecoverable and there is no candidate to diff against. §7.2 pins the FASTA over the
   pre-saturation prefix — **the region the tool does not write** — which is corroboration, not
   proof. This entry stays open permanently; the backstop is that the failure is recoverable from
   backup (§8.5).
   *Note the narrower question, R1(a) "REF-P12 or REF-ASSETS?", is now decidable rather than merely
   mitigated:* §7.2.1's two-reference byte-diff recomputes the written region under both candidates
   and reports the differing-byte count, closing it outright per file where that count is zero. It
   is still a *plan*, not a result — nothing has been run.
2. **No real fragment H5 has been inspected on this machine** (checked `/tmp`, `$HOME`, `tests/`).
   *All* schema knowledge in this doc is read from source. Every structural claim — group layout,
   dataset names, dtypes, chunking, index blocks — is therefore unverified against reality.
   **Stage 1 of §7.4 is also the first time we look at a real file.**
3. **The `round(x, 5)` vectorization is unimplemented and unbenchmarked** (§3.5). We know NumPy's
   and CPython's five-place rounding can disagree and that the disagreement is observable in the
   stored byte; we do not yet know whether the exactly-correct elementwise path is fast enough, and
   the fallback requires an agreement test that has not been written.
4. **HDF5 in-place bloat is unverified at real scale** (§3.2) — one synthetic single-dataset
   experiment. Gated in §7.4 stage 2/4 with a 5% threshold, and `h5repack`'s availability on this
   machine is itself unverified.
5. **S3 `copy_object` + `ChecksumAlgorithm="SHA256"` behaviour is unverified** — specifically
   whether it returns a full-object or composite SHA-256, and whether our IAM policy allows it
   (§6.3). Test on one file first. The fallback (re-download and hash **all** backups) is expensive
   but sound; there is no sampled fallback.
6. **`rebuild_all_frag_h5s.nf`'s FASTA pinning is not verified from source.** The `.nf` lives in the
   separate `omni` monorepo (`omni/pipeline/fragmentomics/specialized_workflows/`) and is not on
   this machine; the claim comes from `AGENT_CONTEXT.md:~668-709`. It is corroborated empirically —
   all three probed `output_rebuild_frag_h5s/` files fingerprint REF-P12, pre- *and* post-fix — but
   corroboration is not verification. This matters only for §7.1's target selection, where the
   fingerprint is direct evidence anyway.
7. **The middle-band change volume is predicted, not measured.** §5b permits `[2**23, 2**24)`
   changes only on N-containing spans. That the permitted set is *sufficient* now has a proof
   (§4.2 Consequence 5: the float32 cumsum error is piecewise constant and changes only at `N`
   positions, so an N-free span has zero error difference), but the proof has not been checked
   against a real file, and neither has the **magnitude** prediction — up to tens of uint8 units on
   a fragment straddling an N-block edge. Stage 1 (clean replay, where the expected middle-band
   change count is zero) and stage 3 both exercise it.
   **Also unmeasured: the float32-accumulator simulation** that §5b uses to classify fragments into
   bands. Classifying on the true float64 cumsum would misplace the `T24` boundary by the
   accumulated N-drift and abort files spuriously (§5b); the simulation is exact in principle but
   has not been written or run.
8. **Whether the per-dataset raw-chunk hashes of §3.1 step 2 can be taken cheaply in h5py.** The
   operand question is settled — hashes are taken **before** mutation, since the backup is in S3 and
   the pre-mutation local bytes do not survive step 7 — but the cost is not.
   `Dataset.id.read_direct_chunk` can read raw compressed chunks without decompressing, which would
   make this nearly free, but I have not confirmed the API surface, that it works for all filter
   configurations here, or how it behaves for contiguous (unchunked) datasets.
9. **Derived H5s built from the 218 are out of scope and unenumerated** (§9.1.1). At least one real
   analysis path (`liftover_porcine_fragments.py:193-197`) rebuilds fragment H5s from a 5-attribute
   allowlist, so any such derivative carries corrupted GC, will not be repaired by this tool, and
   will not advertise that it needs to be. Someone has to inventory them; this doc does not.
10. **Downstream impact of changed GC on already-published CD/UC results** is out of scope for this
    doc but is not out of scope for the project. Noted so it does not fall through.

---

## 13. Self-assessment

**What is solid.** The arithmetic is now settled by measurement rather than argued from assumption.
The REF-P12 sequence alphabet is exactly `{A,C,G,T,N}` (§4.2), so the reachable per-base C+G set is
`{0, 0.5, 1.0}`, the float32 pairwise add at `fragment.py:442` is exact, the float64 cumsum is exact
with a margin of ~4×10^7, and the 2**23 / 2**24 boundaries are derived rather than guessed. The
middle band now has an actual proof behind its rule rather than a plausibility argument: the float32
cumsum error is piecewise constant and changes only at `N` positions (§4.2 Consequence 5), which is
exactly why "change permitted only on N-containing spans" is sound. That argument is conditional on
the reference, and §5 layer 5 turns the condition into a gate that fails closed. Three of the
original residual assumptions (R2, R3, and the three-Stage-A-sites question) are discharged from
source or measurement, with citations.

The safety engineering is also solid. The download-repair-upload shape makes mid-write corruption
in S3 structurally impossible. The per-file backup barrier (§6.4) makes "overwrote before backing
up" structurally impossible; write-once backup keys (§6.5) close the crash-window hole that would
otherwise have let a resumed run destroy the only pristine copy; and §6.6 closes the narrower
version of the same hole, where a resume supplies a different backup prefix or ledger and thereby
addresses a key the write-once check has never seen. `--dry-run` is now genuinely write-free and
also probes the backup keys it will later depend on. The no-op invariant (§5b) means every
production run is preceded by its own dry run on the same bytes.

**The weakest point, first of two.** *Nothing in this document has been checked against a real
fragment H5.* Not one. The entire schema model — that `gc` lives at `data/<contig>/gc`, that it is
`uint8` and parallel to `starts`, that the scaffold contigs are uniformly 255, that in-place rewrite
does not bloat — is read off source code and synthetic experiments. The design is therefore a
well-reasoned structure resting on an unvalidated foundation, and the honest expectation is that
**stage 1 of §7.4 will surface at least one thing this doc gets wrong.** The plan is built to absorb
that (dry-run default, per-file abort, graduated rollout, ledger, backups), but readers should
weight the confidence labels accordingly: they describe confidence in the *reasoning*, not in the
*facts about the files*.

**The weakest point, second of two.** §7.2 was once described as *decisive* for reference identity.
It is not. It pins the FASTA over the pre-saturation prefix, which is **exactly the region the tool
does not write**; every byte the repair emits derives from sequence past ~40.45 Mb, which the oracle
never touches. That limit is structural and no amount of engineering removes it.

What *has* improved is the part of the reference question that was always answerable. §7.2.1 no
longer offers an N-block coverage heuristic — which tested the wrong predicate (N-ness in one FASTA
rather than difference between the two) and could pass on zero fragments — and instead recomputes
the written region under **both** candidate references and reports the differing-byte count. Where
that count is zero, the REF-P12/REF-ASSETS question is *closed* for that file, not narrowed. What
remains open is arbitrary drift of the REF-P12 object itself (§12.1), which has no candidate to diff
against and therefore no possible check. The backstop for that case is that it is recoverable from
backup — which is why §6.5 and §6.6 matter more than they look.

**Grade: B+/A-, and the honest placement is B+.** The correctness argument is genuinely strong —
measured rather than assumed, with its one conditional enforced mechanically, and with the
middle-band rule now proved rather than motivated. The safety design has survived two rounds of
adversarial review and closed everything found: dry-run writing to S3, resume clobbering the backup,
and resume clobbering it via a different prefix. But two things keep it below A-. First, **not one
byte of a real fragment H5 has been read** (§12.2) — every structural claim is inferred from source,
and the honest expectation is still that stage 1 surfaces something wrong. Second, the two
highest-value new mechanisms in this revision — the float32-accumulator simulation (§5b) and the
two-reference byte-diff (§7.2.1) — are *specified but unwritten and unrun*, so this pass converted
"stated wrongly" into "stated correctly and untested", which is progress but not evidence. A
document whose foundation is entirely unvalidated does not get an A-, and claiming one would be the
kind of overclaim this section exists to prevent.

### 13.1 Changes from review

#### Round 1

| # | Finding | How it was addressed |
|---|---|---|
| C1 | Resume can destroy the only backup | New §6.5: backup keys are write-once — `head_object` before `copy_object`, skip if it matches `sha256_local`, **abort and require human review** if it exists and differs. `status:"started"` ledger record added at §3.1 step 9, before the upload. §8.2's false "nothing is left half-done" claim rewritten to name the in-flight file as the dangerous one. |
| C2 | `--dry-run` writes to S3 | §3.1 reordered: recompute and all gates now precede the backup, which moved to step 6, immediately before the write path. Steps 1–5 are provably write-free. §6.4's barrier restated as `verify backup → upload`, still per-file and unconditional. |
| H1 | §4.2 misread `fragment.py:442` | Corrected: `.sum(axis=1)` runs in **float32** before the `.astype(numpy.float64)`; only the cumsum is float64. Discharged explicitly — the pairwise add is `0+0`, `0+1.0`, or `0.25+0.25`, all float32-exact. |
| H2 | The "multiple of 2**-25" lemma is false | Lemma deleted entirely. Replaced by the measured enumeration `{0, 0.5, 1.0}` and the 2**52 float64 bound. The general `sequence.pyx` table is retained in §4.1 solely to motivate the §5 layer-5 gate. |
| H3 | B/V/H/D presence — resolved by measurement | §4.2 now leads with the histogram and its provenance; states that B/V/H/D would have made the usable prefix **zero**, not merely shorter; and converts the check into a mandatory preflight gate (§5 layer 5) that refuses to run on any non-`{A,C,G,T,N}` sequence byte. Moved out of §12. |
| H4 | §7.2 proves prefix identity, not reference identity | R1 row downgraded from "Decisive" to "strong corroboration; does not cover the repaired region". §13 no longer calls it decisive. New §7.2.1 states the limit explicitly and adds the per-contig N-block coverage check across the 18 disagreeing contigs, blocking at §7.4 stage 3. §8.5 and §12.1 updated to match. **Superseded in round 2 by H-A** — the N-block coverage check was unsound and has been replaced by a two-reference byte-diff. |
| H5 | Sampled backup verification | Sampling fallback deleted from §6.3. If the checksum path is unavailable, re-download and hash **all** backups. There is no sampled mode. |
| H6 | Upload has the multipart-ETag problem | §6.4 now mandates single-part `put_object` with `ChecksumAlgorithm="SHA256"`, forbids `upload_file` by name and threshold, requires asserting no `-N` suffix on the returned checksum, and names `ChecksumType=FULL_OBJECT` as the only acceptable multipart form. §3.3 updated. |
| H7 | Vectorization / rounding fidelity | New §3.5: `numpy.round(a,5)` ≠ CPython `round(x,5)`; decision is to reproduce CPython semantics, preferring the elementwise path unless benchmarking forbids it; blocking test asserting elementwise agreement over >=10^7 real fragments plus engineered near-`k+0.5` values. Added to §7.4 stage 0. |
| M1 | Two "unverified risks" are verifiable | R3 discharged with `fragment.py:244-250`, `fragments_h5.py:817`, `:537`. The three-Stage-A-sites question discharged with `fragment.py:496`, `:606`, `:739-743`. Both moved into the §4.3 table (R3, new R7) and deleted from §12 and §13. |
| M2 | Three-region structure never modeled | New table in §2.1 defining `<2**23` / `[2**23,2**24)` / `>=2**24`. §5b replaces the "anything before saturation is a hard error" rule with a per-region rule permitting middle-band changes only on N-containing spans. §7.3 and §7.4 updated; a note explains why `tests/test_gc_cumsum_overflow.py:37`'s `2**24` is correct for that test and should be documented as such. |
| M3 | §12.7's causal claim backwards; §4.2 example wrong | The entire 2**-25 / 2**28 / 7.8% / 2.6× line of reasoning is deleted. No occurrence of those figures remains. |
| M4 | `_repair_history` dropped by real downstream paths | New §9.1.1 citing `fragment_selection_and_build_provenance.md:720-739` and `liftover_porcine_fragments.py:193-197`; notes that derived H5s already built from the 218 are unrepaired and unmarked. Added as §12.9. |
| M5 | TOCTOU between download and backup | §6.3 step 1 records the source `ETag`; step 2 passes `CopySourceIfMatch`, and a mismatch aborts the file with `PreconditionFailed`. Marked mandatory. |
| M6 | `gc is None` → 255 path unspecified | §3.5 states the rule: `length == 0` → force 255, exclude from the divide, exclude from the diff histogram; notes `lengths` is `uint16` so zero must be checked, not assumed absent. §7.3 and §5b carve it out of the `255 → x` hard error. |
| M7 | Resume falls back to the rejected ETag comparison | Fallback deleted from §8.3. Verification is `head_object` with `ChecksumMode="ENABLED"`, or re-download and hash. Noted that `Range: bytes=0-` is a full GET, not free. Resume classification is now a three-row table including the `"started"`-without-`"ok"` crash window. |
| L1 | Contradictory scratch figures | Standardised on "549 GB total / 414 GB free" in §3.1 and §3.4. |
| L2 | No runtime estimate | §3.4 adds an estimate (~5 min cache build; 2–5 min/file; ~2–5 h for the full apply at `--num-processes 4`), labelled an estimate, with the measured 5-min/2-min FASTA-pass datum as its stated basis. §7.4 stage 2 records the first real number. |
| L3 | Citation drift | `fragment.py:435` → `:441` for the `'a'` pad; `fragments_h5.py:566-575` → guard at `:562-564`, remap at `:570-574`. |
| L4 | `_contig_lengths_str` claim overbroad | §5 now scopes it to BAM input and notes `fragments_h5.py:1024-1037` takes contig lengths from the FASTA for TSV/BED. |
| L5 | Orphaned multipart parts | Moot under H6's single-part `put_object`; §3.3 says so explicitly. |
| L6 | Version drift | Header notes `pyproject.toml:7` reads `2.12.0` in this worktree and must be set to `2.13.0` by the implementing PR; `tool_version` and `_repair_history.version` now read `importlib.metadata.version("fragments-h5")` (§8.3, §9.1). |
| scope | Cut `gc_audit`, `_gc_source`, `--restore`, prototype commits, `--yes` | All five adopted: §11 refuses a `gc_audit` module in `src/`; §9.2 refuses `_gc_source`; §8.4 replaces `--restore` with the runbook; §11 commits the TSV **and the measured FASTA histogram** but leaves `/tmp/gcpilot/*` and `fasta_fingerprint.py` as scratch; §10.3 has one prompt and no bypass. |
| also | §7.3 blocking vs advisory | §7.3 split into Blocking (reader opens, dtype/shape, `255↔x`, per-region change rule, non-`gc` byte-identity, oracle still passes, idempotency) and Advisory (change-fraction plausibility, GC unimodal near the measured **0.4102**, file size within 5%). |

#### Round 2

Round 2 diagnosed the round-1 failure mode as *sections fixed locally without a global re-read* —
several findings below were internal contradictions where one section had been corrected and its
restatements elsewhere had not. Each row therefore names the restatement sites it also had to fix.

| # | Finding | How it was addressed |
|---|---|---|
| C-A | Write-once is not bound to a stable backup prefix or ledger — the C1 hole survives in a narrower form | New **§6.6**: write-once is a property of the triple *(source key, backup prefix, ledger)* held constant across resumes. Ledger gains a **header record** pinning `backup_prefix`; a `--backup-prefix` disagreeing with it **refuses the whole run**. Prior records are cross-checked for a differing `backup_uri` root. Populated backup prefix + empty ledger ⇒ refuse to start. All checks run once at startup, never per-file. §6.1 pins the root as a **literal, never date-computed**, and says why. §8.2's overclaim corrected to state the conditional guarantee. Restatements also fixed: §6.5 closing paragraph, §8.5, §10.2 (`--backup-prefix`, new `--ledger`), §11 generality bullet, §13. |
| H-A | §7.2.1's N-block coverage check tested the wrong predicate, could pass on zero evidence, and silently narrowed R1 | Coverage heuristic **deleted**, not caveated, with an explicit note on why it was unsound (tested N-ness in one FASTA rather than difference between two; no minimum fragment count, and reads never align in hard-masked blocks). Replaced by the **primary** check: recompute the written region under **both** references and report the differing-byte count — zero ⇒ question closed for that file, nonzero ⇒ exact blast radius before any write. Supporting per-contig stats (a)/(b)/(c) retained as localizers. R1 split into **(a) two-reference confusion** (now decidable) and **(b) arbitrary REF-P12 drift** (permanently open, corroborated only by the prefix oracle). Restatements also fixed: §7.4 stage 3 gate, §8.5, §12.1, §13's second weakest point. |
| H-B | `_repair_history` falsified five "byte-identical" / "no-op" assertions | `_repair_history` is now carved out **by name** everywhere. §3.1 step 8, §5b invariant, §7.3 blocking all read "every attr except `_repair_history`". Idempotency restated precisely: **the `gc` datasets are idempotent; the file sha256 is not.** §8.3's `"ok"` resume row now compares against the **most recent** record and records `repair_history_len`, so a legitimate re-run cannot brick the run with a hash mismatch it caused itself. §8.3's belt-and-braces sentence corrected. §5b gained an explicit scope note. |
| H-C | §3.1 step 8 was unimplementable — no operand to compare against | Per-dataset raw-chunk hashes are now taken at **step 2, before any mutation**; step 8 compares against those. §12.8 reframed: the operand question is settled, only the cost is open (plus contiguous-dataset behaviour). §3.1's confidence line repointed from step 8 to step 2. |
| H-D | §5 layer 3's override flag contradicted §5b, §10.3 and §11 | Override **deleted**. §5 layer 3 now routes `255 → x` to §5b's hard error. §11 gains an explicit "Refused" bullet naming the deleted override, consistent with the no-bypass stance. |
| M-A | Middle-band mechanism described incorrectly (the rule survives, the rationale was wrong) | §2.1 corrected: in `[2**23, 2**24)` ulp = 1.0, the representable set is the integers, and `x + 0.5` is an **exact tie** — so inside a contiguous `N` run every `+0.5` is **discarded outright** (systematic near-total loss, not slow drift), while isolated `N`s alternate **−0.5 at even `x`, +0.5 at odd `x`**. §5b's "accumulates slowly" and "changes are small" both replaced: the *count* is small, the *magnitude* can be tens of uint8 units where a fragment straddles an N-block edge. New **§4.2 Consequence 5** supplies the missing proof — the cumsum error is **piecewise constant, changing only at `N` positions**, so `err(stop) − err(start) = 0` on any N-free span, which is exactly what licenses §5b's predicate. §12.7 and §13 updated. |
| M-B | Region classification used the true cumsum; the original build saturated on the drifted one | §5b now classifies by a **simulated float32 accumulator** (one extra vectorized pass over the same cached data), not by the float64 cumsum, with the reason stated: the drift is systematically negative and megabase-scale across N blocks, so float64 classification would push already-saturated fragments into the restrictive `[T23, T24)` band and abort files spuriously. Band-widening named and rejected as the inferior alternative. §2.1 and §12.7 updated. |
| M-C | §6.5's table had no row for "backup exists, no `ChecksumSHA256`" — the state §6.3's own fallback creates | Row added: **re-download and hash**, then apply the matching/differing row. Without it the fallback path was unresumable. `head_object` now explicitly specifies **`ChecksumMode="ENABLED"`**, without which the field is simply absent from every response. |
| M-D | No contig-length or coordinate-range preflight | §5 layer 3 now checks, per present contig, `fasta.get_reference_length(c) == h5.contig_lengths[c]` and `max(starts + lengths) <= length`, citing the indexing at `fragment.py:441-442` and naming the failure (`IndexError` at best, negative-index wraparound at worst). |
| M-E | Repaired files are materially larger and the doc never said so | §3.2 states both baselines: **~0.25% vs a from-scratch build, ~2.5× vs the corrupted original** on the `gc` dataset, because zero-filled `gc` gzips to nothing. §3.4's transfer model corrected (upload leg ~1.2–1.5× the download leg). §6.1 adds a revised S3 footprint (**~190–200 GB** additional). §7.3's advisory now names the from-scratch baseline explicitly and requires reporting absolute growth too. |
| M-F | When `_repair_history` is written was unspecified, though `sha256_repaired` depends on it | §3.1 step 7 now writes it, in the same `r+` handle, before close; step 8 hashes the closed file. §9.1 states the ordering and its two consequences. |
| L-a | "24 GB against 96 GB available" vs a 124 GB environment | §3.4 now reads "the **96 GB available** on this 16-CPU host, of 124 GB total". |
| L-b | §9.1.1's cross-references were wrong | Derived-H5 item points at **§12.9**; the allowlist fix is stated as tracked in **neither** §12 entry and belongs as an issue against the architecture doc; §12.10 is identified as the published-results item, a third and separate thing. |
| L-c | Saturation offsets derived from GC-excluding-N, but the accumulator also gains 0.5 per `N` | Accumulator density derived in §4.2: `(632,693,920 + 635,330,310 + 80,665,851.5) / 3,252,208,893 = 0.41470`. Saturation **~40.45 Mb** (was 40.9), half-integer bound **~20.23 Mb** (was 20.5), corrected at all seven sites (§2.1 ×3, §4.2 Consequence 4 ×2, §7.2, §7.2.1 ×2, §7.4, §13). Both are now flagged as **genome-average estimates** varying per contig, with the note that the tool classifies by cumsum value and never by coordinate. |
| L-d | §7.2's "byte-identity below 2**23" and §5b's `< T23` row are the same check under two names | Unified: §5b's region table is the single implementation site; §7.2 explains the bound and names the gate; §3.1 step 4 says so explicitly instead of listing both. |
| L-e | §4.1 quoted Stage B without its `ValueError` range guard | Guard at `fragments_h5.py:842-843` quoted inline, with a note that an out-of-range recomputed `q5` is an internal error rather than something to clamp. |
| L-f | The pre-fix float32-cumsum line was never quoted | §2.1 now quotes `git show 95c76f5^:src/fragments_h5/fragment.py` line 423 (`…sum(axis=1).cumsum()`) beside the fixed HEAD line at `fragment.py:442`. The one load-bearing claim with no quoted evidence now has some. |
| L-g | `--dry-run` never exercised §6.5, so its halt condition was first seen under `--apply` | §10.1 and §3.1 step 5: `--dry-run` performs the read-only `head_object` probe of each backup key and the §6.6 startup cross-checks, and reports what it found. |
| L-h | Duplicate keys in `--target-list` were not forbidden | §10.2 rejects duplicates **at load time, in both modes**, with the failure named: two `_repair_history` timestamps and a ledger record whose `sha256_repaired` §8.3 escalates to "stop the whole run". Cross-referenced from §8.3. |
