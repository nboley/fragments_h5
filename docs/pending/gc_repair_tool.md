# GC Repair Tool

> **Status: PROPOSAL, awaiting approval. Nothing here has been implemented and nothing has been
> run against production data.** Target release: **v2.13.0**. Execution against the 218 production
> files requires a *separate* explicit approval beyond approving this design — see §10.
>
> Version bookkeeping: `main` is tagged v2.12.1 but `pyproject.toml:7` in this worktree reads
> `2.12.0`. The implementing PR must reconcile that and set it to `2.13.0`; the tool reads its own
> version with `importlib.metadata.version("fragments-h5")` and never hardcodes it.

---

## 0. How heavy should this document be? (round-3 decision)

Rounds 1 and 2 calibrated this document for an **irreversible** operation on the only surviving
copies of ~91 GB of scientific data. **That premise is now false**, and the document has been cut
accordingly. The decision and its justification, stated up front because it explains ~250 deleted
lines:

**Fact 1 — the originals are backed up and independently verified.** All 218 have been copied to
`s3://fragmentomics.kariusdx.com/nboley/gc_repair_backup_2026-08-25/`, 218/218 verified by
`ChecksumCRC64NVME` with `ChecksumType=FULL_OBJECT`, and the prefix re-listed afterwards. Ledger:
`~/gc_repair_artifacts/gc_repair_backup_ledger.tsv`. Restore is a mechanical prefix swap
(§8.3). An overwrite is no longer terminal; it is an afternoon.

**Fact 2 — nothing currently reads `gc`.** `return_gc=False` is the default and no code in the
consuming project passes `True`. A wrong repair therefore does not propagate into any result today.

**The simplification these license is asymmetric, and the asymmetry is the point.**

*Safety scaffolding is now dead and has been deleted, not caveated.* Backup orchestration inside the
tool, backup verification, the per-file "verify backup before overwrite" barrier, write-once backup
keys, the pinned `(source key, backup prefix, ledger)` triple, the two-phase `started`/`ok` ledger
protocol, and the startup refusal checks that enforced the triple — all gone. Three separate
Criticals across two review rounds (C1, C-A, and M-C) were spent defending a resumed run against
clobbering the one pristine copy. **The tool no longer has a backup step at all**, so it cannot
address, overwrite, or mis-resolve a backup key. A mechanism that does not exist needs no invariant.
What replaces the entire crash-resume protocol is one property the document already establishes for
another reason: **the repair is idempotent in the datasets** (§5b), so a resumed run may simply
re-download whatever is at the key and repair it again. Idempotency is cheaper than a protocol.

*Correctness machinery is retained in full, and one part of it is strengthened.* The reference-safety
gates (§5), the region invariants (§5b), the two-reference byte-diff (§7.2.1), and the validation
plan (§7) all stay. The justification is that **Fact 2 cuts both ways.** It lowers the urgency of a
bad repair — nothing consumes the result today — and by exactly the same mechanism it removes every
chance of *noticing* one. No test reads `gc`; no downstream run reads `gc`; the repair writes
plausible-looking bytes into a dataset nobody opens. A backup only helps someone who knows they need
it, and detectability is precisely what a backup does not buy. So the checks that make a wrong
repair *loud* are the load-bearing part of this design and are now the *only* load-bearing part.

**Where I disagree with the framing I was given.** Two places.

1. The ledger was proposed for deletion along with the rest. It should not be deleted, only
   collapsed. Its rollback role is dead (the backup ledger already exists and the restore is a prefix
   swap), but its *record-of-what-was-written* role survives and costs one appended JSON line per
   file. It drops from a header record plus two records per file plus three startup cross-checks
   to **one record per file, written after upload, with no startup validation at all** (§8.2).
2. Scale gating is now expressed **once**. §10 previously had both an unbypassable typed-count
   prompt *and* a `--max-files` cap defaulting to 1 — two ceremonies asserting the same fact ("you
   must deliberately opt into scale"), which is the same two-sources-of-truth error §9.2 rejects for
   `_gc_source`. The prompt is deleted; `--max-files` stays.

**Where the document got heavier, deliberately.** The scope grew (§2.2, §3.2): the tool now also
truncates a phantom padding row from **every dataset of every contig** and rebuilds the index and
`fragment_length_counts`. That is a larger blast radius than the original one-dataset repair, and it
brought its own Critical (§2.2.3).

**Net effect on length: roughly flat, and slightly longer.** ~250 lines of safety scaffolding were
deleted; ~370 lines of truncation design, its hazards and its validation gap were added. This
document is not shorter than the one that was graded B+, and it would be dishonest to present the
deletions as a reduction. What changed is *where the weight sits*: it moved off a solved problem
(losing the originals) and onto an unsolved one (proving the repair is right without an oracle that
can reach it). Net effect on risk is roughly flat, which is why the grade in §13 does not move.

---

## 1. Goal

Make the `gc` dataset in 218 corrupted fragment H5 files equal to what a correct build would have
written, remove the phantom padding row those files carry in every dataset, and be able to
demonstrate both.

Secondary goal, stated by the user and binding: the tool must be reusable for future repair/backfill
work on fragment H5s, *without* being over-engineered. §11 resolves that tension concretely, and §0
records how the resolution moved once the backup existed.

## 2. Problem statement

### 2.1 The GC bug

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
minority but they are the ones a `gc == 0` filter misses. The tool recomputes *everything* (§3.1).

### 2.2 The phantom padding row — a second, independent defect

This is new in round 3. It was found by opening a real fragment H5 for the first time, and it is the
single strongest argument in this document for why §7's validation plan cannot be replaced by
reasoning from source.

#### 2.2.1 Mechanism

`mk_dataset` (`fragments_h5.py:882-886`) truncates each per-contig array to the fragment count
before writing it:

```python
def mk_dataset(key, data, dtype):
    f.create_dataset(
        key, data=data[:ff], dtype=dtype,
        compression="gzip", compression_opts=4, chunks=True
    )
```

The buggy form was `data[: ff + 1]` while `ff` was **already** the fragment count, so every dataset
was written one element too long, the extra slot carrying the zero fill from the preallocated
buffer. Introduced `caddb89` (2025-11-19), fixed `778f4d1` (2025-12-17). **No tagged release carries
it** — the earliest tag, v2.2.1, is 2025-12-18. The 218 were built 2025-11-24 from unreleased
`main`, so **all 218 have it**.

#### 2.2.2 It affects every dataset, not just `gc`

`mk_dataset` is one helper called for `starts`, `lengths`, `mapq`, `gc`, `strand`, the methylation
arrays, and `fragment_end_clipped` (`fragments_h5.py:888-901`). Measured on a real file's `chrM`:
every one of those is length **287,534** with a zero/empty final element — `starts=0`, `lengths=0`,
`gc=0`, `mapq=[0,0]`, `strand=b''`.

Note `mapq` is **2-D** (per-mate). Truncation is along **axis 0 only**; nothing may assume 1-D.

#### 2.2.3 Three live consequences, one of which is a query-correctness bug today

1. **`starts` is not sorted.** chrM's tail reads `[..., 16535, 0]`. `read_fragments`
   (`fragments_h5.py:509-510`) runs `numpy.searchsorted` over a sub-array selected by the index; a
   query landing in the final index block scans a range that includes the phantom, and binary search
   over an unsorted array has no defined answer. **This is a live bug in these files right now,
   independent of GC.**
2. **The index counts the phantom.** chrM's index is `[0, 78285, 146135, 248324, 287534]` with
   `index_block_size = 5000` (`fragments_h5.py:172`), and the appended end sentinel
   (`fragments_h5.py:1240`, `numpy.append(index_poss, len(starts))`) equals `n_frags` *including*
   the phantom. Worse, the interior entries were themselves produced by `numpy.searchsorted` over
   the unsorted `starts` (`fragments_h5.py:1235-1237`), so **the trailing entries may be wrong too**
   — the index cannot be repaired by decrementing the sentinel and must be rebuilt from scratch.
3. **`fragment_length_counts` is inflated.** `_add_fragment_length_counts`
   (`fragments_h5.py:723-735`) histograms `lengths` over every contig, so each phantom adds one
   count at bin 0. `FragmentsH5.n_fragments` is `fragment_length_counts.sum()`
   (`fragments_h5.py:321-323`), so **`n_frags` is over-reported by one per contig group** —
   ~189 per file across the 218. This is a second live bug, also independent of GC.

#### 2.2.4 Disposition (user decision, binding)

**Drop the trailing zero-length row from every contig, rebuild the index, and rebuild
`fragment_length_counts`.** The user has explicitly accepted that per-contig fragment counts will
change. Mechanism in §3.2; ordering relative to the GC recompute in §3.3, which is load-bearing.

### 2.3 What is and is not affected

**By the GC bug:** only `gc`. `starts`, `lengths`, `mapq`, `strand`, methylation,
`fragment_end_clipped`, `index/<contig>`, and all file-level attrs are derived from the BAM and
never touch `get_g_or_c_cumsum`. The GC value is computed in `fragment.py:492-496` / `:603-606` /
`:736-743` and consumed only by the quantizer at `fragments_h5.py:841-846`. Confidence: **high** —
`95c76f5` changed only the accumulator dtype and added region support.

**By the padding bug:** every per-contig dataset, `index/<contig>`, and `fragment_length_counts`
(§2.2.2, §2.2.3). File-level attrs are unaffected — and preserving that is what drives the §3.2
mechanism choice.

### 2.4 Blast radius

218 `.h5` objects under `s3://fragmentomics.kariusdx.com/nboley/ibd_v2/build_frag_h5s/`
(prefix total measured: 437 objects / 91.13 GB — 218 `.h5`, 218 zero-byte folder markers, 1 `.gz`),
all from a single build run on 2025-11-24, all fingerprinted to REF-P12 (218/218). They are
consumed by the CD (281) and UC (315) IBD cohorts; 123 CD + 95 UC of the 218 are consumed, and all
218 fall on the *train* side of the split.

**No clean substitute exists.** Every one of the 218 has exactly one alternate copy elsewhere and
every one of those predates the GC feature entirely (no `gc` dataset). Recomputation is the only
option.

**Nothing reads `gc` today.** `return_gc=False` is the default and no code in the consuming project
passes `True`. This bounds the urgency and, more importantly, bounds the *detectability* of a bad
repair — see §0 and §13.

### 2.5 Disposition: overwrite in place, backups already taken

Overwrite the originals in place. The bucket has **no object versioning** — established empirically:
`head_object` on one of the 218 returns no `VersionId`, while the same call against
`karius-biomarker-data-assets` (versioning Enabled) does. `GetBucketVersioning` and
`ListObjectVersions` are both `AccessDenied` for our credentials, so the `head_object` behaviour is
the only evidence — but it is the behaviour that actually matters.

**The backup is already done, out of band, before this tool exists.** 218/218 copied to
`s3://fragmentomics.kariusdx.com/nboley/gc_repair_backup_2026-08-25/` and verified. Two findings
from that run are carried forward because they change what this tool should do:

- **`copy_object(..., IfNoneMatch="*")` gives API-enforced write-once.** Tested: a second copy to
  the same key returns `PreconditionFailed`. This is strictly better than head-then-copy, which
  leaves a TOCTOU window. It is recorded here for future repairs and for §8.3's restore runbook; the
  repair tool itself no longer copies anything.
- **ETag is unusable for verification on these objects.** The sources are 27-part multipart uploads,
  so their ETags are composite `<md5-of-md5s>-27` values; a single-part server-side copy of such a
  source gets a plain-MD5 ETag on the destination. Verifying by ETag produced **218 false
  mismatches**. Use `ChecksumCRC64NVME` via `head_object(ChecksumMode="ENABLED")`. §6.2 adopts
  CRC64NVME for the upload path for the same reason, and because it is the algorithm we have now
  actually exercised end-to-end against these exact objects with these exact credentials.

Because the backup exists and is verified, **this tool has no backup step**. It downloads, repairs,
verifies locally, uploads, and appends one ledger line. Everything that previously guarded the
act of creating a backup is deleted (§0).

---

## 3. Design

### 3.1 What the tool does, end to end

New console script `repair-fragments-h5-gc` (a second `[project.scripts]` entry alongside
`build-fragments-h5 = "fragments_h5.main:main"`; the repo has no `scripts/` dir and no subcommand
structure, so a flat second entry point is the least-surprising shape).

Per target file, strictly in this order. **No step before 7 writes anything to S3**, so `--dry-run`
(the default) provably performs zero writes:

1. **Download** the S3 object to local scratch (`/tmp`, 549 GB total / 414 GB free; largest observed
   file ~740 MB). Compute `sha256_local` of the downloaded bytes.
2. **Preflight / reference safety** (§5) on the local copy. Abort the file on any mismatch. **While
   the file is still pristine, record a `sha256` of the serialized attr set excluding
   `_repair_history`.** Attrs are the only thing this tool promises not to touch (§3.2), so they are
   the only thing whose pre-mutation hash must be captured here.
3. **Detect and truncate the padding row** (§3.2), then **rebuild** `index/<contig>` and
   `fragment_length_counts`. In `--dry-run` this is done on in-memory shapes only and reported, not
   written.
4. **Recompute** all `gc` values, on the **truncated** arrays, for every contig present in both the
   H5 and the FASTA, into in-memory `uint8` arrays (§3.4). The ordering relative to step 3 is
   load-bearing — see §3.3.
5. **Pre-saturation oracle gate** (§7.2) and **diff and gate**: compare recomputed vs stored, per
   contig, per region (§2.1). Report counts. The `< T23` byte-identity requirement of §5b **is** the
   §7.2 oracle — one check, stated once in §5b's region table, not two.
6. If `--dry-run` (the default), **stop here** and emit a report. Nothing local has been modified
   and nothing in S3 has been written.
7. **Write** into the *local* copy, opened `h5py.File(path, "r+")`: the truncated datasets, the
   rebuilt index, the rebuilt `fragment_length_counts`, the new `gc` arrays, **and the
   `_repair_history` element (§9.1), in the same open handle**, before the file is closed. The
   ordering is fixed here because `sha256_repaired` is taken at step 8 and must cover the history
   attr.
8. **Re-verify** the local file: reopen read-only and run the §7.3 blocking checks. Then compute
   `sha256_repaired` over the closed file.
9. **Upload** the local file over the original key with a single-part `put_object` (§6.2).
10. **Verify the upload** (§6.2) and append one ledger record (§8.2).

**Recompute all `gc`, not just the broken values.** Three reasons: (a) spanning fragments are corrupt
without being zero (§2.1); (b) a full recompute makes the clean-file case a *provable `gc` no-op*
(§5b), which is the strongest safety property available here; (c) it removes an entire class of
"which values did we decide to touch" bugs. Cost is bounded and acceptable (§3.6).

**Note what step 2 no longer does.** Rounds 1–2 required a per-dataset raw-chunk hash of every
non-`gc` dataset, taken pre-mutation, to prove byte-identity at step 8. That is now **impossible and
also pointless**: the truncation rewrites every dataset by construction, so there is no dataset whose
bytes are expected to survive. The invariant it enforced is replaced by the shape and content
assertions of §7.3, which are the right checks for a tool that rewrites arrays deliberately. §12's
old entry 8 (can h5py hash raw chunks cheaply?) is therefore **deleted, not resolved** — the question
no longer arises.

### 3.2 Truncation mechanism and index rebuild

#### 3.2.1 Detecting the padding row — fail-closed, and idempotent

`Dataset.resize()` **will not work**: `maxshape == shape` on every dataset in these files, so they
are fixed-size. Shrinking requires delete-and-recreate. That means the detection predicate must be
exact, because a false positive silently deletes a real fragment.

Per contig group, with `n = len(starts)`:

> **Truncate iff both hold:**
> (a) **sortedness violation at the tail**: `n >= 2 and starts[-1] < starts[-2]`;
> (b) **full zero signature**: `starts[-1] == 0`, `lengths[-1] == 0`, and every other present
> per-contig dataset's final element along axis 0 is its zero/empty value (`gc == 0`,
> `mapq == [0, 0]`, `strand == b''`, methyl arrays `0`, `fragment_end_clipped == 0`).
>
> If **exactly one** of (a) and (b) holds, **abort the file** and require human review. Do not
> guess.

Why both, and why this shape:

- (b) alone is not safe. `lengths` is `uint16`, so a genuine zero is representable, and §3.4 already
  has to handle `length == 0` as a legitimate case.
- (a) alone is not safe either, but (a) is the strong half: a correctly built contig's `starts` is
  ascending, so `starts[-1] < starts[-2]` cannot occur in a clean file at all.
- `n >= 2` is always satisfied where the bug is present. `build_sub_fragments_h5` returns early with
  no datasets when `ff == 0` (`fragments_h5.py:873-875`), so any contig group that exists has at
  least one real fragment, and the buggy build therefore wrote at least 2 rows.
- The predicate is **idempotent**: on an already-repaired file the final row is a real fragment,
  neither (a) nor (b) holds, and nothing is truncated. This is what lets §8.1 discard the entire
  crash-resume protocol.

The tool records the per-contig truncation decision in its report and in the ledger. On the 218 the
expected result is "truncated" for **every** contig group in **every** file; a file where some
contigs truncate and others do not is an anomaly and is **reported as blocking**, because the bug is
in a single shared helper and cannot apply selectively.

#### 3.2.2 The mechanism: delete-and-recreate in place

**Decision: delete-and-recreate each dataset in place, in the same `r+` handle**
(`del f[key]` followed by `create_dataset(key, data=arr[:-1], dtype=<same>, compression="gzip",
compression_opts=4, chunks=True)`), reusing `mk_dataset`'s exact parameters so the recreated dataset
matches what a correct build would have produced.

The alternatives and why they lose:

| option | verdict |
|---|---|
| `Dataset.resize()` | **impossible.** `maxshape == shape`; the datasets are not resizable. |
| **delete-and-recreate in place** | **chosen.** Attrs and group structure are never re-serialized. Cost: HDF5 does not return freed space to the filesystem. |
| write a fresh file, copying everything | rejected on attr fidelity — see below. |
| `h5repack` post-pass | retained as a **gated fallback** for size only, not as the primary mechanism. |

**Why in-place, now that we are recreating every dataset anyway.** Rounds 1–2 justified in-place on
the grounds that a full rewrite re-creates every dataset and risks divergence in chunk shape, gzip
level, and attr dtype. Half of that argument is now gone: we *are* recreating every dataset. What
survives is the other half, and it is the decisive half — **file-level and group-level attrs are
never touched.** `_contig_lengths_str` is a `str(dict)` read back with `eval()` at
`fragments_h5.py:305`; a rewrite that turned it from fixed-width `|S` into a variable-length str, or
vice versa, is a silent regression we would have to test for exhaustively. In-place `r+` makes that
class of bug unreachable rather than tested-against. After truncation, attrs are the *only* thing
left that the tool promises not to change, which is exactly why §3.1 step 2 hashes them and §7.3
asserts on them.

**The cost, stated plainly: dead space.** Deleting an HDF5 dataset unlinks it and returns its
allocation to the free-space manager, which may reuse it, but nothing is returned to the filesystem
and the freed extents are not guaranteed to be reused. Per file this is ~189 contigs × 6–10
datasets ≈ **1,200–1,900 delete-and-creates**, plus one B-tree and one object header per recreated
dataset. Bounds:

- **Best case**, and the expected one: each recreated dataset is one element shorter than the
  original, so every chunk but the last is byte-identical in size and lands in the same free-list
  size class. Growth ≈ the metadata leak plus the genuine `gc` growth below.
- **Worst case**: no reuse at all, and the output is roughly **2× a from-scratch build** — the old
  payload stranded alongside the new one.

Separately and inherently, `gc` gets *larger*: a zero-filled `gc` dataset gzips to almost nothing
while real GC bytes do not. Measured on a synthetic chunked+gzip `uint8` dataset (~40% real
beta-distributed GC / ~60% zeroed, rewritten with fully-real data):

| variant | bytes |
|---|---|
| corrupted original | 1,985,045 |
| repaired in place (`r+`) | 4,947,843 |
| built correctly from scratch | 4,935,336 |

Delta vs a from-scratch build: **~0.25%**. Relative to the *corrupted original*: **~2.5×** on the
`gc` dataset. Repaired files will therefore be **materially larger than the objects they replace**,
by roughly the compressed size of one real `gc` dataset each. That is expected, not a defect, and it
must be carried in the upload budget (§3.6), `/tmp` headroom, and the S3 footprint (§6.1). Report
**both** numbers: delta vs from-scratch (the bloat detector) and absolute growth (the budget).

**Evidence level: weak.** That experiment was one write pass on a synthetic *single-dataset* file,
rewriting a dataset **in place at the same shape**. It does not cover 1,200+ delete-and-recreates,
it does not cover HDF5 metadata leak, and it says nothing about free-space reuse across many small
unlinks. **The size question is genuinely open** and is §12.4. The gate: on the first real file
repaired, compare against a from-scratch rebuild expectation; if in-place growth exceeds **5%** of
that, run `h5repack` on the local file before upload. `h5repack` preserves attributes and re-applies
the existing layout by default, so it is the fallback that does *not* resurrect the attr-fidelity
problem a hand-written rewrite would. Its availability on this machine is unverified (§12.4).

Precedent for in-place mutation exists (`build_fragments_h5` reopens its own output `"r+"` at
`fragments_h5.py:1248` to add `fragment_length_counts`; `tests/test_build_provenance.py:109` opens
`h5py.File(..., "r+")`), but **deleting and recreating an existing dataset is a new pattern in this
repo**. Treat it as such: it gets its own test (§7.4 stage 0).

#### 3.2.3 Index rebuild

The index **must be rebuilt from scratch**, not adjusted. §2.2.3 point 2: its interior entries were
produced by `numpy.searchsorted` over an unsorted `starts`, so decrementing the sentinel would
preserve whatever the binary search happened to return near the tail.

Rebuild by **calling the same code path the builder uses**, not by reimplementing it —
`fragments_h5.py:1225-1241`, including both guards:

- skip contigs with `contig_length <= INDEX_BLOCK_SIZE` (`:1228`);
- skip contigs with `len(starts) < MIN_NUM_READS_FOR_INDEX` (`= 100`, `fragments_h5.py:170`, `:1231`);
- `block_indices = range(0, contig_length, INDEX_BLOCK_SIZE)`, `searchsorted(starts, ..., "left")`,
  then `numpy.append(index_poss, len(starts))`.

`INDEX_BLOCK_SIZE` must be read from the file's own `index_block_size` attr (`fragments_h5.py:306`),
**not** from the module constant, which is `5000` at `fragments_h5.py:172` but appears as `10000` at
four other definition sites in the same file (`:31`, `:50`, `:72`, `:95`). A file built with a
different block size must be re-indexed with *its* block size or the index is silently wrong.

**Blocking invariant:** the rebuilt index's *key set* must equal the existing one. A contig that
gains or loses an index entirely means either the guards were evaluated differently at build time or
the truncation crossed the `MIN_NUM_READS_FOR_INDEX` boundary (possible only for a contig with
exactly 100 rows including the phantom). Either way, **abort the file** and require review.

**Blocking invariant:** after truncation `starts` must be non-decreasing on every contig. Assert it
directly — it is one `numpy.diff(...).min() >= 0` per contig and it is the property the index and
`read_fragments` both depend on.

#### 3.2.4 `fragment_length_counts` rebuild

Delete and recreate the root `fragment_length_counts` dataset from the truncated `lengths` arrays.
`_add_fragment_length_counts` (`fragments_h5.py:723-735`) assigns with `self._f[...] = ...`, which
raises if the key exists, so the tool must `del` first — or, preferably, call
`_add_fragment_length_counts` after deleting the key, so there is one implementation of the
histogram and not two.

**Blocking invariant:** `new_counts[0] == old_counts[0] - (number of contigs truncated)` and
`new_counts[k] == old_counts[k]` for every `k >= 1`. This is an exact arithmetic identity — the
phantom rows are the *only* rows removed and every one of them has `length == 0` — so any deviation
means the truncation removed something it should not have. It is the cheapest available end-to-end
check that the truncation did exactly what it claimed, and it costs one array comparison.

Consequence, stated because it is user-visible: **`FragmentsH5.n_fragments` will decrease by the
number of contig groups** (~189) in each repaired file. That is the intended correction of §2.2.3
point 3, and it is the change the user has explicitly accepted.

### 3.3 Ordering: truncate first, then recompute — and why it matters

**Truncation strictly precedes the GC recompute.** This is not a preference; three separate hazards
collapse if and only if the order is this way.

1. **`length == 0` → 0/0 → NaN.** GC is `(cumsum[stop] - cumsum[start]) / length`. On the phantom
   row `length == 0`, so the division is 0/0 → `nan`, and `int(round(nan * 254))` is garbage. This
   would occur on the last row of every contig — roughly **189 contigs × 218 files ≈ 41,000
   occurrences**. Round 1's M6 established the `length == 0` rule without ever noticing that these
   files are *full* of such rows. Truncating first moots all 41,000.
2. **The `255 → x` invariant fires in reverse on scaffolds.** On a contig absent from the FASTA,
   every real row reads `gc = 255` but the phantom reads `gc = 0`. A recompute that ran before
   truncation would move that byte `0 → 255` — an `x → 255` transition, which §5b makes a hard error
   that **aborts the file**. Every one of the 218 would abort, on all ~165 scaffold contigs. Truncate
   first and the transition never exists.
3. **§5 layer 3's all-255 assertion is false on these files as written.** "Contig absent from the
   FASTA → assert its stored `gc` is entirely 255" fails on all 218 today, because of that same
   phantom `0`. §5 layer 3 is therefore restated to apply **after truncation**. This is the kind of
   contradiction a local fix leaves behind; it is fixed at the source and at both restatement sites
   (§5 layer 3, §7.1 pass criterion 2).

**The `length == 0` rule survives truncation anyway**, and is retained deliberately: `lengths` is
`uint16`, a genuine zero is representable, and the tool is meant to be reusable on files this
document has not seen. The rule (§3.4) is unchanged; what changes is that it is now expected to fire
**zero times** on the 218, and a nonzero count there is itself worth reporting.

### 3.4 Recompute: the vectorized path, and the `round(x, 5)` hazard

Step 4 must be vectorized. A Python per-fragment loop over O(10^8) fragments × 218 files is not
viable. But the naive vectorization is **not** guaranteed to reproduce Stage A, and the difference is
observable in the stored byte:

> `numpy.round(a, 5)` and Python's `round(x, 5)` are different functions. NumPy multiplies by 1e5,
> `rint`s, and divides back — a sequence of three inexact float64 operations. CPython's `round(x, 5)`
> uses correctly-rounded decimal conversion. They disagree on some inputs by 1 ulp, and Stage B then
> multiplies by 254, so a 1-ulp disagreement can straddle a `k + 0.5` boundary and change the stored
> `uint8` by one.

Stage B is safe to vectorize: `int(round(y))` (CPython, on a float) and `numpy.rint(y)` both round
half-to-even on the binary value and agree elementwise. Stage A's `round(x, 5)` does not have that
property, so it is the one that must be pinned.

**Decision: reproduce CPython `round(x, 5)` semantics exactly, and prove it by test.** Concretely:

- Compute `q = (cumsum[stops] - cumsum[starts]) / lengths` as a float64 array (numerator exact per
  §4.2), with `lengths == 0` masked out (below).
- Apply a decimal-correct round-to-5-places, `q5`. The implementation is one of:
  (a) an elementwise loop calling CPython's `round` — the semantics are correct by construction;
  benchmark it, and if it fits the §3.6 runtime budget, **use it and stop**;
  (b) a vectorized approximation, permitted **only** if (a) misses the budget, and **only** with the
  agreement test below green.
  Pick (b) on measured need, never on aesthetics.
- `gc_u8 = numpy.rint(q5 * 254).astype(numpy.uint8)`.

**Required test, blocking:** elementwise equality between the chosen path and a plain Python-loop
reference implementing `int(round(round(float(num)/float(den), 5) * 254))`, over **>= 10^7 real
fragments** drawn from an actual fragment H5, *plus* a synthetic set engineered so that `q5 * 254`
lands within 1 ulp of `k + 0.5` for many `k`. "Spot-checked on a few thousand" is not acceptable —
this is the last bit of every stored byte.

**Zero-length fragments (`length == 0`).** Stage A emits `None` when `g_or_c_cumsum is None` **or**
`frag_stop == frag_start` (`fragment.py:492-493`, `:603-604`, `:736-737`), and Stage B maps
`None`/NaN → 255 (`fragments_h5.py:844-846`). The rule:

> `length == 0` → force the recomputed value to **255**, exclude the fragment from the division (no
> divide-by-zero, no NaN entering the quantizer), and exclude it from the diff histogram entirely.

After §3.3's truncation this is expected to fire **zero times on the 218**; it is retained because
`lengths` is `uint16` and because the tool is reusable. A nonzero count is reported, not silent.
Note this is the *only* legitimate source of 255 on a contig that **is** present in the FASTA; every
other 255 means the contig is absent, which §5 layer 3 handles separately.

### 3.5 Where the work happens

Local disk, one file at a time per worker: download → truncate → recompute → verify → upload. Not
streaming, not in-place-over-S3.

This has a load-bearing safety consequence: **an interrupted write can never produce a corrupted
object in S3.** All mutation happens on a throwaway local copy; S3 sees only a complete, verified
`put_object`. The "neither original nor repaired" failure mode is confined to `/tmp`, where the
recovery is `rm` and retry.

**The upload is a single-part `put_object`, never `upload_file`** — see §6.2. `awscli` is not
installed; `boto3` 1.43.74 is. Credentials are the personal IAM user
`arn:aws:iam::573640641260:user/nathan.boley`; `copy_object`, `put_object` and
`head_object(ChecksumMode="ENABLED")` are all confirmed working with them by the backup run (§2.5).
Because nothing here uses multipart, **there is no orphaned-multipart-part liability** and no
lifecycle rule is needed. The tool must **explicitly refuse** any object > 5 GB rather than fail
late; all 218 are far below it.

### 3.6 Parallelism, memory, and runtime

`get_g_or_c_cumsum` materializes the whole contig as a Python `str`, then an `(L, 4)` float32 one-hot
array, before reducing. For chr1 (248,956,422 bp) that is ~4 GB transient for the one-hot plus ~2 GB
for the float64 cumsum: **~6 GB peak per concurrent contig**. There is **no caching** — every call
re-fetches and re-encodes. Naively, 218 files × 25 contigs = 5,450 full-genome encodes.

**Proposal: a one-shot per-contig cumsum cache.** Preflight computes each contig's float64 cumsum
once and writes it to `--cumsum-cache DIR/<sha256-of-fasta>/<contig>.npy`. Workers open them with
`numpy.load(..., mmap_mode="r")`. Cost: ~3.25e9 bases × 8 B ≈ **26 GB on disk** (fine on the 549 GB
`/tmp` overlay, 414 GB free), shared across all workers through the page cache, per-worker resident
memory ≈ 0. Access is near-sequential because fragments are stored sorted by start within a contig.
The cache directory is keyed by the FASTA's sha256, so a different FASTA can never hit a stale cache.

This is ~30 lines and it deletes the entire memory-pressure question, so it clears the "not
over-engineered" bar. If it is rejected, the fallback is: no cache, `--num-processes 4`
(4 × 6 GB = 24 GB peak against the **96 GB available** on this 16-CPU host, of 124 GB total), which
also works but re-encodes the genome 218 times.

Default `--num-processes 4`; one *file* per worker (not one contig per worker) so that the
download/truncate/repair/upload sequence for a file is owned end-to-end by a single process.

**Runtime — an estimate, not a measurement.** Basis: a full gzip-decompress + per-byte histogram
pass over this FASTA took **~5 min** in pure Python and **~2 min** numpy-vectorized (measured, §4.2).
Cache build is therefore ~5 min one-shot. Per file: ~740 MB down + **more than 740 MB up** (§3.2), so
budget the upload leg at ~1.2–1.5× the download leg until stage 2 measures it, plus the truncation
rewrite (~1,200–1,900 delete-and-creates over gzip-compressed data — **the one cost with no
empirical basis at all**, and a plausible new bottleneck) plus the recompute. Call it **3–8
min/file**, wider than round 2's 2–5 because of the rewrite. At `--num-processes 4` and 218 files
that is roughly **3–8 hours wall-clock** for the full apply, plus ~1–2 hours for the §7.4 stage-3
dry-run gate across all 218 (download-bound, no upload). **Every number in this paragraph is an
estimate**; the stage-2 rehearsal produces the first real datum and the runbook should record it.

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
path, and no dependence on array size, so **as written they are not where the risk lives** — the risk
in Stage A is in *reproducing* it, since a vectorized `round(x, 5)` is not the same function (§3.4).
Taking faithful reproduction as given, the remaining question for this section is whether
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

> **Artifact status:** the machine-readable form of this histogram was never committed and `/tmp` has
> since been wiped, so **it must be regenerated** before §5 layer 5 can cache it (~5 min, §3.6). The
> numbers above are the surviving record and they are load-bearing for all of §4.2, §5 layer 5 and
> §7.2. Regenerating is a mandatory stage-0 task (§7.4) and the output is committed under
> `docs/pending/gc_repair/ref_p12_sequence_histogram.json`, keyed by the FASTA's sha256 — at which
> point regeneration also **re-verifies** these numbers rather than merely restoring them.

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

The only B/H/V/D bytes anywhere in the file are inside FASTA *header description text*
(`{'B':27,'D':4,'H':857,'V':33}`), which the encoder never sees. REF-P12 is entirely uppercase: all
of its masking is hard-masking via `N`, and there is no soft-masked sequence at all.

**Consequence 1 — the reachable value set is exactly `{0.0, 0.5, 1.0}`.** `A`/`T` → 0, `N` → 0.5
(from `0.25 + 0.25`), `C`/`G` → 1.0. Every reachable per-base C+G is an exact multiple of 0.5.

**Consequence 2 — the per-base float32 add is exact.** `seq[:, (1,2)].sum(axis=1)` at
`fragment.py:442` runs the pairwise C+G addition in **float32**, *before* the
`.astype(numpy.float64)`; only the `cumsum` is float64. That pairwise add is `0+0`, `0+1.0`, or
`0.25+0.25` — all three exactly representable in float32, all three exactly computed. Obligation
discharged, not glossed.

**Consequence 3 — the float64 cumsum is exact with enormous margin.** A float64 sum of multiples of
0.5 is exact while twice the running total stays below 2**53, i.e. below **2**52 ≈ 4.5e15**. The
largest per-contig G/C cumsum here is ~0.41 × 249 Mb ≈ **1.0e8**. The margin is ~4×10^7. There is no
margin question to argue about.

Given exactness, every partial sum is represented exactly, so `cumsum[stop] - cumsum[start]` is
identical whether the cumsum was built over the whole contig or over any sub-region containing
`[start, stop)`. **Chunk-invariance and process-count-invariance follow.** Stage A and Stage B then
carry it through deterministically (subject to §3.4's `round(x, 5)` care).

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

**Grade: the conclusion holds and the proof now rests on a measurement rather than on an assumption
about the reference.**

### 4.3 Residual assumptions, and how to discharge each

| # | Assumption | Status | How to discharge |
|---|---|---|---|
| R1(a) | The 2025-11-24 build used **REF-P12, not REF-ASSETS** | **Decidable, and the tool decides it.** Both candidates are in hand, so this is not an assumption — it is a computation. | **§7.2.1's two-reference byte-diff**: recompute the written region under both references and count differing bytes. Zero ⇒ reference choice provably does not affect any byte written to that file, and R1(a) is **closed** for it. Nonzero ⇒ exact blast radius reported before any write, human decides at §7.4 stage 3. |
| R1(b) | The REF-P12 object today is **byte-identical** to the REF-P12 object on 2025-11-24 | **Not established, and not establishable.** Fingerprinting is a 5-contig missingness test, not byte equality. The bucket has no versioning and `ListObjectVersions` is `AccessDenied`, so there is no earlier object to diff against. | **The pre-saturation prefix check (§7.2)** is the only evidence: strong corroboration on the real target files, over the region that needs no repair. It is **corroboration, not proof**, and this row stays open permanently (§12.1). The backstop is that the failure is recoverable from the verified backup (§2.5, §8.3). |
| R2 | `one_hot_encode_sequences` emits only float32 values, and the reachable per-base C+G set is `{0, 0.5, 1.0}` | **DISCHARGED for REF-P12 by measurement** (§4.2): sequence bytes are exactly `{A,C,G,T,N}`, zero B/V/H/D, zero K/M/R/Y/S/W, zero lowercase. Table at `sequence.pyx:24-41` verified by reading. | Made a **mandatory preflight gate** (§5 layer 5): the tool rescans the supplied FASTA's sequence-line character inventory and refuses to run if anything outside `{A,C,G,T,N}` (case-insensitive) appears. R2 becomes a per-run fact, not an assumption. |
| R3 | `frag_stop == start + length` | **DISCHARGED from source.** `Fragment.length` → `tlen` → `self.stop - self.start` (`fragment.py:244-250`); the writer stores `lengths_arr[ff] = fragment.length` (`fragments_h5.py:817`); the reader reconstructs `stops = starts + lengths` (`fragments_h5.py:537`). | Nothing further needed; §7.1/§7.2 would also catch an off-by-one loudly. |
| R4 | Stage A/B code paths are unchanged between the 2025-11-24 build and HEAD | **High confidence.** `95c76f5` touched only dtype + region support; the three Stage A sites and the Stage B quantizer are textually unchanged. | Confirmed again by §7.1 passing on a 2026-05-22-built file. |
| R5 | `numpy.cumsum` on float64 is a plain left-to-right sequential scan (no pairwise/SIMD reassociation) | Reassociation would not matter **because the sum is exact** — exact addition is associative. R5 is therefore not load-bearing. | Nothing needed. |
| R6 | The 218 files were built without `--contig-name-map` | Established. | n/a |
| R7 | The three Stage A sites compute the same thing | **DISCHARGED from source.** `fragment.py:496`, `:606`, `:739-743` are arithmetically identical: same expression, same `g_or_c_cumsum is None or frag_stop == frag_start → None` guard, same `gc_offset` subtraction. The only difference is line wrapping. | Nothing further needed. It does not matter which of the three built a given file. |
| R8 | **All 218 carry the padding row, and no clean file does** | **Inferred, not measured across the set.** The mechanism is a shared helper with a known introduce/fix window (`caddb89` 2025-11-19 → `778f4d1` 2025-12-17) and the 218 are a single 2025-11-24 build, so uniformity is near-certain; one real file was opened and confirmed. But 217 have not been checked. | Cheap and mandatory: §3.2.1's detection predicate runs on **every contig of every file** in `--dry-run` and reports the per-contig verdict. §7.4 stage 3 turns 218/218-all-contigs-truncated into a gate. A mixed verdict aborts the file. |

The honest summary: **the arithmetic is settled by measurement** — a measured alphabet, an exact
float32 pairwise add, and a float64 margin of ~4×10^7. **Reference identity is half settled**: R1(a)
is decided per file by §7.2.1's byte-diff, while R1(b), arbitrary drift of the REF-P12 object, has no
possible oracle and remains the real residual risk. **Structure is now partly measured rather than
wholly inferred** — one real file has been opened, which is how R8 exists at all.

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
   the H5. **All `gc` content assertions in this layer are evaluated on the truncated arrays
   (§3.3)** — before truncation the phantom row reads `gc = 0` on scaffolds and would falsify the
   all-255 rule on every one of the 218.
   - contig **present** in the FASTA → **check its geometry, then** recompute its `gc` in full:
     - `fasta.get_reference_length(c) == h5.contig_lengths[c]` — equal lengths, or **abort the
       file**;
     - `max(starts + lengths) <= fasta.get_reference_length(c)` — every stored fragment lies inside
       the contig, or **abort the file**.

       This is not defensive padding. Stage A indexes `g_or_c_cumsum[frag_stop - gc_offset]` into an
       array of length `len(seq) + 1` (`fragment.py:441-442`); a FASTA whose contig is shorter than
       the H5's stored fragments gives an `IndexError` at best and, with negative-index wraparound,
       silently wrong GC at worst. Both checks are O(1) and O(n) respectively over arrays the tool
       already reads.
   - contig **absent** from the FASTA → **assert its post-truncation stored `gc` is entirely 255**;
     if so, leave `gc` completely untouched (the contig is still truncated and re-indexed); if not,
     **abort the file**.

   This rule does the rest of the work. REF-P12 names primaries `chr*` but scaffolds by accession
   while the H5s use UCSC scaffold names, so `chrUn_*` / `chr*_random` legitimately store 255
   throughout — correct, expected behaviour for these files, which the rule reproduces exactly by
   not touching their `gc`. And if someone supplies REF-ASSETS (UCSC-named, 195 contigs), those
   contigs *are* present, so they would be recomputed from 255 to real values — a `255 → x`
   transition, which §5b makes a **hard error that aborts the file**. There is no override flag; see
   §11.
4. **`get_g_or_c_cumsum` returns `(None, 0)` silently when `chrom not in fasta_file`
   (`fragment.py:404`+). The tool must never call it in a way where that silence is possible.**
   Membership is resolved in preflight (layer 3) and the repair path treats a `None` return as an
   internal error, not as "no GC available". Writing 255 over a contig we merely failed to *find* is
   the single easiest way to destroy data here.
5. **Alphabet gate — the tool refuses to run on a FASTA it cannot reason about.** Preflight scans the
   supplied FASTA's **sequence lines** (header lines excluded — REF-P12's headers contain `B`/`D`/
   `H`/`V` in description text that the encoder never sees) and builds a character histogram. If
   **any** byte outside `{A,C,G,T,N}`, case-insensitively, appears, the tool **aborts the entire
   run** with the offending characters and their counts.

   This is not a formality. §4.2's exactness argument — and therefore §7.2's oracle — is valid *only*
   on that alphabet. Had B/V/H/D been present, the usable pre-saturation prefix would have been
   **zero**, not merely shorter: `float32(1/3)` is not on any half-integer grid, so three consecutive
   `H` bases already push the accumulator off the exactly-representable set. K/M/R/Y/S/W would be
   tolerable in principle (they are multiples of 0.5) but they are measured absent, and admitting
   them would mean re-deriving the bounds; the gate keeps one bound for one alphabet.

   The result is cached as `<fasta-sha256>.alphabet.json` alongside the cumsum cache (§3.6), so the
   scan is one pass per FASTA per machine, not one per file. Because the committed histogram was
   lost with `/tmp` (§4.2), the first run of this gate is also what regenerates it.

Confidence: **high** that these layers prevent the wrong-FASTA accident. They do not by themselves
prove the *right* FASTA (that is R1 / §7.2 / §7.2.1).

### 5b. The no-op guard as a safety property

A correct repair of a **clean** file changes nothing but the provenance attr. This is the tool's
central invariant, and it generalises:

> **Invariant.** For any file, `repair(file)` must leave every file-level and group-level attr except
> `_repair_history` byte-identical; must truncate exactly the rows §3.2.1's predicate identifies and
> no others; and the set of `gc` bytes that change must be a subset of what the failure mode
> predicts.

**Scope note, load-bearing, and narrower than in round 2.** Round 2 stated the invariant over
*datasets*. That is no longer available: §3.2's truncation rewrites every dataset in every file that
has a padding row, so "every non-`gc` dataset is byte-identical" is **false by design** on the 218
and has been removed everywhere it appeared. What survives is (i) attr byte-identity, (ii) the exact
arithmetic identity on `fragment_length_counts` (§3.2.4), (iii) the index key-set and sortedness
invariants (§3.2.3), and (iv) the `gc` region rules below. On a **clean** file the truncation
predicate does not fire at all, so on clean files the old, stronger property still holds and §7.1
tests exactly that.

Also unchanged: `--apply` always appends a `_repair_history` element (§9.1), so the **file sha256
always changes**, even on a clean file. `gc` is idempotent; the file hash is not.

Operationally:

- `--dry-run` is the **default**. It performs the full truncation analysis and recompute and reports
  a per-contig diff histogram (`unchanged`, `0 → nonzero`, `nonzero → nonzero` (spanning fragments),
  `255 → x` / `x → 255` (should be **zero** — any occurrence is a hard error)) plus the per-contig
  truncation verdict, without opening the file for writing and without touching S3.
- On a CLEAN file the expected report is *literally all zeros in every change bucket, and zero
  contigs truncated*. A single changed `gc` byte on a known-clean file falsifies the method,
  immediately, with no reference oracle needed.
- On a corrupted file, `255 → x` or `x → 255` transitions must be exactly zero. Round 2 carved out
  `length == 0` fragments as the exception; after §3.3's truncation that exception is expected to be
  **empty** on the 218, and it is retained only for the general case.

**Where changes are allowed: one rule per region (§2.1).** Define the cumsum thresholds
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

The middle-band predicate is cheap: the tool already has the FASTA and `N` is the only non-ACGT code
(§4.2), so a per-contig boolean N-mask plus its prefix-sum answers "does this span contain an N" in
O(1) per fragment. It is also *sound* rather than merely plausible — §4.2 Consequence 5 proves the
float32 cumsum error is piecewise constant, changing only at `N` positions, so an N-free span has
zero error difference and must be byte-identical.

**Do not expect the middle-band changes to be tiny where they occur.** The *count* is small, because
reads essentially never align inside hard-masked blocks. But the *magnitude* is not: where a fragment
straddles an N-block edge the error is up to `0.5 · (N bases in span) / length`, which for a fragment
substantially overlapping an N run approaches 0.5 in GC units — **tens of uint8 units**, not one. A
middle-band diff of 40 is expected behaviour, not a red flag; the flag is a middle-band diff on an
N-free span, of any size.

Because `--apply` is just `--dry-run` plus a write, every production run is preceded by its own dry
run on the same bytes.

---

## 6. S3 mechanics

Round 2 gave this ~175 lines across six subsections, almost all of it about creating and protecting
backups. The backup is done (§2.5) and the tool no longer makes one, so what remains is the upload
and its verification.

### 6.1 Footprint

The 218 originals are 91.13 GB and the verified backup is a second 91.13 GB, already spent. The
repaired originals are larger than what they replace by roughly one compressed real `gc` dataset each
(§3.2.2), and possibly more if the truncation rewrite leaks HDF5 space. Budget **~15–25 GB** of net
S3 growth on the production prefix, plus the 91.13 GB backup, which is released whenever the repair
is accepted (§8.4).

### 6.2 Upload and verification

**The upload is a single-part `put_object`, never `boto3`'s `upload_file`.** `upload_file` switches
to multipart above an 8 MB threshold and every one of the 218 is far above it; a subsequent
`head_object` then returns a **composite** checksum of the form `<base64>-N`, which never equals a
full-object checksum. Single-part `put_object` is available because every file is far below the 5 GB
single-part limit (§3.5). If multipart ever becomes unavoidable, the only acceptable form is
multipart with `ChecksumType=FULL_OBJECT`; that is not needed at these sizes and is not being built.

**Use `ChecksumCRC64NVME`, not SHA-256 and never ETag.** Both parts of that are now grounded in the
backup run rather than proposed:

- **ETag is unusable.** The sources are 27-part multipart uploads, so their ETags are composite
  `<md5-of-md5s>-27` values. Verifying the backup by ETag produced **218 false mismatches** (§2.5).
  The same trap applies one layer up to any upload verification, so there is **no ETag fallback
  anywhere in this document.** `Range: bytes=0-` is a full-object GET, not a free metadata call, so
  it is not a fallback either.
- **CRC64NVME is the algorithm we have actually exercised** end-to-end on these exact objects with
  these exact credentials, and `head_object(ChecksumMode="ENABLED")` is confirmed to return it as a
  `FULL_OBJECT` value. Round 2 specified SHA-256 and listed "does `copy_object` with
  `ChecksumAlgorithm='SHA256'` return a full-object or composite value, and does our IAM allow it?"
  as an open risk. Rather than leave an unverified algorithm in the design when a verified one is in
  hand, switch. The old open item is closed by substitution, not by testing.

Per file:

1. `put_object(..., ChecksumAlgorithm="CRC64NVME")` with the local file's bytes.
2. `head_object(..., ChecksumMode="ENABLED")` on the key; assert `ChecksumCRC64NVME` equals the
   locally computed CRC64NVME of the repaired file, assert `ChecksumType == "FULL_OBJECT"`, assert
   no `-N` suffix, and assert `ContentLength` matches the local file size.
3. Only then append the ledger record (§8.2).

**`copy_object(..., IfNoneMatch="*")` for write-once semantics** is recorded in §2.5 for future
repairs and for the restore runbook. The repair tool does not copy objects and does not need it.

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

**Secondary targets** (same prefix, all fingerprinted REF-P12):
- `.../NC-6724/LB-103927-Lib1/frag_h5s/...` (273,208,801 B, 2026-05-22T17:56:21Z)
- `.../DC4-16772/AC-125674/frag_h5s/AC-125674.fragments.h5` (462,828,323 B, 2026-08-04T22:24:11Z)
- `.../NC-6724/LB-103918-Lib1/frag_h5s/...` (69,303,053 B)
- `.../NC-6724/LB-103913-Lib1/frag_h5s/...` (11,936,498 B — small panel file, probe returned
  INSUFFICIENT_DATA but test_zero_frac = 0.0, directionally clean)

**Pass criteria:**
1. Every recomputed `gc` byte equals the stored byte, on every contig, in every file. Zero tolerance
   — not "99.99%", not "within 1 unit".
2. Every contig absent from REF-P12 (`chrUn_*`, `chr*_random`) is left untouched and was verified
   all-255 by the §5 layer-3 rule. **Scaffold-255 preservation is a pass criterion, not a footnote**
   — a run that "fixes" those has used the wrong reference. On these files the check runs on the
   arrays as stored, because (criterion 4) nothing is truncated.
3. The change histogram (§5b) is all zeros.
4. **Zero contigs truncated**, in every file. These are all post-`778f4d1` builds and have no padding
   row; §3.2.1's predicate must not fire once.

### 7.1.1 The oracle's blind spot, stated rather than assumed

**Clean-file replay cannot validate the truncation path at all, and this is structural.** The
padding bug was fixed on 2025-12-17 (`778f4d1`) and the GC bug on 2026-03-09 (`95c76f5`). A file
that is clean for GC is necessarily built after March 2026 and therefore necessarily has **no
padding row**. Criterion 4 above is a *negative* control — it proves the predicate does not
misfire — and that is all it can be. Every positive assertion in §3.2 is untested by §7.1.

Worse, the usual escape hatch is closed. `output_rebuild_frag_h5s/` holds 489 `.h5` files spanning
2026-01-29 .. 2026-08-04 — 471 pre-GC-fix, 18 post — which round 2 nominated as expendable rehearsal
targets for the apply path. **All of them post-date `778f4d1`, so none has a padding row.** They
remain useful for rehearsing upload/verify/ledger, but they cannot rehearse truncation.

Three responses, in increasing order of what they buy:

1. **A synthetic fixture** built by writing a fragment H5 with a deliberately over-long final row per
   contig (construct the arrays directly with `h5py`, or call the builder with `mk_dataset` patched
   to `data[: ff + 1]`, which reproduces `caddb89` exactly). Cheap, deterministic, belongs in
   `tests/`, and covers every *known* structural difference: detection, truncation, index rebuild,
   `fragment_length_counts` arithmetic, the `n >= 2` edge case, the `MIN_NUM_READS_FOR_INDEX = 100`
   boundary, and a 2-D `mapq`.
2. **A rehearsal on a real affected file** — which now means **a copy of one of the 218**, restored
   from the verified backup into a scratch prefix. This is the only way to exercise the real thing,
   and it is cheap precisely because the backup exists. §7.4 stage 2 does this.
3. **Accepting what neither covers.** A fixture covers differences we thought to model. It cannot
   cover ones we have not discovered — and the padding row itself is the proof of that, since it
   surfaced only when a real file was opened for the first time, after two review rounds that
   reasoned entirely from source. §7.4 stage 3 (dry-run across all 218, reporting per-contig
   truncation verdicts) is the widest net available and is deliberately placed before any write.

### 7.2 Pre-saturation prefix agreement (the only available check on R1(b) — and see its limits)

The corrupted files are corrupt only *past* the saturation point. **Everything before it is
correct**, and it was produced by the exact FASTA bytes used on 2025-11-24. That gives us a per-file
reference oracle on the actual targets:

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
  target files, which the §7.1 clean replay cannot do (its subjects were built 2026-05-22, six months
  on the wrong side);
- simultaneously validates the coordinate convention, the `'a'` pad offset, Stage A rounding and
  Stage B quantization, per file;
- costs nothing extra — the cumsum is already computed.

**It runs on all 218 as a mandatory preflight gate**, and a single mismatched byte in the
pre-saturation region **aborts that file**. It runs on the **truncated** arrays (§3.3); the phantom
row is not a fragment and has no place in an oracle about fragments.

#### 7.2.1 What this check does *not* prove

It proves agreement over the prefix — roughly the first 20.23 Mb of GC-cumsum per contig — which is
precisely **the region that needs no repair**. The bytes the tool actually writes are derived from
sequence past ~40.45 Mb: the untested complement. The two candidate references differ by ~11.9 Mb of
hard-masking (§5), and **that masking could lie entirely past 20.23 Mb**, in which case all 218 files
pass the oracle and the repair still writes 127-over-real-GC across the differing blocks. This is not
hypothetical enough to ignore; the oracle is corroboration, not proof of reference identity.

R1 in §4.3 is actually **two questions**, and only one of them is answerable here:

- **(a) the two-reference confusion** — did the 2025-11-24 build use REF-P12 or REF-ASSETS? Both
  candidates are in hand, so this is decidable, and the check below decides it *outright* rather than
  narrowing it.
- **(b) arbitrary drift of the REF-P12 object itself** — is the object at
  `nboley/resources/GRCh38.p12.genome.fa.gz` today byte-identical to what it was on 2025-11-24? The
  bucket has no versioning, so there is no candidate to diff against and **no check in this document
  closes (b)**. (b) remains open and §12.1 says so.

#### The check that closes (a): recompute under both references and diff

The ~11.9 Mb / 18-contig differential between REF-P12 and REF-ASSETS is a **measured, positionally
known** quantity, and both FASTAs are available. So do not test a proxy — compute the answer:

> For every contig the tool intends to write, recompute the `gc` bytes of the written region under
> **REF-P12** and under **REF-ASSETS**, and report **the number of differing bytes**.

- If the count is **0** for a file, then reference choice **provably does not affect any byte this
  tool writes to that file**, and question (a) is closed for that file. Not narrowed — closed.
- If the count is nonzero, the report names the exact blast radius, per contig, *before anything is
  written*, and the operator decides with a number rather than a hunch.

Cost is one extra float32/float64 cumsum pass over the same cached per-contig data (§3.6) plus one
extra FASTA in the cache — affordable at the scale of §3.6's estimates, and it produces a count, not
a heuristic. This is the **primary form** of the check and is a mandatory preflight report in
`--dry-run`.

**Supporting per-contig statistics**, reported alongside the byte-diff because they localize a
nonzero result:

  (a) differing bases in the pre-saturation prefix that have **nonzero fragment coverage** — i.e. how
      much of the oracle's region actually discriminates between the two references;
  (b) differing bases inside the region the tool actually writes;
  (c) **the number of fragments in the written region that overlap a differing position** — the
      denominator of the exposure.

If (c) == 0 for a file, the byte-diff above is necessarily 0 and the two results corroborate.

**What was here before, and why it is gone.** An earlier draft required that each of the 18
disagreeing contigs have "at least one hard-masked `N` block in its pre-saturation prefix, with the
recomputed `gc` agreeing on fragments overlapping it." That check is unsound in two ways and has been
deleted rather than caveated: it tested `N`-ness in *one* FASTA rather than *difference between the
two*, so a shared telomeric N block at 0–10 kb satisfies it while every differing block sits at
60 Mb; and it imposed no minimum fragment count, so — since reads essentially never align inside
hard-masked blocks — the expected number of overlapping fragments is ~0 and zero comparisons "agree"
vacuously. It could pass on literally no evidence.

Contigs with a nonzero byte-diff are listed by name and are a **blocking** item for the §7.4 stage-3
gate: a human decides whether to proceed, with the exposure quantified.

Confidence: **high** that the prefix check is sound over its region; **high** that the two-reference
byte-diff decides question (a) exactly; **medium** that everything passes cleanly on all 218 first
try; and **explicitly not** a claim that the repaired region is verified against the 2025-11-24
build — question (b) has no oracle.

### 7.3 Post-repair checks on a corrupted file

After writing, for each repaired file. **These are not all the same kind of check** — "the GC
histogram looks unimodal" and "zero `255 → x` transitions" are not both grounds to stop a rollout.

**Blocking — any failure stops the rollout and triggers restore-from-backup (§8.3):**

*GC:*
- **zero** `255 → x` or `x → 255` transitions;
- every change obeys the per-region rule of §5b: none below `2**23`; those in `[2**23, 2**24)` only
  on fragments whose span contains an `N`; unrestricted at or above `2**24`;
- the §7.2 pre-saturation oracle still passes on the repaired file (it must, since the tool does not
  write there — a failure means it wrote somewhere it should not have);
- re-running the tool on the repaired file changes **no `gc` byte** and truncates **no contig**
  (§3.2.1's predicate is idempotent) — the cheapest possible idempotency proof. Note the precise
  claim: **the datasets are idempotent, the file sha256 is not**, because a re-run under `--apply`
  appends a second `_repair_history` element. Run this check in `--dry-run`, which appends nothing.

*Structure (new in round 3):*
- every per-contig dataset is exactly **one element shorter** along axis 0 than before, on every
  contig where §3.2.1 fired, and **unchanged** on every contig where it did not;
- all per-contig datasets within a contig group have **equal length** along axis 0;
- dataset dtypes, chunk shapes, and filter pipelines are unchanged from the originals (`gc` stays
  `uint8`, `mapq` stays 2-D, `strand` stays `|S1`);
- `starts` is non-decreasing on every contig (§3.2.3);
- the rebuilt index's **key set** equals the original's, and every index array's final entry equals
  the new `len(starts)` (§3.2.3);
- `fragment_length_counts` satisfies the exact identity of §3.2.4;
- reopen through `FragmentsH5` succeeds (exercises the reader, `eval()` of `_contig_lengths_str` at
  `fragments_h5.py:305`, the `"gc" not in ...` guard at `:562-564`, and the `/254.0` inversion plus
  255→NaN remap at `:570-574`);
- **a read-back query test**: for a sample of contigs, `read_fragments` over the **final index
  block** returns results consistent with a brute-force scan of the truncated arrays. This is the
  check that the §2.2.3 query-correctness bug is actually fixed, and it is the only one that
  exercises the index end-to-end rather than by invariant.

*Provenance:*
- every attr **except `_repair_history`** hashes to the value recorded at §3.1 step 2;
- `_repair_history` parses as JSON and has gained **exactly one** element.

**Advisory — recorded in the report and reviewed by a human, but do not by themselves halt:**
- `0 → nonzero` fraction consistent with the predicted post-saturation fragment count;
- resulting genome-wide GC distribution is unimodal near the measured **0.4102** (§4.2) with no
  residual spike at 0 — a distribution sanity check, not a correctness proof;
- file size within 5% of a **from-scratch correct build** (§3.2.2) — *not* within 5% of the corrupted
  original, which it will exceed substantially and by design. Report both the from-scratch delta and
  the absolute growth. Advisory because bloat is a storage-cost problem, not a data-integrity
  problem; but exceeding the 5% threshold triggers the `h5repack` decision in §3.2.2 before the next
  stage, so it gates *progression* even though it does not gate *this file*.

### 7.4 Graduated rollout

| Stage | Target | Gate to proceed |
|---|---|---|
| 0 | Unit + regression tests: a delete-and-recreate-dataset test, the **synthetic padding-row fixture** of §7.1.1 with its full check list, the §3.4 rounding-agreement test over >=10^7 fragments, and a §5b float32-accumulator-simulation test asserting the simulated band boundary lands where the pre-fix code saturates. Plus: **regenerate the lost REF-P12 sequence histogram** (§4.2) and commit it. The existing overflow regression test at `8820299` (`tests/test_gc_cumsum_overflow.py`, ~676 MB RSS) must still pass | green CI |
| 1 | `--dry-run` clean replay on the 5 §7.1 targets | §7.1 criteria 1–4, all zeros including **zero contigs truncated** |
| 2 | Full `--apply` rehearsal on **a scratch-prefix copy of one of the 218** (restored from the verified backup), written to the scratch key, not over the original | §7.3 blocking checks all pass, including the read-back query test; file-size advisory reviewed; first real runtime datum recorded (§3.6). **This is the only stage that exercises truncation on real data before production** (§7.1.1) |
| 3 | `--dry-run` across **all 218**: §7.2 pre-saturation gate, §7.2.1 two-reference byte-diff, and the per-contig truncation verdict | 218/218 pass the prefix oracle; **every contig of every file reports "truncate"** (R8) — any mixed verdict is blocking; every file with a nonzero REF-P12/REF-ASSETS byte-diff named, quantified, and accepted by a human |
| 4 | `--apply` on **1** of the 218 | §7.3 all pass; manual review of the diff report |
| 5 | `--apply` on the remaining 217 | — |
| 6 | Post-run: re-probe a sample for `gc == 0`; confirm `n_fragments` dropped by the expected per-file contig count; spot-restore 3 files from backup and confirm they match the backup ledger's recorded checksums | before backups are deleted |

Round 2 had eight stages with a 1 → 10 → 207 ramp. **The 10-file stage is deleted.** Each `--apply`
stage runs the identical gate, so the tenth file tests nothing the first did not, and the failure is
now recoverable from a verified backup in an afternoon. Staging exists to bound the cost of a
mistake; when that cost drops by two orders of magnitude, so should the staging. Stage 3 is kept and
is deliberately placed before *any* production write: it is cheap and it is the check most likely to
surface an unmodelled difference — as §7.1.1 argues, it is the only check with real coverage of the
truncation path across the whole set.

A note on the existing regression test: `OVERFLOW_THRESHOLD = 2**24`
(`tests/test_gc_cumsum_overflow.py:37`) is the right constant *for that test* and needs no change —
its contig is all-`G` with a short `A` patch, so it contains no half-integer bases and the
`[2**23, 2**24)` band of §2.1 simply does not arise. Its docstring should say so, so that a future
reader does not "fix" it to 2**23 or conclude that 2**24 is the universal bound.

---

## 8. Failure, resume and restore

Round 2 gave this five subsections and ~150 lines, of which the majority existed to stop a resumed
run from destroying the only pristine copy. That copy is no longer the only one (§2.5). What remains:

### 8.1 Failure, and why resume is trivial now

**Within one file.** All mutation is on a local scratch copy; the original object is untouched until
a complete, verified file is uploaded (§3.5). An interrupted write leaves a corrupt file in `/tmp`;
recovery is `rm` and retry. Mid-write corruption in S3 is structurally impossible.

**Across the 218.** Each file is independent. A run that dies after 137 files leaves 137 repaired and
81 untouched. The one genuinely awkward state is a crash **between** the upload and the ledger
append: the key holds repaired bytes and advertises nothing.

**That state is now benign, and the reason is idempotency, not protocol.** A resumed run downloads
whatever is at the key and repairs it. If that object is already repaired: §3.2.1's truncation
predicate does not fire (the last row is a real fragment), the `gc` recompute reproduces the same
bytes (§5b), and the only net effect is a second `_repair_history` element recording that it
happened. Nothing is lost and nothing is corrupted. **Rerunning the tool on its own output is
safe by construction, and that single property replaces the entire round-2 apparatus** — the
`started`/`ok` two-phase ledger protocol, write-once backup keys, the pinned
`(source key, backup prefix, ledger)` triple, and the three startup refusal checks that enforced it.
All deleted (§0).

The design obligation this creates is on §3.2.1: the truncation predicate **must** be idempotent and
fail-closed, because it is now the thing standing between a resume and a deleted real fragment. It
is specified to be both, and §7.3's re-run check tests it directly.

### 8.2 Ledger

Append-only JSONL at `--ledger PATH` on **NFS home** (persistent), not `/tmp`. **One record per
file**, appended after upload verification (§6.2 step 3):

```json
{"key": "nboley/ibd_v2/build_frag_h5s/.../X.fragments.h5",
 "status": "ok",
 "sha256_original": "...", "sha256_repaired": "...",
 "crc64nvme_repaired": "...",
 "fasta_sha256": "...", "tool_version": "2.13.0",
 "contigs_repaired": 25, "contigs_skipped_absent_from_fasta": 170,
 "contigs_truncated": 195, "rows_removed": 195,
 "n_fragments_before": 0, "n_fragments_after": 0,
 "gc_bytes_changed": 12345678, "repair_history_len": 1,
 "started_utc": "...", "finished_utc": "..."}
```

`tool_version` is read from `importlib.metadata.version("fragments-h5")` at runtime, never hardcoded
— a hardcoded string is provenance that can silently disagree with the code that wrote it.

**Resume behaviour: skip any key with an `ok` record whose `crc64nvme_repaired` matches the live
object; otherwise process it.** No mismatch escalates to "stop the whole run" any more, because
reprocessing is safe (§8.1). A mismatch is *logged* — it means either a legitimate second repair or
something else wrote the key, and a human should look — but it does not block, and no header record,
prefix cross-check or populated-prefix guard exists to be misconfigured.

Verification uses `head_object(..., ChecksumMode="ENABLED")`. **There is no ETag fallback** (§6.2).

### 8.3 Restore

**There is no `--restore` flag.** Restore is a documented runbook procedure, not a feature: a restore
path that depends on the tool that caused the problem is not a restore path. Once written down it is
ten lines, and shipping a second implementation inside the tool adds a code path exercised only in an
emergency.

The procedure, per key: `copy_object` from
`s3://fragmentomics.kariusdx.com/nboley/gc_repair_backup_2026-08-25/<original key>` back over
`<original key>`, then `head_object(ChecksumMode="ENABLED")` and compare `ChecksumCRC64NVME` against
the value recorded in `~/gc_repair_artifacts/gc_repair_backup_ledger.tsv`. Because backup keys mirror
original keys under one root, this is a mechanical prefix swap, and the backup ledger supplies the
root, the key list and the expected checksums. Do **not** verify by ETag (§2.5).

§7.4 stage 2 exercises the restore leg (it is how the rehearsal target is obtained) and stage 6
spot-checks it on 3 files.

### 8.4 What is and is not recoverable

**Recoverable:** everything the tool can get wrong. Wrong FASTA, wrong GC, an over-eager truncation,
a mis-rebuilt index — all of it is one prefix swap away from the pristine 2025-11-24 bytes. §7.2 is
*weaker* than round 1 assumed (it constrains the reference only over the unrepaired prefix, §7.2.1),
so "the FASTA is wrong in a way §7.2 does not catch" is a live possibility — and it is now a
recoverable one.

**Not recoverable:** losing the backup. The backup prefix is in the same unversioned bucket, and it
is the only rollback. It must not be deleted until §7.4 stage 6 completes and the cohort owners have
re-run at least one downstream analysis. Deletion is a separate, explicitly-approved action. **The
tool never writes to that prefix under any code path**, which is the strongest available guarantee
and is a consequence of deleting the backup step rather than of guarding it.

**Not recoverable by *check*:** arbitrary drift of the REF-P12 object itself (R1(b), §12.1). There is
no oracle for it. The backstop is the backup.

---

## 9. Provenance

The 218 files carry **no build provenance at all** — `_build_version` / `_build_argv` landed in
v2.12.0 (`a3be97b`, 2026-08-21) and these files predate it. The current authority,
`docs/architecture/fragment_selection_and_build_provenance.md`, records those as flat root attrs and
**defines no repair-provenance field**. (Note: `docs/pending/build_provenance_metadata.md` is marked
REJECTED/superseded and should not be used as a reference.)

### 9.1 Proposed attr

A single flat root attr, JSON-encoded, consistent with the existing `_build_argv`-is-JSON convention:

```
_repair_history : str   # JSON array of objects, appended to on each repair
```

Each element:

```json
{"tool": "repair-fragments-h5-gc",
 "version": "2.13.0",
 "argv": ["repair-fragments-h5-gc", "--fasta", "...", "--apply", ...],
 "timestamp_utc": "2026-08-25T17:04:11Z",
 "reason": "gc-float32-cumsum-saturation (fixed by 95c76f5) + trailing padding row (fixed by 778f4d1); see docs/pending/gc_repair_tool.md",
 "datasets": ["data/*/gc", "data/*/* (truncated by 1 row)", "index/*", "fragment_length_counts"],
 "rows_removed_per_contig": 1,
 "fasta_uri": "s3://fragmentomics.kariusdx.com/nboley/resources/GRCh38.p12.genome.fa.gz",
 "fasta_sha256": "...",
 "source_sha256": "...",
 "backup_uri": "s3://fragmentomics.kariusdx.com/nboley/gc_repair_backup_2026-08-25/..."}
```

`version` is `importlib.metadata.version("fragments-h5")`, read at runtime. `backup_uri` names the
out-of-band backup (§2.5) — the tool did not create it, but recording where the pre-repair bytes live
is the single most useful thing this attr can carry.

`datasets` now records the truncation as well as the `gc` rewrite, because a reader who sees only
`data/*/gc` would reasonably assume the fragment counts were preserved. They were not (§3.2.4), and
that is exactly the kind of change provenance exists to announce.

A *list* rather than scalar attrs, because "this file has been repaired more than once" is a real
future state and rewriting scalars would erase history. This is the main piece of forward-looking
generality being built in, and it costs ~5 lines.

**When it is written: §3.1 step 7**, in the same `r+` handle as the datasets, before the file is
closed. The ordering is not incidental — `sha256_repaired` is taken at step 8 over the closed file,
so it must already cover the history attr. Two consequences follow and are stated wherever they bite:
`_repair_history` is **excluded by name** from every "attr unchanged" assertion (§3.1 step 8, §5b,
§7.3), and the file sha256 is **not** idempotent across repeated `--apply` runs even though the
datasets are (§8.1, §8.2).

### 9.1.1 Known limitation: this attr does not survive derivation

`_repair_history` will be **dropped by real downstream paths**, and we know this in advance rather
than discovering it later. `docs/architecture/fragment_selection_and_build_provenance.md:720-739`
documents that two consumers rebuild fragment H5s from a hardcoded 5-attribute allowlist and
therefore "will silently drop `_build_argv` and `_build_version`" — one of them,
`liftover_porcine_fragments.py:193-197`, is a real analysis path, not a test fixture. A sixth attr is
dropped by exactly the same mechanism.

Two consequences, neither quietly absorbed here:

1. **Derived H5s already produced from the 218 are not repaired by this tool and will not advertise
   that they need to be.** They carry corrupted GC *and* the phantom row propagated through
   whatever derivation produced them, and they are indistinguishable from a clean derivation.
   Finding them is a separate task with a separate inventory; the tool's target list covers the 218
   originals only. Tracked as **§12.7**.
2. Any derivation performed *after* the repair loses the record that the repair happened. Widening
   the allowlist is the fix and it is **out of scope for this PR and not tracked in §12 at all** — it
   belongs as a filed issue against `docs/architecture/fragment_selection_and_build_provenance.md`,
   which owns the allowlist.

We are not widening the attr surface to work around this. The honest position is that the provenance
is correct where it is written and lossy where someone else copies it.

### 9.2 Do not backfill `_build_version` / `_build_argv`

**No.** We do not know the original build's argv or version. Writing plausible-looking values would
be **fabricating provenance** — worse than absent provenance, because absent provenance is honestly
absent whereas fabricated provenance is indistinguishable from real. Downstream code that reasons
"`_build_version >= 2.7.2` therefore GC is trustworthy" would then be reasoning from a value we
invented.

The correct signal is the opposite one: `_repair_history` present ⇒ this file's `gc` was produced by
v2.13.0's recompute path and its padding row was removed; `_build_version` absent ⇒ we genuinely do
not know how the rest of it was built. If a downstream consumer needs "is the GC trustworthy", it
should read `_repair_history`, and we should say so in the architecture doc.

**Do not add a `_gc_source` scalar.** An earlier draft proposed one as a cheap alternative for
consumers that do not want to parse JSON. It is two sources of truth for one fact, and the failure
mode is the pair disagreeing after a second repair — precisely the situation `_repair_history` was
made a list to handle. One field, one truth. Consumers parse the JSON.

Confidence: **high** on the do-not-fabricate call; **medium** on the exact attr names, which should
be reviewed against `docs/architecture/fragment_selection_and_build_provenance.md` and updated there
in the same PR.

---

## 10. CLI and the execution gate

Policy is not a mechanism. The mechanism:

1. **`--dry-run` is the default.** Omitting all flags produces a report and touches nothing. There is
   no way to write by accident, and the report covers everything the apply path would do: per-contig
   truncation verdicts, the `gc` diff histogram, the §7.2 oracle result, and the §7.2.1 two-reference
   byte-diff.
2. **`--apply` is required to write**, and additionally requires:
   - `--target-list FILE` — an explicit newline-delimited list of full S3 keys. **Prefix expansion /
     globbing is disabled in `--apply` mode.** You cannot say "everything under this prefix" while
     writing; you must have materialised the list, reviewed it, and checked it in. **Duplicate keys
     are rejected at load time**, in both modes — two workers on one key produce two
     `_repair_history` timestamps and a racing upload, which is pointless even though it is no longer
     dangerous. The list is a set; enforce it where it is read.
   - `--expect-fasta-sha256 <hex>` — pins the reference bytes (§5 layer 2).
   - `--ledger PATH` — where the per-file records go (§8.2).
   - `--max-files N` — hard cap, **default 1**. Scaling from 1 → 218 requires deliberately raising it
     on the command line, which makes §7.4's stages mechanical rather than aspirational.
3. **No production identifiers in the package.** No bucket names, prefixes, FASTA URIs, or file lists
   are hardcoded in `src/`. They live in the target list and the runbook. This is what makes the tool
   "general" in the only sense that matters here — it has no opinion about which files exist.
4. **Human gate.** The target list file and the runbook are committed to the repo and approved in the
   PR. The user's approval to *run* is recorded there, separately from approval of this design.

**Deleted in round 3: the interactive typed-count confirmation.** Round 1 replaced a `--yes` /
`--i-know-this-is-production` pair with a single unbypassable prompt requiring the operator to type
the literal file count. That was the right call *then*. It is now redundant with `--max-files`: both
exist to say "you must deliberately opt into scale", which is two sources of truth for one fact — the
exact error §9.2 rejects for `_gc_source`. `--max-files` is the better of the two because it is
recorded in shell history and in the runbook, whereas a typed prompt leaves no trace of what was
approved. One gate, in the place that logs itself.

**Where the target list comes from.** `/tmp` was wiped, so `ibd_manifest_gc_audit.tsv` — the 763-row
audit that identified the 218 — **no longer exists**. The 218 keys survive inside
`~/gc_repair_artifacts/gc_repair_backup_ledger.tsv`, which is now the only durable record of the
target set. Derive `--target-list` from it and **commit both** under `docs/pending/gc_repair/` before
stage 3. This is a prerequisite, not a nicety: a target list that cannot be reproduced is not
reviewable.

---

## 11. Scope: what we deliberately do not build

The user's instruction — *"a general repair/update script going forward. Don't over-engineer but also
do not treat as a throw-away"* — resolves as: **general in interface, specific in behaviour.** §0
records how the balance moved once the backup existed: the *reusable skeleton* got smaller, the
*correctness checks* did not.

**Generality we are building in** (all cheap, all load-bearing for the current task anyway):
- A flat, config-free CLI with no hardcoded targets (§10.3), so pointing it at a different file set
  is a command-line change.
- The download → mutate → verify → upload → ledger pipeline as the reusable skeleton. The only things
  specific to this incident are two functions — "recompute `gc` for one contig" and "detect and drop
  a trailing padding row" — plus their preflight rules.
- `_repair_history` as an append-only list (§9.1), so a second repair does not erase the first.
- **Idempotency as a design requirement** (§3.2.1, §8.1) rather than a happy accident. This is the
  one property that made it possible to delete the entire crash-resume protocol, and any future
  in-place S3 mutation will want it for the same reason.
- The §5b no-op invariant expressed as `--dry-run`-by-default, which is a property of the skeleton,
  not of GC.

**Refused, explicitly:**
- **In-tool backup orchestration** — backup verification, per-file barriers, write-once backup keys,
  the pinned `(source key, backup prefix, ledger)` triple, and the startup refusal checks. The backup
  is done and verified out of band (§2.5); building a second mechanism to create one would be
  building it for a job that has already been completed. **Deleted in round 3, not gated.** The
  technique that made it safe (`copy_object(IfNoneMatch="*")`) is recorded in §2.5 for reuse.
- **The two-phase `started`/`ok` ledger protocol** — replaced by idempotency (§8.1).
- **A `fragments_h5.gc_audit` module in `src/`.** Detection is a **finished job**. Note the audit TSV
  it produced was lost with `/tmp` and the target set now lives in the backup ledger (§10) — which
  strengthens rather than weakens the case: what has value is the *result*, which must be committed,
  not the *detector*, which will not run again. A future incident will need a different detector.
- **A `--restore` flag** (§8.3) — the runbook is the restore path.
- **A typed-count confirmation prompt** (§10) — redundant with `--max-files`.
- **An override for the `255 → x` hard error** (§5, §5b). That transition means the wrong reference
  was supplied; there is no legitimate use for continuing. Deleted, not gated.
- **A `_gc_source` scalar attr** (§9.2) — one source of truth.
- **A general "shrink any dataset" utility.** §3.2 truncates exactly one trailing row under a
  specific fail-closed predicate. Parameterising the row count, the predicate, or the dataset set
  would be building for a case that does not exist.
- A plugin/registry system for "repair kinds". There are exactly two repairs and they ship together.
  When there is a third, the honest refactor is a third function, not a framework.
- A generic HDF5 schema-migration engine.
- Multipart `put_object` / >5 GB support (§3.5) — refused with a clear error instead.
- Re-deriving anything from BAMs. The whole point of the clean-replay validation is that we never
  need BAM access.
- Bucket scanning / automatic corruption discovery. The tool consumes an explicit list.
- AWS Batch / distributed execution. 91 GB and 16 CPUs is a laptop-scale job.
- Any emulation of object versioning, retention policy automation, or lifecycle rules.
- A progress UI beyond structured log lines.

**Artifact disposition — `/tmp` was wiped and most of the prototype evidence is gone.** The rule
applied: **commit evidence, not scratch code** — which round 2 stated and did not act on in time.

| artifact | status | action |
|---|---|---|
| `~/gc_repair_artifacts/gc_repair_backup_ledger.tsv` | **survives** (NFS home) | **Commit** under `docs/pending/gc_repair/`. It is now the only durable record of the 218 target keys and of their pre-repair checksums. |
| `/tmp/ibd_manifest_gc_audit.tsv` (763 verdicts) | **LOST** | Not regenerated. The 218 are recoverable from the backup ledger; the 545 clean verdicts are not, and are not needed by this tool. If the full audit is wanted again it is a fresh S3 probing job — out of scope here. |
| REF-P12 sequence-byte histogram (§4.2) | **LOST as a file; the numbers survive in §4.2** | **Regenerate** (~5 min, §3.6) as a mandatory §7.4 stage-0 task and commit under `docs/pending/gc_repair/ref_p12_sequence_histogram.json`, keyed by FASTA sha256. Regeneration also re-verifies §4.2's numbers. |
| `/tmp/gcpilot/{probe.py,s3h5.py,fallback.py}`, `/tmp/gcpilot/objects.tsv`, `/tmp/fasta_fingerprint.py` | **LOST** | **Do not rebuild.** The repair tool operates on local files and needs no byte-range S3 reads; the fingerprint script is superseded by the strictly more general §5 layer-3 contig-membership rule. |

The lesson, recorded because it cost real work: **evidence that lives only in `/tmp` is not
evidence.** Round 2 identified two artifacts as "irreplaceable" and scheduled them for commit; both
were lost before the commit happened.

---

## 12. Open questions and unverified claims

Ordered roughly by how much they could hurt. Entries discharged from source or by measurement are not
restated as "resolved risks" — they are gone, and the evidence is in the §4.3 table.

**Closed since round 2, and why** (recorded so a reader does not go looking for them): the
`copy_object` + SHA-256 question is closed **by substitution** — §6.2 now uses CRC64NVME, which the
backup run exercised end-to-end; the raw-chunk-hashing cost question is closed **by deletion** — the
truncation rewrites every dataset, so there is no byte-identity claim left to hash for (§3.1); and
"no real fragment H5 has been inspected" is **closed by inspection**, which is how §2.2 exists.

1. **R1(b): arbitrary drift of the REF-P12 object itself is unprovable, and §7.2 does not close it.**
   The bucket has no versioning and `ListObjectVersions` is `AccessDenied`, so the object's history
   is unrecoverable and there is no candidate to diff against. §7.2 pins the FASTA over the
   pre-saturation prefix — **the region the tool does not write** — which is corroboration, not
   proof. This entry stays open permanently; the backstop is that the failure is recoverable from the
   verified backup (§8.4).
   *The narrower question, R1(a) "REF-P12 or REF-ASSETS?", is decidable:* §7.2.1's two-reference
   byte-diff closes it outright per file where the count is zero. It is still a *plan*, not a result.
2. **The truncation path has no real-data coverage until §7.4 stage 2.** This is now the largest
   untested surface, and it is untestable by the tool's own primary oracle: clean files post-date
   `778f4d1` and have no padding row, and the whole `output_rebuild_frag_h5s/` rehearsal pool does
   too (§7.1.1). Coverage before production is a synthetic fixture plus **one** restored copy of a
   real affected file. A structural difference we have not thought to model would not be caught — and
   the padding row itself is the existence proof for that failure mode, having survived two review
   rounds of source-only reasoning.
3. **R8: 217 of the 218 have not been checked for the padding row.** Uniformity is near-certain
   (single build, single shared helper, known window) but inferred. §7.4 stage 3 makes 218/218
   all-contigs-truncated a gate, so this closes cheaply — but it is open until then, and a *mixed*
   verdict within one file would mean the mechanism is not what §2.2.1 says it is.
4. **HDF5 space behaviour under ~1,200–1,900 delete-and-recreates is unmeasured** (§3.2.2). The one
   synthetic experiment rewrote a *single* dataset *in place at the same shape* and says nothing
   about free-space reuse across many small unlinks or about metadata leak. Worst case is ~2× a
   from-scratch build. Gated at §7.4 stage 2 with a 5% threshold and an `h5repack` fallback whose
   availability on this machine is **itself unverified**. The rewrite is also the largest unmodelled
   term in §3.6's runtime estimate.
5. **The `round(x, 5)` vectorization is unimplemented and unbenchmarked** (§3.4). We know NumPy's and
   CPython's five-place rounding can disagree and that the disagreement is observable in the stored
   byte; we do not yet know whether the exactly-correct elementwise path is fast enough, and the
   fallback requires an agreement test that has not been written.
6. **The middle-band change volume is predicted, not measured**, and so is the float32-accumulator
   simulation §5b uses to classify fragments into bands. The *sufficiency* of the N-containing-span
   predicate has a proof (§4.2 Consequence 5); the proof has not been checked against a real file,
   and neither has the magnitude prediction (tens of uint8 units on a fragment straddling an N-block
   edge). The simulation is exact in principle but has not been written or run.
7. **Derived H5s built from the 218 are out of scope and unenumerated** (§9.1.1). At least one real
   analysis path (`liftover_porcine_fragments.py:193-197`) rebuilds fragment H5s from a 5-attribute
   allowlist, so any such derivative carries corrupted GC **and** an inherited phantom row, will not
   be repaired by this tool, and will not advertise that it needs to be. Someone has to inventory
   them; this doc does not.
8. **`rebuild_all_frag_h5s.nf`'s FASTA pinning is not verified from source.** The `.nf` lives in the
   separate `omni` monorepo and is not on this machine; the claim comes from
   `AGENT_CONTEXT.md:~668-709`. It is corroborated empirically — all three probed
   `output_rebuild_frag_h5s/` files fingerprint REF-P12 — but corroboration is not verification. This
   matters only for §7.1's target selection, where the fingerprint is direct evidence anyway.
9. **Downstream impact of changed GC *and changed fragment counts* on already-published CD/UC
   results** is out of scope for this doc but is not out of scope for the project. The fragment-count
   change is new in round 3 and is the more visible of the two, since `n_fragments` is read by code
   that `gc` is not.

---

## 13. Self-assessment

**What is solid.** The arithmetic is settled by measurement rather than argued from assumption. The
REF-P12 sequence alphabet is exactly `{A,C,G,T,N}` (§4.2), so the reachable per-base C+G set is
`{0, 0.5, 1.0}`, the float32 pairwise add at `fragment.py:442` is exact, the float64 cumsum is exact
with a margin of ~4×10^7, and the 2**23 / 2**24 boundaries are derived rather than guessed. The
middle band has an actual proof behind its rule (§4.2 Consequence 5). That argument is conditional on
the reference, and §5 layer 5 turns the condition into a gate that fails closed. R2, R3 and R7 are
discharged from source or measurement, with citations.

The safety story is now *simpler* and no weaker. The download-repair-upload shape makes mid-write
corruption in S3 structurally impossible. The originals are backed up and verified out of band, and
**the tool has no code path that writes to the backup prefix** — a stronger guarantee than the three
rounds of invariants that previously guarded one. Crash-resume is handled by idempotency rather than
by protocol, which removes roughly 150 lines of mechanism and, with it, the class of bug where the
mechanism is misconfigured.

**Round 3's actual finding, and why it dominates the grade.** Two review rounds reasoned entirely
from source and produced a document that was wrong about the shape of the data. The first time
someone opened a real file, a defect appeared that (a) affects every dataset rather than one,
(b) would have caused ~41,000 NaN divisions had the GC repair run as specified (§3.3), (c) would have
aborted all 218 files on a spurious `x → 255` hard error, and (d) is a **live query-correctness bug
and a live `n_fragments` bug in production today**, independent of GC. None of the three prior
Criticals found it, because none of them looked.

**The weakest point.** The truncation path is the largest thing in this document and it has the least
evidence behind it. Its primary oracle is structurally unavailable — a file that is clean for GC is
necessarily clean for padding, so §7.1 can only prove the detector does not misfire (§7.1.1). Real
coverage before production is a synthetic fixture, which by construction covers only differences we
thought to model, plus exactly **one** restored real file. And §12.2's uncomfortable point stands:
the padding row is itself the proof that "we modelled everything" is not a claim this document has
earned.

The second weakest point is unchanged and structural: §7.2 pins the FASTA over the pre-saturation
prefix, which is **exactly the region the tool does not write**. Every byte the repair emits derives
from sequence past ~40.45 Mb, which the oracle never touches. §7.2.1's two-reference byte-diff closes
the answerable half (REF-P12 vs REF-ASSETS); arbitrary drift of REF-P12 itself has no candidate to
diff against and therefore no possible check.

**Grade: B+, unchanged from round 2, and the stability is the honest signal.**

Three things improved: one real file was finally opened, so the foundation is no longer entirely
inferred; the backup landed and was verified, so the operation is reversible; and ~250 lines of
now-dead mechanism were deleted rather than caveated, which makes the remaining mechanism easier to
review.

Three things offset them exactly. The scope grew from one dataset to **every dataset plus the index
plus `fragment_length_counts`**, so the blast radius of a mistake is larger than what rounds 1–2
graded. The new mechanism is *specified but unwritten and unrun*, which is the same
"stated correctly and untested" position round 2 flagged, now applied to a bigger surface. And the
backup improves **recoverability, not correctness** — a grade on a design should not rise because the
consequences of the design being wrong got cheaper. If anything, §0's Fact 2 (nothing reads `gc`, so
nothing will notice a bad repair) means the detection burden on this document went *up* while its
recovery burden went down.

Resisting the upgrade is the point. The temptation after a round that deleted 250 lines and closed
three open items is to call it an A-. But the single most important event of this round was
discovering that the document's model of the data was wrong, and the correct response to that is not
a higher grade.

### 13.1 Changes from review

> **Round-3 note:** rounds 1 and 2 below are retained as history, but several rows defend mechanisms
> that round 3 **deleted** because the premise that motivated them (irreversibility) is false. Those
> rows are **moot**: C1, H5, M5, M7, C-A, M-C, L-g, and the backup-related half of C2 and of the
> `scope` row. They are not restated as live design anywhere in §§1–12. H-B and H-C are **partly
> moot**: the `_repair_history` carve-out survives (§9.1), but the per-dataset byte-identity claim it
> qualified does not, because §3.2 rewrites every dataset by design.

#### Round 1

| # | Finding | How it was addressed |
|---|---|---|
| C1 | Resume can destroy the only backup | New §6.5: backup keys write-once; `status:"started"` ledger record. **[MOOT — the tool no longer makes backups]** |
| C2 | `--dry-run` writes to S3 | §3.1 reordered so all gates precede any write. The write-free-dry-run half survives; the backup-ordering half is moot. |
| H1 | §4.2 misread `fragment.py:442` | Corrected: `.sum(axis=1)` runs in **float32** before the `.astype(numpy.float64)`. Discharged explicitly — the pairwise add is `0+0`, `0+1.0`, or `0.25+0.25`, all float32-exact. |
| H2 | The "multiple of 2**-25" lemma is false | Lemma deleted. Replaced by the measured enumeration `{0, 0.5, 1.0}` and the 2**52 float64 bound. |
| H3 | B/V/H/D presence — resolved by measurement | §4.2 leads with the histogram; states that B/V/H/D would have made the usable prefix **zero**; converts the check into a mandatory preflight gate (§5 layer 5). |
| H4 | §7.2 proves prefix identity, not reference identity | R1 downgraded from "Decisive"; §7.2.1 states the limit. **Superseded in round 2 by H-A.** |
| H5 | Sampled backup verification | Sampling deleted from §6.3. **[MOOT]** |
| H6 | Upload has the multipart-ETag problem | §6.2 mandates single-part `put_object`, forbids `upload_file` by name and threshold, requires asserting no `-N` suffix. **Strengthened in round 3**: the ETag trap is no longer a prediction — it produced 218 false mismatches in the real backup run (§2.5). |
| H7 | Vectorization / rounding fidelity | §3.4: `numpy.round(a,5)` ≠ CPython `round(x,5)`; blocking agreement test over >=10^7 fragments. |
| M1 | Two "unverified risks" are verifiable | R3 and R7 discharged from source. |
| M2 | Three-region structure never modeled | §2.1's three-region table; §5b's per-region rule. |
| M3 | §12.7's causal claim backwards | The 2**-25 / 2**28 / 7.8% / 2.6× reasoning deleted entirely. |
| M4 | `_repair_history` dropped by real downstream paths | §9.1.1, citing the 5-attribute allowlist. |
| M5 | TOCTOU between download and backup | `CopySourceIfMatch` on the backup copy. **[MOOT — no backup step; and §2.5 records `IfNoneMatch="*"` as the better primitive, now tested]** |
| M6 | `gc is None` → 255 path unspecified | §3.4's `length == 0` rule. **Round 3 found this was not hypothetical**: the phantom row makes it fire ~41,000 times, which is why §3.3 orders truncation first. |
| M7 | Resume falls back to the rejected ETag comparison | Fallback deleted. The no-ETag rule survives (§6.2); the three-row crash-window table is moot. |
| L1 | Contradictory scratch figures | Standardised on "549 GB total / 414 GB free". |
| L2 | No runtime estimate | §3.6 adds one, labelled an estimate. **Widened in round 3** to 3–8 min/file for the truncation rewrite. |
| L3 | Citation drift | `fragment.py:435` → `:441`; `fragments_h5.py:566-575` → guard at `:562-564`, remap at `:570-574`. |
| L4 | `_contig_lengths_str` claim overbroad | §5 scopes it to BAM input, citing `fragments_h5.py:1024-1037` for the TSV/BED path. |
| L5 | Orphaned multipart parts | Moot under single-part `put_object`; §3.5 says so. |
| L6 | Version drift | Header notes `pyproject.toml:7` and the runtime `importlib.metadata` read. |
| scope | Cut `gc_audit`, `_gc_source`, `--restore`, prototype commits, `--yes` | All adopted and all still hold. **Round 3 note:** the "commit the evidence" half was never executed and `/tmp` was wiped — see §11's artifact table. |

#### Round 2

Round 2 diagnosed the round-1 failure mode as *sections fixed locally without a global re-read*.

| # | Finding | How it was addressed |
|---|---|---|
| C-A | Write-once is not bound to a stable backup prefix or ledger | New §6.6: the `(source key, backup prefix, ledger)` triple, ledger header record, three startup refusal checks. **[MOOT — deleted in round 3 with the backup step; §8.1 replaces the whole apparatus with idempotency]** |
| H-A | §7.2.1's N-block coverage check tested the wrong predicate | Coverage heuristic **deleted**; replaced by the two-reference byte-diff. R1 split into (a) decidable and (b) permanently open. **Survives unchanged.** |
| H-B | `_repair_history` falsified five "byte-identical" assertions | Carved out by name everywhere. **Partly moot in round 3**: the carve-out survives, but the per-dataset byte-identity claim it qualified is gone — §3.2 rewrites every dataset, and §5b's invariant is renarrowed to attrs plus the §3.2.3/§3.2.4 arithmetic invariants. |
| H-C | §3.1 step 8 was unimplementable — no operand | Pre-mutation raw-chunk hashes at step 2. **Mostly moot**: only the attr hash survives; the dataset hashes are meaningless once every dataset is rewritten. §12's old entry 8 deleted rather than resolved. |
| M-A | Middle-band mechanism described incorrectly | §2.1's exact-tie analysis; new §4.2 Consequence 5 supplies the proof. **Survives unchanged.** |
| M-B | Region classification used the true cumsum | §5b classifies by a **simulated float32 accumulator**. **Survives unchanged.** |
| M-C | §6.5 had no row for "backup exists, no checksum" | Row added. **[MOOT]** |
| M-D | No contig-length or coordinate-range preflight | §5 layer 3's geometry checks. **Survives**, and round 3 adds the ordering caveat: layer 3's all-255 assertion must run **after** truncation or it fails on all 218 (§3.3). |
| M-E | Repaired files are materially larger and the doc never said so | §3.2.2 states both baselines. **Survives**, and round 3 adds a second, larger and less bounded growth term from the delete-and-recreate rewrite. |
| M-F | When `_repair_history` is written was unspecified | §3.1 step 7, same handle, before close. **Survives.** |
| L-a | "24 GB against 96 GB available" vs a 124 GB environment | Corrected. |
| L-b | §9.1.1's cross-references were wrong | Repointed; round 3 renumbers them again (§12.7, §12.9). |
| L-c | Saturation offsets derived from GC-excluding-N | Accumulator density 0.41470; ~40.45 Mb and ~20.23 Mb, corrected at all sites and flagged as genome-average estimates. |
| L-d | §7.2 and §5b's `< T23` row are the same check | Unified: §5b is the implementation site, §7.2 explains the bound. |
| L-e | §4.1 quoted Stage B without its `ValueError` guard | Guard quoted inline at `fragments_h5.py:842-843`. |
| L-f | The pre-fix float32-cumsum line was never quoted | §2.1 quotes `95c76f5^:fragment.py:423`. |
| L-g | `--dry-run` never exercised §6.5 | Dry-run probes backup keys. **[MOOT]** |
| L-h | Duplicate keys in `--target-list` were not forbidden | Rejected at load time. **Survives**, with the rationale downgraded from "dangerous" to "pointless". |

#### Round 3

| # | Finding | How it was addressed |
|---|---|---|
| W | **Weight recalibration** — the irreversibility premise is false | New **§0** decides the question explicitly and asymmetrically: safety scaffolding deleted, correctness machinery retained in full, with the reason stated (nothing reads `gc`, so nothing will *notice* a bad repair — a backup does not buy detectability). Deleted: the in-tool backup step, backup verification, the per-file barrier, write-once keys, the `(source key, backup prefix, ledger)` triple, the three startup refusal checks, the `started`/`ok` two-phase protocol, the ledger header record, the 10-file rollout stage, and the typed-count prompt. §0 also records **two disagreements** with the framing given: the ledger is collapsed rather than deleted, and scale gating is unified onto `--max-files` alone. |
| C-α | **The padding row** — a second, independent defect affecting every dataset | New **§2.2** (mechanism, three live consequences), **§3.2** (fail-closed idempotent detection predicate, delete-and-recreate mechanism with costs, index rebuild, `fragment_length_counts` rebuild), **§3.3** (ordering). Two of its consequences are **live production bugs today**: unsorted `starts` breaks `searchsorted` in the final index block, and `n_fragments` is over-reported by ~189 per file. |
| C-β | **`length == 0` → 0/0 → NaN**, never mentioned in any prior round | §3.3 hazard 1: it would fire ~41,000 times. Mooted by ordering truncation first; the rule itself is retained for the general case (§3.4) and is now expected to fire **zero** times, with a nonzero count reported. |
| C-γ | The phantom `gc = 0` on scaffolds trips the `255 → x` machinery **in reverse** | §3.3 hazards 2 and 3: a pre-truncation recompute would move `0 → 255` on ~165 scaffold contigs and abort every one of the 218; and §5 layer 3's all-255 assertion is **false as written** on all 218. Both fixed by ordering, at the source and at both restatement sites (§5 layer 3, §7.1 criterion 2). |
| H-α | The validation oracle cannot reach the truncation path | New **§7.1.1**. Clean-for-GC ⇒ clean-for-padding, so §7.1 is a negative control only; and the whole `output_rebuild_frag_h5s/` rehearsal pool post-dates `778f4d1`. §7.4 stage 2's target changed to **a restored copy of one of the 218**. New §7.1 criterion 4 and new §7.3 structural checks incl. a **read-back query test** over the final index block. Stated plainly that a fixture covers only modelled differences. |
| H-β | Backup-run findings supersede two specified mechanisms | §6.2 switches upload verification from SHA-256 to **CRC64NVME** — the algorithm actually exercised end-to-end on these objects with these credentials — closing round 2's §12.5 **by substitution**. §2.5 records that **ETag verification produced 218 false mismatches** (27-part multipart sources) and that `copy_object(IfNoneMatch="*")` is API-enforced write-once, tested. |
| H-γ | The doc treats wiped `/tmp` artifacts as available inputs | §11's artifact table marks `ibd_manifest_gc_audit.tsv`, `gcpilot/*`, `objects.tsv` and `fasta_fingerprint.py` **LOST**; the REF-P12 histogram is **LOST as a file** and its regeneration is a mandatory §7.4 stage-0 task (§4.2 carries a status box). §10 states that the 218 target list now derives from the surviving backup ledger, which must be committed before stage 3. |
| M-α | R8 — 217 files unchecked for the padding row | New §4.3 row R8 and §12.3; §7.4 stage 3 makes 218/218-all-contigs-truncated a gate; a mixed verdict within a file is blocking. |
| M-β | Index rebuild is not a sentinel decrement | §3.2.3: the interior entries were produced by `searchsorted` over unsorted `starts`, so the index must be rebuilt from scratch by **calling** `fragments_h5.py:1225-1241` rather than reimplementing it; `INDEX_BLOCK_SIZE` read from the file's attr, not the module constant (which appears as `10000` at four other sites and `5000` at `:172`); key-set and sortedness invariants added. |
| M-γ | `fragment_length_counts` inflation was undiscovered | §2.2.3 point 3 and §3.2.4, with an **exact arithmetic identity** as the blocking check and the user-visible consequence (`n_fragments` drops by ~189/file) stated. §12.9 flags it as the more visible downstream change, since `n_fragments` is read by code that `gc` is not. |
| M-δ | §5b's invariant was falsified by the truncation | Renarrowed to attrs + §3.2.3/§3.2.4 invariants + the `gc` region rules, with the reason stated. Restatements fixed at §3.1 step 8, §7.3, §8.1, §11 and §13.1's round-2 H-B/H-C rows. |
| L-α | `mapq` is 2-D | §2.2.2: truncation is along **axis 0 only**; nothing may assume 1-D. §7.3 asserts `mapq` stays 2-D. |
| L-β | Runtime and footprint models predate the rewrite | §3.6 widened to 3–8 min/file (3–8 h total) naming the rewrite as the unmodelled term; §6.1's footprint rewritten for a world where the backup is already spent (~15–25 GB net growth, not ~190–200 GB). |
| L-γ | `_repair_history.datasets` understated the change | §9.1 now lists the truncation, the index and `fragment_length_counts`, and records `rows_removed_per_contig` — a reader seeing only `data/*/gc` would wrongly assume fragment counts were preserved. |
