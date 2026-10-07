# Design: a robust single-end fragment store

Status: **DRAFT, implementable.** One open question remains (§7) and it has a stated default.
Scope: `fragments_h5` only. No downstream consumer changes.

Precondition: cutadapt performs adapter trimming only, with no quality trimming and no
fixed-length trimming. §4 depends on this. If upstream trimming changes, revisit §4 first.

## 1. Problem

`single_end_bam_to_fragments` (`fragment.py`) sets the fragment extent from the read's
aligned span:

```python
frag_start = align.pos
frag_stop  = align.aend
```

Two separate defects follow.

**It cannot tell coverage from truncation.** The span is the fragment when the read covered
the whole fragment. It is a lower bound when the fragment extends past the read. Both store
identically. A consumer reading a stored length of 150 cannot know whether that is a 150bp
fragment or the first 150bp of a 400bp one.

**It discards soft-clipped fragment bases.** `align.aend` excludes soft clips. Adapters are
removed by cutadapt *before* alignment, so a surviving soft clip is not adapter — it is real
fragment sequence that failed to align. The span understates the fragment by the clip length.

The CLI help calls `--single-end` "useful for long read technologies". For long reads the
read-covers-fragment assumption holds by construction and neither defect bites. Both bite as
soon as the input is short-read SE data.

Goal: store SE fragments whose extent is correct, and whose end-to-end coverage is **proven**
rather than assumed.

## 2. What already exists — verified, do not rebuild

**The storage field exists.** `Fragment` carries `fragment_end_clipped`, a per-contig `uint8`
dataset with a tri-state encoding already documented in `fetch_array`:

```
0 = False   1 = True   255 = unknown
```

`has_fragment_end_clipped` and `return_fragment_end_clipped` already gate reads on it. The
tri-state is the right shape — "unknown" is a distinct answer from "not clipped".

**GC is already reference-derived.** The SE path computes GC from the FASTA via
`get_g_or_c_cumsum` over `[frag_start, frag_stop)`. It never reads the read's bases. SE GC is
therefore the fragment's genomic GC, and it is free of sequencing error.

**What is missing:** the SE builder passes `fragment_end_clipped=None`, so every SE fragment
stores 255. The field is plumbed end to end and never populated on this path.

## 3. Why the PE detector cannot be reused

The PE builder computes:

```python
if align.has_tag("MC"):
    fragment_end_clipped = not (
        cigar_fragment_end_matches(align.cigarstring, align.is_reverse)
        and cigar_fragment_end_matches(align.get_tag("MC"), align.mate_is_reverse)
    )
```

It needs `MC`, the mate's CIGAR. SE input has no mate.

The geometry is also inverted. `cigar_fragment_end_matches` tests the fragment-**outer** end:
CIGAR starts with `M`/`=` for a forward read, ends with `M`/`=` for a reverse read. Each PE
mate owns one fragment end. An SE read owns both, and the far end is the opposite side:

| read | fragment 5' end | far end (3') |
|---|---|---|
| forward | left, CIGAR start | right, CIGAR end |
| reverse | right, CIGAR end | left, CIGAR start |

SE needs a new predicate, complementary to the existing one. Do not conflate the two; they
answer different questions.

## 4. Coverage criterion — decided

Adapters are trimmed by **cutadapt** before alignment, and cutadapt performs **adapter
trimming only — no quality trimming**. The instrument read length is a **required** CLI
option, not inferred and not defaulted.

Adapter-only trimming is what makes the criterion sound, so state the chain explicitly:
cutadapt found adapter in the read; adapter follows the fragment; therefore the read contains
the whole fragment followed by adapter; therefore the read reached the fragment's far end.
Each step is an implication, not a correlation.

If quality trimming were ever enabled upstream, this chain breaks and the criterion starts
producing false "complete" flags — a read shortened by a poor 3' tail looks identical to one
shortened by adapter. Treat adapter-only trimming as a **precondition of this design**, not an
incidental fact.

The criterion follows directly:

| observation | meaning | `fragment_end_clipped` |
|---|---|---|
| `query_length < read_length` | cutadapt removed adapter, so the read ran past the fragment's far end | `0` — fragment is complete |
| `query_length == read_length` | no adapter was found, so the read did not demonstrably reach the far end | `255` — unknown |

`query_length` is the trimmed read length, because trimming precedes alignment. It includes
soft-clipped bases, so it is the right quantity.

### Why `255` and not `1` for the full-length case
A fragment of exactly `read_length` is covered exactly, yet produces no adapter and so is
indistinguishable from a longer one. Storing `1` would assert "the end is clipped", which is
false for that fragment. `255` is the honest value. Consumers that need provable coverage
filter on `== 0`.

This partitions the store, measured on the libraries in §9: roughly **85–94%** provably
complete, **6–15%** unknown, and the unknown set is exactly the long tail.

### The criterion errs in the safe direction
Worth recording, because it determines how the flag may be used.

**No false "complete".** Given adapter-only trimming, a shortened read provably saw adapter
and therefore the fragment end. A fragment flagged `0` is complete.

**Some false "unknown".** cutadapt needs a minimum adapter overlap to trim — `-O`, default 3.
A fragment 1 or 2 bases shorter than the read length leaves too little adapter to detect, so
the read stays full length and the fragment is flagged `255` despite being complete. With a
150bp read and the default `-O 3`, fragments of 148–150bp fall in this blind spot.

So `255` means "not proven complete", not "incomplete". The `255` bucket slightly over-counts,
by a sliver of mass adjacent to the read length. Consumers may trust `0`. They may not read
`255` as evidence that a fragment is long.

**One configuration would break the comparison.** If cutadapt also applies fixed-length
trimming (`-u` / `--cut`), every read is shortened uniformly and `query_length < read_length`
becomes universally true. The required read-length option must then be given as the
**post-trim** length. Verify no fixed trimming is configured before relying on the flag.

### Why read length must be required
Inferring it as `max(query_length)` fails on any library where no read survives trimming at
full length, and fails silently. The caller knows the value; make them state it. A wrong
value corrupts every flag in the file, so it is not a safe default.

## 5. Soft-clip extension

Because cutadapt removes adapters pre-alignment, a residual soft clip is real fragment
sequence. Extend the span through it so the stored extent is the fragment, not the aligned
subset:

- forward read: `frag_stop = align.aend + (3' soft clip length)`
- reverse read: `frag_start = align.pos - (3' soft clip length)`

Interactions that must hold:

- Extension and the §4 flag are **independent**. Extension corrects the extent. The flag
  records whether the far end was reached. A full-length read with a 3' clip extends and
  still stores `255`.
- Extension must respect `MAX_FRAG_LENGTH` (65535) and `se_max_fragment_length`. Apply the
  checks **after** extension, since extension can only increase the span.
- Extension must not push `frag_start` below 0 or `frag_stop` past the contig length. Clamp,
  and count clamped fragments.
- GC is computed from the span, so extension changes GC. It makes GC more correct, over the
  fragment rather than its aligned subset. This changes output values for existing SE files
  and so is a results-moving change.

**Assumption, stated because it is load-bearing:** extension is correct only if adapter
removal is complete. Residual adapter that cutadapt missed would be extended *into* the
fragment, inflating the length and corrupting the GC. See §7.

## 6. Proposed behaviour

1. Add a required read-length option to the SE build path.
2. Populate `fragment_end_clipped` per §4. Never store `0` without proof.
3. Extend spans per §5.
4. Add a build-time option to **drop** fragments that are not provably complete, default off.
   Consumers needing a clean length distribution opt in; others filter on the flag.
5. Count and report, per contig: fragments flagged unknown, fragments extended, total
   extension in bases, fragments clamped, fragments dropped. A silent filter over 6–15% of
   fragments is the defect this design removes; it must not return as a side effect of the fix.

### Deliberately out of scope
- `se_max_fragment_length` semantics. It caps storage, a different concern.
- Estimating the length of a fragment that is not provably complete. Unknown length stays
  unknown; it is never imputed.
- The PE builder and `cigar_fragment_end_matches`.

## 7. Open questions

1. **Extend both ends, or only the far end?** §5 extends the far end only, which is the
   current recommendation and the stated default. A 5' soft clip is also real sequence under
   the same argument. But the 5' end is the fragment's cut site, and cut sites are a primary
   fragmentomics signal. Moving a cut site on soft-clip evidence is a stronger claim than
   extending the opposite end. If far-end-only stands, state the asymmetry in the code so it
   reads as a decision rather than an oversight.

### Resolved
- **Does cutadapt quality-trim?** No. Adapter trimming only. This was the one question that
  could invalidate §4; it is settled and the criterion stands. Recorded as a precondition in
  §4, because a future change to the upstream trimming step would silently break the flag.

## 8. Validation

Build one **paired-end** BAM twice: once with the PE builder, once with the SE builder reading
R1 only. `TLEN` gives the true fragment length, so the PE result is ground truth for the SE
criterion on identical input.

Report:
- false positive rate: flagged `0` but PE shows the fragment is longer
- false negative rate: flagged `255` but PE shows the read did cover the fragment
- extension accuracy: extended span versus `abs(TLEN)`
- the SE length distribution against the PE distribution, for fragments flagged `0`

Nothing else measures the criterion. Until this exists, "robust" is an assertion, not a
finding.

## 9. Measured context

From 6 IBD round-3 fragment h5 files, whole-genome `fragment_length_counts`, 1.7 billion
fragments: fragments at or below 150bp are **84.72% to 93.90%** across samples, mean 89.66%,
median fragment length 60 to 76.

So at a 150bp read length, 85–94% of fragments are provably complete and 6–15% are unknown.
The unknown set is not random — it is the long tail, and it contains most of the
mononucleosome population (150–167bp is 4.3–4.9% of mass, 167–180bp is 2.5–2.9%).

**Consequence for consumers:** a store filtered to provably-complete fragments has a hard
length ceiling at the read length, and that ceiling cuts through the mononucleosome peak. Such
a store is sound for GC-bias estimation, which needs unbiased GC per fragment. It is censored
for length-distribution work. This is why dropping is opt-in rather than the default.
