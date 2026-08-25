# Fragment Selection Correctness and Build Provenance

Supersedes `docs/pending/build_provenance_metadata.md`. That document proposed a
`requested`/`effective`/`neutralized` provenance schema and was rejected as
over-engineered; do not extend it.

## Goal

Five changes ship together as one release.

Three of them change **which fragments enter an h5**:

1. Secondary alignments are currently accepted and can emit spurious duplicate
   fragments (measured: 3 extra SE fragments, 1 extra PE fragment on synthetic BAMs).
2. Single-end fragment length is unbounded unless the caller passes
   `--se-max-fragment-length`; an over-long reference span reaches a `uint16`
   assignment and fails there.
3. `num_mapped = x.mapped // 2` drops any contig with exactly one mapped
   alignment record — in single-end builds the contig vanishes from the output
   file entirely, with no warning.

Together those three mean that a file built before this release and a file built
after it, from the same BAM with the same arguments, can differ — and today
**nothing in the file records which side of the change it was built on**. That is
what makes change 5 part of the same work rather than a separate feature:

5. Record build provenance — the CLI argument vector and the builder version — so
   the difference is detectable from the file itself. The BAM header is already
   stored (`fragments_h5.py:1149`), so provenance is deliberately just the two
   missing pieces.

The fourth change is adjacent correctness in code already being touched by this bundle:

4. S3 (and any other remote) input is broken by an unconditional
   `os.path.abspath` at `fragments_h5.py:963`, which rewrites
   `s3://bucket/key.bam` to `/<cwd>/s3:/bucket/key.bam`. Two tests in
   `tests/test_s3_bam.py` fail today because of it.

The consistency story is: **1-3 change the output, 5 makes that auditable, 4 is
in the blast radius anyway.** Shipping 1-3 without 5 produces silently
incomparable files. Shipping 5 without 1-3 records provenance for behavior nobody
has a reason to distinguish yet.

A useful property of this bundle, which the design deliberately preserves: because
1-3 are all **unconditional** fixes rather than new flags, the recorded builder
version fully determines them. That is what makes the minimal provenance in
section 5 sufficient (see "Alternatives rejected — argv-only recording").

Recommended release: **2.12.0**. Behavior changes, new attributes, no format
version bump, old files continue to open unchanged.

## Design

### 1. Exclude secondary alignments

**Decision: unconditional. No flag, no opt-out.**

The filters are hand-duplicated, not shared, so the fix lands in two independent
places:

- PE, closure `align_is_invalid` inside `bam_to_align`, `fragment.py:368-384`
- SE, inline in `single_end_bam_to_fragments`, `fragment.py:573-579`

In each, add `align.is_secondary` next to the existing `align.is_supplementary`
term. Nothing else changes. The TSV path (`tsv_to_fragments`) reads plain text
rows with no SAM flags — inapplicable, not a gap.

`is_secondary` appears **zero times** in `src/`, and `git log -S'is_secondary'`
returns zero commits: this was never implemented and never deliberately rejected.
`pysam.fetch()` does not drop secondary alignments, and no other filtering pass
exists.

**Why unconditional rather than opt-out.** A secondary alignment is by SAM
definition a non-representative placement of a read that already has a primary
record elsewhere. There is no fragment-level interpretation under which counting
it is correct: it does not represent an additional molecule. Making it a flag
would (a) require an 18th element in the worker-args tuple at
`fragments_h5.py:735` — the exact coupling this design is required to avoid —
and (b) make the behavior un-determinable from `_build_version` alone, which
undercuts the provenance design in section 5.

**Mechanism, measured, worth stating precisely because the intuitive story is
wrong.** The PE path does **no query-name mate matching** — there is no qname
dict anywhere in `src/`. It selects the fragment-producing record purely by
`0 < align.tlen <= max_tlen`. A secondary record therefore does not "confuse
pairing"; it is evaluated independently and, if it happens to carry a positive
in-range TLEN, silently emits its own duplicate fragment. Measured on synthetic
BAMs with real `build_fragments_h5` runs:

- SE: 5 primary reads produce 5 fragments. The same BAM plus 3 secondary
  alignments produces **8** fragments; the 3 extra starts (5000/5100/5200) are
  exactly the secondary reads' positions.
- PE: 6 pairs produce 6 fragments. The same plus 1 secondary alignment carrying a
  positive TLEN produces **7** fragments (spurious start at 70000).

No crash. Silent output corruption.

**Effect on existing files.** Likely none, but this is a **sampled census plus an
inference from aligner documentation, not a proof**, and the document will not
claim otherwise. Evidence: 0 secondary alignments were observed in ~61,000
sampled reads across local fixtures and two real production BAMs. The aligners in
use are `bwa-mem2 mem -t 30 -k 15 -r 1.25 -T 19 -U 30` and
`bowtie2 --very-sensitive`; neither is given `-a` or `-k`, and neither emits
secondary alignments without them. `biomarker-pipeline`'s microbe minimap2 step
does use `--secondary=yes`, but that BAM does not feed `build_fragments_h5`.

**What if an aligner config later gains `-a`/`-k`?** That is the argument *for*
this change, not against it. Today such a config change would silently inflate
fragment counts in every downstream h5 with no error and no signal. After this
change it is a no-op. The filter is cheap insurance against a change made in a
different repo by someone who has never read this code.

No new worker-tuple element.

### 2. Bound single-end fragment length: raise, do not wrap, do not clamp

**Correction to a widely repeated claim: this does not silently wrap in any
deployed environment.** numpy >= 1.24 raises `OverflowError` on out-of-range
integer assignment; numpy < 1.24 wrapped with a `DeprecationWarning`.
`conda-recipe/recipe.yaml` pins `numpy >=1.26` in both `host:` and `run:`, and
`environment.yml` (used by the Docker image) leaves numpy unpinned so it resolves
to current. `pyproject.toml` also leaves it unpinned. **Silent wraparound is
reachable only through a pip install into a pre-1.24 numpy environment.**

Measured today under numpy 2.2.6, real builds, CIGAR `10M70000N10M` (reference
span 70020), `se_max_fragment_length=None`:

```
File ".../fragments_h5.py", line 1123, in build_fragments_h5
  contig, chunk_start, chunk_stop, sub_h5_path = build_sub_fragments_h5(arg)
File ".../fragments_h5.py", line 801, in build_sub_fragments_h5
  lengths_arr[ff] = fragment.length
OverflowError: Python integer 70020 out of bounds for uint16
```

Boundary sweep, measured: span 65534 OK (`lengths == [20, 65534]`), span
**65535 OK**, span 65536 `OverflowError`, span 65537 `OverflowError`. The
boundary is exactly `uint16` max, which equals `MAX_FRAG_LENGTH`
(`fragments_h5.py:176`).

So the value of this fix is **not** "make it raise" — it already raises. It is:

- **determinism**: the same input must fail the same way regardless of the numpy
  the user happens to have;
- **actionability**: the current error fires inside a multiprocessing worker,
  names no contig, no position, no read, and points at an array assignment that
  means nothing to a user.

**Implementation, corrected.** An earlier draft of this section put the raise
"after the existing `--se-max-fragment-length` filter has had its chance to skip
the read." **That ordering does not exist.** `single_end_bam_to_fragments`
(`fragment.py:548-608`) never receives `se_max_fragment_length` — the call site
(`fragments_h5.py:772-777`) passes it `input_fname, chrom, start, stop, max_tlen,
fasta_file, include_duplicates, fasta_region_start, fasta_region_stop,
fasta_chrom, min_mapq` and nothing else. The filter lives in a different
function, `build_sub_fragments_h5` (`:782`), and only runs on fragments the
generator has *already yielded*.

The fix: add `se_max_fragment_length=None` as a new keyword parameter to
`single_end_bam_to_fragments`. Immediately after computing `frag_start =
align.pos` / `frag_stop = align.aend` (`fragment.py:581-582`), and **only when
`se_max_fragment_length is None`**, raise a `ValueError`:

```
Single-end fragment length 70020 exceeds the maximum storable length 65535
  contig=chr1 start=1000 stop=71020 read=read_a cigar=10M70000N10M
The reference span of a single-end read (including D/N CIGAR operations) is used
as the fragment length. Pass --se-max-fragment-length to skip reads whose span
exceeds a threshold; 65535 is a hard storage limit and cannot be raised.
```

The `is None` guard is load-bearing: an unconditional raise would fire before
`:782` ever gets a chance to skip the read, breaking the existing chunk-level
filter. With the guard, a bounded caller keeps today's skip semantics unchanged
(the generator yields, `:782` discards); an unbounded caller gets the raise.
**`:782` itself is untouched, not redundant, and should not move** — it remains
the only thing that implements the skip when a bound is supplied. A minor
consequence: the check fires on every read `af.fetch()` returns, including ones
outside `[chunk_start, chunk_stop)` that `:779-780` would otherwise discard, so
an oversized read straddling a chunk boundary can raise from two chunks instead
of being rejected once — acceptable, a raise is a raise.

**Wiring needs no 18th tuple element, but not for the reason originally
given.** `se_max_fragment_length` is already worker-tuple element 9, unpacked at
`:735` and in scope in `build_sub_fragments_h5`. The call at `:772-777` passes
one identical kwarg set to whichever of `tsv_to_fragments` /
`single_end_bam_to_fragments` / `bam_to_fragments` was selected;
`tsv_to_fragments` absorbs unknown kwargs via `**kwargs` (`fragment.py:645`) but
`bam_to_fragments` does not (`:444-456`), so `se_max_fragment_length` must be
added only when `single_end` is true — a one-line conditional, not a tuple
change.

**Resolving a circular import.** The raise needs `MAX_FRAG_LENGTH`
(`fragments_h5.py:176`), but `fragment.py` cannot import it from there:
`fragments_h5.py:153` imports `fragment.py` *before* `:176` defines the
constant, so `from fragments_h5.fragments_h5 import MAX_FRAG_LENGTH` inside
`fragment.py` is a circular, partially-initialized-module import. The
dependency direction is one-way (`fragments_h5.py` depends on `fragment.py`,
never the reverse; `fragment.py`'s own imports touch only
`fragments_h5.sequence`) — so **move `MAX_FRAG_LENGTH = 65535` into
`fragment.py` and have `fragments_h5.py:153` import and re-export it** under the
same name, leaving its existing use sites (`:765`, `:774`, `:873`, `:1148`)
unchanged. A function-local import was rejected — it would run per yielded read.

Raising from a `Pool` worker is acceptable: multiprocessing re-raises the
exception in the parent with its message intact.

**Exact boundary semantics: `>`, not `>=`.** A span of exactly 65535 remains
legal. This matches (a) the measured `uint16` limit, (b) the PE comparison
`abs(align_.tlen) > max_tlen` at `fragment.py:372`, and (c) the TSV comparison
`if frag_stop - frag_start > max_tlen: continue` at `fragment.py:707`. All three
paths share the identical predicate *shape* (`> bound`), but not an identical
constant in general: `max_tlen` is a **parameter** in PE and TSV, defaulting to
`1000` (`fragment.py:307`, `:450`, `:641`), equal to `MAX_FRAG_LENGTH` only
because `fragments_h5.py:774` passes it as `max_tlen`. A library caller doing
`bam_to_fragments(..., max_tlen=500)` gets a 500 PE bound while SE stays
hardcoded at `65535` — the SE bound is a constant, not a parameter.

**Reconciling SE-raises against PE/TSV-silently-skips.** This asymmetry is
deliberate and it turns on what the quantity means, not on which code path it is.

- In PE, `tlen` is a routine aligner-reported quantity. Absurdly large TLENs are
  expected and common (chimeric pairs, translocation-spanning pairs, mismapped
  mates). `max_tlen` is threaded in from the caller as a *filter parameter*
  (`fragments_h5.py:774` passes `MAX_FRAG_LENGTH`), and skipping is its
  documented, intended job. Same for TSV.
- In SE there is no TLEN. The "length" is `aend - pos`, the reference footprint
  of a single read. A span above 65535 means a D or N CIGAR operation spanning an
  intron-scale gap — that is not a long fragment, it is a read whose reference
  span is not a fragment length at all. It signals that the input or the
  parameters are wrong.

The reconciliation, matching the corrected implementation above: **the skip
behavior already exists for SE, is opt-in, and is untouched by this change.**
`--se-max-fragment-length` is a chunk-level filter at `fragments_h5.py:782`,
measured to correctly exclude the 70020-span read when set to 1000; the
generator-level raise is guarded on `se_max_fragment_length is None` precisely so
it never fires while that filter is active:

- caller sets `--se-max-fragment-length` -> over-long SE reads are skipped,
  exactly like PE/TSV, unchanged from today;
- caller sets nothing -> hard error naming the read, rather than a
  numpy-version-dependent outcome.

**The CLI can never trigger the raise.** `main.py:82-84` caps
`--se-max-fragment-length` at `1`-`65535`, and `main.py:92-94` makes it mandatory
with `--single-end` for BAM input, so every CLI-driven SE build is already
bounded before `build_fragments_h5` is called. Honest framing: **this raise
protects library callers; CLI callers are already bounded** — which makes the
SE/PE asymmetry easier to accept, since it is "a caller who bypasses the CLI's
own validation gets a named error instead of an unnamed crash," not "SE users get
a worse experience than PE users."

`se_max_fragment_length` defaults to `None` in the library signature
(`fragments_h5.py:927`); only the CLI enforces the bound above, so **library
callers get no bound at all today.** One of the two production SE call sites
(`biomarker-projects/projects/lupus/scripts/build_frag_h5s.py:42`) passes
`--single-end` with no `--se-max-fragment-length`, via `os.system()` with the
return code unchecked — but it also passes `--reference hg38 --sample-id
{sample_id}`, neither of which `main.py`'s parser defines (`main.py:15-61`), so
that invocation already exits at argparse today. The conclusion (errors there
are invisible, since the return code is unchecked) still holds; the premise that
it is a currently-working SE build does not. A different repo's problem, out of
scope here.

**Dead parameter — dead, but NOT removable. Do not touch it.**
`single_end_bam_to_fragments` declares `max_tlen=1000` (`fragment.py:554`) and the
name never appears again in the body. It is dead and actively misleading — a reader
reasonably concludes SE is bounded at 1000.

An earlier revision of this document prescribed removing it. **That was wrong, and
it was proven wrong by execution**, not by reading. The hazard is not external
callers; it is the internal call site described at line 153 of this document. The
dispatch at `fragments_h5.py:737-742` selects a *function object*, and `:772-777`
is a **single shared call expression** that passes `max_tlen=MAX_FRAG_LENGTH`
unconditionally to whichever generator was selected. Deleting the parameter
therefore raises

    TypeError: single_end_bam_to_fragments() got an unexpected keyword argument 'max_tlen'

on **every** single-end build — including
`biomarker-projects/scripts/build_se_fragment_h5s.nf:140`, the one live SE
production pipeline. A reviewer applied this document's own prescription verbatim
and both of its named SE scenarios died at that `TypeError`.

Removing it safely would require moving `max_tlen` out of the shared call into the
same conditional-kwargs mechanism used for `se_max_fragment_length`, while still
delivering it to `tsv_to_fragments` (which genuinely uses it at `fragment.py:707`)
and to `bam_to_fragments`. That is a real design decision, not a cleanup, and it is
**out of scope here**. Leave the parameter alone and document it as dead-by-design.

**DECIDED (user, 2026-08-21): add `numpy>=1.24` to `pyproject.toml` dependencies.**
Previously flagged as optional; it is now in scope. The explicit check makes the
numpy version irrelevant for this specific field, but a floor closes the last path
by which any out-of-range assignment in this codebase could wrap instead of raise.
Note what the floor actually repairs: `conda-recipe/recipe.yaml` already pins
`numpy >=1.26` in both `host:` and `run:`, while `pyproject.toml` and
`environment.yml` pin nothing. That gap is the entire reason the uint16 failure
mode is environment-dependent — a hard crash in every deployed environment, and a
silent wrap only for someone who pip-installs into a numpy <1.24 env. The floor
makes the two agree.

No new worker-tuple element: `se_max_fragment_length` is already tuple element 9
(see "Implementation" above); the hard `65535` ceiling is a relocated module
constant, not a parameter.

### 3. Contigs with exactly one mapped alignment

`fragments_h5.py:1029-1033`:

```python
num_mapped = {
    x.contig: x.mapped // 2
    for x in bam_fp.get_index_statistics()
    if allowed_contigs is None or x.contig in allowed_contigs
}
```

`x.mapped` comes from the BAM index and counts mapped **alignment records** (each
mate of a pair counted separately). The `// 2` is a paired-end assumption applied
unconditionally.

The only consumer is the skip at `fragments_h5.py:1062-1065`:

```python
for bam_contig in contig_lengths.keys():
    # skip contigs with zero mapped reads (BAM only; TSV tabix only lists contigs with data)
    if not is_tsv_input and num_mapped[bam_contig] == 0:
        continue
```

Grep-confirmed: `num_mapped` is referenced at **exactly this one site**. It does
not feed chunking, work splitting, or progress (`total_bases` is computed
independently at `fragments_h5.py:1086`). The fix is fully localized.

**Replacement: drop the `// 2` and rename.**

```python
num_mapped_alignments = {
    x.contig: x.mapped
    for x in bam_fp.get_index_statistics()
    if allowed_contigs is None or x.contig in allowed_contigs
}
```

The predicate stays `== 0`, which is now what it always claimed to be: skip
contigs with literally zero mapped alignment records. Since the dict was only
ever used as a zero test, `// 2` bought nothing and cost a threshold artifact.
The rename is included because "num_mapped" was ambiguous between reads and
records and that ambiguity is the bug.

**Consequence today, measured** (contigA fixed at 5 reads, contigB varied, SE):

| contigB reads | output contigs | total fragments | contigB |
|---|---|---|---|
| 1 | `['contigA']` | 5 | **absent**, no warning on stderr |
| 2 | `['contigA','contigB']` | 7 | starts `[5000, 5100]` |
| 3 | `['contigA','contigB']` | 8 | starts `[5000, 5100, 5200]` |

The `continue` at `:1065` precedes chunk-arg construction (`:1062-1083`), so the
contig is never dispatched, never lands in `completed_results`, and the merge
step (`:1153-1186`) — which builds `f['data/<contig>']` groups strictly from
`completed_results` — creates no group for it. The contig is completely absent
and indistinguishable from one that was never in the BAM. Note also that `// 2`
only gates *inclusion*: once a contig is included, all of its reads' fragments
appear. This is a threshold artifact, not undercounting.

**Does this affect paired-end?** In principle yes — a PE contig with exactly one
mapped record is dropped by the same expression. In observable output, almost
certainly not: the PE filter (`fragment.py:368-384`) requires
`not align_.mate_is_unmapped` and `0 < tlen <= max_tlen`, so a lone mapped record
whose mate is unmapped is filtered anyway, and a record whose mate is on another
contig conventionally carries TLEN=0 and is filtered by `align_.tlen == 0`.
**Confidence: high on the filter logic (read directly from source); moderate on
the cross-contig TLEN=0 convention, which is an inference from SAM/samtools
behavior and was not measured here.** The PE case was not directly testable with
the synthetic reads used, because those carry no TLEN and PE therefore yields
zero fragments and raises `ValueError: No fragments were extracted...` for all
three counts.

**Correction: PE output does not observably change.** An earlier claim here — that
such a contig gets an empty `data/<contig>` group after this fix — is false:
`fragments_h5.py:858-859` returns `None` in place of a sub-h5 whenever a chunk
yields zero fragments, and the merge loop only appends to `completed_results`
when that path is not `None` (`:1111`, `:1125`), so a zero-fragment chunk
produces **no group at all** — byte-identical to being skipped. Measured:
contigB with 2 duplicate-flagged mapped records passes the gate, is dispatched,
yields 0 fragments, and is absent from `/data`. Item 3's PE story is therefore
**no observable change** — the accurate claim, stated plainly. Writing an empty
group instead would be a legitimate future change but is out of scope here.

**Add a log line.** The current silent skip is half the harm. Log **one
aggregated line**, not one per contig — a full GRCh38+ALT+decoy reference has
~3,366 sequences and most will be empty in a targeted BAM:

```
INFO skipping 3312 contigs with zero mapped alignments (e.g. chr1_KI270706v1_random, ...)
```

**Not a risk: the merge path never sees an empty array.** Whether
`:1153-1186` handles zero-fragment chunks is not a risk to verify — it is
unreachable, per the correction above: `:858-859` returns before a sub-h5 for a
zero-fragment chunk is ever created, so `numpy.concatenate` at `:1180` only ever
receives arrays from chunks with `ff > 0`. Same for the whole-file `ValueError:
No fragments were extracted...` guard at `:1129-1135` — a single empty contig
among otherwise non-empty ones contributes nothing to `completed_results` either
way and cannot trigger it. The real residual question is the one already
answered above: item 3 changes no PE output at all.

No new worker-tuple element — this is main-process only.

### 4. Remote (S3) input paths

`fragments_h5.py:963` applies `os.path.abspath(input_fname)` unconditionally, and
`:958` does the same for `fasta_filename` when one is given. Measured:
`os.path.abspath('s3://bucket/key.bam')` -> `'/<cwd>/s3:/bucket/key.bam'`;
`posixpath.normpath` collapses the `//` and prepends cwd.

Measured failure, `python -m pytest tests/test_s3_bam.py -v` from `cwd=/tmp`:

```
test_pysam_opens_s3_bam                            PASSED
test_build_fragments_h5_from_s3_bam                FAILED
test_build_fragments_h5_from_s3_bam_with_s3_fasta  FAILED
2 failed, 1 passed in 1.64s
```

```
FileNotFoundError: [Errno 2] Could not open alignment file: No such file or directory:
  '/tmp/s3:/fragmentomics.kariusdx.com/nboley/fragments_h5_test_data/small.chr6.bam'
```

Raised at `fragments_h5.py:1023` (`pysam.AlignmentFile(input_fname)` in the
**main** process, extracting contig lengths) — before any worker or chdir is
reached, immediately after the `:963` `abspath`. Not a credentials problem:
`pysam.AlignmentFile('s3://.../small.chr6.bam')` opens successfully with
`references == ('chr6',)` and `has_index() == True`, and
`test_pysam_opens_s3_bam` — which never calls `abspath` — passes.

**Why `abspath` is there, and why the fix must preserve it.**
`git log -S'abspath' -- src/` surfaces two commits from one debugging session:
`ee2c87d` ("Fix FASTA file path resolution bug causing OSError on AWS Batch",
adds `:958`) and `7b2cd99` ("Fix BAM file path resolution bug causing
FileNotFoundError on AWS Batch", adds `:963`), whose message reads: *"Same root
cause as the FASTA fix in ee2c87d: `build_sub_fragments_h5()` changes cwd via
`_temporary_working_directory()`, breaking relative paths."*

`_temporary_working_directory()` (`fragments_h5.py:892-904`, own docstring at
`:894-897`) `os.chdir()`s into a fresh tempdir and restores cwd on exit; it is
entered inside the per-chunk worker at `:744`. The docstring of the *caller*,
`build_sub_fragments_h5` (`:730-731`), says the temp-dir-per-worker trick exists
specifically *"so S3/remote BAM index/cache files do not overlap when multiple
workers run."* So the mechanism was built **for** S3, and the fix for the side
effect it caused is what now breaks S3.

**The key observation that makes the fix trivial:** `abspath` exists to make paths
cwd-independent before workers chdir. A remote URL is *already* cwd-independent.
So the correct rule is "absolutize local paths only" — it preserves the entire
reason `abspath` was added, with no loss.

**Fix.** Add to `fragments_h5.py`, next to the path handling:

```python
_REMOTE_URL_RE = re.compile(r"^[a-z][a-z0-9+.\-]*://", re.IGNORECASE)

def _is_remote_url(path):
    """True for anything htslib will treat as a URL (s3://, gs://, http(s)://, ftp://)."""
    return bool(_REMOTE_URL_RE.match(str(path)))

def _resolve_input_path(path):
    """Absolutize local paths so workers can chdir; leave remote URLs untouched."""
    if path is None or _is_remote_url(path):
        return path
    return os.path.abspath(path)
```

Use `_resolve_input_path` at `:958` and `:963`.

**Avoiding a duplicate predicate.** `_is_remote_url` currently exists only at
`main.py:66-68` and is not imported by `fragments_h5.py`. Given this repo's
documented history of duplicated logic drifting (the two hand-copied fragment
filters in item 1 are the same pathology), the definition moves to
`fragments_h5.py` and `main.py` imports it. **One definition, two call sites.**
A new `paths.py` module for a single four-line predicate would be
over-abstraction; rejected.

**`gs://`: yes, via scheme matching rather than a prefix list.** The current
`main.py` predicate matches `s3://`, `http://`, `https://` and misses `gs://`. A
generic scheme regex is correct for any scheme htslib supports now or gains
later, and removes the need to revisit the list. Note carefully what this does
and does not claim: it does **not** add `gs://` support. htslib may or may not be
built with GCS support in a given image. What it does is stop *corrupting* the
path, so the user gets htslib's real error rather than
`No such file or directory: /cwd/gs:/bucket/...`. Windows drive letters are not a
concern (POSIX-only project; `C:\` has no `//`).

**`main.py` needs no other change.** Its index handling (`:105-127`) and FASTA
validation (`:130-142`) already run on the un-mangled `args.input_file` /
`args.fasta` and already pass; the corruption happens later, inside
`build_fragments_h5`. `main.py` never downloads — it hands the raw URL to pysam,
and htslib streams S3 natively. Widening the predicate to `gs://` does change one
`main.py` behavior: a `gs://` BAM missing an index now takes the remote branch and
raises `SystemExit` advertising that the index must pre-exist, instead of trying
to `samtools index` a nonexistent local file. That is the desired behavior.

**Remote output is out of scope.** `ofname` is never `abspath`'d, never checked
against `_is_remote_url`, and is opened directly with `h5py.File(ofname, "x")` in
the main process. S3 output is not implemented and not advertised. No change.

No new worker-tuple element: `input_fname` and `fasta_filename` are already
elements 1 and 7 of the tuple at `:735`; only their values change.

### 5. Build provenance

The user's framing, verbatim: *"why can't we just save the command arguments in
the h5? Maybe the bam header too?"* The BAM header is already stored
(`fragments_h5.py:1149`, from `str(bam_fp.header)` at `:1024`). So this is only
the two missing pieces.

**Two new root attributes. Nothing else.**

| Attribute | Type | Written when | Content |
|---|---|---|---|
| `_build_argv` | `str` (JSON array of strings) | only when the caller supplies argv | `json.dumps(sys.argv)` as captured in `main.py` |
| `_build_version` | `str` | always, when determinable | `importlib.metadata.version("fragments-h5")` |

Written in the file-creation block alongside the existing five attributes
(`fragments_h5.py:1147-1151`). Attributes, not datasets, matching the house rule
visible there: datasets for array-like things, attributes for small
scalars/strings. Underscore prefix matches `_bam_header` / `_source_format` /
`_contig_lengths_str` and signals a forensic field.

**Size is not a constraint.** Measured twice independently: attribute strings of
1KB through 1,000,000 chars write and read back correctly as `str`, with file
size tracking content linearly. The old ~64KB compact-attribute limit does not
apply; HDF5 transparently promotes to dense storage. For reference, real
`str(header)` sizes measured: 5,622 bytes for a 194-contig BAM, 487 for
`small.chr6.bam`, 299 for `test_duplicates.bam` (~29 bytes per `@SQ` line). An
argv list is orders of magnitude smaller than the header already stored.

**How argv gets in — `sys.argv` is never read inside `build_fragments_h5`.**
This is a hard constraint: a library or test caller would record
`['pytest', '-v', ...]`, which is garbage that looks authoritative. Add a
keyword-only parameter to the signature at `fragments_h5.py:916-934`:

```python
def build_fragments_h5(input_fname, ofname, ..., min_mapq=None, *, build_argv=None):
```

`main.py` passes `build_argv=sys.argv` at its call site (`main.py:166-182`, which
already passes everything but the first two arguments by keyword). Defaults to
`None`, keyword-only, so no existing caller changes behavior and no library
caller can record argv by accident. All direct callers outside `main.py` are
tests and tooling — `tests/test_fragments_h5.py` (13 sites),
`tests/test_cli_validation.py` (7), `tests/test_s3_bam.py` (2),
`tests/specialized/compare_chunked_vs_unchunked.py` (2) — and
`pyproject.toml:33-34` declares exactly one console script,
`build-fragments-h5 = "fragments_h5.main:main"`. Nothing else needs updating.

**Library callers record `_build_version` but not `_build_argv`.** There is no
argv to record and inventing one would be a lie. The version still identifies the
code, which is the part that determines items 1-3.

**JSON, not `str()` + `eval()` — a deliberate departure from house precedent.**
`fragments_h5.py:293` reads `self.contig_lengths = eval(self._f.attrs["_contig_lengths_str"])`,
and there is no `json.dumps` anywhere in the writer. That precedent is unsafe:
`eval()` on content read from an h5 downloaded from S3 is arbitrary code
execution, and provenance — a field whose whole purpose is to carry
externally-supplied strings — is exactly where that lands. JSON is also readable
by `h5dump` and by non-Python tools. Changing `_contig_lengths_str` itself is out
of scope (it would break every existing file's reader); noted as a follow-up.

**Reading side.** Follow the `_source_format` precedent at
`fragments_h5.py:296`, which is the established backward-compatible pattern
(`.attrs.get(...)` with a default) versus the required reads at `:293-295`. In
`FragmentsH5.__init__` (`:288-304`):

```python
_argv = self._f.attrs.get("_build_argv")
self.build_argv = json.loads(_argv) if _argv is not None else None
self.build_version = self._f.attrs.get("_build_version")
```

Both are `None` for files built before this release. h5py attrs are an open dict
and there is no validator rejecting unknown keys, so adding attributes is safe
for old readers; new code only breaks old files if it does `attrs["new_key"]`
instead of `.get()`, which the above avoids. No schema-version attribute exists
and none is added.

**`.get()` is also required for a second, more immediate reason, not just old
files.** `build_fragments_h5` closes the file at `:1210` and reopens it as
`FragmentsH5(ofname, "r+")` at `:1214` to add `fragment_length_counts` — the
modified `__init__` above runs *inside the builder itself*, on a file that, for
a library caller (`build_argv=None`), has `_build_version` but no
`_build_argv`. `.attrs["_build_argv"]` there would `KeyError` on every
library-caller build. Provenance does survive the reopen: attrs are written at
`:1147-1151`, before `f.close()`, and `_add_fragment_length_counts` (`:708-720`)
only adds a dataset, never touching `.attrs`.

**Full local paths are recorded verbatim.** `input_file`, `output_frags_h5`,
`--fasta` and `--contig-name-map` are path-bearing and will embed directory
layout, and on shared systems a username. Judgment: **acceptable, record them.**
These are internal scientific data files, not user secrets; the paths are the
single most useful forensic field in the record; and any redaction scheme is
exactly the over-engineering that got the previous design rejected. This also
does not create a new exposure class: `_bam_header`, already stored today
(`fragments_h5.py:1149`), already carries the aligner's `@PG` header line, and
SAM `@PG` records include a `CL:` (command line) field that routinely embeds the
same kind of local input/reference paths. `_build_argv` widens an existing
exposure rather than introducing one. Stated explicitly here so the user can
overrule it — the alternative would be a one-line change to strip directories, at
the cost of most of the value.

**Known limitation: `--contig-name-map` records the path, not the contents.** If
the mapping file changes between two builds, argv alone will not detect it.
Recording a hash or embedding the contents was considered and rejected below.

**Version caveat that must be stated.** `importlib.metadata.version` reads
INSTALLED dist-info, not the working tree. Verified empirically in the dev env:
both `version("fragments-h5")` and `version("fragments_h5")` return `'2.11.0'`,
matching `pyproject.toml:7`. But an editable install with stale metadata reports a
stale version — this exact repo previously had dist-info reporting 2.9.1 while
the tree said 2.11.0. So `_build_version` pins the *installed distribution*,
which is exactly right for containers and can mislead in a dev checkout. Related:
the `Makefile` `docker-build` target now gates on tree cleanliness, the
`docker-build` counterpart of the `tag` target's `pyproject.toml` check: it
refuses to build when there are tracked changes (staged or unstaged) or any
untracked files. This closes the gap that previously let a container be built
from a dirty tree with a stale version.

If `importlib.metadata.PackageNotFoundError` is raised (uninstalled tree), omit
the attribute rather than writing a sentinel like `"unknown"`. Absent means
unknown; there is no second thing to disambiguate.

**No new worker-tuple element** — provenance is written in the main process at
file creation.

### Code revision is library-side, unlike argv

`_build_code_revision` is determined in `fragments_h5.py` via
`_resolve_build_code_revision()`, not at the CLI boundary. This is a deliberate
departure from `_build_argv`, which is passed in from `main.py`. The reasoning:
argv is a property of the *invocation*, known only to the CLI entry point, while
code revision is a property of the *installed package*, derivable from
`__file__`. Placing it in the library means the ~30 library and test callers of
`build_fragments_h5()` also get a revision stamp, rather than only CLI builds.

The resolver returns a self-labeling string with one of four prefixes (`git:`,
`baked:`, `dist:`, `dist-editable:`) or `None`. The prefix names both the value
and the oracle that produced it, so consumers know what guarantee applies.

### Worker-args tuple: nothing is added

Explicit answer to the standing question. The 17-element positional tuple
unpacked at `fragments_h5.py:735` is **unchanged by all five items**:

| Item | Needs a tuple element? | Why not |
|---|---|---|
| 1 secondary | No | unconditional; no parameter exists |
| 2 SE bound | No | `se_max_fragment_length` is already element 9; bound is a constant |
| 3 `num_mapped` | No | main process only |
| 4 remote paths | No | mutates existing elements 1 and 7 |
| 5 provenance | No | written at file creation, main process |

This was a design constraint, not a coincidence: it is one of the reasons item 1
is unconditional rather than opt-out. The tuple remains the acknowledged root
cause of this bug class and has its own deferred design; it is **out of scope
here** and is not touched.

## Rebuild guidance

Position per item. None of these mandate a blanket rebuild.

**Item 1 — no rebuild.** The sampled census (0 secondary alignments in ~61,000
sampled reads across fixtures and two production BAMs) plus the aligner-flag
inference (`bwa-mem2`/`bowtie2` without `-a`/`-k`) makes it very likely that no
existing h5 contains a secondary-derived fragment. This is a **census and an
inference, not a proof**, and a file cannot be checked from its own contents.
The check is on the source BAM, not the h5:

```
samtools view -c -f 256 <input.bam>
```

Nonzero means that BAM's h5 should be rebuilt. Zero means it is provably
unaffected.

**Item 2 — no silently-bad files exist in any deployed environment.** On conda
(`numpy>=1.26`) and in the Docker image, an over-long SE span raises
`OverflowError` — the failure mode is *absent output*, not corrupt output. Only a
pip install into a pre-1.24 numpy environment could have produced wrapped
lengths. There is no reliable way to detect wraparound from the file: a wrapped
70020 becomes 4484, which is indistinguishable from a genuine 4484 in
`fragment_length_counts`. So the check is on the **build environment**: if a
single-end h5 was built without `--se-max-fragment-length` by a pip install whose
numpy version you cannot rule out as pre-1.24, rebuild it. Of the two nominal SE
call sites, `biomarker-projects/scripts/build_se_fragment_h5s.nf:139-143` passes
`--se-max-fragment-length 120` (unaffected either way, including under the
corrected implementation in section 2); the other,
`biomarker-projects/projects/lupus/scripts/build_frag_h5s.py:42`, passes none —
but also passes flags `main.py`'s parser does not define (see section 2), so it
is not a currently-working call site regardless.

**Item 3 — this is the one that could genuinely have dropped data, in SE builds.**
Check without rebuilding, by comparing the h5's contigs against the source BAM
index:

```python
import pysam, h5py
h5_contigs = set(h5py.File(path)["data"].keys())
bam_contigs = {x.contig for x in pysam.AlignmentFile(bam).get_index_statistics() if x.mapped > 0}
missing = bam_contigs - h5_contigs   # non-empty => affected
```

For a `--contigs`-restricted build, intersect with the requested list first. If
`missing` is non-empty, the dropped contigs each had exactly one mapped record;
decide from the biology whether one read on that contig matters. Paired-end
builds are very likely unaffected in output (see item 3 above, moderate
confidence).

**Nothing reaches production until container pins are bumped.** Observed pins:
`build_fragment_h5s.config` and `catlas_demux_fragments.config` -> 2.9.1;
`build_se_fragment_h5s.config` -> 2.10.1;
`biomarker-pipeline/containers/container-versions.config:48` -> 2.8.1;
`biomarker-pipeline/tests/specialized/test_build_fragments_h5_batch.nf:33`
hardcodes 2.8.0, bypassing that shared config. **No consumer references 2.11.0**,
so no production consumer has adopted even the currently-released SE/MAPQ flags.
Bumping pins is a separate, uncoordinated, per-repo action and is out of scope.

**Format compatibility.** No format version bump. Old files open unchanged; the
two new attributes are read with `.get()` and default to `None`. New files open in
old readers, because h5py attrs are an open dict and nothing validates the key
set.

## Known limitation: derived h5s drop provenance

Two verified downstream consumers construct a fragments h5 by writing a
**hardcoded 5-attribute allowlist** into a fresh `h5py.File`, dropping everything
else:

1. `biomarker-projects/tests/tf_footprinting/fixtures/create_test_h5s.py:28-32` —
   a test-fixture generator (also blanks `_bam_header` unconditionally).
2. `biomarker-projects/tf_binding_site_classification/projects/xeno_cardiac/scripts/liftover_porcine_fragments.py:193-197`
   (`:187-188` are the reads from the source h5) — a *different* 5-attribute
   allowlist, and a **real analysis path**, not a
   fixture. It passes `_bam_header` and `_contig_lengths_str` through and
   overwrites `_source_format` with `'pyliftover_from_porcine_h5'`.

Both will silently drop `_build_argv` and `_build_version`. Repo-scoped searches
found no other h5-creating/copying code or `.attrs` reader in
`biomarker-projects`, `biomarker-pipeline`, or `fragmentomics_tools`.

So provenance survives only along paths that copy attributes generically, and a
derived h5 may carry no provenance while looking like a normal file.

**Is absence-of-provenance distinguishable from built-before-this-change?**
Partially, and the design does not pretend otherwise:

- `_build_version` present -> built by >= 2.12.0 from an installed distribution.
- `_build_version` absent -> pre-2.12.0, **or** built from an uninstalled tree,
  **or** derived by one of the two allowlist scripts above.

Those three are not distinguishable from the file alone. Adding a sentinel or a
format version to disambiguate them was considered and rejected (below): it buys
a distinction between three cases that all reduce to the same user action —
consult the pipeline that produced the file. **Cross-repo follow-up: extend the
two allowlists.** Per scope constraints, no changes to those repos are proposed
here.

## Files to modify

**Line numbers below are addresses on the pristine tree and DRIFT as you edit.**
The very first prescribed edit — removing the `MAX_FRAG_LENGTH` definition at
`fragments_h5.py:174-176` — deletes three lines, so every address after it shifts
by three. A reviewer observed the traceback line move from `:1023` to `:1021` after
patching. Anchor on the quoted code, not the number, from step two onward.

`src/fragments_h5/fragment.py`
- `:368-384` add `align_.is_secondary` to the PE filter
- `:573-579` add `align.is_secondary` to the SE filter
- `:554` **leave `max_tlen=1000` exactly as it is.** It is dead, and removing it
  breaks every single-end build — see "Dead parameter" in section 2.
- move `MAX_FRAG_LENGTH = 65535` here from `fragments_h5.py:176` (resolves
  the circular import in section 2)
- `:548-560` add `se_max_fragment_length=None` to the signature
- after `:581-582` raise `ValueError` when `se_max_fragment_length is None and
  frag_stop - frag_start > MAX_FRAG_LENGTH`

`src/fragments_h5/fragments_h5.py`
- `:140-153` add `import json`, `import re`, `import importlib.metadata` —
  none of the three is currently imported, and section 4's `_REMOTE_URL_RE` and
  section 5's `json.dumps`/`json.loads`/`version()` need them
- `:153` also import (re-export) `MAX_FRAG_LENGTH` from `fragments_h5.fragment`;
  remove the `:176` definition
- add `_is_remote_url` and `_resolve_input_path` helpers
- `:958`, `:963` use `_resolve_input_path` instead of `os.path.abspath`
- `:1029-1033` drop `// 2`, rename to `num_mapped_alignments`
- `:1062-1065` predicate unchanged; add one aggregated log line for skipped contigs
- `:772-777` pass `se_max_fragment_length=` to the call only when `single_end`
  is true (`bam_to_fragments` has no such parameter and no `**kwargs`).

  **This is not a one-line edit, and it is the only non-obvious mechanical step in
  the change.** `:737-742` assigns a *function object*, and `:772-777` is a single
  shared call expression — there is no per-branch `if single_end: ... else: ...`
  call to add a kwarg to. Use a conditional kwargs dict and a `**` splat.
  The mechanism below was implemented in a throwaway tree and executed against
  synthetic BAMs across all three dispatch branches — SE BAM receives the kwarg,
  PE BAM correctly does not, and a direct `bam_to_fragments(..., se_max_fragment_length=…)`
  raises `TypeError` as expected:

      _se_kwargs = {"se_max_fragment_length": se_max_fragment_length} if single_end else {}
      for fragment in input_to_fragments(
          input_fname, chrom=bam_contig, start=chunk_start, stop=chunk_stop,
          max_tlen=MAX_FRAG_LENGTH, fasta_file=fasta_file,
          include_duplicates=include_duplicates, min_mapq=_min_mapq, **_se_kwargs
      ):

  Take the exact argument list from the existing call at `:772-777` rather than
  from this snippet — the point being illustrated is the `_se_kwargs` splat, and
  the surrounding arguments are reproduced only so the shape is clear.

  **Dispatch precedence:** `is_fragment_file(...)` is tested BEFORE `single_end`
  (`:737-742`), so a single-end *TSV* build routes to `tsv_to_fragments` and will
  receive `se_max_fragment_length` in this splat. That is harmless —
  `tsv_to_fragments` absorbs unknown keywords via `**kwargs` (`fragment.py:645`) —
  but it is worth knowing, since this document elsewhere stresses that
  `bam_to_fragments` has no `**kwargs`, which invites the wrong conclusion about
  the third branch.
- `:916-934` add keyword-only `build_argv=None`
- `:1147-1151` write `_build_argv` (when provided) and `_build_version`
- `:288-304` read both back via `.attrs.get()`

`src/fragments_h5/main.py`
- `:66-68` delete the local `_is_remote_url`, import it from `fragments_h5`
- `:166-182` pass `build_argv=sys.argv`

`pyproject.toml` — add the `numpy>=1.24` floor (now decided, in scope).
**Do NOT bump the version.** 2.12.0 is bumped and tagged together from `main`
after this lands — see Open Question 6 for why splitting them is how this repo
produced three mislabeled tags.

Not modified: the worker-args tuple at `:735`; `_contig_lengths_str` encoding;
output path handling; any other repo.

## Testing

Synthetic BAMs are built with `pysam.AlignmentFile(..., 'wb', header=...)` plus
`pysam.index` — the same construction used to produce every measured result cited
in this document, so these tests are known-buildable. There is no `conftest.py`
in `tests/`; keep fixtures local to each module, matching the existing layout.

**`tests/test_fragment_selection.py` (new)**

- `test_secondary_alignments_excluded_single_end` — 5 primary reads plus 3
  secondary alignments at 5000/5100/5200. Assert 5 fragments and that no fragment
  starts at 5000/5100/5200. *(Measured pre-change: 8 fragments.)*
- `test_secondary_alignments_excluded_paired_end` — 6 proper pairs plus 1
  secondary record carrying a positive in-range TLEN at 70000. Assert 6 fragments
  and no start at 70000. *(Measured pre-change: 7 fragments.)*
- `test_se_oversized_span_raises` — one read with CIGAR `10M70000N10M`
  (reference span 70020) plus a 20bp companion, `se_max_fragment_length=None`.
  Assert `ValueError` whose message contains the contig, the start position, the
  read name and `70020`. *(Measured pre-change: `OverflowError` from
  `fragments_h5.py:801`, naming none of those.)*
- `test_se_span_boundary` — spans 65535 and 65536. Assert 65535 builds
  successfully and 65536 raises, pinning `>` rather than `>=`. *(Measured
  pre-change: 65535 OK, 65536 `OverflowError` — same boundary, different
  exception and message.)*
- `test_se_max_fragment_length_still_skips` — the 70020-span read with
  `se_max_fragment_length=1000`. Assert no exception and that only the 20bp
  fragment survives. This guards the SE-raises / SE-skips reconciliation.
  *(Measured: this is today's behavior and must not regress.)*
- `test_contig_with_single_mapped_read_included` — contigA 5 reads, contigB 1
  read, single-end. Assert both contigs present, 6 fragments total, contigB start
  5000. *(Measured pre-change: contigB absent, 5 fragments, no warning.)*
- `test_contig_with_zero_mapped_reads_skipped` — contigA 5 reads, contigB 0.
  Assert contigB absent and that the aggregated skip message is logged (via
  `caplog`; `pytest.ini` sets `log_level=INFO`).

**`tests/test_s3_bam.py` (existing)**

- `test_build_fragments_h5_from_s3_bam` and
  `test_build_fragments_h5_from_s3_bam_with_s3_fasta`: currently FAILED, expected
  to PASS with no test-side changes.
- `test_resolve_input_path` (new, pure unit, no network) — assert
  `_resolve_input_path` leaves `s3://`, `gs://`, `https://` unchanged and makes a
  relative local path absolute. Run it with `os.chdir` into a tempdir to prove
  cwd-independence, which is the property `abspath` was added for.

**`tests/test_build_provenance.py` (new)**

- `test_provenance_absent_for_library_caller` — direct `build_fragments_h5(...)`.
  Assert `_build_argv` not in `f.attrs` and `_build_version` equals
  `importlib.metadata.version("fragments-h5")`.
- `test_provenance_recorded_when_argv_passed` — pass an explicit
  `build_argv=['build-fragments-h5', 'in.bam', 'out.h5', '--single-end']`.
  Assert `json.loads(f.attrs['_build_argv'])` round-trips exactly, including
  paths.
- `test_provenance_accessors` — `FragmentsH5(path).build_argv` /
  `.build_version` return the parsed list and the version string.
- `test_file_without_provenance_opens` — build a file, delete both attributes
  with raw h5py, reopen with `FragmentsH5`, assert both accessors return `None`.
  This is the backward-compatibility guarantee for every file currently under
  `s3://karius-biomarker-data-assets/`.

**Expected suite baseline.** Measured today from `cwd=/tmp`:

```
2 failed, 47 passed, 6 errors in 273.44s
```

The 2 failures are item 4. The 6 errors are all
`subprocess.CalledProcessError` from `/bin/sh: 1: build-fragments-h5: not found`
— the console script is not on PATH in that environment. They are an environment
artifact, not a code defect.

**Expected after this change: `0 failed, 47 + N passed, 6 errors`**, where N is
the count of new tests above. The 2 S3 failures become passes; **the 6 errors
remain** and must not be treated as a regression. No existing test asserts the
full attribute set (`attrs.keys()`, `set(f.attrs)`, `len(f.attrs)` all return zero
grep hits), so adding attributes breaks nothing; `test_docker_build.py:92-99`
already compares attributes defensively via `.attrs.get(attr)`.

## Alternatives rejected

**Secondary-alignment opt-out flag.** Would need an 18th worker-tuple element and
would make the behavior undeterminable from `_build_version`, undermining
section 5. No use case was identified for counting secondary alignments as
fragments.

**Clamping SE length to 65535.** Explicitly excluded by the user, and correctly:
a clamped value is a plausible-looking wrong number that flows into
`fragment_length_counts` and every downstream analysis with no signal.

**Widening the length dtype past `uint16`.** Doubles the size of the largest
array in the file to accommodate spans that are not fragment lengths in any
meaningful sense. `max_fragment_length` is also a public attribute
(`fragments_h5.py:1148`) that consumers read.

**Making SE silently skip like PE/TSV.** Rejected because the skip behavior
already exists behind `--se-max-fragment-length` and because, for SE, an
over-long span signals a parameter error rather than a normal biological outlier.
See the reconciliation in section 2.

**Downloading remote files locally instead of fixing `abspath`.** `main.py`
deliberately never downloads; htslib streams S3 natively, and
`_temporary_working_directory()` already exists specifically to keep per-worker
remote index caches from colliding, per `build_sub_fragments_h5`'s docstring
(`fragments_h5.py:730-731`). Downloading would
discard a working design to avoid a four-line predicate.

**A `paths.py` module for the remote-URL predicate.** Over-abstraction for one
function; `fragments_h5.py` with an import in `main.py` gives the same
single-definition guarantee.

**A prefix list (`s3://`, `http://`, `https://`, `gs://`) instead of a scheme
regex.** The existing list already drifted — it is missing `gs://`. A scheme
match never needs revisiting.

**`requested`/`effective`/`neutralized` provenance schema.** The previous
design's proposal. Explicitly rejected as over-engineered and not resurrected in
any form.

**Argv-only recording, as rejected by the previous design.** That document
rejected argv-only on the grounds that it "records what was requested, never what
was effective." Confronting it directly: that objection was correct *for a design
in which behaviors are configurable and can be silently neutralized*. It is not
correct for this bundle, because items 1-3 are all unconditional. Secondary
exclusion has no flag; the SE bound raises rather than varying; the `num_mapped`
predicate is fixed. Every behavior argv does not capture is a constant of the
code, and `_build_version` pins the code. **argv plus version is fully
determining for CLI builds** — which is precisely why bundling 1-3 with 5 works,
and why making 1-3 configurable would have broken it.

**~35 flat typed attributes, one per parameter.** No clean `None` encoding for
the 12 nullable keys.

**`str(dict)` + `eval()`, matching `_contig_lengths_str`.** Arbitrary code
execution from content read off S3. JSON instead; see section 5.

**A sidecar `.provenance.json`.** Separable from the file it describes; loses
self-description, which is the entire point.

**Recording `contig_name_map` as a boolean or a content hash.** Destroys the only
useful information. The path is recorded verbatim as part of argv; the residual
gap (contents can change under a stable path) is documented rather than papered
over with a hash that tells the user nothing actionable.

**A new HDF5 `/provenance` group.** Two small strings do not need a group; every
other scalar metadata item in this file is a root attribute.

**A build timestamp.** The user's framing was deliberately minimal. Filesystem
and S3 object metadata already carry mtime, and a timestamp adds no
reproducibility information beyond `_build_version`.

**Bumping a format version and refusing old files.** Would break everything under
`s3://karius-biomarker-data-assets/`.

**Adding a `format_version` attribute, or repurposing `max_fragment_length`.**
No reader needs to branch on version: the two new attributes are additive and
read with `.get()`.

## Self-assessment

**Grade: A-.**

**Review history — the most useful calibration signal here.** Independently
reviewed once before this grade and graded **B-**. That review found two
critical defects, both fixed: item 2's prescribed raise location contradicted
the document's own reconciliation and its own named test, and — measured against
a real downstream `.nf` pipeline — would have turned a silent skip into a
production regression; fixed by threading `se_max_fragment_length` into
`single_end_bam_to_fragments`, guarded on `is None`. Separately, the prescribed
implementation could not be written as specified — importing `MAX_FRAG_LENGTH`
into `fragment.py` is a circular import — fixed by relocating the constant to
`fragment.py` and re-exporting it. Two high-severity findings on item 3's PE
story (a false claimed benefit; a test covering a structurally unreachable path)
and four medium findings (an overstated "identical constant" claim, an
unaddressed reopen-safety reason for `.get()`, missing imports, a
self-contradicting size number) were also fixed.

This is a design document, not a verified implementation. Every output block,
table of measurements, traceback, and test count quoted above is reproduced
verbatim from an actual run reported in the research findings. Where something is
an inference or a sample, it is labeled inline.

**A previous revision of this paragraph claimed that "every correction here was
re-verified line-by-line against `fragment.py` and `fragments_h5.py`." That claim
was false, and it should not be reinstated.** A third round-2 Critical
(the `max_tlen` removal, now retracted in section 2) survived that supposed
re-verification, and the fact refuting it was already stated in section 2 of this
same document, hundreds of lines above the prescription that contradicted it. Reading a file is not verifying a prescription against it. The round-2
reviewer caught it by *patching a throwaway copy and running it* — which is the
only method that has reliably worked on this document across three rounds.
Treat any assurance of thoroughness here as weaker evidence than one executed test.

**Confident about:**

- Item 1's mechanism. The PE path does no qname mate matching; a secondary record
  with positive in-range TLEN emits its own fragment. Measured 5->8 (SE) and
  6->7 (PE).
- Item 2's actual failure mode. It raises `OverflowError` on every deployed
  environment; "silently wraps" is wrong and is not repeated here. The `>` vs
  `>=` boundary is measured, not reasoned (65535 OK, 65536 raises).
- Item 3's blast radius. `num_mapped` has exactly one consumer; the contig is
  fully absent, not empty; measured at 1/2/3 reads.
- Item 4's root cause and the fix's safety. The bug fires at
  `fragments_h5.py:1023` in the main process, before any chdir, and remote URLs
  are already cwd-independent so nothing `abspath` was added for is lost.
- Item 5's compatibility. No test asserts the attribute set; `.attrs.get()` is
  the established precedent at `:296` and is load-bearing for the new-build
  reopen at `:1214` too, not just old files; attribute size is not a constraint
  (measured up to 1,000,000 chars, ~1.1MB on disk).
- The worker tuple genuinely does not need a new element.
- Item 3's PE story is a genuine no-op, traced to the `ff == 0` early return at
  `:858-859`, which also makes the previously-flagged empty-array merge risk
  structurally unreachable.

**Weakest points:**

1. **The SE-raises / PE-skips reconciliation is the load-bearing judgment call
   and it is a judgment, not a fact.** The argument — that TLEN is a routine
   aligner quantity while an SE reference span above 65535 signals a parameter
   error — is defensible but a reasonable reviewer could prefer symmetry
   (skip everywhere, warn once). If the user prefers that, section 2 changes but
   nothing else in the document does.
2. **The PE analysis in item 3 is partly inference.** "A cross-contig mate
   carries TLEN=0" is SAM/samtools convention, not something measured in this
   round. If it is wrong, the PE case has an observable output change, not just
   the (now-corrected) claim that it produces an empty group.
3. **RETRACTED — this entry named the wrong risk, and that is instructive.** It
   previously read: *"Removing the dead `max_tlen=1000` parameter assumes no
   external callers of `single_end_bam_to_fragments`."* The removal is now dropped
   entirely (section 2), because the hazard was never external callers — it was the
   **internal** shared call site at `fragments_h5.py:772-777`, described at line 153
   of this document. Naming a plausible-but-wrong risk is worse than naming none:
   it signals the question was considered and closed, so a reader stops looking.
   Every prescription in this document that has failed did so in a section that
   read as already-checked.
4. **Item 1's no-op claim rests on a sampled census** (~61,000 reads) plus an
   aligner-flag inference. It is labeled as such everywhere it appears, but it is
   still the shakiest empirical claim in the rebuild guidance.
5. **The `gs://` widening changes `main.py` index-validation behavior** as a side
   effect. Believed to be an improvement, argued as such, but it is a behavior
   change that was not requested.
6. **Section 2's original implementation placement was wrong in an earlier draft
   of this document** (see review history above), and the corrected version has
   not itself been executed, only re-checked line-by-line against source. The
   same caution that caught the first two defects should apply to the
   correction.

**Open questions for the user:**

1. **SE raise vs. skip** — confirm the asymmetry in section 2, or ask for
   symmetry with PE/TSV.
2. **Full local paths in `_build_argv`** — recorded verbatim, including
   `/home/<user>/...`, in files that may be published to
   `s3://karius-biomarker-data-assets/`. Acceptable, or strip directories?
3. ~~**The optional `numpy>=1.24` floor**~~ — **CLOSED: IN.** Decided by the user
   2026-08-21. See section 2.
4. ~~**Removing the dead `max_tlen` parameter**~~ — **CLOSED, answer is "leave it
   alone."** Not a judgement call: removing it raises `TypeError` on every
   single-end build, proven by execution (section 2).
5. ~~**Widening the remote predicate to any URL scheme (including `gs://`)**~~ —
   **CLOSED: YES, and it MUST ship with a test.** Decided by the user 2026-08-21,
   who was explicit that this one needs test coverage.

   **The test that matters is not the obvious one.** A test asserting that
   `gs://bucket/x.bam` is classified as remote will pass trivially and proves
   almost nothing. The reason this was an open question is that widening the
   predicate changes `main.py`'s **index-validation branch** as a side effect: a
   `gs://` input newly takes the remote path at `:105-127` and therefore skips
   local index handling. That behavioral side effect is what needs pinning, so a
   future change cannot silently revert it. Test the branch, not the predicate.
6. **Release number 2.12.0** — confirm, and confirm that bumping the four
   downstream container pins (2.8.0 / 2.8.1 / 2.9.1 / 2.10.1) is out of scope
   here.

   **CLOSED 2026-08-21: 2.12.0, but the bump happens AFTER this work lands, not as
   part of it.** Remove the version bump from "Files to modify" — the implementer
   must NOT touch `pyproject.toml`'s version. Bump and tag together, from `main`,
   once merged.

   This is not bureaucratic sequencing. This repo has produced mislabeled tags
   three times — `v2.5.0`→2.3.0, `v2.6.0`→2.5.0, `v2.10.1`→2.10.0 — and the
   mechanism each time was a version living in one commit while the tag labelled
   another (the `Makefile` `tag:` target reads VERSION from the WORKING TREE while
   tagging HEAD). A version bumped during implementation and tagged later is that
   exact setup. Verify the tag afterwards by reading pyproject **from the tagged
   commit**: `git show v2.12.0:pyproject.toml | grep '^version'`.
   Downstream container pins remain out of scope, confirmed.
7. **Item 3 rebuild scope** — run the contig-comparison check against existing
   single-end h5s, or accept the exposure?
