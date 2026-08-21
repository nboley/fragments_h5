# Build Provenance Metadata for fragments h5 files

**Status:** proposal (not implemented). **Target:** `fragments_h5` @ `main` / `v2.11.0` (`aa753c7`).

A fragments h5 records nothing about how it was built. Several build-time filters are *destructive* — filtered fragments never enter the file and no reader can recover them — yet two h5s built from the same BAM with different filter settings carry byte-identical metadata. This document proposes recording the effective build parameters, the builder version, and the CLI invocation as a single JSON-encoded root attribute `_build_provenance`, plus a `has_build_provenance` reader property so that *absence* is a first-class, detectable state rather than silently masquerading as defaults. The design follows the existing `_source_format` backward-compatibility precedent (`fragments_h5.py:296`), requires no change to the multiprocessing worker interface, and deliberately fixes one thing it inherits: it uses `json.loads` rather than the `eval()` the `_contig_lengths_str` attribute uses today (`fragments_h5.py:293`).

**Scope limit, stated up front:** this design is *prospective only*. It makes files built by the release that ships it self-describing. It does nothing for the files that already exist — for those, the answer to "what filters were applied?" is permanently "unknown", and §8 is the only available recourse. See §1.3.

---

## 1. Problem statement

### 1.1 Build-time filters are destructive and permanent

There are two categories of filtering in this stack, and they are not symmetric:

- **Read-time filters** are applied by the consumer at analysis time and can be varied freely. Example: `fragmentomics_tools`' `RegionFragmentArray.from_fragments_h5` calls `fetch_array` and then applies its *own* client-side MAPQ mask (`fragment_array/fragment_array.py:1754-1775`). A different analysis can choose a different threshold from the same file.
- **Build-time filters** are applied in `build-fragments-h5` and drop reads before they are ever written. There is no undo. A file built with `--min-mapq 30` cannot be turned back into one built with `--min-mapq 0`; the reads are gone.

A destructive transform that is *also* unrecorded is permanent and invisible. That combination is the problem.

### 1.2 Nothing in the file says how it was made

The complete set of root attributes written today is five, all at `fragments_h5.py:1146-1151`:

```python
f.attrs["index_block_size"] = INDEX_BLOCK_SIZE          # int 5000     (fragments_h5.py:169)
f.attrs["max_fragment_length"] = MAX_FRAG_LENGTH        # int 65535    (fragments_h5.py:176) — hardcoded
f.attrs["_bam_header"] = bam_header_str                 # str(bam_fp.header), "" for TSV
f.attrs["_source_format"] = "TSV" if is_tsv_input else "BAM"
f.attrs["_contig_lengths_str"] = contig_lengths_str      # str(dict), read back with eval()
```

Of these, exactly one — `_source_format` — carries any information about the build. `max_fragment_length` looks like it does but does not: it is the module constant, written unconditionally, unrelated to any CLI flag (see §7.3). There is **no schema or format version attribute anywhere** in the file and no version check on read (`fragments_h5.py:266-304`).

`_bam_header` may contain an `@PG` line from the aligner, which describes how the BAM was made. It says nothing about how the h5 was made.

### 1.3 Goal, and the scope it is achievable within

**Notation: the first release carrying the record is written `2.12.0` throughout this document** — in the §4.1.4 error string, in both worked examples' `builder.version`, and in the table below. `2.12.0` is the next minor after the current `v2.11.0` and is therefore the presumptive target, but **the number is provisional and must be set to the actual release version at implementation time**, in all of those places at once. It is bound to a concrete value here only so that the error message and the examples are not written against an unresolved placeholder.

**Goal (scoped):** *For files built by fragments-h5 2.12.0 or later*, a consumer should be able to determine what was done to the data before it reached them by inspecting the file, without locating the pipeline invocation that produced it. Filtered and unfiltered files must be distinguishable by inspection.

**What files built before 2.12.0 get: nothing.** This is not a limitation that better design removes — the information was never written and cannot be reconstructed (§10, "Retro-labeling existing S3 files"). For the existing corpus in `s3://karius-biomarker-data-assets/`, this design delivers exactly one thing: `has_build_provenance` returns `False`, which converts "silently assumed to be defaults" into "explicitly unknown". That is a real improvement over the status quo but it is not the goal above, and this document should not be read as claiming otherwise. The only tool that says anything about older files is the distributional forensics in §8, and that is one-directional (it can disprove a filter, not confirm one).

Concretely: the two blocks below both hold, and both matter.

| file built by | `has_build_provenance` | what a consumer can conclude |
|---|---|---|
| >= 2.12.0 (this design shipped) | `True` | the full `requested` / `effective` / `neutralized` record, per §4 |
| < 2.12.0 | `False` | nothing about filters. Not "defaults were used". Nothing. §8 forensics may *disprove* specific filters. |

### 1.4 The motivating incident — corrected

The incident that prompted this work was originally described as: *`build_se_fragment_h5s.nf` passed `--se-max-fragment-length` to a container that did not accept it, failing silently under `errorStrategy 'ignore'`, so h5s exist that were supposed to carry the cap and do not.*

**That description is wrong in its conclusion, and the correction matters for the design.** Verified facts:

- `--se-max-fragment-length` and `--min-mapq` were both introduced in the *same* commit, `f45b6f6` ("Expose `--se-max-fragment-length` and `--min-mapq` via the CLI; release 2.11.0"), `v2.11.0` → `aa753c7` (2026-08-21). Note that `aa753c7` is the **merge commit** ("Merge resolve-uncommitted-cli-flags: ...", parents `c667fd7` and `f4f5084`) that brought that branch onto `main`, *not* `f45b6f6` itself; `f45b6f6` is an ancestor of the tag, not the tagged object. This matters only for reproducing the checkout — `git show v2.11.0:...` resolves through the merge and gives the released source either way. The previous tag `v2.10.1` → `dbed0ae` (2026-06-08); `git show v2.10.1:src/fragments_h5/main.py` contains neither flag and `git show v2.10.1:src/fragments_h5/fragments_h5.py` contains no `se_max` or `min_mapq` identifier at all.
- The pipeline's `build_h5` process (`/home/nathanboley/src/biomarker-projects/scripts/build_se_fragment_h5s.nf:139-147`) invokes `build-fragments-h5 --single-end --verbose --se-max-fragment-length 120 --num-processes N --include-duplicates --min-mapq 30 --fasta F --contig-name-map M --contigs ... -- BAM OUT`, with `publishDir` to `s3://karius-biomarker-data-assets/projects/ibd/frag_h5s_se/` (`.nf:122`) and `output: path "${result_id}.fragments.h5"` (`.nf:131`).
- The container is pinned to `ghcr.io/nboley/fragments-h5:2.10.1` (`build_se_fragment_h5s.config:29,69`). A stale comment at `.nf:113` says "v2.10.0".
- **`errorStrategy` is profile-dependent, and only one profile ignores.** `.nf:123` declares `errorStrategy 'ignore'` in the script, but `build_se_fragment_h5s.config` sets `process.errorStrategy` in *both* profiles, and Nextflow config takes precedence over in-script directives: `standard` sets `errorStrategy = 'terminate'` (`.config:21`), `remote` sets a closure that retries up to `task.attempt <= 2` and only then returns `'ignore'` (`.config:58-61`, `maxRetries = 2` at `:62`). So under `-profile standard` the whole run **terminates loudly** on the first failed `build_h5`; only under `-profile remote` is the failure absorbed per-sample. An earlier version of this document cited `.nf:123` as operative — that was wrong, and it matters directly to the claim "the failure is silent".
- **Exactly one commit has ever touched that `.nf`: `df0aff3`.** `git log --oneline -- scripts/build_se_fragment_h5s.nf` returns `df0aff3` and nothing else. `8b0879b` ("pin the fragmentomics-tools container instead of tracking :latest") and `9d363c2` ("repin containers 1.2.0 -> 1.3.0...") touch `build_se_fragment_h5s.config` only, and repin the **fragmentomics-tools** container, not fragments-h5. The fragments-h5 pin is `2.10.1` at all three revisions. An earlier version of this document asserted all three revisions touched the `.nf`; that came from a `git log` run against both paths at once and is incorrect. The conclusion below does not depend on it — a single revision that has always passed both flags supports it equally.

**Corrected failure mode: no path through the committed pipeline can produce an h5 under container 2.10.1.** Four independent reasons, any one of which is sufficient:

1. `v2.10.1`'s `main.py` has no `--se-max-fragment-length` and no `--min-mapq`, and does not call `parse_known_args`, so argparse exits 2 on the first unrecognized option — before anything else in `main()` runs.
2. Even if argparse were bypassed, `main.py:98-103` exits 1 if the output path already exists, so no partial output is ever overwritten in place.
3. The output file is created by `h5py.File(ofname, "x")` at `fragments_h5.py:1139`, which runs **after** every worker has finished; workers write only into a `tempfile.TemporaryDirectory` (`:1059`) at paths under `tmp_dir_name` (`:862`). An abort at any earlier point leaves no file at `ofname`.
4. `output: path "${result_id}.fragments.h5"` (`.nf:131`) makes Nextflow fail the task if the file is missing, regardless of exit code, so nothing reaches `publishDir`.

The failure is *absence* of output, not silently-uncapped output. Under `-profile standard` that absence is loud (run terminates); under `-profile remote` it is quiet (task ignored after 2 retries) and only shows up as a short sample list downstream.

This strengthens the case for provenance but shifts the argument:

1. The S3 prefix can accumulate h5s from heterogeneous sources — hand-run containers, uncommitted working copies of the `.nf`, reruns after the container is repinned — with nothing in the files to tell them apart.
2. A build that produced nothing and a build that was never attempted are equally invisible downstream.
3. Once 2.11.0 is deployed, newly-capped h5s will land in the *same* prefix alongside whatever is already there.

**Open question we cannot close from git:** an uncommitted local working copy of the `.nf`, or a manual container invocation, may have produced uncapped h5s in that prefix. Git history rules out the *committed* pipeline as a source of uncapped output; it says nothing about anything run by hand. This is stated as unresolved rather than asserted either way. Resolving it requires forensics on the S3 objects themselves (§8).

---

## 2. Current state and precedent

`_source_format` was introduced in commit `2ed461c` ("Add TSV/BED fragment file input support..."), whose message explicitly says *"store `_source_format` attribute for provenance"*. This proposal is a continuation of that idea, not a new mechanism.

Existing backward-compatibility precedents in the reader, in the order of preference this design follows:

1. **`.attrs.get(name, default)` with a default that was implicitly true of every older file.** `self.source_format = self._f.attrs.get("_source_format", "BAM")` (`fragments_h5.py:296`). Correct there because every pre-`2ed461c` file *was* built from a BAM.
2. **`has_X` boolean properties doing dataset-presence checks**: `has_methyl` (`:328-331`), `has_gc` (`:345-349`), `has_fragment_end_clipped` (`:351-356`).
3. **Two-level capability checks for files that are permanently wrong in the field**: `has_strand` (`:333-343`) checks both presence *and* `len(shape) == 1`, to distinguish an old 2-bit encoding from the current 1-byte one, with a narrative comment noting the old "small frag" assay h5s were never rebuilt.
4. **Hard failure with a clear message** for optional per-fragment data: `fetch_array` raises `ValueError("The referenced h5 file does not contain gc info.")` and siblings (`:548-549, 563-564, 573-576, 585-588`).

Two gaps worth naming: no test anywhere constructs an h5 *missing* an attribute to assert fallback behavior (not even for `_source_format`), and there is no `CHANGELOG.md`. `docs/` currently contains `plan_chunk_based_parallelization.md` and `pending/` (which holds this document).

---

## 3. Scope: what gets recorded

Build-time behavior splits into three groups. They need different treatment.

### 3.1 Group A — parameterized by CLI (recorded as values)

Column 2 gives the **CLI** default (`main.py`) and, where it differs, the **library** default (`build_fragments_h5`, `fragments_h5.py:916-934`). They are not always the same, and the library parameter name is not always the flag name — the record uses the *library* parameter names (§4.1).

| CLI flag | library param | CLI default | library default | applied at | destructive |
|---|---|---|---|---|---|
| `--min-mapq` | `min_mapq` | `None` (validated `>= 0`, `main.py:85-86`) | `None` | PE `fragment.py:379-382`; SE `fragment.py:578` | yes |
| `--se-max-fragment-length` | `se_max_fragment_length` | `None` (validated 1..65535, `main.py:82-84`) | `None` | `fragments_h5.py:782` | yes |
| `--include-duplicates` | `include_duplicates` | `False` (dups excluded) | `False` | PE `fragment.py:375`; SE `fragment.py:576` | yes |
| `--single-end` | `single_end` | `False` | `False` | mode switch; changes extent computation | changes the meaning of every length |
| `--contigs` | `allowed_contigs` | `None` | `None` | `fragments_h5.py:1008-1018` / `1026-1028`, loop `:1062` | yes, by omission |
| `--contig-name-map` | `contig_name_map` | `None` (parsed to `dict`, `main.py:150-164`) | `None` | `fragments_h5.py:1039-1048` | not to fragments; **lossy on contig identifiers** — see below |
| `--set-mapq-255-to-none` | `set_mapq_255_to_none` | `False` | `False` | `fragments_h5.py:805-820` | no (value substitution), but aborts the build if unset and 255 is seen |
| `--fasta` | `fasta_filename` | `None` | `None` | enables `gc` (`:957-961`) | no |
| `--exclude-strand` | `read_strand` (inverted) | flag `False` → `read_strand=True` | `read_strand=True` | disables `strand` | no |
| `--read-methyl` | `read_methyl` | `False` | `False` | enables methyl fields | no |
| `--no-store-fragment-end-clipped` | `store_fragment_end_clipped` | flag absent → `True` | `True` | disables dataset | no |
| `--num-processes` | `num_processes` | string `'1'`, parsed at `main.py:145-148`; `'all'` → `None` | `None` | `fragments_h5.py:1088` | no (but see §4.8) |
| `--skip-chunking` | `skip_chunking` | `False` | `False` | `fragments_h5.py:1073` | no |

All of Group A is recorded. The parallelism flags are recorded too — they should not affect output, and recording them is exactly what makes that claim checkable if a chunk-boundary bug is ever suspected (`fragments_h5.py:779-780`).

**`--contig-name-map` must be recorded as its contents, not as a boolean.** An earlier draft recorded `contig_name_map_given: bool`. That throws away the only part that matters. `fragments_h5.py:1044-1045` rebuilds `_contig_lengths_str` from the **output** names only (`contig_lengths_output = {_map_name(c): length ...}`); the input names are not written anywhere in the file. In the motivating pipeline (`build_se_fragment_h5s.nf:145`, `grch38_refseq_to_ucsc.tsv`) the map converts the GenBank accessions listed at `.nf:29-35` to UCSC names for every contig, so a consumer holding the h5 cannot tell whether `chr1` was `chr1` or `CM000663.2` in the source BAM — i.e. cannot map the file back to its source. **Decision: record the full mapping dict** (`requested.contig_name_map`, `null` when not given). Size is not an argument against it: the map is bounded above by the contig count of the source, and `_contig_lengths_str` already enumerates every contig with its length, so this adds no new order of magnitude. A hash + entry count was considered and rejected — a hash lets you check two files used the *same* map but never lets you recover a name, which is the actual use.

`--fasta`, `--read-methyl`, `--exclude-strand` and `--no-store-fragment-end-clipped` are already inferable from dataset presence via `has_gc` / `has_methyl` / `has_strand` / `has_fragment_end_clipped`. Recording them anyway is cheap and removes the inference step.

**Correction on which side wins when the record and the datasets disagree.** An earlier draft said "the datasets are authoritative". That is wrong for strand, and strand is the one case where it has actually happened. `has_strand` (`fragments_h5.py:333-343`) returns `False` in two distinct situations: (a) the `strand` dataset is absent because `--exclude-strand` was passed, and (b) the dataset is *present* but 2-D, the old "small frag" 2-bit encoding that was never rebuilt. In case (b) a provenance record saying `read_strand: true` is **correct** and `has_strand == False` is the misleading signal. Corrected rule:

- The `has_*` properties answer *"can I read this data through the current API?"*
- The provenance record answers *"was this data requested at build time?"*
- These are different questions. A disagreement is a signal to investigate, not a contradiction to be resolved by fiat. For strand specifically, `has_strand == False` with `read_strand: true` recorded means "old encoding, or zero fragments" — not "the record is lying".
- (Moot in practice for strand: any file carrying a provenance record was built by a version that writes 1-byte strand. The rule is stated because it generalizes to the next encoding change.)

### 3.2 Group B — hardcoded, no flag, determined solely by the code version (recorded by *version*, not enumerated)

- **PE positive-TLEN-only selection.** `fragment.py:388`: `if include_neg_tlen or 0 < align.tlen <= max_tlen:`. `include_neg_tlen` defaults to `False` (`fragment.py:308`) and is never exposed or passed, so it is always `False`. This keeps one representative mate per template. Destructive, unconditional, unrecorded. It also makes the `align.tlen == 0` / `abs(tlen) > max_tlen` checks in `align_is_invalid` (`fragment.py:371-372`) and the `raise AssertionError("properly paired reads should not have a tlen of 0")` (`fragment.py:289-296`) dead code in the current call graph.
- **Hard 65535 length ceiling.** Signatures default `max_tlen=1000` (`fragment.py:307,450,554,641`) but the only caller passes `max_tlen=MAX_FRAG_LENGTH` = 65535 (`fragments_h5.py:774`); the 1000 default is unreachable. SE mode never references `max_tlen` at all — `single_end_bam_to_fragments` (`fragment.py:548-608`) declares and ignores it, so the only SE length bound is `--se-max-fragment-length`. The TSV path does apply it (`fragment.py:707-708`).
- **Read-flag filtering.** PE excludes `is_qcfail`, `is_supplementary`, `is_unmapped`, `mate_is_unmapped` (`fragment.py:373-377`); SE excludes `is_qcfail`, `is_supplementary`, `is_unmapped` (`fragment.py:574-577`). **`is_secondary` is never tested in either path** — secondary alignments are *not* excluded. Destructive, unconditional, unrecorded.
- **Zero-mapped-read contig skip, with an SE-specific off-by-half.** `fragments_h5.py:1064`: `if not is_tsv_input and num_mapped[bam_contig] == 0: continue`. Destructive by omission. `num_mapped` is built at `:1029-1030` as `{x.contig: x.mapped // 2 for x in bam_fp.get_index_statistics() ...}` — the `// 2` assumes every mapped read is half of a pair. **In single-end mode that assumption is false**, so a contig with exactly 1 mapped read gets `num_mapped == 0` and is dropped entirely. The motivating incident (§1.4) is an SE pipeline, which is why this belongs here rather than in a general-cleanup list. It affects only contigs with a single mapped read, so it is small — but it is silent, and it is one more thing `builder.version` has to carry.
- **PE-only per-chunk overlap drop.** `fragment.py:482-485`: `if chrom is not None and not intervals_intersect_with_none(frag_start, frag_stop, start, stop): continue`, inside `bam_to_fragments`. There is **no SE counterpart** — `single_end_bam_to_fragments` (`fragment.py:548-608`) relies solely on `af.fetch(chrom, start, stop)` plus the caller's start-position check at `fragments_h5.py:779-780`. Whether the two paths therefore admit exactly the same fragment set is a question this design does not answer; recording `single_end` in `effective` at least makes it possible to ask.
- **TSV row validation.** `fragment.py:676-704` drops rows with `ncols < 3`, wrong `ncols`, non-integer or negative coordinates, or `frag_stop <= frag_start`, each logged as a warning; if skipped/total > 0.5 for a contig it raises `ValueError` (`fragment.py:739-743`).
- **Fixed-width storage with no bound check.** Length is stored `uint16` (`fragments_h5.py:751`, dataset at `:875`), guarded only by `assert MAX_FRAG_LENGTH <= 2**16 - 1` at `:873` — which checks the constant, not the data. This is safe on the PE and TSV paths because both apply `max_tlen = MAX_FRAG_LENGTH = 65535`. **It is not obviously safe on the SE path**, which ignores `max_tlen` entirely (see the previous bullet) and computes its extent from `align.aend`, a CIGAR-derived reference span that D and N ops inflate above the sequenced length (§7.2). An SE alignment whose reference extent exceeds 65535 silently wraps in `uint16`. `--se-max-fragment-length` prevents this when it is set — which is one more reason it is worth recording — but it is optional for library callers, and `main.py:93-94` only requires it for BAM input in `--single-end` mode. Methyl counts are `uint8` (`:758`, written `:836-843`) with no bound check either; the justification in the comment at `:756-757` rests on a 150 bp read-length assumption that nothing in the code enforces.

**Decision: do not enumerate Group B as data in the provenance record.** Record `builder.version` and rely on it to pin these behaviors. Rationale:

- These are not parameters. There is no value to record — only "the code at version X did what version X does". A hand-written enumeration (`"excludes_qcfail": true, "excludes_secondary": false, ...`) is a second source of truth that will silently drift from the code the first time someone changes a predicate without remembering to update the list. A drifted provenance record is worse than none, because it is trusted.
- A version string is exact and externally verifiable: given `2.11.0` you can `git show v2.11.0:src/fragments_h5/fragment.py` and read the actual predicates.
- The cost is that the reader must consult source to interpret Group B. That is acceptable: Group B is *constant across all files built by a given version*, so it is a one-time lookup per version, not a per-file question.

Two exceptions, recorded as values because they are cheap and disambiguate something confusing:

- `effective.max_tlen` (= 65535, the value actually passed at `fragments_h5.py:774`, or `null` in SE mode where it is never consulted), because the root `max_fragment_length` attribute is misleading (§7.3) and because SE ignores `max_tlen` entirely — recording it under `effective` alongside `single_end` makes the SE/PE asymmetry legible.
- `builder.python` and the `pysam` version, since CIGAR interpretation (`aend`, §7.2) is pysam behavior.

### 3.3 Group C — requested is not effective (both recorded)

**All ten TSV/BED divergence paths, exhaustively, from `fragments_h5.py:971-1002` plus the SE/methyl gate at `:954-955`.** An earlier draft listed six and had no way to express "warned but not forced", which is a real and materially different outcome.

| line(s) | parameter | what the code does | disposition |
|---|---|---|---|
| `:954-955` | `single_end` + `read_methyl` | `raise NotImplementedError` (**not** TSV-specific — applies to BAM too) | aborts |
| `:972-973` | `fasta_filename is None` | `raise ValueError("--fasta is required for TSV/BED input ...")` | aborts |
| `:974-975` | `read_methyl` | `raise ValueError("--read-methyl is not supported for TSV/BED input ...")` | aborts |
| `:976-978` | `single_end` | warn, then `single_end = False` | `forced` |
| `:981-983` | `se_max_fragment_length` | warn, then `se_max_fragment_length = None` | `forced` |
| `:984-986` | `min_mapq` | warn, then `min_mapq = None` | `forced` |
| `:987-988` | `include_duplicates` | warn **only** — the variable is *not* reassigned | `warned_only` |
| `:989-990` | `set_mapq_255_to_none` | warn **only** — the variable is *not* reassigned | `warned_only` |
| `:991-993` | `store_fragment_end_clipped` | warn, then `store_fragment_end_clipped = False` | `forced` |
| `:998-1002` | `read_strand` when `num_columns < 6` | warn, then `read_strand = False` | `forced` |

The two `warned_only` rows are not equivalent to the `forced` rows and must not be recorded as if they were. `include_duplicates=True` and `set_mapq_255_to_none=True` are still passed to the worker in the tuple (`:1076-1083`, elements 14 and 13) and still reach `build_sub_fragments_h5`. They have no *effect* on the TSV path — `tsv_to_fragments` swallows `include_duplicates` via `**kwargs` (`fragment.py:645`), and `set_mapq_255_to_none` is consulted at `fragments_h5.py:805/813` but can never fire because every TSV fragment is constructed with `mapq1=None, mapq2=None` (`fragment.py:730-731`). But the *reason* they are inert is "the data has no such field", not "the value was overwritten", and a debugger tracing the worker will see `True`, not `False`. The record should say so.

**On the `disposition` enum — a partial disagreement with the review, stated rather than silently applied.** The review asked for three values: `forced` / `warned_only` / `rejected`. I am specifying only two. The three `aborts` rows above all `raise` before `h5py.File(ofname, "x")` at `:1139`, so no output file exists and no provenance record is ever written. A `"rejected"` value would therefore be unreachable in every conforming record — an enum member that can never be observed is a schema lie of its own, and a reader written to branch on it would have dead code. The three abort paths are documented in the table above because a *human* reading this needs them; they are not part of the on-disk vocabulary. If the schema later grows a way to record a failed build (it should not — a failed build has no file to attach to), `rejected` can be added under the additive-only rule in §5.3.

CLI validation (`main.py:92-96`) requires `--single-end` and `--se-max-fragment-length` together for BAM input only; TSV skips that pairing check. `main.py:82-86` range-validates both `--se-max-fragment-length` (1..65535) and `--min-mapq` (`>= 0`) before either check. `tests/test_cli_validation.py` covers this.

**Requirement: the record must state what was EFFECTIVELY applied, not what was requested.** A record that says `min_mapq: 30` for a TSV-built file is an actively harmful lie. But recording only the effective value loses the fact that the operator *asked* for 30 and did not get it — which is precisely the class of bug that motivated this document. So record both, in separate blocks, plus an explicit `neutralized` list so that the discrepancies the build *warned about* are greppable rather than requiring a diff of two objects.

**`neutralized` is not a complete index of `requested`/`effective` divergences, and must not be sold as one.** Its entries are generated from the warning sites at `:976-1002`, two of which (`:987`, `:989`) are guarded by `if <param>:` and so stay silent when the parameter is already `False` — while `effective` for those same parameters is `null` on every TSV build. A consumer that wants "where did requested and effective differ?" must diff the two blocks; `neutralized` answers the narrower "what did the build announce it was ignoring?". §4.1.1 states this normatively and §4.1.6 shows a conforming record that exhibits it.

**Implementation trap: the neutralized parameters are rebound in place.** `single_end` (`:978`), `se_max_fragment_length` (`:983`), `min_mapq` (`:986`), `store_fragment_end_clipped` (`:993`) and `read_strand` (`:1002`) are all reassigned *to the same local names* that hold the requested values. By the time control reaches `:1139`, where the record is written, the requested values are **gone**. `requested` must therefore be snapshotted at the top of `build_fragments_h5`, before `:954`. See §4.1's capture-point column; this is the single most likely way to implement this design incorrectly.

### 3.4 The "version pins Group B" argument has a verifiable hole in this repo's own release path

§3.2 rests entirely on `builder.version` being a faithful label for the source that ran. It is not guaranteed to be, and the counterexample is not hypothetical — it is producible today with the `Makefile` at `aa753c7`:

- `Makefile:114-130` (`tag`) refuses to tag only when **`pyproject.toml`** has uncommitted changes (`git diff --quiet HEAD -- pyproject.toml`). It does not look at `src/` at all.
- `Makefile:94-97` (`docker-build`) checks nothing whatsoever and runs `docker build -t $(IMAGE_NAME):$(VERSION) ... .` against the **working tree**, where `VERSION` is `grep`ped out of `pyproject.toml` (`:19`).
- `Makefile:46` chains them: `all: login tag conda docker clean`.

So: commit the version bump to `pyproject.toml`, leave a modified `fragment.py` in the working tree, run `make all`, and you get an image tagged `2.11.0`, containing a dist that reports `fragments-h5 2.11.0` to `importlib.metadata`, built from source that is not `v2.11.0`. Every h5 it produces would carry `"version": "2.11.0"` in its provenance record and that string would be false. This is exactly the scenario §3.2's argument cannot survive, and nothing currently prevents it. Editable/dev installs (`pip install -e .`) have the same property with less ceremony.

`RELEASE.md:148` states: *"**Version mismatch**: Not possible by construction — `pyproject.toml` is the only place the version is declared."* That is true of the *declaration* and false of the *correspondence between the declared version and the code*. This document contradicts a documented invariant and says so rather than quietly working around it; `RELEASE.md:148` should be reworded.

**The git-SHA alternative is concrete, not speculative.** Today there is no SHA available in-image: `.dockerignore:2` excludes `.git`, and `Dockerfile:12` is `COPY --chown=... . /tmp/fragments_h5`, so the build context has no git metadata. Making one available is small and fully specified:

1. Add `ARG GIT_SHA=unknown` / `ENV FRAGMENTS_H5_GIT_SHA=$GIT_SHA` to the `Dockerfile`.
2. In `Makefile:94-97`, compute the SHA plus a dirty marker and pass `--build-arg GIT_SHA=...`.
3. In the provenance record, emit `builder.git_sha` from `os.environ.get("FRAGMENTS_H5_GIT_SHA")`, `null` when unset.

**Decision: defer, explicitly.** It is out of scope for *this* change because it touches the release tooling and container build rather than the library, and because `builder.git_sha: null` is a perfectly honest value for a non-container build. But it is not "worth considering" — it is the single highest-value follow-up to this document, because §3.2 is load-bearing and `builder.version` alone does not carry it. Recorded in §10 as a named follow-up with an owner-shaped description, not as a vague aspiration. The schema (§4.1) reserves `builder.git_sha` as a nullable key from `schema_version: 1` so that adding it later is not a schema change.

---

## 4. Representation

### 4.1 Decision: one root attribute, `_build_provenance`, containing a JSON object

The attribute is a single UTF-8 `str`, written in the main process immediately after the existing block at `fragments_h5.py:1146-1151`. Everything in this section is **normative**: an implementation that differs from it is not implementing this design.

#### 4.1.1 Normative definition of `requested`, `effective` and `neutralized`

These three words carry the whole design and were used inconsistently in earlier drafts. Exactly one definition each:

> **`requested.<param>`** — the value of the `build_fragments_h5` parameter *as the caller supplied it*, snapshotted at the top of the function **before line 954**, before any validation, defaulting or neutralization. Its `null` means "the caller did not supply this parameter"; it never means anything else.
>
> **`effective.<param>`** — the value that **governed the behavior as actually applied to the data**, replicated in the main process. Its `null` means **"this parameter governed nothing — no executed code path consulted it for this build."** It does *not* mean "unset", and it is *not* necessarily the value that appears in the worker tuple.
>
> **`neutralized`** — a list of records, one per parameter for which the build **emitted a warning** at `fragments_h5.py:976-1002`. Empty list, never `null`, when no such warning fired. Sorted by `param` ascending (§4.1.2).

**Normative, and the thing a consumer is most likely to get wrong: `neutralized` is NOT the set of `requested`/`effective` divergences.** It is the set of *warned* divergences, and the two are not the same set. `requested` and `effective` may legitimately differ with **no** `neutralized` entry, whenever the parameter was already sitting at a value the code had no reason to warn about. The concrete cases, all on the TSV path:

- `:987` is `if include_duplicates:` and `:989` is `if set_mapq_255_to_none:`. A TSV build that left either at its default `False` emits no warning and therefore produces no `neutralized` entry — but `effective.include_duplicates` and `effective.set_mapq_255_to_none` are `null` for **every** TSV build regardless (§4.1.2), because TSV data has no duplicate marks and no MAPQ field. So `requested: false` / `effective: null` with an empty `neutralized` is a correct, conforming record. §4.1.6 shows exactly this for `set_mapq_255_to_none`.
- The same shape arises for `effective.max_tlen` and `effective.se_max_fragment_length`, which are `null` in PE and SE mode respectively by the "governed nothing" rule, with nothing warned and nothing to warn about.

A consumer that wants *"did requested and effective differ?"* must compare the two blocks. `neutralized` answers the narrower and more actionable question *"did the build tell the operator it was ignoring them?"* — which is the one tied to the class of bug in §1.4. Do not use one as a proxy for the other.

The `effective` definition has one consequence the implementation must get right, and it is the single most subtle point in this document:

**`effective` is not "the worker tuple".** An earlier draft defined `effective` as "as actually handed to the workers" and then asserted `effective.min_mapq == 0` for a BAM build with no `--min-mapq`. Those two statements contradict each other. What is handed to the worker is element 17 of the tuple at `fragments_h5.py:1082`, which is the raw `min_mapq` — i.e. `None`. The collapse to `DEFAULT_MIN_MAPQ` happens *inside* the worker at `:771` (`_min_mapq = min_mapq if min_mapq is not None else DEFAULT_MIN_MAPQ`), and it is `_min_mapq` that is passed to `bam_to_fragments` / `single_end_bam_to_fragments` at `:776` and compared against `align.mapq` at `fragment.py:379-382` / `:578`.

So the main process **must replicate the `:771` collapse** when computing `effective.min_mapq` for BAM input. Recording the raw tuple value would make `effective.min_mapq` identical to `requested.min_mapq` in every case, which is both useless and false.

And for TSV input the correct value is neither: `min_mapq` is forced to `None` at `:986`, and even if it were not, `tsv_to_fragments` accepts it into `**kwargs` (`fragment.py:645`) and never reads it. Recording `0` there — as a naive "apply the `:771` rule everywhere" implementation would — would be an outright lie, asserting that a MAPQ filter of 0 was applied to data that has no MAPQ field at all.

Worked values, normative:

| build | `requested.min_mapq` | `effective.min_mapq` | `neutralized` entry |
|---|---|---|---|
| BAM, no `--min-mapq` | `null` | `0` | none |
| BAM, `--min-mapq 0` | `0` | `0` | none |
| BAM, `--min-mapq 30` | `30` | `30` | none |
| TSV, no `--min-mapq` | `null` | `null` | none |
| TSV, `--min-mapq 30` | `30` | `null` | `{"param":"min_mapq","requested":30,"effective":null,"disposition":"forced","code":"tsv_no_mapq_field"}` |

Rows 1 and 2 are deliberately identical in `effective` and deliberately different in `requested`. That is the point: the data is genuinely indistinguishable (MAPQ is unsigned, so a threshold of 0 filters nothing), and the record says both "nothing was asked for" and "nothing was done" without pretending they are the same statement. See §4.4.

#### 4.1.2 Normative key table

`source expression` is written against the names in `build_fragments_h5`'s real signature (`fragments_h5.py:916-934`) — note `fasta_filename`, not `fasta`, and `allowed_contigs`, not `contigs`. `req_*` denotes the snapshot taken at the capture point.

**Capture points.** `A` = top of `build_fragments_h5`, before `:954`. `B` = at the write site, `:1146-1151`. Every `requested.*` key is capture point `A` and **must** be, because `single_end`, `se_max_fragment_length`, `min_mapq`, `store_fragment_end_clipped` and `read_strand` are rebound in place at `:976-1002` (§3.3).

| key | type | nullable | source expression | capture |
|---|---|---|---|---|
| `schema_version` | int | no | literal `1` | — |
| `builder.package` | str | no | literal `"fragments-h5"` | — |
| `builder.version` | str | no | `importlib.metadata.version("fragments-h5")` (§3.4 caveat) | — |
| `builder.git_sha` | str | **yes** | `os.environ.get("FRAGMENTS_H5_GIT_SHA")` — always `null` until §3.4's follow-up ships; reserved now so adding it is not a schema change | — |
| `builder.python` | str | no | `platform.python_version()` | — |
| `builder.pysam` | str | no | `pysam.__version__` | — |
| `builder.optimization_level` | int | no | `sys.flags.optimize` (0/1/2) — see §7.1 | — |
| `invocation` | list[str] | **yes** | the `invocation` parameter; `null` for library callers (§4.6) | A |
| `source_format` | str | no | `"TSV" if is_tsv_input else "BAM"` — duplicates the root `_source_format` attribute deliberately, so the record is self-contained | B |
| `requested.min_mapq` | int | yes | `req_min_mapq` | A |
| `requested.se_max_fragment_length` | int | yes | `req_se_max_fragment_length` | A |
| `requested.include_duplicates` | bool | no | `req_include_duplicates` | A |
| `requested.single_end` | bool | no | `req_single_end` | A |
| `requested.set_mapq_255_to_none` | bool | no | `req_set_mapq_255_to_none` | A |
| `requested.read_strand` | bool | no | `req_read_strand` | A |
| `requested.read_methyl` | bool | no | `req_read_methyl` | A |
| `requested.store_fragment_end_clipped` | bool | no | `req_store_fragment_end_clipped` | A |
| `requested.skip_chunking` | bool | no | `req_skip_chunking` | A |
| `requested.num_processes` | int | yes | `req_num_processes` — `null` covers *both* `--num-processes all` and the library default; see §4.8 | A |
| `requested.allowed_contigs` | list[str] | yes | `None if req_allowed_contigs is None else [str(c) for c in req_allowed_contigs]` | A |
| `requested.contig_name_map` | dict[str,str] | yes | `req_contig_name_map` (the parsed dict, `null` when not given) — §3.1 | A |
| `requested.fasta_given` | bool | no | `req_fasta_filename is not None` — the *path* is not recorded; it is machine-local and appears in `invocation` when relevant | A |
| `effective.min_mapq` | int | yes | `null` if TSV; else `min_mapq if min_mapq is not None else DEFAULT_MIN_MAPQ` (replicates `:771`) | B |
| `effective.se_max_fragment_length` | int | yes | `se_max_fragment_length if single_end else null` — `:782` gates on `single_end`, so in PE mode the value governs nothing even if a library caller passed one | B |
| `effective.max_tlen` | int | yes | `null` if `single_end` (SE never reads `max_tlen`); else `MAX_FRAG_LENGTH` (65535, `:774`; TSV applies it at `fragment.py:707-708`) | B |
| `effective.include_duplicates` | bool | yes | `null` if TSV (no duplicate marks; `tsv_to_fragments` swallows it via `**kwargs`); else `include_duplicates` | B |
| `effective.single_end` | bool | no | `single_end` (post-`:978`) | B |
| `effective.set_mapq_255_to_none` | bool | yes | `null` if TSV (`:805`/`:813` are executed but can never fire — every TSV fragment has `mapq1=mapq2=None`, `fragment.py:730-731`); else `set_mapq_255_to_none` | B |
| `effective.read_gc` | bool | no | `read_gc` (`:959`/`:961`) — named `read_gc`, not `fasta`, because that is what the worker consumes (tuple element 10) | B |
| `effective.read_strand` | bool | no | `read_strand` (post-`:1002`) | B |
| `effective.read_methyl` | bool | no | `read_methyl` | B |
| `effective.store_fragment_end_clipped` | bool | no | `store_fragment_end_clipped` (post-`:993`) | B |
| `effective.contigs` | list[str] | no | `sorted(str(c) for c in contig_lengths)` — the resolved **input-name** contig set, after `--contigs`/TSV-scan resolution but before the zero-mapped skip at `:1064`. Never `null`. | B |
| `effective.contigs_source` | str | no | `"cli"` if `req_allowed_contigs is not None` else `"source_default"` | B (depends on A)* |
| `effective.parallelism` | object | no | `{"multiprocess": bool, "num_processes": int\|null}` — see §4.8 | B |
| `effective.skip_chunking` | bool | no | `skip_chunking` | B |
| `neutralized` | list[object] | no (may be `[]`) | see below | B |

*`effective.contigs_source` is written at capture point B, but its source expression tests `req_allowed_contigs` — the capture-point-A snapshot — not the live `allowed_contigs` variable. This matters because `allowed_contigs` is rebound in place at `:1009` (TSV) and `:1027` (BAM) before B is reached. An implementer who reads the live `allowed_contigs` at B instead of the `req_allowed_contigs` snapshot from A would emit `"cli"` for every build, including `source_default` ones, because by B the variable is never `None`.

**Why `effective.contigs_source` exists.** `allowed_contigs` is rebound at `:1009` (TSV) and `:1027` (BAM, to `bam_fp.references`). Recording the post-default value alone would make "the operator asked for these 24 contigs" and "the operator asked for nothing and the BAM happened to contain these 24 contigs" indistinguishable — a smaller version of the same `None`-vs-default confusion this whole document is about. `requested.allowed_contigs` (`null` or the CLI list) plus `effective.contigs_source` resolves it. *Minor correction to the review, which called `bam_fp.references` "not JSON-serializable": pysam returns a `tuple` of `str`, and `json.dumps` encodes tuples as arrays, so it would serialize. The defect is the conflation, not the type. The explicit `[str(c) for c in ...]` coercion in the table is still required, for determinism and to avoid depending on that pysam detail.*

**`effective.contigs` has a real size cost in the `source_default` case, and it is accepted rather than bounded.** §3.1's size argument covers `contig_name_map` only ("bounded by the contig count of the source, and `_contig_lengths_str` already enumerates every contig"). It does not cover this key. For a BAM built **without** `--contigs`, `:1027` sets `allowed_contigs = bam_fp.references`, so `effective.contigs` is the *entire* reference list — for a full hg38 with alts, HLA and decoys that is on the order of 3,000 strings, tens of kilobytes, and it substantially duplicates `_contig_lengths_str` (which holds the same set under output names, with lengths). In the `cli` case it is bounded by what the operator typed (~24 in the motivating pipeline) and is a strict duplicate of `requested.allowed_contigs` (§4.1.5).

**Decision: accept the cost, do not truncate.** Reasons, in order: (a) the file it is attached to is typically hundreds of megabytes, so tens of kilobytes of ASCII in one attribute is not a meaningful fraction; (b) `_contig_lengths_str` already pays essentially this cost today, uncompressed, for the same corpus, so this is not a new order of magnitude for the format; (c) every bounding scheme considered — a count, a hash, "first N plus an ellipsis" — reintroduces exactly the §3.1 objection to hashing `contig_name_map`: it lets you check that two files match but never lets you recover a name, and recovering the input-name set is the entire purpose of the key. A reader that finds the size objectionable should note that HDF5 does not compress attributes, so the only real remedy would be dropping the key in the `source_default` case — which is the one case where it carries information the rest of the record does not (§4.1.5). That trade is not worth making.

**`neutralized` element order is normative: sorted by `param` ascending.** `json.dumps(sort_keys=True)` sorts *object* keys; `neutralized` is a JSON **array**, and `json.dumps` preserves list order verbatim. Without a stated rule, two conforming implementations could emit byte-different records for the same build, which would defeat §4.1.3's byte-identity requirement and §9 test 8. Sorting by `param` — not by `code`, not by source line — is chosen because `param` is the key a reader greps for and because the `code` vocabulary is expected to grow (§5.3) while the parameter names are fixed by `build_fragments_h5`'s signature.

**Why there is no `effective.contig_name_map`.** Every other Group A parameter appears in both blocks, and this one deliberately does not. The map has no neutralization path: `:1039-1040` only replaces `None` with `{}`, and `_map_name` (`:1041`) is then applied unconditionally at `:1044` and `:1067`. So `effective.contig_name_map` would restate `requested.contig_name_map` in every conforming record, differing only in spelling `null` as `{}` — a key that can never carry information, which is the same objection §3.3 raises against the unreachable `"rejected"` disposition. The map's *effect* is recoverable without it: `requested.contig_name_map` gives the mapping and the root `_contig_lengths_str` gives the post-map output names. (Contrast `effective.contigs`, which is retained because it *does* differ from `requested.allowed_contigs` in the `source_default` case — see §4.1.5.)

**`neutralized` record shape** — exactly five keys, all required:

```json
{"param": "min_mapq", "requested": 30, "effective": null,
 "disposition": "forced", "code": "tsv_no_mapq_field"}
```

- `param` — the `build_fragments_h5` parameter name.
- `requested` / `effective` — the same values as in the two blocks above, duplicated here so the entry is readable standalone.
- `disposition` — `"forced"` or `"warned_only"`, exactly as in §3.3's table. No other values; see §3.3 for why `"rejected"` is deliberately absent.
- `code` — a **stable identifier**, not free text and not a `file.py:line` citation. An earlier draft used `"ignored for TSV/BED input (fragment.py:645)"`; line numbers embedded in on-disk data are wrong by the next release and cannot be corrected retroactively. The `schema_version: 1` vocabulary is fixed at: `tsv_no_mapq_field`, `tsv_not_single_end`, `tsv_no_se_max_fragment_length`, `tsv_no_duplicate_marks`, `tsv_no_mapq_255`, `tsv_no_fragment_end_clipped`, `tsv_too_few_columns_for_strand`. Adding a code is additive (§5.3); readers must tolerate an unrecognized `code`.

#### 4.1.3 Serialization

```python
# PROPOSED — fragments_h5.py, adjacent to the existing attrs block at :1146-1151
f.attrs["_build_provenance"] = json.dumps(
    record,
    sort_keys=True,        # key order independent of insertion order
    separators=(",", ":"), # no insignificant whitespace
    ensure_ascii=True,     # pure-ASCII output regardless of contig/path encoding
    allow_nan=False,       # NaN/Infinity are not JSON; fail loudly rather than emit them
)
```

All four arguments are normative. Together with §4.7 (no timestamp, hostname or username) they make the attribute **byte-identical for two builds that differ only in when and where they ran**, given the same `invocation`. `sort_keys` alone does not achieve that.

h5py stores this as a variable-length UTF-8 string; `ensure_ascii=True` means the payload is ASCII regardless, so no encoding question arises on read.

**New imports required in `fragments_h5.py`.** The current import block is `:140-157` and contains **none** of the modules the §4.1.2 source expressions and the §4.1.3 / §4.1.4 sketches assume. Four must be added:

| module | needed for |
|---|---|
| `json` | `json.dumps` (§4.1.3) and `json.loads` / `json.JSONDecodeError` (§4.1.4) |
| `sys` | `sys.flags.optimize` → `builder.optimization_level` (§7.1) |
| `platform` | `platform.python_version()` → `builder.python` |
| `importlib.metadata` | `.version("fragments-h5")` → `builder.version` |

Already available at `:140-157` and needing nothing new: `os` (`:140`) for `builder.git_sha`'s `os.environ.get`, `pysam` (`:149`) for `builder.pysam`, and — the one that would otherwise be easy to miss — **`DEFAULT_MIN_MAPQ`, which is already imported at `:156`**, so replicating the `:771` collapse for `effective.min_mapq` (§4.1.1) requires no new import and no duplicated literal.

#### 4.1.4 Reader side

```python
# PROPOSED — fragments_h5.py, FragmentsH5
@property
def has_build_provenance(self) -> bool:
    return "_build_provenance" in self._f.attrs

@property
def build_provenance(self) -> dict:
    raw = self._f.attrs.get("_build_provenance")
    if raw is None:
        raise ValueError(
            f"'{self._f_fname}' does not record build provenance. It was built by "
            f"fragments-h5 < 2.12.0, or by a caller that did not pass one. The "
            f"build parameters are unknown; do not assume defaults. "
            f"Use `has_build_provenance` to test for this case."
        )
    try:
        return json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"'{self._f_fname}' has a '_build_provenance' attribute that is not "
            f"valid JSON: {exc}. The record is corrupt; this is NOT the same as "
            f"the record being absent, and must not be handled as such."
        ) from exc
```

Both branches name the file, matching what §6 and §9 assert about them. The `try/except` is not optional — without it a corrupt record surfaces as a bare `json.JSONDecodeError`, which is a `ValueError` subclass and so would satisfy a naive `pytest.raises(ValueError)` while giving the caller no filename and no indication that "corrupt" and "absent" are different states.

#### 4.1.5 Worked example A — the motivating SE BAM build

For `build-fragments-h5 --single-end --verbose --se-max-fragment-length 120 --num-processes 4 --include-duplicates --min-mapq 30 --fasta hg38.fa --contig-name-map map.tsv --contigs CM000663.2 CM000664.2 -- in.bam out.h5` (abbreviated to two contigs; the real pipeline passes ~24). Shown pretty-printed for readability — on disk it is minified per §4.1.3.

**`--contigs` takes INPUT (BAM) names, not output names.** `fragments_h5.py:1025` builds `all_contig_lengths` from `bam_fp.references`, and `:1028` does `{c: all_contig_lengths[c] for c in allowed_contigs}` — so passing an output name that is not also a BAM contig name raises `KeyError` and the build dies long before the record is written at `:1139`. The motivating pipeline passes GenBank accessions (`build_se_fragment_h5s.nf:29-35`, `CM000663.2 … CM000686.2`) and lets `--contig-name-map` rename them on the way out. This example does the same.

```json
{
  "builder": {"git_sha": null, "optimization_level": 0, "package": "fragments-h5",
              "pysam": "0.22.1", "python": "3.10.14", "version": "2.12.0"},
  "effective": {"contigs": ["CM000663.2", "CM000664.2"],
                "contigs_source": "cli",
                "include_duplicates": true,
                "max_tlen": null,
                "min_mapq": 30,
                "parallelism": {"multiprocess": true, "num_processes": 4},
                "read_gc": true, "read_methyl": false, "read_strand": true,
                "se_max_fragment_length": 120,
                "set_mapq_255_to_none": false,
                "single_end": true,
                "skip_chunking": false,
                "store_fragment_end_clipped": true},
  "invocation": ["--single-end", "--verbose", "--se-max-fragment-length", "120",
                 "--num-processes", "4", "--include-duplicates", "--min-mapq", "30",
                 "--fasta", "hg38.fa", "--contig-name-map", "map.tsv",
                 "--contigs", "CM000663.2", "CM000664.2", "--", "in.bam", "out.h5"],
  "neutralized": [],
  "requested": {"allowed_contigs": ["CM000663.2", "CM000664.2"],
                "contig_name_map": {"CM000663.2": "chr1", "CM000664.2": "chr2"},
                "fasta_given": true, "include_duplicates": true, "min_mapq": 30,
                "num_processes": 4, "read_methyl": false, "read_strand": true,
                "se_max_fragment_length": 120, "set_mapq_255_to_none": false,
                "single_end": true, "skip_chunking": false,
                "store_fragment_end_clipped": true},
  "schema_version": 1,
  "source_format": "BAM"
}
```

Two things to read off it that no existing attribute can tell you: `max_tlen` is `null`, so the 65535 ceiling did **not** apply and the only length bound was 120; and `contig_name_map` lets you take `chr1` in this file back to `CM000663.2` in the source BAM, which `_contig_lengths_str` cannot, because `:1044-1045` writes output names only.

**What `effective.contigs` does and does not add here.** When `--contigs` is given, `contig_lengths` is built from exactly those names (`:1028`, and `contig_lengths = eval(contig_lengths_str)` at `:1035`), so `effective.contigs == sorted(requested.allowed_contigs)` and carries **no new information**. Its only distinct value is the `contigs_source: "source_default"` case, where `allowed_contigs` was `None` and `:1027` filled it from `bam_fp.references` — there, and only there, does `effective.contigs` tell you something `requested` cannot. The contig map's effect is *not* visible in this pair at all; it is visible only by comparing `requested.contig_name_map` against the root `_contig_lengths_str` attribute, which holds the post-map output names.

(A related asymmetry, a real property of the code rather than of this record: `allowed_contigs` is matched against BAM contig names at `:1026-1028`, while the FASTA is looked up under *output* names at `:770` (`fasta_chrom = output_contig if output_contig != bam_contig else None`). So the contig map must line the BAM up with the FASTA. This design records the values; it does not reconcile them.)

#### 4.1.6 Worked example B — a TSV build with six neutralizations

For `build-fragments-h5 --min-mapq 30 --single-end --se-max-fragment-length 120 --include-duplicates --fasta hg38.fa -- frags.bed.gz out.h5` where `frags.bed.gz` has 4 columns:

```json
{
  "builder": {"git_sha": null, "optimization_level": 0, "package": "fragments-h5",
              "pysam": "0.22.1", "python": "3.10.14", "version": "2.12.0"},
  "effective": {"contigs": ["chr1", "chr2"],
                "contigs_source": "source_default",
                "include_duplicates": null,
                "max_tlen": 65535,
                "min_mapq": null,
                "parallelism": {"multiprocess": false, "num_processes": null},
                "read_gc": true, "read_methyl": false, "read_strand": false,
                "se_max_fragment_length": null,
                "set_mapq_255_to_none": null,
                "single_end": false,
                "skip_chunking": false,
                "store_fragment_end_clipped": false},
  "invocation": ["--min-mapq", "30", "--single-end", "--se-max-fragment-length", "120",
                 "--include-duplicates", "--fasta", "hg38.fa", "--", "frags.bed.gz", "out.h5"],
  "neutralized": [
    {"code": "tsv_no_duplicate_marks", "disposition": "warned_only",
     "effective": null, "param": "include_duplicates", "requested": true},
    {"code": "tsv_no_mapq_field", "disposition": "forced",
     "effective": null, "param": "min_mapq", "requested": 30},
    {"code": "tsv_too_few_columns_for_strand", "disposition": "forced",
     "effective": false, "param": "read_strand", "requested": true},
    {"code": "tsv_no_se_max_fragment_length", "disposition": "forced",
     "effective": null, "param": "se_max_fragment_length", "requested": 120},
    {"code": "tsv_not_single_end", "disposition": "forced",
     "effective": false, "param": "single_end", "requested": true},
    {"code": "tsv_no_fragment_end_clipped", "disposition": "forced",
     "effective": false, "param": "store_fragment_end_clipped", "requested": true}
  ],
  "requested": {"allowed_contigs": null, "contig_name_map": null, "fasta_given": true,
                "include_duplicates": true, "min_mapq": 30, "num_processes": 1,
                "read_methyl": false, "read_strand": true,
                "se_max_fragment_length": 120, "set_mapq_255_to_none": false,
                "single_end": true, "skip_chunking": false,
                "store_fragment_end_clipped": true},
  "schema_version": 1,
  "source_format": "TSV"
}
```

Compare `effective.min_mapq` here (`null` — no MAPQ filter was applied, because the data has no MAPQ) with `30` in example A. This is the case the earlier draft got backwards in both directions: §4.4 would have written `0` (asserting a filter that was never applied to a field that does not exist) and §9.4 would have asserted `null` for a reason that did not match §4.1's own definition. One definition, applied consistently, produces `null` here and `0` for a BAM build with no flag — and those are different facts.

Note also `effective.max_tlen: 65535` here versus `null` in example A: TSV *does* apply `max_tlen` (`fragment.py:707-708`) while SE does not. The `null` convention makes that visible without the reader having to know either fact in advance.

Three further details in this record are easy to get wrong and are called out explicitly:

- **`effective.parallelism.num_processes` is `null`, not `1`.** No `--num-processes` was passed, so `main.py:51`'s default `'1'` parses to the int `1` (`main.py:148`), and `fragments_h5.py:1088` (`num_processes is not None and num_processes != 1`) is **False** — the build takes the serial path at `:1117-1126`. Nothing was ever handed to `ctx.Pool(processes=...)`, so per §4.1.1's definition of `effective` (`null` = "this parameter governed nothing") the value is `null`. `requested.num_processes` is `1`, which is where the raw value is preserved. §4.8 and §9 test 9 use the same convention.
- **`set_mapq_255_to_none` has `effective: null` but no `neutralized` entry.** The flag was not passed, so it was already `False`, so `:989`'s `if set_mapq_255_to_none:` never fired and no warning was emitted. `effective` is still `null` because `:805`/`:813` cannot fire on TSV data. This is the `requested`/`effective` divergence-without-a-`neutralized`-entry case; see §4.1.1.
- **The `neutralized` array is sorted by `param` ascending**, per §4.1.2. `sort_keys=True` orders object keys but not array elements, so the order is fixed here rather than left to the implementation (§4.1.3's byte-identity requirement depends on it).

### 4.2 Why a JSON blob rather than ~35 individual typed attributes

Individual typed attributes match the surface style of `index_block_size` / `max_fragment_length`. They are rejected because:

- **`None` has no encoding.** h5py attributes have no null, and the §4.1.2 key table has **twelve** nullable top-level keys within `requested`/`effective`/`builder` — five in `requested` (`min_mapq`, `se_max_fragment_length`, `num_processes`, `allowed_contigs`, `contig_name_map`), five in `effective` (`min_mapq`, `se_max_fragment_length`, `max_tlen`, `include_duplicates`, `set_mapq_255_to_none`), plus `invocation` and `builder.git_sha`. (`effective.parallelism.num_processes` is additionally nullable as a nested value inside the non-nullable `parallelism` object — see §4.1.6/§4.8 — and is not counted among the twelve since it is not itself a top-level key.) For most of them `null` is not "unset" but a distinct semantic state (§4.1.1). Encoding that in typed attrs requires either a sentinel per key (a `max_fragment_length`-style trap, twelve times over) or a parallel `_min_mapq_is_set` boolean per nullable key. JSON has native `null` and native `true`/`false`.
- **Requested-vs-effective plus the nested `builder` and `parallelism` objects** puts the flat-name count around 35, plus the `neutralized` list, which is a list of five-key records and does not flatten at all.
- **Every future flag adds a new attribute name**, and every reader that wants to be robust must `.get()` each one with the right per-attribute default. One attribute with an internal `schema_version` extends without touching the reader's presence check.
- **There is already precedent for a serialized structure in a string attribute**: `_contig_lengths_str` (`fragments_h5.py:1151`). This proposal follows that shape while fixing its deserialization.

**JSON, not `repr(dict)` + `eval`.** `fragments_h5.py:293` currently does `eval(self._f.attrs["_contig_lengths_str"])` — arbitrary code execution driven by file content, on files routinely fetched from S3 by automated pipelines. Do not propagate that pattern. `json.loads` is the same amount of code and cannot execute anything. The existing `eval` at `:293` is pre-existing debt; see §10 for why this design does not fix it.

### 4.3 Naming

`_build_provenance`, leading underscore, matching the existing convention for machine-oriented attributes (`_bam_header`, `_source_format`, `_contig_lengths_str`) as distinct from the two "public" numeric ones (`index_block_size`, `max_fragment_length`).

### 4.4 `min_mapq`: what the requested/effective split actually buys, stated precisely

This section was the source of the worst error in the previous draft and is now subordinate to the normative definitions in §4.1.1. It exists to explain *why* those definitions are shaped the way they are.

**The mechanism.** `fragments_h5.py:771` does `_min_mapq = min_mapq if min_mapq is not None else DEFAULT_MIN_MAPQ`, with `DEFAULT_MIN_MAPQ = 0` (`fragment.py:15`). `_min_mapq` is then passed at `:776` and compared at `fragment.py:379-382` (PE) / `:578` (SE). Since MAPQ is unsigned, a threshold of 0 excludes nothing. **For BAM input, "no `--min-mapq`" and "`--min-mapq 0`" are behaviorally identical and produce byte-identical data.**

**Where line 771 lives, and why it matters to the implementation.** `:771` is inside `build_sub_fragments_h5` — the *worker*. The main process passes the raw `min_mapq` as tuple element 17 (`:1082`), so the value crossing the process boundary is `None`, not `0`. An implementation that records "what was handed to the workers" gets `None` and has recorded nothing at all. **The main process must replicate the `:771` expression** to compute `effective.min_mapq` (and must not apply it for TSV input, where the parameter is forced to `None` at `:986` and discarded by `tsv_to_fragments`'s `**kwargs` at `fragment.py:645` regardless). §4.1.1 gives the five worked cases.

**What the split buys.** For a BAM build:

- `requested.min_mapq` distinguishes `null` from `0` — i.e. it records operator *intent*, which the data cannot.
- `effective.min_mapq` is `0` in both cases — which is the truth, and is *why* the data cannot distinguish them.

For a TSV build with `--min-mapq 30`, `requested` is `30`, `effective` is `null`, and `neutralized` carries the discrepancy explicitly. That is the case that motivated the whole split: without it, either the record lies about a filter that was never applied, or it silently drops the fact that one was asked for.

**What the split does not buy.** It does not make a filtered file recoverable, and it does not let you tell "no filter" from "threshold 0" by looking at the data — only by looking at the record. If the record is absent (§5.1), neither question is answerable.

### 4.5 No worker-tuple change is required

`build_sub_fragments_h5` (`fragments_h5.py:723`) writes only *data* into per-chunk temp h5s; it writes **no attributes**. All five attributes today are written by the main process in `build_fragments_h5` (`:916`, attrs at `:1146-1151`), which already has every CLI parameter in scope as named arguments. Recording provenance therefore requires **no change to the 17-element positional worker tuple** (unpacked at `:735`, constructed at `:1076-1083`). This is the main reason this design is cheap. The worker tuple is a real maintenance hazard and may deserve its own design; it is not coupled to this one and is not addressed here.

### 4.6 `invocation`: fully specified, and re-justified

`build_fragments_h5` is an importable function. Reading `sys.argv` inside it would record whatever launched the *process* — `pytest`, a notebook kernel, an unrelated wrapper — and label it as the build command. That is a misleading record, which is worse than a missing one. So it is passed in, never read from `sys.argv` inside the library. That much was already decided; what follows is the part the previous draft left unspecified, and it is where the risk actually is.

**Normative:**

- Signature addition: `build_fragments_h5(..., invocation: list[str] | None = None)`. **The default is `None`**, so every existing library call site keeps working unchanged and records `"invocation": null`.
- `main()` passes **`sys.argv[1:]`**, not `sys.argv`. `argv[0]` is an install-dependent absolute path (`/opt/conda/bin/build-fragments-h5` in the container, something else under a dev install) that varies between environments for reasons unrelated to the build, and which carries no information not already in `builder.package` / `builder.version`. Excluding it is also what makes the attribute byte-identical across a container-vs-local comparison of the same logical build (§4.7).
- Recorded **verbatim**, as a JSON array of strings. No shell re-quoting, no path canonicalization, no truncation.
- The list is not parsed by anything. It is for a human to read and grep.

**Honest limits, which a reader of the record must be told:**

- Under Nextflow, the recorded paths are **per-task work-directory paths**, and the exact mechanism differs by profile — but the conclusion does not. `build_h5` declares its inputs as `path(bam)`, `path fasta`, `path contig_map` (`.nf:126-128`) and interpolates the staged names into the command (`.nf:144-147`), so what `build-fragments-h5` sees is a bare basename resolved against a work directory that exists only for the duration of the task. Under `-profile standard` those staged entries are symlinks into the local work tree and the task additionally runs in a scratch directory (`scratch = true`, `.config:23`, **`standard` only**). Under `-profile remote` — the profile that produced the S3 corpus — there is no `scratch` setting; inputs are **copied**, not linked (`stageInMode = 'copy'`, `.config:56`), into an AWS Batch task directory backed by `workDir = "s3://fragmentomics.kariusdx.com/nboley/nextflow-work/"` (`.config:39`). Either way the recorded path does not resolve anywhere after the run and is near-useless as provenance about *which* input was used. `_bam_header` is the better handle on input identity.
- In the motivating pipeline, `--contigs !{params.contigs.join(' ')}` (`.nf:146`) expands to roughly 24 tokens, so `invocation` is dominated by contig names that appear verbatim in both `requested.allowed_contigs` and `effective.contigs`. It is redundant with both. (Not with `_contig_lengths_str`, which that pipeline's `--contig-name-map` rewrites to UCSC names — §4.1.5.)
- It records *requested* only, never *effective*. A TSV build's `invocation` says `--min-mapq 30` and nothing applied it (§4.1.6).

**Does `invocation` earn its slot? Yes — but only just, and for one reason.** The review is right that the previous justification (an alleged `test_docker_build.py` constraint) was fictitious, so the question deserves a real answer. The structured `requested` block strictly dominates `invocation` for every *machine* use: it is typed, complete, and does not require reimplementing argparse. The one thing `invocation` has that `requested` does not is **fidelity to what the operator actually typed, including things this design does not model** — `--verbose` and future flags added before the key table is updated, argument order, and the raw `--fasta` path (which `requested.fasta_given` reduces to a boolean). When someone is trying to reproduce a build, the first thing they want is the command line, and a design that made them reconstruct it from a key table would be worse. Kept, with the limits above documented in this section so nobody mistakes it for authoritative.

*(Rejected alternative: recording only `invocation` and dropping `requested`. See §11 — it cannot express `effective` at all.)*

### 4.7 No timestamp, no hostname, no username

Deliberately excluded, for two reasons that stand on their own:

1. **Reproducibility.** Two builds of the same input with the same flags should produce the same file. A wall-clock timestamp guarantees they never do, which forecloses byte-comparison as a debugging tool for exactly the class of question ("did these two containers behave identically?") that this document is otherwise trying to make answerable. §4.1.3's `sort_keys` / `separators` / `ensure_ascii` / `allow_nan` settings exist to make the attribute byte-stable; adding a timestamp would defeat all four. Build time is recoverable from the S3 object's `LastModified` anyway, at no cost to the file.
2. **PII.** Hostname and username identify a person and a machine, and these files land in a shared bucket that is read by people and pipelines with no relationship to whoever ran the build. There is no analysis question they answer.

**Correction: the previous draft justified this with a test constraint that does not exist.** It claimed `tests/test_docker_build.py:91-99` "diffs two h5s attribute-by-attribute using `f.attrs.get(attr)`" and would therefore break. The actual code at `tests/test_docker_build.py:91` is:

```python
for attr in ['index_block_size', 'max_fragment_length']:
```

— a hardcoded two-element allow-list. It never enumerates `f.attrs` and would never observe `_build_provenance`. `tests/specialized/compare_chunked_vs_unchunked.py:32` uses the identical fixed tuple `("index_block_size", "max_fragment_length")`. Those two files are the only places under `tests/` that touch `.attrs` at all (the remaining hits, `test_docker_build.py:98-99`, read `_contig_lengths_str` by name). And `test_docker_build.py:5` states it "is NOT run by pytest automatically".

**Consequence: no existing test compares the full attribute set, so adding `_build_provenance` breaks nothing.** No test change is required by this design — the previous §9.8 mandated one and has been deleted. The reproducibility and PII arguments above are sufficient on their own and are the actual reasons.

(`invocation` does still contain paths that differ between a local and a containerized build of the same input. Since nothing compares the attribute, that is a property of the record, not a problem to solve. It is documented in §4.6.)

### 4.8 `num_processes` is recorded as resolved parallelism, because the raw value encodes a bug

The raw parameter cannot be recorded as-is without propagating a defect. `main.py:145-146` maps `--num-processes all` to `None`. `fragments_h5.py:1088` then tests `if num_processes is not None and num_processes != 1:` — so `None` takes the **single-process branch** at `:1117-1126`. Therefore:

- `--num-processes all` runs **serially**, which is the opposite of what it says.
- `--num-processes all` and the library default (`num_processes=None`, `:928`) are the same value, so a single `null` in the record would cover both an explicit request for maximum parallelism and no request at all.

**Decision: record both, and do not paper over it.**

- `requested.num_processes` — the raw value, `null` for both cases above. Faithful to what was passed.
- `effective.parallelism` — `{"multiprocess": <bool>, "num_processes": <int|null>}`, where `multiprocess` replicates `:1088`'s predicate exactly and `num_processes` is the value handed to `ctx.Pool(processes=...)` at `:1102`, or **`null` whenever `multiprocess` is `false`**, because on the serial path (`:1117-1126`) no `Pool` is constructed and nothing governed a process count. `num_processes` is therefore non-`null` if and only if `multiprocess` is `true`.

  The three cases, normative (and pinned by §9 test 9):

  | invocation | `requested.num_processes` | `effective.parallelism` |
  |---|---|---|
  | `--num-processes all` | `null` | `{"multiprocess": false, "num_processes": null}` |
  | none (CLI default `'1'`) | `1` | `{"multiprocess": false, "num_processes": null}` |
  | `--num-processes 4` | `4` | `{"multiprocess": true, "num_processes": 4}` |

  Row 1 is the truth about the `all` bug and makes it visible in any file built that way. **Row 2 is the one an implementer is most likely to get wrong**: the value is `null`, not `1`, even though `1` is what was passed and what `requested` records. `:1088` short-circuits on `num_processes != 1`, so a literal `1` takes the serial branch just as `None` does, and `effective`'s "what governed" definition (§4.1.1) makes both `null`. Recording `1` there would assert that a one-worker `Pool` existed, which is false — worker code ran in the main process via a direct call at `:1123`.

Fixing `:1088`/`main.py:145-146` so that `all` means `os.cpu_count()` is a real bug fix and is **out of scope here** (§10): it changes runtime behavior, and this document's job is to record behavior, not change it. Recording it accurately is the cheapest way to make it discoverable.

---

## 5. Backward compatibility

### 5.1 The core issue: absent is not the same as default

Precedent #1 (`.attrs.get(name, default)`) works for `_source_format` because the default `"BAM"` was *provably true of every file built before the attribute existed* (`fragments_h5.py:296`). That property does not hold here. An h5 built by 2.10.1 might have been built with any combination of flags that version supported; defaulting its `include_duplicates` to `False` would be a guess presented as a fact.

So a plain `.get(default)` per parameter is **rejected**. With it, absent and present-but-default are indistinguishable, and every consumer that reads `prov["effective"]["min_mapq"] == 0` would conclude "no MAPQ filter" for a file that might have had one.

**Decision:** absence is a first-class state, detected by `has_build_provenance` (§4.1), following the `has_X` convention (`fragments_h5.py:328-331, 345-349, 351-356`). Nothing in the provenance record has a fallback default. The record is present in full or absent in full; there is no partial-credit path.

This is the same reasoning as `has_strand`'s two-level check (`:333-343`), which exists because old files in the field are permanently wrong and the reader must cope rather than assume.

### 5.2 Existing readers keep working — and the search that establishes it

Adding a root attribute changes nothing for existing readers. `has_build_provenance` is a sufficient guard, and **no existing consumer needs any change**. That is a claim about code that is not in this repo, so here is the search that backs it and its exact scope.

**Search performed** (at `fragments_h5` `aa753c7`, `fragmentomics_tools` and `biomarker-pipeline` at their current checkouts):

1. `FragmentsH5.__init__` (`fragments_h5.py:266-304`) read line by line. It reads **four** attributes by name — `_contig_lengths_str` (`:293`), `index_block_size` (`:294`), `max_fragment_length` (`:295`), `_source_format` (`:296`) — and one dataset by name (`:297-298`). It never iterates `self._f.attrs`, never calls `.attrs.keys()`, and has no `**` expansion over attributes. An unknown sixth attribute is invisible to it.
2. `grep -rn "\.attrs" fragmentomics_tools --include=*.py` → **zero matches.** So does grepping for `max_fragment_length`, `_source_format` and `_contig_lengths_str`. `fragmentomics_tools` reads no h5 attribute of any kind.
3. `RegionFragmentArray.from_fragments_h5` (`fragmentomics_tools/.../fragment_array/fragment_array.py:1708-1847`) reaches the file only through `FragmentsH5.fetch_array` (called at `:1758-1763`). Its capability checks at `:1740-1743` are `has_methyl` (`:1740`), `has_strand` (`:1741`) and `has_gc` (`:1743`), which are dataset-presence checks on `self.data` (`fragments_h5.py:328-356`), not attribute reads.
4. Within `fragments_h5` itself, the only code that enumerates attributes is in `tests/` — `tests/test_docker_build.py:91` and `tests/specialized/compare_chunked_vs_unchunked.py:32` — and both use fixed two-element allow-lists (§4.7), so neither sees a new attribute either.

**Scope of the claim.** This covers the three repositories readable here. It does **not** cover ad-hoc notebooks, scripts outside these repos, or any consumer that opens the h5 with raw `h5py` and iterates `f.attrs` — such a consumer would newly see a `_build_provenance` key. That is a benign addition (h5py attribute iteration is order-agnostic and a new key breaks nothing that was not already fragile), but it is the one category the search cannot rule out, and it is stated rather than glossed.

Old readers opening new files ignore `_build_provenance`. New readers opening old files get `has_build_provenance == False`. No existing dataset, dtype, or attribute changes.

### 5.3 Should a file-format schema version attribute be introduced now?

**Decision: no.** Add `schema_version` *inside* the provenance JSON, scoped to the provenance record only. Reasons:

1. A root `format_version` would be absent from every existing file, so it needs its own fallback — presumably `0` — which conveys nothing that the existing per-feature presence checks do not already convey.
2. This codebase's evolution model is per-feature capability detection (`has_methyl`, `has_gc`, `has_fragment_end_clipped`, and `has_strand`'s two-level check). That model has already absorbed a genuine on-disk encoding change (strand, 2-bit → 1-byte, `:333-343`) without a version number. Introducing a version now would give two competing mechanisms for the same question.
3. A version integer invites `if version >= N:` gates, which is exactly the brittle pattern that breaks against the heterogeneous files actually in the field.

**Why adding `schema_version` at record scope does not re-import objection 3.** Objection 3 is about *file* scope: a root `format_version` would be the gate on whether the reader can interpret the file at all, and getting that gate wrong makes production data unreadable (§11's "bump a format version" row). At record scope the blast radius is one optional attribute — a reader that mishandles `schema_version` loses provenance and nothing else. And the record, unlike the file, has no legacy: every record that exists was written by code that knew the schema. The scope change defuses the objection because it changes what a mistake costs, not because it makes mistakes less likely.

**Evolution rule for the record — normative, and the thing the previous draft omitted:**

- **`schema_version: 1` is additive-only.** Within version 1, future releases may add keys to any object and may add members to the `neutralized` `code` vocabulary. They may **not** remove a key, change a key's type, change a key's nullability, or change what a key means.
- **Readers MUST ignore unknown keys** and MUST NOT validate the key set exhaustively. A reader that rejects a record because it contains a key the reader does not recognize is non-conforming.
- **Readers MUST NOT reject a record because `schema_version` is higher than they know**, unless they depend on a key that a later version removed — which, under the additive-only rule, can only happen at version 2 or above. The correct default behavior on `schema_version > 1` is to proceed and ignore what you do not recognize.
- **`schema_version` increments only on a breaking change** — removing a key, changing a type or nullability, or changing a meaning. That is the only situation in which a reader is entitled to branch on it, and at that point a `>=` gate is warranted precisely because there is a real incompatibility to gate on.
- **The version does not track the builder version.** `builder.version` moves every release; `schema_version` is expected to stay at `1` indefinitely. A record with `builder.version: "3.4.0"` and `schema_version: 1` is normal, not stale.

---

## 6. Reader-side API

Two members on `FragmentsH5`, no more:

| member | type | absent | present and valid | present and unparseable |
|---|---|---|---|---|
| `has_build_provenance` | `bool` property | `False` | `True` | **`True`** — it is a presence check, not a validity check |
| `build_provenance` | `dict` property | raises `ValueError` naming the file, the likely cause, and `has_build_provenance` | returns the parsed dict | raises `ValueError` naming the file and stating that corrupt ≠ absent |

Both `ValueError` messages are given verbatim in §4.1.4, including the `try/except json.JSONDecodeError` that produces the third column. That code is the specification; this table is a summary of it. (An earlier draft described the corrupt-record behavior here in §6 but omitted both the filename and the `try/except` from the §4.1 sketch, so the sketch would have surfaced a bare `json.JSONDecodeError` with no filename. They now agree.)

**Why `build_provenance` raises rather than returning `None`:** the hard-fail-with-clear-message precedent in `fetch_array` (`:548-549` and siblings) applies directly, and there is a specific hazard here. If it returned `None`, the natural defensive idiom is `(prov or {}).get("effective", {}).get("min_mapq")` → `None`, which is *the same value* the record uses for "no filter requested" (§4.4). Returning `None` on absence would collapse "unknown" into "unfiltered" — the exact confusion this design exists to eliminate. Raising makes the caller decide.

Callers that want tolerance write the two-line form:

```python
# PROPOSED — consumer pattern
if not h5.has_build_provenance:
    log.warning("%s predates build provenance; filters unknown", path)
else:
    eff = h5.build_provenance["effective"]
```

**Not proposed:** convenience accessors like `h5.effective_min_mapq`. They would each need their own absent-behavior decision, and each is a place where "unknown" can leak back into a value-shaped return. One dict, one presence flag.

`__init__` is **not** modified to parse the JSON eagerly. Parsing is deferred to the property so that opening a malformed file still works for data access, and so `__init__` (`:290-304`) stays as cheap as it is today. If the attribute is present but unparseable, the property raises `ValueError` naming the file (§4.1.4) — a corrupt record must not be silently treated as absence, and `has_build_provenance` deliberately still returns `True` for it, because the honest answer to "does this file have a provenance record?" is yes, and the honest answer to "can I read it?" is the exception.

---

## 7. Verified findings that constrain the design

### 7.1 The MAPQ 255 sentinel is structural

MAPQ is stored `uint8`, shape `(N, 2)`. `fragments_h5.py:822-823` writes `mapq_arr[ff, 0] = 255 if fragment.mapq1 is None else fragment.mapq1`. On read, `fetch_array` (`:543-545`) casts to `int32` and maps `255 → -1`; `fetch()`'s `_negative_1_to_none` (`:626-631`, used at `:687-688`) maps `-1 → None`. **A genuine MAPQ of 255 is structurally unrepresentable.** Recording `set_mapq_255_to_none` in provenance is therefore meaningful: it distinguishes "255s were seen and deliberately mapped to None" from "255s were never seen".

**The `python -O` hazard, and the conclusion the previous draft declined to draw.** The `assert False` at `:809-812` / `:817-820` is a sentinel-collision guard, not a QC check: it is the only thing standing between a real MAPQ 255 and its silent conversion to `None`. `assert` statements are removed entirely by the compiler under `python -O` / `PYTHONOPTIMIZE`. So if optimization is ever enabled, the guard vanishes, `mapq_arr[ff, 0] = 255 if fragment.mapq1 is None else fragment.mapq1` at `:822` stores a real 255 identically to a `None`, and the file silently contains `None` MAPQs that were never `None`. Worse for this design specifically: the record would say `set_mapq_255_to_none: false`, which a reader would correctly interpret as "no 255-to-None substitution was requested" and incorrectly conclude as "therefore any `None` MAPQ in this file is a genuine missing value". **`set_mapq_255_to_none: false` is actively misleading under `-O`, and this is a case where the provenance record could make things worse rather than better.**

**Decision: record `builder.optimization_level` = `sys.flags.optimize` (§4.1.2).** One integer, no runtime cost, evaluated in the main process. `0` means the guard was live and `set_mapq_255_to_none: false` is trustworthy; `1` or `2` means it was compiled out and the flag's `false` value says nothing about whether 255s were encountered. This is cheaper than any alternative (removing the `assert`, or converting it to a raised exception, are both algorithmic changes requiring approval — see §10), and it is the minimum needed to keep the record from lying.

Caveat on scope: `sys.flags.optimize` is read in the main process. Workers are created with `multiprocessing.get_context('fork')` (`:1098`), which inherits the parent's interpreter flags, so the value is correct for the workers too. It would not be under a `spawn` context with different flags — not a configuration this code can currently produce, but the reason the value is meaningful is worth stating rather than assuming.

### 7.2 SE extent is the reference span

`fragment.py:581-582`: `frag_start = align.pos`, `frag_stop = align.aend`. `aend` is pysam's `reference_end`, CIGAR-derived, counting reference-consuming ops (M/D/N/=/X) and excluding I and soft-clips. Deletions and N-skips inflate the extent relative to sequenced bases. This is why `--se-max-fragment-length 120` is an artifact guard, not a read-length cap — and why recording it matters: you cannot infer it from read length.

### 7.3 `max_fragment_length` is a misleading constant and is load-bearing

`fragments_h5.py:1148` writes the module constant 65535 unconditionally. It is not the observed maximum and is unrelated to `--se-max-fragment-length` (there is no `--max-fragment-length` flag at all). On read (`:295`) it becomes `self.max_fragment_length`, used as (a) the default backward search window in `fetch_array` (`:449`, `max_frag_len = self.max_fragment_length`) and (b) the size of the `fragment_length_counts` histogram (`:713`, `numpy.zeros(self.max_fragment_length + 1)`). `README.md:133` documents it as "the maximum fragment length stored", which is wrong.

**Design implication: it is load-bearing in two places, not one, and the second is worse.** Besides the backward search window at `:449`, `max_frag_len` is also a **filter term in the overlap mask** at `:524`:

```python
mask = (starts < region_stop) & (stops > region_start) & (lengths <= max_frag_len)
```

So lowering the stored constant would not merely cause `fetch_array` to *miss* fragments outside the search window — it would **actively drop** every fragment longer than the new value from every read, on every query, silently. A file whose `max_fragment_length` was "corrected" to an observed maximum of, say, 180 would return no fragment longer than 180 even for fragments physically present in the `lengths` dataset. Do **not** repurpose it. This design adds `effective.max_tlen` under a different name and leaves `max_fragment_length` untouched (§10).

### 7.4 GC is a lossy 8-bit quantization, and can be silently absent per-contig

Two distinct facts, neither recorded today and neither fixable by this design — documented so that `effective.read_gc: true` is not over-read.

**Quantization.** GC is computed as a float, rounded to 5 decimal places (`fragment.py:491` for PE, `:587` for SE, `:719-723` for TSV), then stored as `uint8` via `int(round(fragment.gc * 254))` (`fragments_h5.py:828-830`) with 255 reserved as the null sentinel, and read back as `_tmp_gc.astype("float32") / 254.0` (`:557`). The round trip is therefore lossy to roughly 1/508 ≈ 0.002 in GC fraction, and the 5-decimal rounding upstream is irrelevant next to it. This is the same class of structural sentinel-plus-quantization decision that MAPQ 255 gets a whole section for (§7.1), and it deserves naming for the same reason: a consumer comparing GC values across files, or thresholding near a boundary, is working with 254 buckets, not a float.

**Per-contig silent absence.** `get_g_or_c_cumsum` returns `(None, 0)` when the contig is not in the FASTA (`fragment.py:419-421`, `if chrom not in fasta_file`). Every fragment on that contig then takes the `gc = None` branch (`fragment.py:487-488` / `:584-585`), is stored as 255, and reads back as `NaN` (`fragments_h5.py:555-559`). **`has_gc` stays `True`** throughout, because it is a dataset-presence check (`:345-349`) and the `gc` dataset is created for every contig when `read_gc` is set (`:877-878`). So a file can be GC'd on 20 contigs and all-NaN on 4, and nothing in the file distinguishes that from a build with no FASTA at all except counting NaNs.

**Not addressed by this design, deliberately.** `effective.read_gc` records whether GC was *requested*, which is all the main process knows: whether a given contig actually resolved in the FASTA is determined inside the worker and is not returned to `build_fragments_h5` (the worker's return tuple is `(output_contig, chunk_start, chunk_stop, ofname)`, `:889`). Recording per-contig GC coverage would require a worker-interface change, which §4.5 is explicitly built to avoid. A consumer who needs this can compute it: `numpy.isnan(gc).all()` per contig. Noted here so the gap is deliberate rather than overlooked.

### 7.5 No downstream consumer reads `max_fragment_length`

Only `FragmentsH5` itself does. `fragmentomics_tools` reads no h5 attributes at all — `fragment_array/fragment_array.py:19` imports `FragmentsH5` and `RegionFragmentArray.from_fragments_h5` (`:1708-1847`) goes through `fetch_array` (`:1758-1763`). This means the new attribute has no existing consumer to break and no existing consumer to automatically benefit; adoption is opt-in. (The full search backing the "no consumer breaks" half of that claim is in §5.2, which is where it belongs — this section is about `max_fragment_length` specifically.)

---

## 8. Identifying already-built files

**This is remediation, not design, and it is out of scope for the implementation** — but it belongs in this document because it determines whether remediation is even possible, and because the answer is "partially, already, without any new metadata."

Existing files are already partially forensically classifiable:

- **`fragment_length_counts`** is a top-level dataset already written for every file. `_add_fragment_length_counts` (`fragments_h5.py:708-720`, writing the dataset at `:720`) is invoked unconditionally at `:1214-1216`, after the output file is closed and reopened `"r+"`, on every successful build. It is the observed length distribution. An SE file built with `--se-max-fragment-length 120` has zero counts above 120.
- **The `mapq` dataset** reveals the observed minimum MAPQ. A file built with `--min-mapq 30` contains no fragment below 30.
- **`_contig_lengths_str`** reflects the *final* contig set, so `--contigs` and `--contig-name-map` are already partially visible in output names.
- **`_source_format`** distinguishes BAM from TSV input for files built since `2ed461c`.

**Limits — this is evidence, not proof:**

- Absence of fragments longer than 120 is not proof a cap was applied. A shallow, short-insert, or heavily-degraded library can produce the same distribution. Conversely, presence of a fragment at 121 *is* proof the 120 cap was **not** applied — the negative direction is sound, the positive is not.
- A file could coincidentally contain no fragment below MAPQ 30, particularly for small contig subsets or high-quality libraries. Again asymmetric: observing MAPQ 12 disproves `--min-mapq 30`; not observing anything below 30 does not prove it.
- Group B behaviors (§3.2) leave no distributional signature at all and are unrecoverable by inspection.
- `_bam_header` may carry an aligner `@PG` line but records nothing about the `build-fragments-h5` invocation.

So for the specific open question in §1.4 — whether uncapped SE h5s exist in `s3://karius-biomarker-data-assets/projects/ibd/frag_h5s_se/` — the *disproving* direction is available today: scan each object's `fragment_length_counts` for any nonzero bin above 120, and scan `mapq` for any value below 30. Any hit is conclusive evidence of an unfiltered build. No hits is suggestive but not conclusive. That scan needs no code from this design and should be run independently of it.

---

## 9. Testing strategy

New test module `tests/test_build_provenance.py`. **No changes to any existing test are required** — see §4.7; the two places under `tests/` that touch `.attrs` both use fixed two-element allow-lists and cannot observe a new attribute.

1. **Round-trip.** Build a small BAM-derived h5 with a known flag set; assert `has_build_provenance` and that every `effective` key matches the §4.1.2 source expression for those flags.
2. **Absent-metadata fallback** — *the class of test that does not exist today, for any attribute*. Build an h5, then `del f.attrs["_build_provenance"]` with h5py, reopen, and assert `has_build_provenance is False` and that `build_provenance` raises `ValueError` whose message contains the filename. Add the equivalent test for `_source_format` (`:296`) while there — currently nothing verifies that fallback either.
3. **`min_mapq` None vs 0 (BAM).** Two builds, one with no `--min-mapq`, one with `--min-mapq 0`. Assert identical fragment data, `requested.min_mapq` `null` vs `0`, and `effective.min_mapq == 0` in **both** — because for BAM input the main process replicates the `:771` collapse (§4.1.1, rows 1-2). A test that expected `null` here would be asserting the previous draft's bug.
4. **TSV neutralization, exhaustive.** Build from a 4-column TSV with `--fasta REF --min-mapq 30 --single-end --se-max-fragment-length 120 --include-duplicates --set-mapq-255-to-none`. **`--fasta` is mandatory here, not incidental:** `fragments_h5.py:972-973` raises `ValueError("--fasta is required for TSV/BED input ...")` before any of the neutralization sites at `:976-1002` are reached, so an invocation without it aborts and produces no file and no record to assert against. Assert:
   - `requested` holds all five as passed;
   - `effective.min_mapq is None` (**not** `0` — TSV has no MAPQ field, so no threshold governed anything; §4.1.1 row 5), `effective.single_end is False`, `effective.se_max_fragment_length is None`, `effective.include_duplicates is None`, `effective.set_mapq_255_to_none is None`, `effective.read_strand is False`, `effective.store_fragment_end_clipped is False`, `effective.max_tlen == 65535`;
   - `neutralized` contains one entry per §3.3 row that fired, with `disposition == "forced"` for `min_mapq` / `single_end` / `se_max_fragment_length` / `store_fragment_end_clipped` / `read_strand` and `disposition == "warned_only"` for `include_duplicates` / `set_mapq_255_to_none`;
   - `[e["param"] for e in neutralized]` equals its own `sorted()` — this is what pins the §4.1.2 ordering rule. Test 8's byte-identity check cannot pin it, since that test compares one implementation against itself and would pass just as well on stable source-line order.
   Mirror the cases in `tests/test_cli_validation.py`.
5. **Requested survives rebinding.** The regression test for the §3.3 trap. Build from TSV with `--fasta REF --single-end --se-max-fragment-length 120 --min-mapq 30` (again `--fasta` is required, per test 4) and assert `requested.single_end is True`, `requested.se_max_fragment_length == 120`, `requested.min_mapq == 30`. An implementation that snapshots after `:971` instead of before it passes tests 1 and 3 and fails only this one.
6. **Library invocation.** Call `build_fragments_h5(...)` directly with no `invocation=`; assert `record["invocation"] is None` and that the rest of the record is fully populated.
7. **Corrupt record.** Set `_build_provenance` to `"{"`; assert `has_build_provenance is True` and `build_provenance` raises `ValueError` naming the file. Assert the message is *not* the absent-record message — corrupt and absent must be distinguishable from the exception alone.
8. **No `eval`, and byte-stability.** Assert the attribute is a `str` and `json.loads` succeeds. Then build the same input twice to different output paths, passing the *same* `invocation`, and assert the two `_build_provenance` attribute strings are **byte-identical** — the test that §4.1.3's `sort_keys` / `separators` / `ensure_ascii` and §4.7's no-timestamp rule actually hold.
9. **`num_processes`, both serial spellings and the parallel one.** Three builds, pinning §4.8:
   - `--num-processes all` → `requested.num_processes is None`, `effective.parallelism == {"multiprocess": False, "num_processes": None}`. This pins the bug-recording behavior so that a later fix to `:1088` has to update the test deliberately.
   - **no `--num-processes`** (CLI default `'1'` → int `1`, `main.py:51`/`:148`) → `requested.num_processes == 1` but `effective.parallelism == {"multiprocess": False, "num_processes": None}` — **`None`, not `1`**, because `:1088` is False and nothing reached `ctx.Pool` (`:1102`). This is the case §4.1.6 shows, and asserting `1` here would be asserting a draft error.
   - `--num-processes 2` → `requested.num_processes == 2`, `effective.parallelism == {"multiprocess": True, "num_processes": 2}`.
10. **Contig keys.** The keys with the highest error rate in this document's own drafting history, and previously untested. Build a small BAM whose references are input-style names (e.g. `CM000663.2`, `CM000664.2`, plus a third the build is not asked for) with **both** `--contigs CM000663.2 CM000664.2` and a `--contig-name-map` file mapping them to `chr1` / `chr2`. Assert:
    - `requested.contig_name_map == {"CM000663.2": "chr1", "CM000664.2": "chr2"}` — the full dict, not a boolean (§3.1);
    - `effective.contigs == ["CM000663.2", "CM000664.2"]` — **input** names, sorted, and specifically *not* the mapped output names;
    - `effective.contigs_source == "cli"`, and the excluded third reference is absent from `effective.contigs`;
    - `_contig_lengths_str` parses to output names (`chr1`, `chr2`), i.e. the record and the legacy attribute disagree by exactly the map — which is the only place the map's effect is observable (§4.1.5).
    Then repeat the build **without** `--contigs` and assert `effective.contigs_source == "source_default"` and that `effective.contigs` now contains all three references while `requested.allowed_contigs is None`. Also assert that passing an *output* name to `--contigs` (`--contigs chr1`) fails the build rather than producing a record — `:1028` raises `KeyError` — so the "input names only" rule is pinned rather than merely documented.
11. **Schema version and reader tolerance.** Pins §5.3's normative reader rules, which are otherwise untested. Assert `record["schema_version"] == 1` on a freshly built file. Then, on a copy, rewrite `_build_provenance` with (a) an added unrecognized key at top level and inside `effective`, and (b) `schema_version` set to `2`, and assert `build_provenance` returns the parsed dict without raising in both cases — a reader that rejects unknown keys or a higher `schema_version` is non-conforming.

---

## 10. Explicitly out of scope

| Not doing | Why |
|---|---|
| Fixing `max_fragment_length` semantics | Load-bearing in **two** places: `fetch_array`'s backward search window (`:449`) *and* the `lengths <= max_frag_len` term of the overlap mask (`:524`), plus the histogram size (`:713`). Changing it is a silent-correctness hazard — too low actively **drops** fragments on every read (§7.3), no error. A separate design with its own test plan. |
| Fixing the stale `README.md` attribute block (`:126-134`) | The block is not merely imprecise at `:133` — it is broadly wrong. It documents `ref` (`:128`) and `sample_id` (`:129`) as `FragmentsH5` attributes; **neither exists** (the class defines `name` at `:241`, `filename` at `:246`, and the four `has_*` properties at `:328-356`). It also omits `has_gc`, and describes `max_fragment_length` as "the maximum fragment length stored" (§7.3). Deferring this is a judgment call, not an oversight: correcting it means deciding what the documented public surface *is*, and this design adds two members to that surface. The right sequencing is to fix the block in the same change that adds `has_build_provenance` / `build_provenance` to it — i.e. this deferral should be revisited at implementation time rather than treated as settled. |
| Threading a git SHA into the container (§3.4) | The highest-value follow-up, and fully specified in §3.4 (three steps: `ARG GIT_SHA` in `Dockerfile`, `--build-arg` in `Makefile:94-97`, read the env var in the record). Out of scope here because it changes release tooling rather than the library, and because `builder.git_sha: null` is honest in the meantime. The key is reserved in `schema_version: 1` so shipping it later is not a schema change. **`RELEASE.md:148` should be reworded in the same change**, since it currently asserts an invariant that the `Makefile` does not enforce. |
| Fixing `--num-processes all` (`main.py:145-146`, `fragments_h5.py:1088`) | `all` currently resolves to `None`, which takes the single-process branch — i.e. it runs serially (§4.8). A real bug and a small fix, but it changes runtime behavior, and this document's job is to record behavior, not change it. §4.8 records it accurately so it becomes discoverable; test 9 in §9 pins the current behavior so a fix has to be deliberate. |
| Fixing the SE `num_mapped // 2` contig drop (`fragments_h5.py:1029-1030`, `:1064`) | An SE contig with exactly one mapped read is silently dropped (§3.2). Algorithmic — it changes which data enters the file — so it requires approval and its own evaluation. Named here so it is not rediscovered as a mystery. |
| Refactoring the 17-element worker tuple (`:735`, `:1076-1083`) | Genuinely bad, genuinely unrelated. Attributes are written in the main process (§4.5), so this design does not touch it. Coupling them would make a cheap change expensive. |
| Retro-labeling existing S3 files | Cannot be done honestly — the information does not exist (§8). Writing a *guessed* provenance record into an old file is strictly worse than leaving it absent, because it would be trusted. |
| Removing the `eval` at `:293` | Real security debt (arbitrary code execution from S3-sourced file content) but a separate, riskier change: it must handle files whose `_contig_lengths_str` is a `repr` that is not valid JSON, i.e. a compatibility shim of its own. Flagged as a **follow-up worth doing**. This design at least stops propagating the pattern. |
| Fixing the `is_secondary` gap (`fragment.py:373-377`, `574-577`) | An algorithmic change to which fragments enter the file. Requires explicit approval and its own before/after evaluation. Documented here so the decision is deliberate rather than inherited. |
| Removing dead code in `fragment.py` (`:289-296`, `:371-372`; the unreachable `max_tlen=1000` defaults) | Unrelated cleanup. Noted so a reader of `fragment.py` is not misled into thinking `include_neg_tlen` or `max_tlen=1000` are live. |
| Making `--min-mapq` None-vs-0 meaningful at the filter level (`:771`) | It is already a structural no-op and correct as such: MAPQ is unsigned, so 0 filters nothing. The distinction only needs to survive into *metadata*, which §4.1.1 and §4.4 achieve. Changing filter semantics is an algorithmic change requiring approval. |
| Replacing the MAPQ-255 `assert False` (`:809-812`, `:817-820`) with a raised exception | The `assert` is a correctness guard that `python -O` deletes (§7.1). Converting it to a `raise` is the right fix and is an algorithmic change (it changes which builds succeed). This design instead records `builder.optimization_level` so the record cannot lie about it — a mitigation, not the fix. |
| Recording per-contig GC coverage | A contig missing from the FASTA yields all-NaN GC while `has_gc` stays `True` (§7.4). Detecting it requires information that lives only inside the worker and is not in its return tuple (`:889`), so recording it would require the worker-interface change §4.5 exists to avoid. |
| A `CHANGELOG.md` | The repo has none. Worth adding, unrelated to this. |

---

## 11. Alternatives considered and rejected

| Alternative | Rejected because |
|---|---|
| **~35 individual typed root attributes** (matching `index_block_size` style) | No encoding for `None`, and the key table has twelve nullable keys where `null` is a distinct semantic state (§4.1.1); requested-vs-effective plus the nested `builder`/`parallelism` objects roughly triples the name count; `neutralized` is a list of five-key records and does not flatten; every new flag adds a name every robust reader must `.get()` individually. §4.2. |
| **`str(dict)` + `eval()`**, following the `_contig_lengths_str` precedent (`:1151`, `:293`) | Arbitrary code execution driven by file content, on files routinely pulled from S3 by automated pipelines. `json.loads` is identical effort and cannot execute. Following this precedent would double the blast radius of a bad pattern. |
| **Sidecar `.provenance.json` next to the h5** | Trivially separable from the file it describes. `publishDir` copies, S3 `cp` of a single key, and hand-moved files all lose it. The whole requirement is that the file be self-describing under inspection. Also introduces a "which one wins when they disagree" question. |
| **Record only the full argv** | Records *requested*, never *effective* — exactly the failure mode of §3.3, where a TSV build's argv says `--min-mapq 30` and nothing applied it. Also unavailable for library callers (§4.6), and unparseable without reimplementing argparse and all ten neutralization paths of `fragments_h5.py:954-1002` in the reader. Recorded *in addition to* the structured values, not instead — with its limits (non-resolving per-task work-directory paths, redundancy with `requested`) documented in §4.6 rather than assumed away. |
| **Record `contig_name_map` as a boolean or a hash** | A boolean discards the mapping entirely; `:1044-1045` writes only output names, so the input names exist nowhere else in the file and a consumer cannot get from `chr1` back to `CM000663.2` (§3.1). A hash lets you check two files used the same map but never lets you recover a name, which is the actual use. The map is bounded by the source contig count, which `_contig_lengths_str` already enumerates in full, so recording it costs nothing new. |
| **A new HDF5 group (`/provenance`) with one dataset per field** | Heavier: group + N datasets vs one attribute. Attributes are the established home for file-level metadata here (`:1146-1151`). Datasets invite chunking/compression/dtype questions for scalar strings. Presence-checking a group is no cheaper than presence-checking an attribute. No upside. |
| **Bump a format version and refuse to read old files** | Every production file in `s3://karius-biomarker-data-assets/` becomes unreadable, and `biomarker-pipeline` is pinned to `fragments_h5 2.8.1` (`containers/container-versions.config:48`) while `build_se_fragment_h5s.config:29,69` pins 2.10.1. The hard constraint is backward compatibility. Non-starter. |
| **Reuse / repurpose `max_fragment_length`** to mean the effective cap | Silent-correctness hazard (§7.3): it is the `fetch_array` backward search window (`:449`). Setting it to 120 for an SE-capped file would *appear* to work and would break any file where the cap was later raised or where the histogram size assumption (`:713`) matters. Never overload a load-bearing constant with a semantic it did not have. |
| **A root `format_version` attribute in addition** | §5.3: absent from all existing files, duplicates the existing `has_X` capability-detection model, and invites brittle `>=` gates. Versioning is scoped inside the provenance record instead. |
| **`.attrs.get()` with per-parameter defaults** (strict precedent #1) | The `_source_format` default was *provably* true of all older files; no such default exists here. Would make absent indistinguishable from present-but-default and would present a guess as a fact. §5.1. |

---

## 12. Self-assessment

### 12.1 Round 1 — what the first self-assessment got wrong

The previous version of this section graded the document **B+** and listed four weaknesses. That assessment was not defensible, and the way it failed is more informative than the grade:

- Three of its four listed weaknesses were real but **secondary** (version-string fidelity, Group B-by-version being a judgment call, whether a property should raise). None of them would have caused a wrong implementation.
- The fourth — "the §1.4 open question is genuinely open" — was listed as a weakness but is a **strength**. Refusing to assert that the S3 prefix is clean, when git cannot establish it, is the correct epistemic move. Grading yourself down for honesty inflates the apparent seriousness of the list while costing nothing.
- It **missed both defects that would actually have produced a wrong implementation**: the three mutually contradictory definitions of `effective` (§4.1.1 now), and the incomplete TSV neutralization enumeration with no way to express "warned but not forced" (§3.3 now). Neither was hedged, flagged, or listed. Both were stated confidently and both were wrong.
- Its single most confident-sounding paragraph — the `test_docker_build.py` constraint that supposedly justified §4.7 and mandated a §9.8 test change — was **fabricated from two cited line numbers that were never read in context**. The actual code is a hardcoded two-element allow-list (§4.7). The self-assessment even said "I have not read that test's full comparison logic, only the cited lines", and then let a normative requirement ("This is a required change, not optional") stand on top of that admission. Flagging an unverified premise is not a substitute for verifying it before building on it.

The pattern: the previous self-assessment audited the document's *judgment calls* and did not audit its *factual claims*. Judgment calls were where it felt uncertain; factual claims were where it was wrong.

### 12.2 Round 2 — the same failure recurred, in the material written to fix round 1

§12.1 named that pattern explicitly. The revision that added §12.1 then reproduced it. A second independent review of that revision confirmed all five round-1 fixes were genuine, and found **two new blocking defects, both factual, both in the worked examples added in round 1**:

- **§4.1.5's worked example described a build that cannot execute.** It passed `--contigs chr1 chr2` — output names — against a BAM whose references were RefSeq/GenBank accessions. `fragments_h5.py:1025` keys `all_contig_lengths` by `bam_fp.references`, and `:1028` does `{c: all_contig_lengths[c] for c in allowed_contigs}`, so the build raises `KeyError` and dies ~110 lines before the record is written at `:1139`. The example then explained the resulting `requested`/`effective` mismatch with a sentence — *"the map is what connects them"* — that was false, and that the **very next paragraph of the same section contradicted**, correctly stating that `allowed_contigs` is matched against BAM names at `:1026-1028`. Two adjacent paragraphs asserting incompatible things is not a subtle error; it is what happens when an example is composed rather than traced.
- **§4.1.6's worked example contradicted the normative key table on `effective.parallelism.num_processes`.** §4.1.2/§4.8 said `null` on the serial path; the example said `1` for a build that provably takes the serial path (`main.py:51` default `'1'` → `:148` int `1` → `:1088` `1 != 1` is False → `:1117-1126`); §9's test said `null`. An implementer reading the table and an implementer reading the flagship example would have written different code, with no way to tell which was intended.

Two secondary defects in the same material are worth recording for the same reason: §4.1.6 was headed "three neutralizations" while listing six correct ones, and the `neutralized` array's element order was left unspecified while §4.1.3 and §9 test 8 made byte-identity normative — `sort_keys=True` orders object keys, not array elements, so two conforming implementations could have emitted different bytes.

**The conclusion this forces.** Round 1's diagnosis was "judgment calls were audited, factual claims were not". Round 2 shows that naming a failure mode does not prevent it, and it localizes *where* it survives: **the worked examples are the highest-risk material in this document.** The prose is hedged, cited, and argued; the examples are not argued at all — they are transcribed, and then copied. That asymmetry is exactly backwards from how they will be used. §4 positions §4.1.5-4.1.6 as the normative reference an implementer copies from, so an error there propagates directly into code, while an error in a paragraph of §3.2 gets caught by whoever reads the cited line. Every round-1 and round-2 blocking defect that survived review was in material that *looked* concrete enough not to need checking.

**What an implementer must re-verify in the examples before copying them**, in priority order:

1. **Every contig name in §4.1.5 is an input (BAM) name.** Re-derive `effective.contigs` from `contig_lengths` at `:1035`/`:1018` — not from `contig_lengths_output` (`:1044`), which is what `_contig_lengths_str` holds. If a name in `effective.contigs` also appears as a *value* in `requested.contig_name_map`, that is the round-2 bug returning.
2. **`effective.parallelism` in §4.1.6 against the actual branch taken at `:1088`.** `num_processes` is non-`null` if and only if `multiprocess` is `true`. Both `None` and `1` take the serial branch.
3. **The `neutralized` list in §4.1.6 against the warning sites at `:976-1002`, one by one**, including the two `if <param>:` guards at `:987`/`:989` that do *not* fire on a default-`False` parameter — and confirm the array is sorted by `param`.
4. **Every `effective.*` value against its §4.1.2 source expression**, evaluated in the branch the example's invocation actually takes, rather than against the example's own prose summary.
5. **That each example's invocation would pass `main.py`'s validation and reach `:1139` at all** — the §9 test 4 defect (a TSV invocation missing the `--fasta` that `:972-973` makes mandatory) was the same class of error in the test plan rather than the examples.

The general rule extracted: in this document, a concrete artifact — an example record, a test invocation, a line-number citation — is *less* trustworthy than a paragraph of argument, and should be checked first, not last.

### 12.3 Grade for the current version: B

Unchanged from the previous revision's self-grade, deliberately. Fixing defects that this document introduced does not earn a higher grade than the revision that introduced them; at best it returns the document to where it claimed to already be. Two consecutive review rounds have each found blocking factual defects, and both times the self-assessment had rated the document ready. That is the strongest available evidence about this document's reliability, and it is evidence against, not for.

What the B rests on: the schema is genuinely specified. §4.1.1 gives one normative definition of `requested` / `effective` / `neutralized` with worked `min_mapq` cases and an explicit statement of what `neutralized` does *not* cover; §4.1.2 gives a complete key table with types, nullability, source expressions and capture points; §3.3 enumerates all ten divergence paths with a `disposition` field; §5.3 gives the evolution rules; §9 now covers the contig keys and the reader-tolerance rules that previously had no test. §8 (forensic classification of the existing corpus, actionable today with no code from this design) and §5.1 (absence as a first-class state) remain the strongest sections, and §3.4 remains the most useful thing the revisions added.

**What the grade is contingent on.** It assumes the round-3 corrections in §4.1.5, §4.1.6 and §9 are themselves correct, and they were produced by the same process that produced two rounds of errors in that material — reading source and transcribing, without executing anything. The grade should be revised down if either worked example fails when actually run. **The cheap way to retire this risk, and the thing that should happen before implementation: build the two example files with the real CLI and diff the emitted records against §4.1.5 and §4.1.6.** That converts the highest-risk material in the document from transcription into output, and it is the only step here that would have caught all four blocking defects found so far.

### 12.4 Real residual risks, in descending order

0. **The worked examples may still be wrong.** Two review rounds, two sets of factual defects, both in examples, both stated confidently. Round 3 corrected them by the same means that produced them. This is listed first and separately from the judgment-shaped risks below because §12.1 diagnosed exactly this and §12.2 shows the diagnosis alone did not help. Mitigation is mechanical, not editorial: run the examples (§12.3).
1. **The design's central value proposition is unmeasured.** Nothing in this document establishes that anyone will *read* `_build_provenance`. §5.2's search found zero existing attribute consumers — that is offered there as evidence of safety, but it is equally evidence that nobody in this stack currently reads h5 attributes for anything. Adoption is opt-in (§7.5) and no consumer change is proposed. The realistic failure mode is not a wrong record; it is a correct record that no analysis ever consults, while the operators who would have benefited keep reading the `.nf`.
2. **`builder.version` can be false, and §3.2 depends on it entirely.** §3.4 documents the exact `Makefile` path that produces a mislabeled image and specifies the git-SHA fix, but the fix is deferred (§10). Until it ships, "the version pins Group B" is an argument with a known, reproducible hole. This is the single largest correctness risk in the design as specified, and it is not mitigated — only documented.
3. **The `effective` `null` convention is subtle and will be implemented wrong at least once.** "This parameter governed nothing" is not an obvious reading of `null`, and it differs from `requested`'s `null` ("not supplied") in the same document. The §4.1.2 table gives per-key source expressions precisely because the convention cannot be applied by intuition, and §9's tests 3, 4 and 5 exist to catch specific wrong implementations. It is still a place where a plausible-looking implementation produces a subtly false record — which is the failure mode this document elsewhere calls "worse than none".
4. **`effective.set_mapq_255_to_none: null` for TSV rests on a reachability argument, not a code path.** `:805`/`:813` *are* executed on the TSV path; they cannot fire only because `fragment.py:730-731` constructs every TSV fragment with `mapq1=mapq2=None`. If that ever changes — a TSV format with a MAPQ column, say — the `null` becomes wrong silently, and nothing tests the linkage. This is the one place where the record encodes a conclusion rather than a value.
5. **The record cannot describe what the worker actually did.** §4.5's "no worker-tuple change" is the reason this design is cheap, and the cost is that everything the worker learns is unrecoverable: per-contig GC coverage (§7.4), whether a MAPQ 255 was ever encountered (§7.1), how many TSV rows were dropped (`fragment.py:676-704`). All three are cases where the record says what was *requested* and cannot say what was *observed*. That is an honest limit, but it means `effective` is "the parameters that governed", not "what happened".
6. **`invocation` remains the weakest-justified field.** §4.6 now argues for it explicitly rather than resting on a fictitious test constraint, but the argument is thin — fidelity to untyped operator input, for a human reader — and its recorded paths are non-resolving per-task work-directory paths under the one pipeline that motivated this document (§4.6). A reviewer who wanted to cut it would have a case. It is kept, not defended enthusiastically.
7. **Group B by version is still a judgment call**, unchanged from the previous version and correctly identified there. A hand-maintained enumeration of unconditional filters is a defensible opposite conclusion; the argument against it (a drifted record is trusted and therefore worse than none) is one I hold but not overwhelmingly. Note that §3.4 makes this weaker than it was, since the version string it delegates to is itself not guaranteed.
8. **`README.md`'s attribute block is wrong and this design touches that surface.** §10 defers fixing it while adding two members to the public API that the block purports to document. That deferral is flagged in §10 as one to revisit at implementation time rather than settled, but as written the document leaves the documentation more inconsistent than it found it.

### 12.5 What a reviewer should verify before implementation

Beyond the mechanical check in §12.3 (run both worked examples and diff), the remaining unverified surfaces are all judgment calls rather than facts — which, per §12.2, is the *less* dangerous category and should be looked at second:

- **The `neutralized` `code` vocabulary** (§4.1.2) is invented by this document. It should be sanity-checked against how anyone would actually query it before it is frozen into `schema_version: 1`, because under the additive-only rule (§5.3) the strings cannot be changed afterward.
- **`effective.contigs` before-or-after the zero-mapped skip** (§4.1.2 specifies *before*, at `:1062`). The contigs actually written to `data/` are the post-skip set and are already recoverable from `_contig_lengths_str` — but `_contig_lengths_str` holds *output* names while `effective.contigs` holds *input* names, so the two are not directly comparable. Someone should confirm the pre-skip input-name set is the more useful of the two, since it is not obvious.
- **Whether `requested.fasta_given: bool` is the right reduction** (§4.1.2). It drops the FASTA path on the grounds that it is machine-local and appears in `invocation`. For a library caller `invocation` is `null`, so for that caller the FASTA identity is unrecorded entirely. Reference identity is exactly the kind of thing provenance is for; this may be a reduction too far.
