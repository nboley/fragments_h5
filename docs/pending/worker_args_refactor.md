# Worker Argument Structure Refactor

**Status:** PENDING — proposal. **Scope:** replace the 17-element positional tuple consumed by `build_sub_fragments_h5` with a module-scope, keyword-constructed dataclass. Verified against tag `v2.12.1` (`b44bd7d`); all line numbers below are pinned to that tag, not to a moving HEAD.

Every empirical claim is marked **EXECUTED** (observed in this session unless stated) or **INFERRED**. Unrun code is labelled `# PROPOSED`.

---

## 1. Problem statement

`build_sub_fragments_h5(args)` (`src/fragments_h5/fragments_h5.py:738`) takes one positional parameter and destructures it on a single line, `:750`, into 17 names. The tuple is built at exactly one pack site, `:1102-1109`, inside a contig × chunk loop in `build_fragments_h5`.

**There is a third positional accessor, and it is the one that has actually broken.** At `:1120`:

```python
total_bases = sum(a[4] - a[3] for a in args)  # chunk_stop - chunk_start
```

This is a raw-index read of the same tuples, in the same function, 370 lines away from the unpack (`:1120` − `:750`). Any design that only reconciles pack-vs-unpack does not address the site with the demonstrated failure history.

### 1.1 Why positional is unsafe here (EXECUTED — read of `:1102-1109`)

Positional order and type: 1 `input_fname` str · 2 `bam_contig` str · 3 `output_contig` str · 4 `chunk_start` int · 5 `chunk_stop` int · 6 `contig_length` int · 7 `fasta_filename` str|None · 8 `single_end` bool · 9 `se_max_fragment_length` int|None · 10 `read_gc` bool · 11 `read_strand` bool · 12 `read_methyl` bool · 13 `set_mapq_255_to_none` bool · 14 `include_duplicates` bool · 15 `store_fragment_end_clipped` bool · 16 `tmp_dir_name` str · 17 `min_mapq` int|None.

Positions 10–15 are **six consecutive bools**; 1–3 three consecutive strs; 4–6 three consecutive ints. Any transposition inside a run type-checks and cannot raise at the boundary. It surfaces, if at all, as wrong data in the output h5.

### 1.2 The released defect (EXECUTED — the central evidence)

The tuple grew by insertion, and the derived accessor at `:1120` drifted:

| Commit | Date | Change | `total_bases` |
|---|---|---|---|
| `9d4f72a` | 2025-11-21 | 8-element tuple introduced via `pool.map` (same day as `e24ff97`, which had 7 *separate* parameters — no tuple, no `Pool` — one of which, `input_to_fragments`, never made it into the tuple) | — |
| `95c76f5` | 2026-03-09 (v2.8.0) | inserted `chunk_start, chunk_stop, contig_len` → 14 | correctly updated to `a[3] - a[2]` |
| `5c84100` | 2026-06-08 | inserted `output_contig` at index 2 → 15 | **not updated; still `a[3] - a[2]`** |
| `e876cfa` | 2026-07-09 | appended `se_max_fragment_length`, `min_mapq` → 17 | fixed |

At `5c84100`, `a[2]` is `output_contig` (str) and `a[3]` is `chunk_start` (int). Checked out `5c84100` into `/tmp` (copying the prebuilt `sequence*.so`, unchanged across the range) and ran `build_fragments_h5('tests/data/small.chr6.bam', out, num_processes=N)`:

- `num_processes=2` → `TypeError: unsupported operand type(s) for -: 'int' and 'str'`
- `num_processes=1` → same `TypeError`

**Every code path was fatally broken**, not an edge case, for 5 commits spanning 2026-06-08 to 2026-07-09.

**It was released.** `git merge-base --is-ancestor 5c84100 v2.10.1` → yes. `git show v2.10.1:src/fragments_h5/fragments_h5.py` line 1063 → `a[3] - a[2]`. Checking out tag `v2.10.1` into `/tmp` and running the same build → `TypeError` at both `num_processes=1` and `num_processes=4`. **v2.10.1 is the tag pinned by the production Nextflow config.**

**Nuance, stated to avoid overclaiming (INFERRED — from a prior session's verification, not re-executed here):** a `docker pull ghcr.io/nboley/fragments-h5:2.10.1` was found to have been built from an ahead-of-tag working tree, ~3 minutes after the tag, and does work; a corpus of 48 SE h5s was built with it. So production escaped the defect *by accident of the image not matching its label*. That is not a mitigation — it is this repo's signature failure mode (a label disagreeing with its contents) appearing a second time in the same incident.

Two conclusions follow, and they set the design's acceptance criteria:

1. **Pack and unpack were never the problem.** Both were correctly updated at every single insertion, including `5c84100`. The defect came from a *third* accessor that nobody was looking at. A fix must eliminate positional access at every site — the three in `src/` and the two in `tests/` enumerated in §2 — and must make a further site impossible to add silently.
2. **"The tests pass" is weak evidence in this repo.** A totally-fatal defect reached a released, production-pinned tag.

### 1.3 This work is already anticipated

`docs/pending/build_provenance_metadata.md:518-520` — "The worker tuple is a real maintenance hazard and may deserve its own design; it is not coupled to this one." `docs/architecture/fragment_selection_and_build_provenance.md:637-653` — "The tuple remains the acknowledged root cause of this bug class and has its own deferred design; it is out of scope here." `COORDINATION.resolve-uncommitted-cli-flags.md` records the explicit prohibition given to a prior agent ("no worker-tuple refactor") and an EM verification note at `a3be97b` confirming it was honoured. The deferral is documented fact, not folklore.

---

## 2. Scope / Out of scope

**In scope:** the shape of the data passed from `build_fragments_h5` to `build_sub_fragments_h5`; the positional read/write sites — pack `:1102-1109`, unpack `:750`, derived accessor `:1120` — plus the two test monkeypatch spies that read the packed tuple positionally, `tests/test_fragment_selection.py:353` and `:453` (required changes, not optional — see §4 step 6). `:1157` is a pass-through consumer (it receives one packed object and forwards it, reading nothing positionally) and needs confirming but is not itself a positional-access site.

**Out of scope, explicitly:**

- **`max_tlen` for single-end.** Declared at `fragment.py:559`, never referenced in the SE body — dead, but the shared call at `:788-793` passes it unconditionally and `single_end_bam_to_fragments` has no `**kwargs`, so removing it raises `TypeError` on every SE build. This was proven by execution previously and the removal retracted. **Hazard note only: do not touch it.**
- **The 4-element return tuple** `(output_contig, chunk_start, chunk_stop, ofname_or_None)` (`:875`/`:905`, unpacked at `:1143` and `:1157`). It is a genuine second positional shape and deserves the same treatment. It is excluded here because it is only 4 elements, is unpacked immediately at both consumers with no third derived accessor, and has no defect history. Folding it in doubles the diff for the weaker half of the problem. If pursued, it must be a **separate later commit**, not a stage of this one.
- `se_max_fragment_length` semantics (live at `fragment.py:591` gating a hard-limit `ValueError`; the actual SE length filter is separately at `fragments_h5.py:798` — an older note calling it dead was wrong).
- The S3 `abspath` bug; the missing `--contig-name-map` test coverage; annotating the rest of `fragments_h5.py`; the `build-fragments-h5` PATH errors in the test env; `pytest-timeout` not being installed; any cross-repo change.

**API surface (EXECUTED — grep):** `build_sub_fragments_h5` appears only at `:738` (def), `:1139`, `:1157`, plus monkeypatch spies in `tests/test_fragment_selection.py:353` and `:453` (which unconditionally read the packed argument positionally, `args[1]  # bam_contig` — must be updated, see §4 step 6) and a stale line-number mention in a docstring at `tests/test_cli_validation.py:248`. It is not in `__init__.py` (`__all__ = ["FragmentsH5"]`) and not in `main.py`. **Zero external callers — the signature change breaks no public API.** The test sites above are therefore part of this change's scope, alongside `src/`. Precedent for an internal-only change of this kind (`a5fe651`/`452b214`): the change commit itself carries no version bump; `pyproject.toml` stayed at `2.12.0` through both, and a separate `Release 2.12.1` commit (`b44bd7d`) bumped it afterward. See §4.

---

## 3. Representation

**Decision: a module-scope `@dataclass(frozen=True, slots=True)` named `SubBuildArgs`, constructed with keywords only.**

### 3.1 Why not `NamedTuple`

A `NamedTuple` is a tuple subclass. Under it, `a[4] - a[3]` at `:1120` and the 17-name unpack at `:750` **both keep working unchanged**. For a refactor whose entire purpose is eliminating positional access, that backwards compatibility is a liability: it silently preserves the exact accessor that produced the released v2.10.1 defect, and it lets a future insertion re-introduce the identical bug with the identical mechanism. A `NamedTuple` would make the code *look* fixed while leaving the failure reachable.

A dataclass makes `args[3]` a `TypeError` at first execution. Every positional site *must* be converted for the code to run at all. Omission becomes impossible rather than merely unlikely — and this is not hypothetical: it is exactly how the two test spies at `tests/test_fragment_selection.py:353`/`:453` were found (§4 step 6), by the prototype failing loudly rather than by anyone grepping for them.

Pickle size: **EXECUTED**, 17 fields, module-scope classes, protocol 5. Report *deltas against a plain tuple*, not absolutes — absolute sizes scale with the payload and are not reproducible across measurements (two independent passes on short vs. realistic path values gave tuple=48 B and tuple=128 B). The deltas are payload-invariant:

| Representation | Δ vs plain tuple |
|---|---|
| `NamedTuple` | +20 B |
| `dataclass` / `dataclass(frozen=True)` | ≈ +301 B |
| `dataclass(slots=True)` | ≈ +305 B |
| **`dataclass(frozen=True, slots=True)`** (the shape proposed in §3.3) | **+25 B** |

`frozen=True, slots=True` generates a `__getstate__`/`__setstate__` pair that pickles a plain tuple of values (confirmed: `"__getstate__" in DCFS.__dict__` is `True`, `DCS` is `False`), so it lands 5 bytes above a `NamedTuple` — not the ~300-byte cost that applies to the unslotted or unfrozen forms. There is no pickle-size cost to this decision.

### 3.2 Constraints the representation must satisfy (all EXECUTED)

- **Pickling through `imap_unordered`, tested under both `fork` and `spawn`:** plain tuple OK; module-scope `NamedTuple` OK; `@dataclass` OK; `frozen=True` OK; `slots=True` OK. **Function-local `NamedTuple` and function-local dataclass BOTH FAIL** under both start methods: `AttributeError: Can't pickle local object '...<locals>.X'`. → **the struct must be defined at module scope**, not inside `build_fragments_h5`.
- **Python floor:** `pyproject.toml` declares `requires-python = ">=3.10"`; `biomarker_env` is 3.10.19. `slots=True` is available. Dependencies are only `numpy, h5py, pysam, tqdm` (dev `pytest, pytest-timeout`) — **no pydantic, no attrs**, and none should be added.
- **Fork safety:** start method is an explicit `multiprocessing.get_context('fork')` at `:1132` (the `forkserver`→`fork` switch landed in `ebfb7fc`, whose commit message reads "Added s3 support"). The parent's `pysam` handles at `:1023`/`:1041` are inside `with` blocks that close before `args` and the Pool are built at `:1083+`; the worker opens its own FASTA at `:761`; `_pool_worker_init` (`:923`) only ignores SIGINT; module globals are constants and are never mutated (one, `methyl_keys` at `:168`, is a mutable list but is never written). **No hidden fork dependency — a plain-data struct does not touch this boundary.**
- **`chdir`:** the worker body runs inside `with _temporary_working_directory():` (`:759`, defined `:908-920`), which `os.chdir`s into a fresh temp dir. All paths must already be absolute at pack time (`tmp_dir_name` carries the comment "absolute path from main process" at `:877`). The struct changes nothing here, and **the design must not move any path resolution into the worker.**

### 3.3 Shape

```python
# PROPOSED — but this exact shape has been prototyped and executed; see §7
@dataclass(frozen=True, slots=True)
class SubBuildArgs:
    input_fname: str
    bam_contig: str
    output_contig: str
    chunk_start: int
    chunk_stop: int
    contig_length: int
    fasta_filename: Optional[str]
    single_end: bool
    se_max_fragment_length: Optional[int]
    read_gc: bool
    read_strand: bool
    read_methyl: bool
    set_mapq_255_to_none: bool
    include_duplicates: bool
    store_fragment_end_clipped: bool
    tmp_dir_name: str
    min_mapq: Optional[int]
```

`frozen=True` because the worker never mutates any of these and freezing removes a class of aliasing question across the fork; `slots=True` because it is free at this Python floor.

Field order mirrors the current tuple exactly. This is deliberate: with keyword-only construction the order is semantically irrelevant, and preserving it keeps the pack-site diff a pure `(` → `SubBuildArgs(` + `name=` annotation, reviewable line-by-line against the old tuple.

**Style note, stated plainly (corrected — the class names were transposed in an earlier draft):** `fragment.py` already uses `@dataclass(slots=True, frozen=True)` for `MethylCounts` (`:26`) — genuinely frozen and slotted, confirmed by inspection (`frozen=True slots_attr=True dict=False`) — so `frozen=True, slots=True` is idiomatic **for the package**. `Fragment` (`:154`) uses plain `@dataclass` with a manually declared `__slots__` list (`:155`), not `slots=True`, and is **not** frozen (`frozen=False`) — it cannot be, since the worker mutates it at `:823` (`fragment.mapq1 = None`), under the file's own TODO at `:819-820`. But `fragments_h5.py` — where all of this code lives — has essentially no type annotations and no dataclass/`NamedTuple` usage at all. This is **novel for this file**. Do not claim otherwise in review.

---

## 4. Staging

### Stage 1 — the whole positional-access removal, one commit

The commit must remove *all* positional access, because a partial conversion leaves exactly the split-brain state that caused the v2.10.1 defect.

1. Define `SubBuildArgs` at module scope in `fragments_h5.py`, above `build_sub_fragments_h5`.
2. Pack site `:1102-1109`: `args.append(SubBuildArgs(input_fname=..., bam_contig=..., ...))` — all 17 keyword.
3. Worker `:750`: replace the one 17-name unpack with **17 one-to-one rebinds**:
   ```python
   # PROPOSED — but prototyped and executed; see §7
   input_fname = args.input_fname
   bam_contig = args.bam_contig
   ...
   ```
   The remaining ~150 lines of the worker body are then **untouched**, including the dispatch at `:752-757`, the shared generator call at `:788-793`, the `_se_kwargs` conditional at `:787`, and the `min_mapq` normalisation at `:786`. The diff is mechanical and locally verifiable.
4. `:1120`: `total_bases = sum(a.chunk_stop - a.chunk_start for a in args)`.
5. Serial path: the `else` branch at `:1151-1160` (guarded by `if num_processes is not None and num_processes != 1:` at `:1122`) calls `build_sub_fragments_h5(arg)` at `:1157`. No signature change is needed under this design, so this is a no-op — **but it must be explicitly confirmed in review**, because the second consumer of the packed objects is easy to forget and it is the path most tests actually exercise (see §5.2).
6. **Required, not conditional:** update the two monkeypatch spies in `tests/test_fragment_selection.py:353` and `:453` (`called_contigs.append(args[1])  # bam_contig` → `args.bam_contig`). Both inspect the packed argument positionally, unconditionally. Verified against a prototype: before the fix, `test_contig_with_zero_mapped_reads_skipped` and `test_pe_contig_with_single_mapped_read_skipped` fail with `TypeError: 'SubBuildArgs' object is not subscriptable` (2 failed, 6 passed); after `args[1]` → `args.bam_contig` at both lines, the full suite matches baseline exactly.

### Stage 2 — inlining `args.field` and deleting the 17 rebinds: **dropped, not deferred**

Replacing every use of `bam_contig` with `args.bam_contig` through the 150-line body touches dozens of lines for zero safety gain — after Stage 1 the positional hazard is already gone, and the rebind block *is* the readable name→field mapping. It converts a mechanically reviewable diff into a large one in the exact function with a released-defect history. Recommend dropping it. If someone wants it later it is pure cosmetics and should be argued on its own.

**Version:** patch bump (2.12.1 → 2.12.2), consistent with the `a5fe651`/`452b214` precedent for internal-only changes — but note what that precedent actually shows: neither commit touched `pyproject.toml` (it stayed at `2.12.0` through both); the bump to `2.12.1` happened in a separate, later `Release` commit (`b44bd7d`). So this change's own commit should not bump the version either — that happens in its own `Release` commit afterward. `make tag` (Makefile:114-133) refuses to tag with `pyproject.toml` uncommitted; `pyproject.toml` is the single source of version truth. No CHANGELOG.md — add a dated section to `RELEASE.md` and `AGENT_CONTEXT.md`.

---

## 5. Testing strategy

"Tests pass" is explicitly insufficient here (§1.2). Three independent checks, in increasing order of strength.

### 5.1 Byte-exact output comparison (the strong check)

**EXECUTED — byte-for-byte comparison is valid in this repo.** Built `small.chr6.bam` + the chr6 FASTA three times through the real API: twice at `num_processes=1`, once at `num_processes=4`. All three were **md5-identical and `cmp`-identical**, `a2b806bdcc52e8c5ef9361c38d326d80`, with zero attribute differences. Verified mechanism: `h5py.File(...).libver` is `('earliest', 'v114')`, and HDF5's earliest object-header format (v1) has no ctime/mtime field at all — timestamps appear only with v2 headers / `libver='latest'`.

**Caveat, must be stated in the PR:** this holds only while `libver` stays at the default and no differing provenance kwargs are passed. If either changes, this check silently degrades to a no-op.

So: **a before/after md5 matrix**, built at `v2.12.1` and at the refactor commit, over the available fixtures:

| Fixture | Config |
|---|---|
| `tests/data/small.chr6.bam` (198736 B, 2600 mapped, 604 dup) | np=1, np=4; GC off |
| same + `GRCh38.p12.genome.chr6_99110000_99130000.fa.gz` | np=1, np=4; GC on |
| `scATAC_breast_v1_chr6_99118615_99121634.hg38.bam` (36K) | np=1, np=4 |
| synthetic SE BAM built in-test with pysam | np=1, np=4 |
| `test_duplicates.bam` (354B), `include_duplicates` both ways | np=1 |

Every cell must match its `v2.12.1` counterpart exactly.

**Honest coverage gaps (EXECUTED — fixture and test survey):**
- There is **no stored SE fixture**; `--single-end` is only exercised via in-test pysam-built BAMs.
- **`--contig-name-map` has ZERO test coverage anywhere.** It is the flag that makes `output_contig != bam_contig`. Since `output_contig` is one of the 17 fields — and is the very field whose insertion caused the v2.10.1 defect — the refactor touches genuinely untested territory. The md5 matrix cannot close this. Reviewers should treat `output_contig` handling as verified only by inspection plus §5.3.
- `--contigs` appears inside the `target_h5_path` fixture (the 6 erroring tests) and 4× in `tests/specialized/compare_chunked_vs_unchunked.py` (standalone, not pytest-collected). It has never successfully executed **in the collected suite**.

### 5.2 Argument-capture diff (proves the mapping independently of output) — EXECUTED

**Method:** two `/tmp` copies of `src` (never the repo) — one unmodified `v2.12.1`, one with `SubBuildArgs` exactly as §3.3/§4 specify — each instrumented at the pack site to dump every packed argument set as a sorted `name → repr(value)` mapping (`tmp_dir_name` redacted). Ran `build_fragments_h5` over `small.chr6.bam` via `PYTHONPATH=/tmp/<copy>/src` (verified to take precedence over the editable `.pth`) at `num_processes=1` and `4`, both with and without the chr6 FASTA — 4 configs × 2 trees, 8 runs total.

**Result:** all 8 runs succeeded; each config produced 18 packed argument sets (18 genomic chunks). `diff` between the old-tree and new-tree dumps was **empty on all 4 configs** (np=1/no-FASTA, np=1/FASTA, np=4/no-FASTA, np=4/FASTA) — the name→value mapping is byte-identical before and after the refactor.

### 5.3 Mutation battery (proves the defect class is unreachable)

**Baseline (EXECUTED):** `65 passed, 0 failed, 6 errors`. All 6 errors are the module-scoped `target_h5_path` fixture calling `subprocess.run("build-fragments-h5 ...", shell=True, check=True)` and getting `/bin/sh: 1: build-fragments-h5: not found`, exit 127 — an environment PATH problem, verified by reading real stderr, not a code defect.

**Baseline mutation (EXECUTED):** in a `/tmp` copy only, swapped the *order* of `include_duplicates` and `store_fragment_end_clipped` at the pack site without touching the unpack. Full suite: **`4 failed, 61 passed, 6 errors`** (same 6 PATH errors). Caught by `test_read_all`, `test_fetch`, `test_include_duplicates`, `test_fragment_end_clipped_storage_and_read`, all in `test_fragments_h5.py`. Mechanism: the swap feeds `store_fragment_end_clipped`'s default `True` into the `include_duplicates` slot, so duplicates get included, cascading into fragment-identity mismatches even in tests that never mention `include_duplicates`. **The multiprocessing-specific tests did not catch it.** So the suite catches this transposition *incidentally*, via unrelated fixture defaults — not because any test verifies argument order.

**Wall-clock, not reproducible:** an earlier pass reported `574.91s`/`492.59s`; independent re-measurement on the same host/env gave `388.51s`/`384.62s`. Pass/fail counts reproduced exactly; the timings did not. Treat the full-suite runtime as "roughly six minutes," not a fixed figure.

Post-refactor battery, predicted at draft time, `# PROPOSED`; rows marked **verified** were confirmed against an independent reviewer's prototype (§7):

| Mutation | Expected |
|---|---|
| Reorder two fields in the `SubBuildArgs` class body | **No behavioural change** — keyword construction is order-independent; output md5 unchanged. Demonstrates the insertion-order hazard is gone. **Verified.** |
| Delete one `name=` from the pack site | `TypeError: missing required argument` at first pack. Immediate, total, unmissable. **Verified.** |
| Restore `total_bases = sum(a[4] - a[3] ...)` | `TypeError: 'SubBuildArgs' object is not subscriptable`. **This is the v2.10.1 defect made unreachable-by-construction. Verified.** |
| Insert a new required field mid-class (no default) and not update any accessor | `TypeError: SubBuildArgs.__init__() missing 1 required positional argument` at first pack — loud and immediate. **Verified.** |
| Same, but the new field has a default | `OK`, no failure. This is the only silent case, and it requires a default (a defaulted field must also trail the non-default ones, or `TypeError: non-default argument follows default argument` at class-definition time). Strictly better than the tuple: `5c84100`'s failure mode (silent, positional) requires no default and no keyword to reproduce; here it requires deliberately opting into both. |
| Swap the *values* on two same-typed keywords, e.g. `include_duplicates=store_fragment_end_clipped` | Still catchable only by output tests; expect the same 4 failures as the baseline mutation. |

**Be precise about what remains possible.** Keyword construction eliminates order-dependence and makes omission fatal. It does **not** prevent binding the wrong *value* to a correctly-named field. That residual risk is bounded by §5.1/§5.2 and by the fact that the pack site is a single location with 17 lines of `name=value`, where the name and the source expression are adjacent and usually identical.

**Note (EXECUTED):** `pytest-timeout` is **not installed** in this environment (`Unknown config option: timeout`), so every `@pytest.mark.timeout(N)` on the parallel/chunking tests is silently inert. A refactor of the parallel path is therefore not protected from a hang by those markers — run the np=4 cases with an external wall-clock bound. Neither `mutmut` nor `cosmic-ray` is installed; the battery above is manual by necessity.

### 5.4 Where to run it

**EXECUTED:** the editable install is a plain `.pth` — `site-packages/__editable__.fragments_h5-2.11.0.pth` containing the single line `/home/nathanboley/src/fragments_h5/src`. `PYTHONPATH=/tmp/<copy>/src` takes precedence over it (proved with a marker print that appeared against the copy and 0 times against the repo). So all mutation work happens in `/tmp` copies with **zero repo modification**. (The dist-info is stale: it reports 2.11.0 while `pyproject.toml` declares 2.12.1 — cosmetic, not a functional issue.)

Existing diff tooling is **not** sufficient for §5.1 and should not be reused as-is: `tests/specialized/compare_chunked_vs_unchunked.py` (standalone, not pytest-collected) compares only two named root attrs (`index_block_size`, `max_fragment_length`), not all attrs; `tests/test_docker_build.py` has a near-duplicate `compare_h5_files()`. Plain `md5sum` is strictly stronger than both and needs no new code.

---

## 6. Alternatives considered and rejected

**`NamedTuple`.** Rejected — §3.1. The pickle-size argument for it **does not exist**: the proposed `frozen=True, slots=True` dataclass pickles 5 bytes above a `NamedTuple`, payload-independent (measured). Its only remaining property is that it preserves the exact accessor pattern that shipped a fatal defect — a cost, not a benefit.

**Keep the tuple, add a comment / assert `len(args) == 17`.** Rejected. `5c84100` kept pack and unpack in sync; the length assert would have passed. It does not address positional *meaning*, only positional *count*.

**Multi-argument worker via `pool.starmap` / `starmap_async`.** Rejected on measured grounds. **EXECUTED**, 5 tasks with delays `[0.6, 0.1, 0.1, 0.1, 0.1]` submitted in index order: `imap_unordered` with a single tuple arg returned `1,2,3,4,0` (lazy, unordered); `starmap` blocked the full 0.618 s and returned in submission order (eager, ordered); `starmap_async().get()` returned from submission in 0.012 s but `.get()` blocked 0.619 s — all-or-nothing, no partial delivery. The full public method list of `multiprocessing.pool.Pool` is `apply, apply_async, close, imap, imap_unordered, join, map, map_async, starmap, starmap_async, terminate` — **there is no lazy + unordered multi-argument method.** The tqdm bar at `:1141`/`:1155` consumes `total=total_bases` incrementally, which is exactly what streaming buys. The one-argument shape is genuinely forced.

**`functools.partial`-bound invariant config + a small per-chunk struct.** Real and tempting: ~12 of the 17 fields are invariant across every task, and only 5 vary per chunk (`bam_contig`, `output_contig`, `chunk_start`, `chunk_stop`, `contig_length`). **EXECUTED:** `imap_unordered` + `functools.partial` works, with the same streaming/unordered behaviour. Rejected anyway, on three grounds: (a) it changes the worker's *signature*, which forces a corresponding change on the serial path at `:1157` too, widening the diff in the riskiest function; (b) it creates **two** structures where one suffices, and a future field then needs a judgement call about which one it belongs in — a new way to get it wrong, in a refactor whose purpose is removing ways to get it wrong; (c) the only benefit is pickle bytes, already measured as negligible. Against the repo's stated preference for minimum complexity, one flat struct wins.

**Fold in the 4-element return tuple.** Rejected for *this* commit — §2. Doubles the diff for the weaker half of the problem. Note that the deferred merge of worker outputs into the output h5 (comment at `:1078-1080`, done only after all workers finish, for HDF5 fork-safety) is unaffected either way.

**`pydantic` / `attrs` validation model.** Rejected: a new runtime dependency for a 17-field internal plain-data carrier, in a package whose entire dependency set is `numpy, h5py, pysam, tqdm`. `frozen=True` + type annotations give the useful part for free.

---

## 7. Self-assessment

**Grade: A-** (revised from the first draft's B+, credited to independent review, not to this doc's own execution.)

The B+ rating's stated reason — "the design itself has not been executed" — no longer holds. An independent reviewer prototyped `SubBuildArgs` exactly as specified (module scope, `frozen=True, slots=True`, 17 keywords at the pack site, 17 rebinds at `:750`, `a.chunk_stop - a.chunk_start` at `:1120`) and got builds succeeding at np=1/2/4 that were **md5-identical to v2.12.1 on all six §5.1 matrix cells**, and the full suite at `65 passed, 6 errors` — identical to baseline — once the two test spies (§4 step 6) were fixed. All of §5.3's post-refactor mutation rows with a concrete predicted error message (delete a `name=`, restore the `:1120` subscript, insert a required field with no default) hit that message verbatim. §5.2 (argument-capture diff) was also run, in this pass — see above — and came out clean.

Not raised to A: the same independent review found one High-severity defect in the doc's own measurements (§3.1/§6, a pickle-size number for a representation the design does not propose), a second High (§2/§4, a mandatory test fix mis-labelled optional, with wrong line numbers), and five Medium defects (swapped class names, a fabricated version-bump precedent, an internally inconsistent site count, an imprecise mutation-battery row, non-reproducible timings presented as measured). None affect the design's correctness, but see the new Weakest-points item below.

**Confident about:**
- The problem is real, and the specific mechanism is proven: three positional accessors, and the *third* one — not the pack/unpack pair — produced a `TypeError` on every code path in released tag `v2.10.1`, which is production-pinned. Executed at both `5c84100` and `v2.10.1`.
- The dataclass-over-`NamedTuple` choice, on both grounds now: `NamedTuple` keeps `a[4] - a[3]` working (that accessor is the defect), and it no longer even has a pickle-size advantage to trade off against.
- The multiprocessing constraints. Single-argument callable is forced (measured across all `Pool` methods); the struct must be module scope (function-local pickling fails under both `fork` and `spawn`); no fork-safety interaction.
- Byte-exact h5 comparison is a valid verification tool here, with a verified mechanism (`libver=('earliest','v114')`, v1 object headers carry no timestamp) and a stated caveat.
- The design as specified builds successfully and reproduces `v2.12.1` output exactly — independently prototyped and confirmed, not merely argued.

**Weakest points:**
- **`--contig-name-map` has zero coverage anywhere**, so `output_contig != bam_contig` is untested — and `output_contig` is precisely the field whose insertion caused the released defect. The md5 matrix cannot close this gap; §5.2's argument-capture diff (now executed) used the default `contig_name_map=None` and so does not exercise this case either. This is the single largest residual risk, and it is still open.
- The first draft of this document had four claims labelled EXECUTED that measured the wrong thing (pickle numbers for an unproposed representation, a fabricated version-bump precedent, transposed class names, an internally inconsistent accessor count). Worse, the pickle measurement was then re-run twice more and produced *three different sets of absolute byte counts* — all correct, all EXECUTED, all incomparable, because the size scales with the payload. §3.1 now reports payload-invariant deltas instead. The lesson generalises: an EXECUTED label certifies that a command was run, not that it measured the quantity the argument needs. Re-verify independently rather than trusting the label — which is the same failure mode §1.2 documents in the v2.10.1 incident, reproduced inside this document.
- The claim that the deployed `2.10.1` image escaped the defect because it was built from an ahead-of-tag tree is **INFERRED from a prior session's verification, not re-executed here.** If it matters to a decision, re-pull and re-check.
- The Stage-2 "drop the inlining" recommendation is a judgement call, not an empirical result. A reviewer who values uniform access style over diff size could reasonably disagree.
- The excluded return tuple leaves a second positional shape in the file. Defensible, but it means the file is not positional-free after this work.
- `pytest-timeout` being inert means the np=4 verification runs unprotected against hangs; the external wall-clock bound in §5.3 is a workaround, not a fix.

**What a reviewer should verify:**
1. That the pack site after conversion has exactly 17 `name=` keywords and that each name's source expression is the same expression that occupied that tuple position at `v2.12.1` — line-by-line against `git show v2.12.1:...` lines `1102-1109`. This is the one check no automated test performs.
2. That `:1120` reads `a.chunk_stop - a.chunk_start` and that **no other subscript of any element of `args`** remains anywhere in `src/` **or `tests/`** — `grep` for `a\[` and `args\[` in both, not just inside `build_fragments_h5`. This is what catches the two test spies at `tests/test_fragment_selection.py:353`/`:453`.
3. That the serial branch at `:1151-1160` was actually exercised, not just read. Most of the test suite runs `num_processes=1`, so a break there is loud — but confirm it ran.
4. That `SubBuildArgs` is at module scope, not nested inside `build_fragments_h5`. Nested definitions pickle-fail under both start methods, and the failure would only appear on the `num_processes > 1` path.
5. That `max_tlen` is still passed unconditionally in the shared generator call at `:788-793`. Its removal is a known-fatal SE regression that has already been made and retracted once.
6. That no path resolution moved into the worker, which runs under `os.chdir` into a temp dir.
7. The md5 matrix in §5.1 — reproduced, not quoted. And confirm `libver` is still `('earliest', 'v114')` at the time it is run, or the check is silently vacuous.
