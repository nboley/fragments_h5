"""Console script repair-fragments-h5-gc.

Repairs two defects in fragment H5 files:
1. GC float32 cumsum saturation (commit 95c76f5, fixed 2026-03-09)
2. Trailing phantom padding row (commit 778f4d1, fixed 2025-12-17)

See docs/pending/gc_repair_tool.md for the full design document.
"""

import argparse
import hashlib
import importlib.metadata
import json
import logging
import os
import struct
import sys
import tempfile
from datetime import datetime, timezone

import h5py
import numpy
import pysam

from fragments_h5.fragments_h5 import (
    FragmentsH5,
    INDEX_BLOCK_SIZE,
    MIN_NUM_READS_FOR_INDEX,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

T23 = 2**23  # 8,388,608 — float32 exact to half-integers below this
T24 = 2**24  # 16,777,216 — float32 fully saturated above this

ALLOWED_SEQ_BYTES = frozenset(ord(c) for c in "ACGTNacgtn")

# Max single-part S3 upload size (5 GB)
MAX_UPLOAD_BYTES = 5 * 1024**3


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------

class RepairAbort(Exception):
    """A file-level abort — the file cannot be repaired safely."""
    pass


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def compute_file_sha256(path):
    """SHA-256 hex digest of a file's bytes."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_attrs_excluding(h5_file, exclude=("_repair_history",)):
    """Deterministic SHA-256 of all root-level HDF5 attrs, excluding named keys.

    Used to verify attrs are unchanged after repair (design §3.1 step 2).
    """
    h = hashlib.sha256()
    for name in sorted(h5_file.attrs.keys()):
        if name in exclude:
            continue
        val = h5_file.attrs[name]
        h.update(name.encode("utf-8"))
        h.update(b"\x00")
        if isinstance(val, numpy.ndarray):
            h.update(val.dtype.str.encode("utf-8"))
            h.update(val.tobytes())
        elif isinstance(val, (numpy.integer, int)):
            h.update(struct.pack("<q", int(val)))
        elif isinstance(val, (numpy.floating, float)):
            h.update(struct.pack("<d", float(val)))
        elif isinstance(val, (bytes, numpy.bytes_)):
            h.update(bytes(val))
        else:
            h.update(str(val).encode("utf-8"))
        h.update(b"\x01")
    return h.hexdigest()


def check_gc_presence(h5_data_group):
    """Check GC dataset presence across all contig groups.

    Returns 'all', 'none', or 'partial'.
    Partial is a blocking abort per §3.2.5.
    """
    has = [
        "gc" in h5_data_group[contig]
        for contig in h5_data_group
    ]
    if all(has):
        return "all"
    elif not any(has):
        return "none"
    else:
        return "partial"


# ---------------------------------------------------------------------------
# FASTA preflight
# ---------------------------------------------------------------------------

def scan_fasta_alphabet(fasta_path):
    """Scan FASTA sequence bytes for the §5 layer-5 alphabet gate.

    Returns a dict mapping uppercase character -> count.
    Raises ValueError if any byte outside {A,C,G,T,N} (case-insensitive) is found.
    """
    total_counts = numpy.zeros(256, dtype=numpy.int64)
    with pysam.FastaFile(fasta_path) as fasta:
        for contig in fasta.references:
            seq = fasta.fetch(contig)
            arr = numpy.frombuffer(seq.encode("ascii"), dtype=numpy.uint8)
            total_counts += numpy.bincount(arr, minlength=256)

    result = {}
    disallowed = {}
    for i in range(256):
        if total_counts[i] > 0:
            c = chr(i)
            cu = c.upper()
            result[cu] = result.get(cu, 0) + int(total_counts[i])
            if i not in ALLOWED_SEQ_BYTES:
                disallowed[c] = int(total_counts[i])

    if disallowed:
        raise ValueError(
            f"FASTA contains characters outside {{A,C,G,T,N}} — aborting entire run: {disallowed}"
        )
    return result


def build_cumsum_cache(fasta_path, cache_dir, fasta_sha256):
    """Build per-contig float64 cumsum arrays and cache them.

    Each file is written atomically via os.replace so a crash leaves no half-file.
    Also caches the alphabet histogram as <sha256>/alphabet.json.
    """
    contig_dir = os.path.join(cache_dir, fasta_sha256)
    os.makedirs(contig_dir, exist_ok=True)

    # Alphabet gate — scans all sequence bytes
    alphabet_path = os.path.join(contig_dir, "alphabet.json")
    if not os.path.exists(alphabet_path):
        logger.info("Scanning FASTA alphabet (§5 layer 5)...")
        histogram = scan_fasta_alphabet(fasta_path)
        tmp = alphabet_path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(histogram, f, indent=2)
        os.replace(tmp, alphabet_path)
        logger.info("Alphabet: %s", histogram)
    else:
        with open(alphabet_path) as f:
            histogram = json.load(f)
        # Still validate
        disallowed = {c: n for c, n in histogram.items() if c.upper() not in "ACGTN"}
        if disallowed:
            raise ValueError(f"Cached alphabet contains disallowed characters: {disallowed}")

    # Build per-contig cumsum arrays
    from fragments_h5.fragment import get_g_or_c_cumsum

    with pysam.FastaFile(fasta_path) as fasta:
        for contig in fasta.references:
            npy_path = os.path.join(contig_dir, f"{contig}.npy")
            if os.path.exists(npy_path):
                continue
            logger.info("Building cumsum cache for %s...", contig)
            cumsum, offset = get_g_or_c_cumsum(fasta_path, contig)
            assert offset == 0, f"Whole-contig fetch returned offset={offset}"
            tmp = npy_path + ".tmp"
            with open(tmp, "wb") as tmpf:
                numpy.save(tmpf, cumsum)
            os.replace(tmp, npy_path)

    # Report cache size
    total_bytes = sum(
        os.path.getsize(os.path.join(contig_dir, f))
        for f in os.listdir(contig_dir)
    )
    logger.info("Cumsum cache at %s: %.1f GB", contig_dir, total_bytes / 1e9)
    return histogram


def load_cached_cumsum(cache_dir, fasta_sha256, contig):
    """Load a cached float64 cumsum array via mmap (read-only)."""
    npy_path = os.path.join(cache_dir, fasta_sha256, f"{contig}.npy")
    return numpy.load(npy_path, mmap_mode="r")


# ---------------------------------------------------------------------------
# Truncation detection and execution (§3.2)
# ---------------------------------------------------------------------------

def detect_padding_row(f, contig):
    """Detect the phantom trailing padding row per §3.2.1.

    Returns 'truncate', 'clean', or raises RepairAbort.
    """
    grp = f["data"][contig]
    starts = grp["starts"]
    n = len(starts)
    if n < 2:
        return "clean"

    # Condition (a): sortedness violation at the tail
    last_start = int(starts[-1])
    penult_start = int(starts[-2])
    cond_a = last_start < penult_start

    # Condition (b): full zero signature
    cond_b_parts = [
        last_start == 0,
        int(grp["lengths"][-1]) == 0,
    ]
    if "gc" in grp:
        cond_b_parts.append(int(grp["gc"][-1]) == 0)
    if "mapq" in grp:
        mapq_last = grp["mapq"][-1]
        cond_b_parts.append(numpy.all(mapq_last == 0))
    if "strand" in grp:
        cond_b_parts.append(grp["strand"][-1] == b"")
    for ds_name in grp:
        if ds_name in ("starts", "lengths", "gc", "mapq", "strand"):
            continue
        # methyl arrays, fragment_end_clipped: check final element is 0
        cond_b_parts.append(int(grp[ds_name][-1]) == 0)
    cond_b = all(cond_b_parts)

    if cond_a and cond_b:
        return "truncate"
    elif not cond_a and not cond_b:
        return "clean"
    else:
        raise RepairAbort(
            f"Contig {contig}: exactly one of (sortedness violation, zero signature) holds — "
            f"cond_a={cond_a}, cond_b={cond_b}. Requires human review."
        )


def truncate_contig_datasets(f, contig):
    """Delete and recreate all datasets for a contig, dropping the trailing row.

    Uses mk_dataset's exact parameters: gzip, compression_opts=4, chunks=True (§3.2.2).
    Handles 2-D datasets (mapq) by truncating along axis 0 only (§2.2.2).
    """
    grp = f["data"][contig]
    dataset_names = list(grp.keys())
    truncated = {}

    for ds_name in dataset_names:
        ds = grp[ds_name]
        arr = ds[:]
        dtype = ds.dtype
        # Truncate: drop last element along axis 0
        truncated[ds_name] = (arr[:-1], dtype)

    # Delete and recreate
    for ds_name in dataset_names:
        del grp[ds_name]

    for ds_name, (arr, dtype) in truncated.items():
        grp.create_dataset(
            ds_name,
            data=arr,
            dtype=dtype,
            compression="gzip",
            compression_opts=4,
            chunks=True,
        )


def rebuild_contig_index(f, contig, contig_length, index_block_size):
    """Rebuild index for a contig from the (truncated) starts array (§3.2.3).

    Replicates fragments_h5.py:1225-1241 including both skip guards.
    """
    starts = f[f"data/{contig}/starts"][:]
    n_frags = len(starts)
    index_key = f"index/{contig}"

    # Guard 1: skip contigs shorter than index block size
    if contig_length <= index_block_size:
        if index_key in f:
            del f[index_key]
        return False

    # Guard 2: skip contigs with too few reads
    if n_frags < MIN_NUM_READS_FOR_INDEX:
        if index_key in f:
            del f[index_key]
        return False

    block_indices = numpy.arange(0, contig_length, index_block_size)
    index_poss = numpy.searchsorted(starts, block_indices, side="left")
    index_poss = numpy.append(index_poss, n_frags)

    if index_key in f:
        del f[index_key]
    f[index_key] = index_poss
    return True


# ---------------------------------------------------------------------------
# GC recomputation (§3.4)
# ---------------------------------------------------------------------------

def recompute_gc_for_contig(starts, lengths, cumsum):
    """Recompute GC uint8 values for one contig.

    Uses CPython's round(x, 5) semantics per §3.4 to guarantee exact reproduction
    of the original Stage A computation.

    Returns a uint8 array of the same length as starts.
    """
    n = len(starts)
    gc_u8 = numpy.full(n, 255, dtype=numpy.uint8)

    # Compute stops
    stops = starts.astype(numpy.int64) + lengths.astype(numpy.int64)

    # Mask for valid (nonzero-length) fragments
    valid = lengths > 0
    valid_idx = numpy.where(valid)[0]

    if len(valid_idx) == 0:
        return gc_u8

    # Vectorized numerator computation (exact in float64)
    nums = cumsum[stops[valid_idx]] - cumsum[starts[valid_idx]]
    denoms = lengths[valid_idx].astype(numpy.float64)
    q = nums / denoms

    # Validate range — original build would have raised ValueError
    if numpy.any((q < -1e-10) | (q > 1.0 + 1e-10)):
        out_of_range = q[(q < -1e-10) | (q > 1.0 + 1e-10)]
        raise RepairAbort(
            f"Recomputed GC fraction out of [0, 1]: {out_of_range[:5]}"
        )

    # Stage A: CPython round(x, 5) — element-wise loop for exact semantics
    q5 = numpy.empty(len(valid_idx), dtype=numpy.float64)
    for i in range(len(valid_idx)):
        q5[i] = round(float(q[i]), 5)

    # Stage B: int(round(gc * 254)) — numpy.rint is equivalent to CPython round()
    gc_u8[valid_idx] = numpy.rint(q5 * 254.0).astype(numpy.uint8)

    return gc_u8


# ---------------------------------------------------------------------------
# Float32 accumulator simulation (§5b)
# ---------------------------------------------------------------------------

def simulate_float32_cumsum(float64_cumsum):
    """Simulate the pre-fix float32 accumulator from the cached float64 cumsum.

    The pre-fix code did: seq[:, (1,2)].sum(axis=1).cumsum() in float32.
    We reconstruct the per-base contributions by diff'ing the exact float64 cumsum,
    cast to float32, and re-accumulate in float32 to reproduce the saturation behavior.
    """
    # diff gives per-base C+G contributions (exact {0.0, 0.5, 1.0} in float64)
    per_base = numpy.diff(float64_cumsum).astype(numpy.float32)
    # Prepend the pad contribution (0.0)
    all_per_base = numpy.empty(len(float64_cumsum), dtype=numpy.float32)
    all_per_base[0] = numpy.float32(0.0)
    all_per_base[1:] = per_base
    return numpy.cumsum(all_per_base)


def build_n_prefix_sum(float64_cumsum):
    """Build a prefix sum of N-base positions from the float64 cumsum.

    N contributes exactly 0.5 to the cumsum; A/T contribute 0; C/G contribute 1.0.
    So position i (after pad) is N iff diff(cumsum)[i] == 0.5.
    """
    per_base = numpy.diff(float64_cumsum)
    is_n = numpy.isclose(per_base, 0.5, atol=1e-12)
    prefix = numpy.zeros(len(per_base) + 1, dtype=numpy.int64)
    prefix[1:] = numpy.cumsum(is_n)
    return prefix


def classify_gc_changes(old_gc, new_gc, starts, lengths, f32_cumsum, n_prefix):
    """Classify and validate GC changes per §5b's region rules.

    Returns a dict with:
      - region counts and change statistics
      - any violations (blocking errors)
    """
    n = len(starts)
    stops = starts.astype(numpy.int64) + lengths.astype(numpy.int64)

    # Classify by simulated float32 cumsum at frag_stop
    f32_at_stop = f32_cumsum[stops]

    below_t23 = f32_at_stop < T23
    between_t23_t24 = (f32_at_stop >= T23) & (f32_at_stop < T24)
    above_t24 = f32_at_stop >= T24

    changed = old_gc != new_gc

    violations = []
    stats = {
        "total_frags": n,
        "below_T23": int(numpy.sum(below_t23)),
        "between_T23_T24": int(numpy.sum(between_t23_t24)),
        "above_T24": int(numpy.sum(above_t24)),
        "changed_total": int(numpy.sum(changed)),
        "changed_below_T23": 0,
        "changed_between_T23_T24": 0,
        "changed_between_T23_T24_no_N": 0,
        "changed_above_T24": 0,
        "zero_to_nonzero": int(numpy.sum((old_gc == 0) & (new_gc != 0) & changed)),
        "nonzero_to_nonzero": int(numpy.sum((old_gc != 0) & (new_gc != 0) & changed)),
        "x255_to_other": int(numpy.sum((old_gc == 255) & (new_gc != 255) & changed)),
        "other_to_x255": int(numpy.sum((old_gc != 255) & (new_gc == 255) & changed)),
    }

    # Check transitions: 255→x and x→255 must be zero
    if stats["x255_to_other"] > 0:
        violations.append(
            f"{stats['x255_to_other']} 255→x transitions (wrong reference?)"
        )
    if stats["other_to_x255"] > 0:
        violations.append(
            f"{stats['other_to_x255']} x→255 transitions (wrong reference?)"
        )

    # Region: below T23 — byte-identical required
    changed_below = changed & below_t23
    stats["changed_below_T23"] = int(numpy.sum(changed_below))
    if stats["changed_below_T23"] > 0:
        violations.append(
            f"{stats['changed_below_T23']} changes in pre-saturation region (< T23) — abort"
        )

    # Region: between T23 and T24 — change only if span contains N
    changed_middle = changed & between_t23_t24
    stats["changed_between_T23_T24"] = int(numpy.sum(changed_middle))
    if stats["changed_between_T23_T24"] > 0:
        # Check which of these DON'T contain an N
        middle_changed_idx = numpy.where(changed_middle)[0]
        for idx in middle_changed_idx:
            s, e = int(starts[idx]), int(stops[idx])
            has_n = (n_prefix[e] - n_prefix[s]) > 0
            if not has_n:
                stats["changed_between_T23_T24_no_N"] += 1
        if stats["changed_between_T23_T24_no_N"] > 0:
            violations.append(
                f"{stats['changed_between_T23_T24_no_N']} middle-band changes on N-free "
                f"spans — abort"
            )

    # Region: above T24 — change freely
    stats["changed_above_T24"] = int(numpy.sum(changed & above_t24))

    return stats, violations


# ---------------------------------------------------------------------------
# Provenance (§9.1)
# ---------------------------------------------------------------------------

def build_repair_history_entry(
    fasta_uri, fasta_sha256, source_sha256, backup_uri,
    contigs_truncated, has_gc, argv,
):
    """Build one _repair_history JSON element per §9.1."""
    try:
        version = importlib.metadata.version("fragments-h5")
    except importlib.metadata.PackageNotFoundError:
        version = "unknown"

    datasets = []
    if contigs_truncated > 0:
        datasets.append("data/*/* (truncated by 1 row)")
        datasets.append("index/*")
        datasets.append("fragment_length_counts")
    if has_gc:
        datasets.append("data/*/gc")

    entry = {
        "tool": "repair-fragments-h5-gc",
        "version": version,
        "argv": argv,
        "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "reason": (
            "gc-float32-cumsum-saturation (fixed by 95c76f5) + trailing padding row "
            "(fixed by 778f4d1); see docs/pending/gc_repair_tool.md"
        ),
        "datasets": datasets,
        "rows_removed_per_contig": 1 if contigs_truncated > 0 else 0,
        "fasta_uri": fasta_uri,
        "fasta_sha256": fasta_sha256,
        "source_sha256": source_sha256,
        "backup_uri": backup_uri,
    }
    return entry


def write_repair_history(f, entry):
    """Append a repair history entry to the _repair_history root attr."""
    raw = f.attrs.get("_repair_history", "[]")
    history = json.loads(raw)
    history.append(entry)
    f.attrs["_repair_history"] = json.dumps(history)


# ---------------------------------------------------------------------------
# S3 operations (§6.2)
# ---------------------------------------------------------------------------

def download_from_s3(bucket, key, local_path):
    """Download an S3 object to a local file."""
    import boto3
    s3 = boto3.client("s3")
    logger.info("Downloading s3://%s/%s -> %s", bucket, key, local_path)
    s3.download_file(bucket, key, local_path)
    size = os.path.getsize(local_path)
    logger.info("Downloaded %.1f MB", size / 1e6)
    return size


def upload_to_s3(bucket, key, local_path):
    """Upload a local file to S3 with CRC64NVME verification (§6.2).

    Uses single-part put_object — refuses files > 5 GB.
    Returns the CRC64NVME value from the response.
    """
    import boto3
    size = os.path.getsize(local_path)
    if size > MAX_UPLOAD_BYTES:
        raise RepairAbort(f"File {local_path} is {size} bytes (> 5 GB limit)")

    s3 = boto3.client("s3")
    logger.info("Uploading %s (%.1f MB) -> s3://%s/%s", local_path, size / 1e6, bucket, key)
    with open(local_path, "rb") as f:
        response = s3.put_object(
            Bucket=bucket,
            Key=key,
            Body=f,
            ChecksumAlgorithm="CRC64NVME",
        )
    crc = response.get("ChecksumCRC64NVME")
    logger.info("Uploaded, CRC64NVME=%s", crc)
    return crc


def verify_s3_upload(bucket, key, expected_crc, expected_size):
    """Verify the uploaded object matches expectations (§6.2)."""
    import boto3
    s3 = boto3.client("s3")
    head = s3.head_object(Bucket=bucket, Key=key, ChecksumMode="ENABLED")

    actual_crc = head.get("ChecksumCRC64NVME")
    actual_size = head["ContentLength"]

    errors = []
    if actual_crc is None:
        errors.append("No CRC64NVME in head_object response")
    elif actual_crc != expected_crc:
        errors.append(f"CRC64NVME mismatch: uploaded={expected_crc}, head={actual_crc}")
    if "-" in (actual_crc or ""):
        errors.append(f"CRC64NVME has composite suffix: {actual_crc}")
    if actual_size != expected_size:
        errors.append(f"Size mismatch: local={expected_size}, remote={actual_size}")
    if errors:
        raise RepairAbort(f"Upload verification failed for s3://{bucket}/{key}: {errors}")

    logger.info("Upload verified: CRC=%s, size=%d", actual_crc, actual_size)


# ---------------------------------------------------------------------------
# Ledger (§8.2)
# ---------------------------------------------------------------------------

def append_ledger_record(ledger_path, record):
    """Append one JSON line to the ledger file."""
    with open(ledger_path, "a") as f:
        f.write(json.dumps(record) + "\n")
    logger.info("Ledger record appended to %s", ledger_path)


def load_ledger_ok_keys(ledger_path):
    """Load all keys with status='ok' from the ledger, keyed by (key, crc64nvme)."""
    completed = {}
    if not os.path.exists(ledger_path):
        return completed
    with open(ledger_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            if rec.get("status") in ("ok", "ok_no_gc"):
                completed[rec["key"]] = rec.get("crc64nvme_repaired")
    return completed


# ---------------------------------------------------------------------------
# Core repair pipeline (§3.1)
# ---------------------------------------------------------------------------

def repair_local_file(
    local_path,
    fasta_path=None,
    cumsum_cache_dir=None,
    fasta_sha256=None,
    dry_run=True,
    fasta_uri=None,
    backup_uri=None,
    argv=None,
):
    """Core repair logic for one local H5 file.

    This is the per-file pipeline from §3.1:
    1. Preflight and reference safety
    2. Detect and truncate padding row
    3. Recompute GC
    4. Diff and validate
    5. Write (if not dry_run)

    Returns a report dict with all per-file metrics and outcomes.
    """
    report = {
        "path": local_path,
        "dry_run": dry_run,
        "status": "pending",
    }

    # Step 1: sha256 of the original file
    sha256_original = compute_file_sha256(local_path)
    report["sha256_original"] = sha256_original

    # Read the file and run all analysis
    with h5py.File(local_path, "r") as f:
        # Step 2: Preflight
        gc_presence = check_gc_presence(f["data"])
        report["gc_presence"] = gc_presence

        if gc_presence == "partial":
            raise RepairAbort(
                f"Partial gc dataset — some contigs have gc, others don't. Abort."
            )

        has_gc = gc_presence == "all"

        if has_gc and fasta_path is None:
            raise RepairAbort("File has gc dataset but --fasta not provided")

        # Record attr hash (excluding _repair_history)
        attr_hash = hash_attrs_excluding(f)
        report["attr_hash"] = attr_hash

        # Read metadata
        contig_lengths = eval(f.attrs["_contig_lengths_str"])
        index_block_size = int(f.attrs["index_block_size"])

        contigs = list(f["data"].keys())
        report["contigs"] = contigs

        # FASTA reference info (if gc is present)
        fasta_ref = None
        if has_gc:
            fasta_ref = pysam.FastaFile(fasta_path)

        # ── Step 3: Detect padding per contig ──
        truncation_verdicts = {}
        contig_arrays = {}  # contig -> {dataset_name: array}

        for contig in contigs:
            verdict = detect_padding_row(f, contig)
            truncation_verdicts[contig] = verdict

            # Read all arrays into memory
            grp = f["data"][contig]
            arrays = {}
            for ds_name in grp:
                arrays[ds_name] = grp[ds_name][:]
            contig_arrays[contig] = arrays

        # Check: all contigs must agree on truncation
        truncated_contigs = [c for c, v in truncation_verdicts.items() if v == "truncate"]
        clean_contigs = [c for c, v in truncation_verdicts.items() if v == "clean"]

        # A file where some contigs truncate and others don't is blocking
        if truncated_contigs and clean_contigs:
            raise RepairAbort(
                f"Mixed truncation verdicts — {len(truncated_contigs)} truncated, "
                f"{len(clean_contigs)} clean. Requires human review."
            )

        report["truncation_verdicts"] = truncation_verdicts
        report["contigs_truncated"] = len(truncated_contigs)

        # Truncate in memory
        for contig in truncated_contigs:
            for ds_name, arr in contig_arrays[contig].items():
                contig_arrays[contig][ds_name] = arr[:-1]

        # Validate: starts must be non-decreasing after truncation (§3.2.3)
        for contig in contigs:
            starts = contig_arrays[contig]["starts"]
            if len(starts) >= 2:
                if numpy.diff(starts).min() < 0:
                    raise RepairAbort(
                        f"Contig {contig}: starts not non-decreasing after truncation"
                    )

        # Record existing index key set
        existing_index_contigs = set(f["index"].keys()) if "index" in f else set()
        report["existing_index_contigs"] = sorted(existing_index_contigs)

        # Record old fragment_length_counts
        old_flc = f["fragment_length_counts"][:] if "fragment_length_counts" in f else None
        report["old_n_fragments"] = int(old_flc.sum()) if old_flc is not None else None

        # ── Step 4: Recompute GC (on truncated arrays) ──
        gc_results = {}

        if has_gc:
            for contig in contigs:
                starts = contig_arrays[contig]["starts"]
                lengths = contig_arrays[contig]["lengths"]
                old_gc = contig_arrays[contig]["gc"]

                # §5 layer 3: check contig against FASTA
                if contig in fasta_ref.references:
                    # Geometry checks
                    fasta_len = fasta_ref.get_reference_length(contig)
                    h5_len = contig_lengths.get(contig)
                    if h5_len is not None and fasta_len != h5_len:
                        raise RepairAbort(
                            f"Contig {contig}: FASTA length {fasta_len} != H5 contig_lengths {h5_len}"
                        )
                    stops = starts.astype(numpy.int64) + lengths.astype(numpy.int64)
                    if len(stops) > 0 and int(stops.max()) > fasta_len:
                        raise RepairAbort(
                            f"Contig {contig}: fragment extends past FASTA contig "
                            f"(max stop={int(stops.max())}, contig_len={fasta_len})"
                        )

                    # Load cumsum and recompute
                    cumsum = load_cached_cumsum(cumsum_cache_dir, fasta_sha256, contig)
                    new_gc = recompute_gc_for_contig(starts, lengths, cumsum)

                    # Simulate float32 accumulator for region classification
                    f32_cumsum = simulate_float32_cumsum(cumsum)
                    n_prefix = build_n_prefix_sum(cumsum)

                    # Classify and validate changes
                    stats, violations = classify_gc_changes(
                        old_gc, new_gc, starts, lengths, f32_cumsum, n_prefix
                    )
                    gc_results[contig] = {
                        "action": "recomputed",
                        "new_gc": new_gc,
                        "stats": stats,
                        "violations": violations,
                    }

                    if violations:
                        raise RepairAbort(
                            f"Contig {contig}: GC validation violations: {violations}"
                        )
                else:
                    # Contig absent from FASTA — assert gc is all 255
                    if not numpy.all(old_gc == 255):
                        non_255 = numpy.sum(old_gc != 255)
                        raise RepairAbort(
                            f"Contig {contig}: absent from FASTA but gc has "
                            f"{non_255} non-255 values — abort"
                        )
                    gc_results[contig] = {
                        "action": "skipped_absent",
                        "stats": {"total_frags": len(starts), "all_255": True},
                        "violations": [],
                    }

        report["gc_results"] = {
            contig: {k: v for k, v in res.items() if k != "new_gc"}
            for contig, res in gc_results.items()
        }

        if fasta_ref is not None:
            fasta_ref.close()

    # ── Step 6: If dry-run, stop and report ──
    if dry_run:
        report["status"] = "dry_run_ok"
        _log_report(report)
        return report

    # ── Step 7: Write into the local copy ──
    logger.info("Writing repairs to %s", local_path)
    with h5py.File(local_path, "r+") as f:
        # Truncate datasets and write gc
        for contig in contigs:
            if truncation_verdicts[contig] == "truncate":
                truncate_contig_datasets(f, contig)

            # Write new gc values
            if has_gc and contig in gc_results and gc_results[contig]["action"] == "recomputed":
                new_gc = gc_results[contig]["new_gc"]
                gc_key = f"data/{contig}/gc"
                if gc_key in f:
                    del f[gc_key]
                f.create_dataset(
                    gc_key,
                    data=new_gc,
                    dtype="uint8",
                    compression="gzip",
                    compression_opts=4,
                    chunks=True,
                )

        # Rebuild index for all contigs
        new_index_contigs = set()
        for contig in contigs:
            contig_len = contig_lengths.get(contig, 0)
            if rebuild_contig_index(f, contig, contig_len, index_block_size):
                new_index_contigs.add(contig)

        # Verify index key set matches
        if new_index_contigs != existing_index_contigs:
            added = new_index_contigs - existing_index_contigs
            removed = existing_index_contigs - new_index_contigs
            raise RepairAbort(
                f"Index key set changed: added={added}, removed={removed}"
            )

        # Write repair history
        history_entry = build_repair_history_entry(
            fasta_uri=fasta_uri,
            fasta_sha256=fasta_sha256,
            source_sha256=sha256_original,
            backup_uri=backup_uri,
            contigs_truncated=len(truncated_contigs),
            has_gc=has_gc,
            argv=argv or [],
        )
        write_repair_history(f, history_entry)

    # ── Step: rebuild fragment_length_counts (§3.2.4 close-and-reopen) ──
    fm_h5 = FragmentsH5(local_path, "r+")
    try:
        old_counts = fm_h5._f["fragment_length_counts"][:]
        del fm_h5._f["fragment_length_counts"]
        fm_h5._add_fragment_length_counts()

        # Verify the exact arithmetic identity (§3.2.4)
        new_counts = fm_h5._f["fragment_length_counts"][:]
        n_truncated = len(truncated_contigs)
        if n_truncated > 0:
            if int(new_counts[0]) != int(old_counts[0]) - n_truncated:
                raise RepairAbort(
                    f"fragment_length_counts[0]: expected {int(old_counts[0]) - n_truncated}, "
                    f"got {int(new_counts[0])}"
                )
        if not numpy.array_equal(new_counts[1:], old_counts[1:]):
            diff_bins = numpy.where(new_counts[1:] != old_counts[1:])[0] + 1
            raise RepairAbort(
                f"fragment_length_counts changed at bins {diff_bins[:10]}"
            )

        report["new_n_fragments"] = int(new_counts.sum())
    finally:
        fm_h5.close()

    # ── Step 8: Re-verify ──
    with h5py.File(local_path, "r") as f:
        # Verify attrs unchanged
        new_attr_hash = hash_attrs_excluding(f)
        if new_attr_hash != attr_hash:
            raise RepairAbort(
                f"Attr hash changed: {attr_hash} -> {new_attr_hash}"
            )

        # Verify _repair_history is valid JSON with exactly one new element
        raw = f.attrs.get("_repair_history")
        if raw is None:
            raise RepairAbort("_repair_history attr missing after write")
        history = json.loads(raw)
        if not isinstance(history, list) or len(history) < 1:
            raise RepairAbort(f"_repair_history malformed: {raw[:200]}")

        # Verify structural invariants (§7.3)
        for contig in contigs:
            grp = f["data"][contig]
            ds_lengths = {ds: len(grp[ds]) for ds in grp}
            # All datasets same length along axis 0
            lengths_set = set()
            for ds, length in ds_lengths.items():
                if grp[ds].ndim > 1:
                    lengths_set.add(grp[ds].shape[0])
                else:
                    lengths_set.add(length)
            if len(lengths_set) > 1:
                raise RepairAbort(
                    f"Contig {contig}: datasets have different lengths: {ds_lengths}"
                )

            # starts non-decreasing
            s = grp["starts"][:]
            if len(s) >= 2 and numpy.diff(s).min() < 0:
                raise RepairAbort(f"Contig {contig}: starts not non-decreasing after write")

            # Index final entry equals len(starts)
            if contig in f.get("index", {}):
                idx = f[f"index/{contig}"][:]
                if idx[-1] != len(s):
                    raise RepairAbort(
                        f"Contig {contig}: index sentinel {idx[-1]} != len(starts) {len(s)}"
                    )

        # Verify FragmentsH5 can open and read the file
        fm = FragmentsH5(local_path, "r")
        try:
            _ = fm.n_fragments
            _ = fm.contig_lengths
        finally:
            fm.close()

    sha256_repaired = compute_file_sha256(local_path)
    report["sha256_repaired"] = sha256_repaired
    report["status"] = "ok" if has_gc else "ok_no_gc"

    _log_report(report)
    return report


def _log_report(report):
    """Log a human-readable summary of the repair report."""
    logger.info("=== Repair report for %s ===", report.get("path", "?"))
    logger.info("  Status: %s", report.get("status"))
    logger.info("  GC presence: %s", report.get("gc_presence"))
    logger.info("  Contigs truncated: %s", report.get("contigs_truncated"))
    if report.get("old_n_fragments") is not None:
        logger.info("  n_fragments before: %s", report.get("old_n_fragments"))
    if report.get("new_n_fragments") is not None:
        logger.info("  n_fragments after: %s", report.get("new_n_fragments"))

    gc_results = report.get("gc_results", {})
    total_changed = 0
    for contig, res in gc_results.items():
        stats = res.get("stats", {})
        changed = stats.get("changed_total", 0)
        if changed > 0:
            total_changed += changed
            logger.info("  Contig %s: %d gc changes", contig, changed)
    if gc_results:
        logger.info("  Total gc bytes changed: %d", total_changed)


# ---------------------------------------------------------------------------
# Orchestration — S3-aware wrapper
# ---------------------------------------------------------------------------

def process_one_target(key, bucket, args, fasta_sha256, ledger_ok_keys):
    """Download, repair, upload, ledger — one S3 key end-to-end."""
    # Resume: skip if already done
    if key in ledger_ok_keys:
        import boto3
        s3 = boto3.client("s3")
        head = s3.head_object(Bucket=bucket, Key=key, ChecksumMode="ENABLED")
        remote_crc = head.get("ChecksumCRC64NVME")
        if remote_crc == ledger_ok_keys[key]:
            logger.info("Skipping %s — already in ledger with matching CRC", key)
            return {"key": key, "status": "skipped_already_done"}
        else:
            logger.warning(
                "Key %s in ledger but CRC mismatch (ledger=%s, remote=%s) — reprocessing",
                key, ledger_ok_keys[key], remote_crc,
            )

    with tempfile.TemporaryDirectory(prefix="repair_h5_") as tmpdir:
        local_path = os.path.join(tmpdir, os.path.basename(key))

        # Download
        download_from_s3(bucket, key, local_path)

        # Repair
        report = repair_local_file(
            local_path=local_path,
            fasta_path=args.fasta,
            cumsum_cache_dir=args.cumsum_cache,
            fasta_sha256=fasta_sha256,
            dry_run=args.dry_run,
            fasta_uri=args.fasta_uri,
            backup_uri=args.backup_uri,
            argv=sys.argv,
        )

        if args.dry_run:
            return report

        # Upload
        crc = upload_to_s3(bucket, key, local_path)
        verify_s3_upload(bucket, key, crc, os.path.getsize(local_path))

        # Ledger
        gc_results = report.get("gc_results", {})
        total_gc_changed = sum(
            res.get("stats", {}).get("changed_total", 0)
            for res in gc_results.values()
        )
        contigs_repaired = sum(
            1 for res in gc_results.values() if res.get("action") == "recomputed"
        )
        contigs_skipped = sum(
            1 for res in gc_results.values() if res.get("action") == "skipped_absent"
        )

        ledger_record = {
            "key": key,
            "status": report["status"],
            "sha256_original": report["sha256_original"],
            "sha256_repaired": report["sha256_repaired"],
            "crc64nvme_repaired": crc,
            "fasta_sha256": fasta_sha256,
            "tool_version": _get_version(),
            "contigs_repaired": contigs_repaired,
            "contigs_skipped_absent_from_fasta": contigs_skipped,
            "contigs_truncated": report["contigs_truncated"],
            "rows_removed": report["contigs_truncated"],
            "n_fragments_before": report.get("old_n_fragments"),
            "n_fragments_after": report.get("new_n_fragments"),
            "gc_bytes_changed": total_gc_changed if gc_results else None,
            "repair_history_len": 1,
            "started_utc": report.get("started_utc"),
            "finished_utc": datetime.now(timezone.utc).isoformat(),
        }
        append_ledger_record(args.ledger, ledger_record)
        report["ledger_record"] = ledger_record
        return report


def _get_version():
    try:
        return importlib.metadata.version("fragments-h5")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="repair-fragments-h5-gc",
        description="Repair GC float32 saturation and phantom padding rows in fragment H5 files.",
    )
    parser.add_argument(
        "--fasta",
        help="Path to the reference FASTA file. Required when any target has a gc dataset.",
    )
    parser.add_argument(
        "--target-list",
        help="Newline-delimited file of S3 keys to process.",
    )
    parser.add_argument(
        "--bucket",
        help="S3 bucket name (no s3:// prefix).",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--dry-run",
        action="store_true",
        default=True,
        help="Analyze and report without writing (default).",
    )
    mode.add_argument(
        "--apply",
        action="store_true",
        help="Actually write repairs to the files.",
    )
    parser.add_argument(
        "--expect-fasta-sha256",
        help="Expected SHA-256 of the FASTA file (§5 layer 2). "
             "Required in --apply mode whenever --fasta is given.",
    )
    parser.add_argument(
        "--ledger",
        help="Path for the append-only JSONL ledger file. Required in --apply mode.",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=1,
        help="Maximum number of files to process (default: 1).",
    )
    parser.add_argument(
        "--cumsum-cache",
        default=os.path.join(tempfile.gettempdir(), "repair_h5_cumsum_cache"),
        help="Directory for per-contig cumsum cache files.",
    )
    parser.add_argument(
        "--rebuild-cache",
        action="store_true",
        help="Delete and rebuild the cumsum cache for the current FASTA.",
    )
    parser.add_argument(
        "--fasta-uri",
        help="S3 URI of the FASTA file (recorded in provenance, not downloaded).",
    )
    parser.add_argument(
        "--backup-uri",
        help="S3 URI of the backup prefix (recorded in provenance).",
    )
    parser.add_argument(
        "--local-file",
        action="append",
        default=[],
        help="Process a local file directly (for testing). Can be repeated.",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output.",
    )
    args = parser.parse_args(argv)

    # --apply overrides --dry-run
    if args.apply:
        args.dry_run = False

    # Validation for --apply mode
    if not args.dry_run:
        if args.target_list and not args.ledger:
            parser.error("--ledger is required in --apply mode")
        if args.target_list and not args.bucket:
            parser.error("--bucket is required when using --target-list")
        # §5 layer 2: the reference bytes must be pinned before any write. --fasta
        # is itself required for any target that has a gc dataset (repair_local_file
        # aborts otherwise), so gating on --fasta covers every gc-bearing target.
        if args.fasta and not args.expect_fasta_sha256:
            parser.error("--expect-fasta-sha256 is required in --apply mode when --fasta is given")

    return args


def main(argv=None):
    args = parse_args(argv)

    # Configure logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        level=level,
    )

    # Load target list
    targets = []
    if args.target_list:
        with open(args.target_list) as f:
            raw_keys = [line.strip() for line in f if line.strip()]
        # Reject duplicates
        if len(raw_keys) != len(set(raw_keys)):
            dupes = [k for k in raw_keys if raw_keys.count(k) > 1]
            logger.error("Duplicate keys in target list: %s", dupes[:5])
            sys.exit(1)
        targets = raw_keys[:args.max_files]
    elif args.local_file:
        targets = args.local_file[:args.max_files]

    if not targets:
        logger.error("No targets specified. Use --target-list or --local-file.")
        sys.exit(1)

    logger.info("Processing %d target(s) (max-files=%d, dry-run=%s)",
                len(targets), args.max_files, args.dry_run)

    # FASTA preflight
    fasta_sha256 = None
    if args.fasta:
        logger.info("Computing FASTA SHA-256...")
        fasta_sha256 = compute_file_sha256(args.fasta)
        logger.info("FASTA SHA-256: %s", fasta_sha256)

        if args.expect_fasta_sha256 and args.expect_fasta_sha256 != fasta_sha256:
            logger.error(
                "FASTA SHA-256 mismatch: expected %s, got %s",
                args.expect_fasta_sha256, fasta_sha256,
            )
            sys.exit(1)

        # Rebuild cache if requested
        if args.rebuild_cache:
            cache_dir = os.path.join(args.cumsum_cache, fasta_sha256)
            if os.path.exists(cache_dir):
                import shutil
                shutil.rmtree(cache_dir)

        # Build cumsum cache
        build_cumsum_cache(args.fasta, args.cumsum_cache, fasta_sha256)

    # Ledger for resume
    ledger_ok_keys = {}
    if args.ledger and os.path.exists(args.ledger):
        ledger_ok_keys = load_ledger_ok_keys(args.ledger)
        logger.info("Loaded %d completed keys from ledger", len(ledger_ok_keys))

    # Process targets
    results = []
    is_local = bool(args.local_file)

    for target in targets:
        try:
            if is_local:
                report = repair_local_file(
                    local_path=target,
                    fasta_path=args.fasta,
                    cumsum_cache_dir=args.cumsum_cache,
                    fasta_sha256=fasta_sha256,
                    dry_run=args.dry_run,
                    fasta_uri=args.fasta_uri,
                    backup_uri=args.backup_uri,
                    argv=sys.argv,
                )
            else:
                report = process_one_target(
                    key=target,
                    bucket=args.bucket,
                    args=args,
                    fasta_sha256=fasta_sha256,
                    ledger_ok_keys=ledger_ok_keys,
                )
            results.append(report)
        except RepairAbort as e:
            logger.error("ABORT on %s: %s", target, e)
            results.append({"target": target, "status": "aborted", "error": str(e)})
        except Exception as e:
            logger.error("ERROR on %s: %s", target, e, exc_info=True)
            results.append({"target": target, "status": "error", "error": str(e)})

    # Summary
    ok = sum(1 for r in results if r.get("status", "").startswith("ok") or r.get("status") == "dry_run_ok")
    aborted = sum(1 for r in results if r.get("status") == "aborted")
    errored = sum(1 for r in results if r.get("status") == "error")
    skipped = sum(1 for r in results if r.get("status") == "skipped_already_done")

    logger.info(
        "Summary: %d ok, %d aborted, %d errors, %d skipped (of %d targets)",
        ok, aborted, errored, skipped, len(targets),
    )

    if args.fasta:
        cache_path = os.path.join(args.cumsum_cache, fasta_sha256)
        if os.path.exists(cache_path):
            cache_size = sum(
                os.path.getsize(os.path.join(cache_path, f))
                for f in os.listdir(cache_path)
            )
            logger.info("Cumsum cache: %s (%.1f GB)", cache_path, cache_size / 1e9)

    if aborted or errored:
        sys.exit(1)


if __name__ == "__main__":
    main()
