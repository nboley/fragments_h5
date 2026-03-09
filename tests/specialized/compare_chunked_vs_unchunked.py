"""Compare build-fragments-h5 output with chunking enabled vs disabled,
and optionally against a previous Docker image version.

Usage:
    # Compare chunked vs unchunked (current code)
    python tests/specialized/compare_chunked_vs_unchunked.py \
        /tmp/debugging_frag_h5/mapped.sorted.bam \
        /tmp/debugging_frag_h5/sequence.fa.bgz

    # Also compare against Docker image 2.7.2 on chr1
    python tests/specialized/compare_chunked_vs_unchunked.py \
        /tmp/debugging_frag_h5/mapped.sorted.bam \
        /tmp/debugging_frag_h5/sequence.fa.bgz \
        --compare-docker 2.7.2 --contigs chr1
"""
import argparse
import subprocess
import sys
import tempfile
import os

import h5py
import numpy

from fragments_h5.fragments_h5 import build_fragments_h5


def compare_h5_files(path_a, path_b, label_a="A", label_b="B"):
    errors = []
    with h5py.File(path_a, "r") as a, h5py.File(path_b, "r") as b:
        # Compare attributes
        for attr in ("index_block_size", "max_fragment_length"):
            va, vb = a.attrs[attr], b.attrs[attr]
            if va != vb:
                errors.append(f"Attribute {attr} differs: {va} vs {vb}")

        # Compare data contigs
        contigs_a = set(a["data"].keys())
        contigs_b = set(b["data"].keys())
        if contigs_a != contigs_b:
            errors.append(f"Data contigs differ: only in {label_a}={contigs_a - contigs_b}, only in {label_b}={contigs_b - contigs_a}")
            contigs = contigs_a & contigs_b
        else:
            contigs = contigs_a

        total_frags = 0
        for contig in sorted(contigs):
            ds_a = set(a[f"data/{contig}"].keys())
            ds_b = set(b[f"data/{contig}"].keys())
            if ds_a != ds_b:
                errors.append(f"{contig}: dataset keys differ: {ds_a} vs {ds_b}")
                continue

            for ds_name in sorted(ds_a):
                arr_a = a[f"data/{contig}/{ds_name}"][:]
                arr_b = b[f"data/{contig}/{ds_name}"][:]

                if arr_a.shape != arr_b.shape:
                    errors.append(f"{contig}/{ds_name}: shape {arr_a.shape} vs {arr_b.shape}")
                    continue

                if not numpy.array_equal(arr_a, arr_b):
                    mismatches = numpy.where(arr_a != arr_b)
                    n_mismatch = len(mismatches[0])
                    first_idx = mismatches[0][0]
                    errors.append(
                        f"{contig}/{ds_name}: {n_mismatch}/{len(arr_a)} values differ "
                        f"(first at index {first_idx}: {arr_a[first_idx]} vs {arr_b[first_idx]})"
                    )

            if f"data/{contig}/starts" in a:
                total_frags += len(a[f"data/{contig}/starts"])

        # Compare index
        idx_a = set(a["index"].keys()) if "index" in a else set()
        idx_b = set(b["index"].keys()) if "index" in b else set()
        if idx_a != idx_b:
            errors.append(f"Index contigs differ: only in {label_a}={idx_a - idx_b}, only in {label_b}={idx_b - idx_a}")
        for contig in sorted(idx_a & idx_b):
            ia = a[f"index/{contig}"][:]
            ib = b[f"index/{contig}"][:]
            if not numpy.array_equal(ia, ib):
                errors.append(f"index/{contig}: arrays differ")

        # Compare fragment length counts
        if "fragment_length_counts" in a and "fragment_length_counts" in b:
            flc_a = a["fragment_length_counts"][:]
            flc_b = b["fragment_length_counts"][:]
            if not numpy.array_equal(flc_a, flc_b):
                errors.append("fragment_length_counts differ")

        print(f"Compared {len(contigs)} contigs, {total_frags} total fragments")

    return errors


def build_with_docker(image_tag, bam_path, fasta_path, output_path, contigs=None, num_processes=4):
    """Build a fragments H5 using a Docker image of build-fragments-h5."""
    bam_path = os.path.abspath(bam_path)
    fasta_path = os.path.abspath(fasta_path)
    output_path = os.path.abspath(output_path)

    # Ensure the output directory is writable by the Docker container
    output_dir = os.path.dirname(output_path)
    os.chmod(output_dir, 0o777)

    # Collect all unique parent directories that need mounting
    dirs = set()
    for p in [bam_path, fasta_path, output_path]:
        dirs.add(os.path.dirname(p))

    cmd = ["docker", "run", "--rm"]
    for d in sorted(dirs):
        cmd += ["-v", f"{d}:{d}"]

    cmd += [
        f"ghcr.io/nboley/fragments-h5:{image_tag}",
        "build-fragments-h5",
        "--verbose",
        "--fasta", fasta_path,
        "--num-processes", str(num_processes),
        "--include-duplicates",
        bam_path, output_path,
    ]
    # --contigs uses nargs="+", so it must come after the positional args
    # to avoid consuming them as contig names
    if contigs:
        cmd += ["--contigs"] + contigs

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"STDOUT:\n{result.stdout}")
        print(f"STDERR:\n{result.stderr}")
        raise RuntimeError(f"Docker build failed with exit code {result.returncode}")


def run_comparison(label_a, path_a, label_b, path_b):
    """Run comparison and report results. Returns True if passed."""
    print(f"\n{'='*60}")
    print(f"Comparing: {label_a} vs {label_b}")
    print(f"{'='*60}")
    errors = compare_h5_files(path_a, path_b, label_a=label_a, label_b=label_b)
    if errors:
        print(f"\nFAILED — {len(errors)} difference(s):")
        for e in errors:
            print(f"  - {e}")
        return False
    else:
        print("PASSED — outputs are identical")
        return True


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bam", help="Input BAM file")
    parser.add_argument("fasta", help="Reference FASTA file")
    parser.add_argument("--num-processes", type=int, default=4)
    parser.add_argument("--include-duplicates", action="store_true", default=True)
    parser.add_argument("--contigs", nargs="+", default=None, help="Restrict to specific contigs")
    parser.add_argument("--compare-docker", metavar="TAG", default=None,
                        help="Also compare against a Docker image version (e.g. 2.7.2)")
    args = parser.parse_args()

    all_passed = True

    with tempfile.TemporaryDirectory(prefix="chunk_compare_") as tmpdir:
        chunked_path = os.path.join(tmpdir, "chunked.fragments.h5")
        unchunked_path = os.path.join(tmpdir, "unchunked.fragments.h5")

        print(f"Building chunked H5 (v2.8.0, num_processes={args.num_processes})...")
        build_fragments_h5(
            args.bam, chunked_path,
            fasta_filename=args.fasta,
            num_processes=args.num_processes,
            include_duplicates=args.include_duplicates,
            allowed_contigs=args.contigs,
            skip_chunking=False,
        )

        print(f"\nBuilding unchunked H5 (v2.8.0, num_processes={args.num_processes})...")
        build_fragments_h5(
            args.bam, unchunked_path,
            fasta_filename=args.fasta,
            num_processes=args.num_processes,
            include_duplicates=args.include_duplicates,
            allowed_contigs=args.contigs,
            skip_chunking=True,
        )

        if not run_comparison("v2.8.0-chunked", chunked_path, "v2.8.0-unchunked", unchunked_path):
            all_passed = False

        if args.compare_docker:
            docker_path = os.path.join(tmpdir, f"docker_{args.compare_docker}.fragments.h5")
            print(f"\nBuilding H5 with Docker image v{args.compare_docker}...")
            build_with_docker(
                args.compare_docker, args.bam, args.fasta, docker_path,
                contigs=args.contigs, num_processes=args.num_processes,
            )

            if not run_comparison(
                f"v2.8.0-chunked", chunked_path,
                f"v{args.compare_docker}-docker", docker_path,
            ):
                all_passed = False

    print(f"\n{'='*60}")
    if all_passed:
        print("ALL COMPARISONS PASSED")
    else:
        print("SOME COMPARISONS FAILED")
        sys.exit(1)


if __name__ == "__main__":
    main()
