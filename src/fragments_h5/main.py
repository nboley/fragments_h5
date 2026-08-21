import argparse
import os
import sys

import pysam

from fragments_h5.fragments_h5 import build_fragments_h5
from fragments_h5.fragment import is_fragment_file
import fragments_h5._logging as logging

logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(parents=[logging.build_log_parser()])
    parser.add_argument("input_file", metavar="INPUT_FILE", help="Input BAM or bgzipped TSV/BED fragment file")
    parser.add_argument("output_frags_h5", help="where to write the new fragments h5")

    parser.add_argument(
        "--fasta", default=None, help="Path to a fasta file containing the reference genome.",
    )
    parser.add_argument(
        "--contigs", default=None, nargs="+", help="Restrict building the fragment h5 over these contigs.",
    )
    parser.add_argument(
        "--contig-name-map", default=None,
        help="TSV file mapping input contig names to output names (two columns: input_name, output_name). "
             "Contigs not in the map are kept as-is. FASTA contig names must match the OUTPUT names.",
    )

    parser.add_argument("--set-mapq-255-to-none", action="store_true", help="set mapqs of 255 to None")
    parser.add_argument("--exclude-strand", default=False, action="store_true", help="Exclude strand info")
    parser.add_argument(
        "--read-methyl", default=False, action="store_true", help="Parse cpg's and converted cpg's from YN tag"
    )
    parser.add_argument(
        "--single-end", default=False, action="store_true", help="Sequencing is single ended (useful for long read technologies)"
    )
    parser.add_argument(
        "--se-max-fragment-length", type=int, default=None,
        help="Maximum fragment length to include in single-end mode (required with --single-end for BAM input). "
             "Fragments longer than this are excluded. Must be between 1 and 65535.",
    )
    parser.add_argument(
        "--min-mapq", type=int, default=None,
        help="Minimum mapping quality to include a fragment (default: 0, i.e. keep all). Must be >= 0.",
    )
    parser.add_argument(
        "--include-duplicates", default=False, action="store_true", help="Include duplicate-marked fragments in the output (default: exclude duplicates)"
    )
    parser.add_argument(
        "--num-processes", default='1', help="Num of processes to use (defaults to 1 -- use 'all' for all cores)"
    )
    parser.add_argument(
        "--no-store-fragment-end-clipped", dest="store_fragment_end_clipped",
        action="store_false", default=True,
        help="Do not store fragment_end_clipped flag (default: store it)",
    )
    parser.add_argument(
        "--skip-chunking", default=False, action="store_true",
        help="Disable chunk-based parallelization and process each contig as a whole",
    )

    return parser, parser.parse_args()


def _is_remote_url(path: str) -> bool:
    """True if path looks like a remote URL (S3, HTTP) that should not be indexed via samtools."""
    return path.startswith("s3://") or path.startswith("http://") or path.startswith("https://")


def main():
    parser, args = parse_args()

    logging.configure_root_logger_from_args(args)

    # Argument-consistency checks run before the input-file indexing below,
    # which can shell out to samtools/tabix and write a .bai/.tbi next to the
    # input — that must not happen for an invocation we are about to reject.
    input_is_tsv = is_fragment_file(args.input_file)

    # Range validation
    if args.se_max_fragment_length is not None:
        if args.se_max_fragment_length < 1 or args.se_max_fragment_length > 65535:
            parser.error("--se-max-fragment-length must be between 1 and 65535")
    if args.min_mapq is not None and args.min_mapq < 0:
        parser.error("--min-mapq must be >= 0")

    # For BAM input, --se-max-fragment-length is required with --single-end
    # and --se-max-fragment-length requires --single-end.
    # For TSV/BED input, --single-end and --se-max-fragment-length are ignored
    # (the library warns and neutralizes them).
    if not input_is_tsv:
        if args.single_end and args.se_max_fragment_length is None:
            parser.error("--se-max-fragment-length is required when using --single-end")
        if args.se_max_fragment_length is not None and not args.single_end:
            parser.error("--se-max-fragment-length can only be used with --single-end")

    if os.path.exists(args.output_frags_h5):
        logger.error(
            f"Output file '{args.output_frags_h5}' already exists. "
            f"Please remove it or specify a different output path."
        )
        sys.exit(1)

    if not is_fragment_file(args.input_file):
        with pysam.AlignmentFile(args.input_file) as bam:
            if not bam.has_index():
                if _is_remote_url(args.input_file):
                    raise SystemExit(
                        "Remote BAM (s3:// or http(s)://) has no index. "
                        "build-fragments-h5 requires an indexed BAM; ensure the .bai exists "
                        "(e.g. s3://bucket/file.bam.bai for s3://bucket/file.bam)."
                    )
                import subprocess
                subprocess.run(["samtools", "index", args.input_file], check=True)
    else:
        try:
            with pysam.TabixFile(args.input_file) as _tbx:
                pass
        except (OSError, ValueError):
            if _is_remote_url(args.input_file):
                raise SystemExit(
                    "Remote TSV/BED file has no tabix index (.tbi). "
                    "Create one with 'tabix -p bed <file>' before uploading."
                )
            import subprocess
            subprocess.run(["tabix", "-p", "bed", args.input_file], check=True)

    # Validate FASTA file if provided
    if args.fasta:
        if _is_remote_url(args.fasta):
            # For remote FASTA files, verify they can be opened (pysam will check for index)
            try:
                with pysam.FastaFile(args.fasta) as fasta:
                    _ = fasta.references  # Access to verify it's readable
            except Exception as e:
                raise SystemExit(
                    f"Failed to open remote FASTA file '{args.fasta}'. "
                    f"Ensure the index file (.fai) exists at the same S3 path. "
                    f"For compressed FASTA (.gz), the .gzi index is also required. "
                    f"Error: {e}"
                )

    # TODO: Add validation for num_processes with user-friendly error message
    if args.num_processes.lower() == 'all':
        num_processes = None
    else:
        num_processes = int(args.num_processes)

    # Parse contig name map if provided
    contig_name_map = None
    if args.contig_name_map:
        contig_name_map = {}
        with open(args.contig_name_map) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) != 2:
                    raise SystemExit(
                        f"Invalid contig name map line (expected 2 tab-separated columns): {line!r}"
                    )
                contig_name_map[parts[0]] = parts[1]

    build_fragments_h5(
        args.input_file,
        args.output_frags_h5,
        fasta_filename=args.fasta,
        allowed_contigs=args.contigs,
        set_mapq_255_to_none=args.set_mapq_255_to_none,
        read_strand=not args.exclude_strand,
        read_methyl=args.read_methyl,
        single_end=args.single_end,
        se_max_fragment_length=args.se_max_fragment_length,
        num_processes=num_processes,
        include_duplicates=args.include_duplicates,
        store_fragment_end_clipped=args.store_fragment_end_clipped,
        skip_chunking=args.skip_chunking,
        contig_name_map=contig_name_map,
        min_mapq=args.min_mapq,
    )


if __name__ == "__main__":
    main()