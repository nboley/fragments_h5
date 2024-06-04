import argparse
import os.path


from fragments_h5.fragments_h5 import build_fragments_h5
import fragments_h5._logging as logging


REFERENCE_ANNOTATIONS = ["hg16", "hg17", "hg18", "hg19", "hg38"]


def parse_args():
    parser = argparse.ArgumentParser(parents=[logging.build_log_parser()])
    parser.add_argument("input_bam", help="bam file to read fragments from")
    parser.add_argument("output_frags_h5", help="where to write the new fragments h5")

    parser.add_argument(
        "--reference",
        choices=sorted(REFERENCE_ANNOTATIONS),
        required=True,
        help="The reference genome of input_bam (e.g. hg38).",
    )
    parser.add_argument(
        "--sample-id",
        required=True,
        help="Unique identifier (typically a sample id or cell type).",
    )
    parser.add_argument(
        "--fasta", default=None, help="Path to a fasta file containing the reference genome.",
    )
    parser.add_argument(
        "--contigs", default=None, nargs="+", help="Restrict building the fragment h5 over these contigs.",
    )

    parser.add_argument("--set-mapq-255-to-none", action="store_true", help="set mapqs of 255 to None")
    parser.add_argument("--exclude-strand", default=False, action="store_true", help="Exclude strand info")
    parser.add_argument(
        "--read-methyl", default=False, action="store_true", help="Parse cpg's and converted cpg's from YN and YC tags"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    logging.configure_root_logger_from_args(args)

    if args.input_bam.endswith(".bam") and (
            not os.path.isfile(args.input_bam + ".bai")
    ):
        import subprocess
        subprocess.run(f"samtools index {args.input_bam}", shell=True, check=True)

    build_fragments_h5(
        args.input_bam,
        args.output_frags_h5,
        args.sample_id,
        reference=args.reference,
        fasta_file=args.fasta,
        allowed_contigs=args.contigs,
        set_mapq_255_to_none=args.set_mapq_255_to_none,
        read_strand=not args.exclude_strand,
        read_methyl=args.read_methyl,
    )


if __name__ == "__main__":
    main()
