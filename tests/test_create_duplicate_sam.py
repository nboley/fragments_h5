"""Test script to create a SAM file with duplicate-marked reads.

This script generates a SAM file containing paired-end reads where one fragment
is marked as a duplicate. It is used to create test data for the 
test_include_duplicates unit test, which verifies that the include_duplicates
parameter correctly includes/excludes duplicate-marked fragments.

The reads are placed in the non-N region of chr6 (99110000-99130000) as indicated
by the FASTA filename GRCh38.p12.genome.chr6_99110000_99130000.fa.gz, and the
sequences match the actual FASTA sequence at those positions.

Usage:
    # Generate SAM file
    python tests/test_create_duplicate_sam.py > test_duplicates.sam
    
    # Convert to BAM, sort, and index
    samtools view -bS test_duplicates.sam > test_duplicates.bam
    samtools sort -o test_duplicates.sorted.bam test_duplicates.bam
    samtools index test_duplicates.sorted.bam
    
    # Copy to test data directory
    cp test_duplicates.sorted.bam tests/data/test_duplicates.bam
    cp test_duplicates.sorted.bam.bai tests/data/test_duplicates.bam.bai

    # View with info (shows flag explanations)
    python tests/test_create_duplicate_sam.py --show-info
"""

def create_duplicate_sam_string():
    """Create a SAM file string with duplicate-marked reads.
    
    Reads are placed in the non-N region of chr6 (99110000-99130000) as indicated
    by the FASTA filename GRCh38.p12.genome.chr6_99110000_99130000.fa.gz
    """
    # Non-N region from FASTA filename: chr6:99110000-99130000
    # Place reads within this region
    read1_start = 99110001  # Start of non-N region + 1 (1-based SAM coordinates)
    read1_mate_start = 99110001 + 116  # Fragment length 116bp
    
    # Sequences from FASTA at these positions
    read1_seq = "AACAAGAGCTTGCCGA"  # Forward strand at 99110001-99110016
    read2_seq = "GCCTGCTTTACTTTAA"  # Reverse complement of 99110117-99110132 (reverse strand)
    qual = "IIIIIIIIIIIIIIII"  # Quality scores
    
    sam_lines = [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr6\tLN:170805979",
        # First fragment (non-duplicate): read1
        f"read1\t99\tchr6\t{read1_start}\t60\t16M\t=\t{read1_mate_start}\t116\t{read1_seq}\t{qual}\tMQ:i:60",
        # First fragment (non-duplicate): read2 (mate)
        f"read1\t147\tchr6\t{read1_mate_start}\t60\t16M\t=\t{read1_start}\t-116\t{read2_seq}\t{qual}\tMQ:i:60",
        # Second fragment (duplicate): read1 with duplicate flag (0x400 = 1024)
        f"read2\t1123\tchr6\t{read1_start}\t60\t16M\t=\t{read1_mate_start}\t116\t{read1_seq}\t{qual}\tMQ:i:60",
        # Second fragment (duplicate): read2 (mate) with duplicate flag
        f"read2\t1651\tchr6\t{read1_mate_start}\t60\t16M\t=\t{read1_start}\t-116\t{read2_seq}\t{qual}\tMQ:i:60",
    ]
    return "\n".join(sam_lines) + "\n"


if __name__ == "__main__":
    import sys
    sam_content = create_duplicate_sam_string()
    
    if len(sys.argv) > 1 and sys.argv[1] == "--show-info":
        print("SAM file content:")
        print("=" * 80)
        print(sam_content)
        print("=" * 80)
        # Flag explanations:
        print("\nFlag explanations:")
        print("99  = 0x63 = paired + read1 + proper_pair + forward")
        print("147 = 0x93 = paired + read2 + proper_pair + reverse")
        print("1123 = 0x463 = 99 + 0x400 (duplicate) = paired + read1 + proper_pair + forward + duplicate")
        print("1651 = 0x673 = 147 + 0x400 (duplicate) = paired + read2 + proper_pair + reverse + duplicate")
    else:
        # Just output SAM content for piping
        print(sam_content, end='')
