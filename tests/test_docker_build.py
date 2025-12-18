#!/usr/bin/env python
"""
Integration test to verify Docker builds produce identical results to local builds.

This test is NOT run by pytest automatically. Run it manually with:
    python tests/test_docker_build.py

Requirements:
    - Docker must be installed and running
    - The fragments_h5 package must be installed locally
    - Test data must be available in tests/data/
"""
import os
import sys
import subprocess
import tempfile
import hashlib
import h5py
import numpy as np
from pathlib import Path


# Paths
SCRIPT_DIR = Path(__file__).parent.absolute()
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_DIR = SCRIPT_DIR / "data"
TEST_BAM = DATA_DIR / "small.chr6.bam"
TEST_FASTA = DATA_DIR / "GRCh38.p12.genome.chr6_99110000_99130000.fa.gz"

DOCKER_IMAGE = "fragments-h5:test"


def run_command(cmd, check=True, capture_output=True):
    """Run a shell command and return the result."""
    print(f"Running: {cmd}")
    result = subprocess.run(
        cmd, shell=True, check=check, capture_output=capture_output, text=True
    )
    if result.returncode != 0 and not check:
        print(f"Command failed with code {result.returncode}")
        print(f"stdout: {result.stdout}")
        print(f"stderr: {result.stderr}")
    return result


def build_docker_image():
    """Build the Docker image."""
    print("\n=== Building Docker image ===")
    result = run_command(
        f"docker build -t {DOCKER_IMAGE} {PROJECT_ROOT}",
        capture_output=False
    )
    return result.returncode == 0


def build_local_h5(output_path):
    """Build fragments H5 using local installation."""
    print("\n=== Building H5 with local installation ===")
    cmd = f"build-fragments-h5 {TEST_BAM} {output_path} --fasta {TEST_FASTA}"
    run_command(cmd)


def build_docker_h5(output_path):
    """Build fragments H5 using Docker container."""
    print("\n=== Building H5 with Docker ===")
    output_dir = Path(output_path).parent
    output_name = Path(output_path).name
    
    # Mount the data directory and output directory
    cmd = (
        f"docker run --rm "
        f"-v {DATA_DIR}:/data/input:ro "
        f"-v {output_dir}:/data/output "
        f"{DOCKER_IMAGE} "
        f"/data/input/{TEST_BAM.name} /data/output/{output_name} "
        f"--fasta /data/input/{TEST_FASTA.name}"
    )
    run_command(cmd)


def compare_h5_files(file1, file2):
    """Compare two H5 files and return whether they are equivalent."""
    print(f"\n=== Comparing H5 files ===")
    print(f"  File 1: {file1}")
    print(f"  File 2: {file2}")
    
    differences = []
    
    with h5py.File(file1, 'r') as f1, h5py.File(file2, 'r') as f2:
        # Compare attributes
        for attr in ['index_block_size', 'max_fragment_length']:
            v1 = f1.attrs.get(attr)
            v2 = f2.attrs.get(attr)
            if v1 != v2:
                differences.append(f"Attribute '{attr}' differs: {v1} vs {v2}")
        
        # Compare contig lengths (stored as string, so compare parsed)
        cl1 = eval(f1.attrs.get('_contig_lengths_str', '{}'))
        cl2 = eval(f2.attrs.get('_contig_lengths_str', '{}'))
        if cl1 != cl2:
            differences.append(f"Contig lengths differ: {cl1} vs {cl2}")
        
        # Compare fragment length counts
        if 'fragment_length_counts' in f1 and 'fragment_length_counts' in f2:
            flc1 = f1['fragment_length_counts'][:]
            flc2 = f2['fragment_length_counts'][:]
            if not np.array_equal(flc1, flc2):
                differences.append("Fragment length counts differ")
        
        # Compare data for each contig
        contigs1 = set(f1['data'].keys()) if 'data' in f1 else set()
        contigs2 = set(f2['data'].keys()) if 'data' in f2 else set()
        
        if contigs1 != contigs2:
            differences.append(f"Different contigs: {contigs1} vs {contigs2}")
        
        for contig in contigs1 & contigs2:
            for dataset in ['starts', 'lengths', 'mapq', 'gc', 'strand']:
                path = f'data/{contig}/{dataset}'
                if path in f1 and path in f2:
                    d1 = f1[path][:]
                    d2 = f2[path][:]
                    if not np.array_equal(d1, d2):
                        differences.append(f"Dataset '{path}' differs")
                elif path in f1 or path in f2:
                    differences.append(f"Dataset '{path}' exists in only one file")
        
        # Compare index
        if 'index' in f1 and 'index' in f2:
            idx_contigs1 = set(f1['index'].keys())
            idx_contigs2 = set(f2['index'].keys())
            if idx_contigs1 != idx_contigs2:
                differences.append(f"Different indexed contigs: {idx_contigs1} vs {idx_contigs2}")
            
            for contig in idx_contigs1 & idx_contigs2:
                i1 = f1[f'index/{contig}'][:]
                i2 = f2[f'index/{contig}'][:]
                if not np.array_equal(i1, i2):
                    differences.append(f"Index for '{contig}' differs")
    
    return differences


def main():
    """Run the Docker comparison test."""
    print("=" * 60)
    print("Docker Build Comparison Test")
    print("=" * 60)
    
    # Check prerequisites
    if not TEST_BAM.exists():
        print(f"ERROR: Test BAM not found: {TEST_BAM}")
        sys.exit(1)
    
    if not TEST_FASTA.exists():
        print(f"ERROR: Test FASTA not found: {TEST_FASTA}")
        sys.exit(1)
    
    # Check Docker is available
    result = run_command("docker --version", check=False)
    if result.returncode != 0:
        print("ERROR: Docker is not installed or not running")
        sys.exit(1)
    
    # Check local installation
    result = run_command("build-fragments-h5 --help", check=False)
    if result.returncode != 0:
        print("ERROR: fragments_h5 is not installed locally")
        print("Install with: pip install -e .")
        sys.exit(1)
    
    # Build Docker image
    if not build_docker_image():
        print("ERROR: Failed to build Docker image")
        sys.exit(1)
    
    # Create temp directory for outputs
    with tempfile.TemporaryDirectory() as tmpdir:
        local_h5 = Path(tmpdir) / "local.fragments.h5"
        docker_h5 = Path(tmpdir) / "docker.fragments.h5"
        
        # Build with both methods
        build_local_h5(local_h5)
        build_docker_h5(docker_h5)
        
        # Compare results
        differences = compare_h5_files(local_h5, docker_h5)
        
        if differences:
            print("\n" + "=" * 60)
            print("FAILED: Files differ!")
            print("=" * 60)
            for diff in differences:
                print(f"  - {diff}")
            sys.exit(1)
        else:
            print("\n" + "=" * 60)
            print("PASSED: Docker and local builds produce identical results!")
            print("=" * 60)
    
    # Cleanup: remove test image
    print("\n=== Cleaning up ===")
    run_command(f"docker rmi {DOCKER_IMAGE}", check=False)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

