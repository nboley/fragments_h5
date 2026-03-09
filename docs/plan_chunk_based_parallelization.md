# Chunk-Based Parallelization Refactor

## Goal

Replace the current per-contig multiprocessing strategy with fixed-size genomic chunks distributed via a work queue. This balances worker load — currently, workers assigned large contigs (e.g., chr1 at ~249M bases) bottleneck the entire build while workers for small contigs (e.g., chrM at ~16K bases) finish instantly.

## Design

### Chunking Strategy

- **Chunk size:** 10M bases (constant, configurable).
- **Chunk creation:** For each contig, divide `[0, contig_length)` into consecutive `[chunk_start, chunk_stop)` intervals of 10M bases. The final chunk may be shorter.
- **Small contigs:** Contigs shorter than the chunk size remain as a single chunk (no splitting).
- **`--contigs` filter:** Only create chunks for user-specified contigs when this flag is set.
- **Work queue:** All chunks across all contigs are placed into a single flat list and distributed to workers via `imap_unordered`.

### Fragment Ownership

- A fragment belongs to the chunk containing its **start position**.
- Since fragments are sorted by start in the BAM, and chunks are non-overlapping contiguous intervals, each fragment is assigned to exactly one chunk.
- Concatenating chunk results in chunk order produces a globally sorted array per contig (no re-sort needed).

### BAM Reading Per Chunk

- Use `pysam.fetch(contig, chunk_start, chunk_stop)` to get reads overlapping the chunk region.
- This returns a superset — filter to fragments where `fragment_start >= chunk_start` and `fragment_start < chunk_stop`.
- Since we only keep reads with `tlen > 0`, `read.pos` equals `fragment_start`, so all relevant reads are captured by the fetch region.

### FASTA / GC Handling Per Chunk

- Fetch only the needed FASTA region: `fasta.fetch(contig, chunk_start, min(chunk_stop + MAX_FRAG_LENGTH, contig_length))`.
- The extra `MAX_FRAG_LENGTH` (65535) buffer accounts for fragments whose start is near `chunk_stop` but whose end extends beyond it.
- Build the GC cumulative sum array from this sub-sequence. The array is zero-indexed relative to `chunk_start`.
- Offset all genomic coordinates: `gc = (cumsum[frag_stop - chunk_start] - cumsum[frag_start - chunk_start]) / length`.

### Worker Task Signature

Each worker task becomes a tuple like:

`(contig, chunk_start, chunk_stop, bam_path, fasta_path, tmp_dir, ...options)`

Replacing the current per-contig signature which has no start/stop.

### Temp File Naming

Change from `tmp.fragment_h5.{contig}.h5` to `tmp.fragment_h5.{contig}.{chunk_start}.h5` to avoid collisions between chunks of the same contig.

### Merge Process

The merge step becomes more complex since multiple chunks may contribute to the same contig.

1. Group completed chunk results by contig.
2. Sort each contig's chunks by `chunk_start`.
3. For each contig, **concatenate** the arrays (starts, lengths, mapq, gc, strand, etc.) from all chunks in order.
4. Write the concatenated arrays into `/data/{contig}/` in the output HDF5.
5. Build the block-based index on the merged per-contig data (unchanged logic).
6. Add fragment length distribution (unchanged).

### Progress Bar

- tqdm updates weighted by chunk size (genomic bases), same as current contig-length weighting.
- Since chunks are mostly equal size (~10M), progress will advance more smoothly than the current per-contig approach.

## Files to Modify

- **`src/fragments_h5/fragments_h5.py`** — Chunk generation, updated worker function, new merge logic (concatenation), updated multiprocessing dispatch.
- **`src/fragments_h5/fragment.py`** — Update `get_g_or_c_cumsum()` and `bam_to_fragments()` to accept region bounds for FASTA sub-fetching and coordinate offsetting.
- **`src/fragments_h5/main.py`** — No changes expected (CLI interface unchanged).

## Testing

- **Merge correctness test:** Build the same BAM with chunk-based parallelization (multiple chunks per contig) and with single-process mode. Compare all datasets in the output HDF5 files — they must be identical.
- **Boundary fragment test:** Verify fragments near chunk boundaries are assigned to exactly one chunk and appear correctly in the merged output.
- **GC accuracy test:** Verify GC values from chunk-based build match a single-process build (validates the FASTA offset logic).
- **Small contig test:** Verify contigs smaller than chunk size are handled correctly as single chunks.
- **Edge cases:** Empty chunks (genomic regions with no mapped reads), single-chunk contigs, contigs with reads only at boundaries.
