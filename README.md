# FLoRe: Fingerprint-based Long-Read Overlap Reconstruction

FLoRe is an open-source tool for fast and accurate detection of overlaps between long DNA sequencing reads (common in third-generation sequencing data). It introduces a novel fingerprinting approach to efficiently identify candidate overlapping reads, followed by multiple verification strategies to confirm true overlaps even in the presence of high error rates. FLoRe encodes each k-mer (substring of length k) as a unique integer (bijective k-mer encoding) and compresses these fingerprints via run-length encoding to eliminate local redundancies. A compact inverted index then maps each fingerprint value to the list of reads containing it, enabling rapid pre-filtering of promising read pairs. For overlap validation, FLoRe employs a hierarchy of seed-and-extend algorithms – FGOE, AOE, PSH, and CORE – with a suffix-automaton-based fallback to ensure high accuracy and robustness against sequencing errors. The tool is implemented in modern C++ with multi-threading support and includes an integrated profiler to measure execution time and peak memory usage.

## Key Features and Algorithms

### Bijective k-mer Fingerprinting
Each DNA k-mer is encoded into a unique integer ID (using a 2-bit per base representation), ensuring collision-free fingerprints and preserving lexicographic order. This provides an efficient numerical representation of reads for quick comparison.

### Run-Length Compression
The sequence of k-mer fingerprints in each read is compressed by collapsing consecutive identical IDs, which removes local repeats (e.g. homopolymers or low-complexity regions). This compressed fingerprint list is much shorter than the raw sequence and retains the essential pattern of each read.

### Inverted Index
FLoRe builds an inverted index mapping each fingerprint ID to all reads (and orientations) in which it appears. This index allows fast lookup of candidate overlaps – read pairs that share common fingerprints – drastically reducing the number of pairs to examine in detail.

### Candidate Filtering
Using the index, the algorithm quickly intersects fingerprint sets to find read pairs with a significant number of shared fingerprints. Efficient set-intersection methods (e.g. hybrid sorted intersection) are used to handle these comparisons at scale, filtering out most non-overlapping pairs before expensive alignment steps.

### Overlap Validation Strategies
FLoRe validates each candidate pair through a series of seed-and-extend methods:

**FGOE (Fingerprint-Guided Overlap Extension):** a fast, greedy extension using matching fingerprint seeds at the read ends to catch short overlaps quickly (very stringent).

**AOE (Adaptive Overlap Extension):** a more flexible overlap search that scans the fingerprint lists for longer common segments, allowing detection of larger overlaps (though with higher worst-case cost).

**PSH (Progressive Seed Hopping):** an approach that iteratively tries smaller k-mer seeds to bridge gaps between matches, improving sensitivity for gapped overlaps.

**CORE (Combined Overlap Refinement Engine):** integrates results from FGOE, AOE, and PSH – taking the best candidate alignment “seeds” from each – and refines them using a lightweight dynamic programming step. CORE effectively stitches together partial overlaps into a consistent alignment.

If all the above methods fail to confirm an overlap, a suffix automaton fallback is used to find the Longest Common Substring between the two reads in linear time. This guarantees that any true overlap (common segment) will be found, at the cost of a heavier computation in these rare cases.

**High Error Tolerance:** The multi-strategy validation ensures that true overlaps are distinguished from spurious matches due to sequencing errors or repetitive sequences. By progressing from quick, specific checks (FGOE) to more exhaustive methods (suffix automaton) only as needed, FLoRe achieves both speed and high sensitivity.

**Performance and Profiling:** The C++ implementation is optimized with multithreading (overlap computations are done in parallel). FLoRe logs execution time and peak memory usage for each run, helping users evaluate performance on different datasets. It has been shown to outperform baseline tools in speed while maintaining or improving overlap detection accuracy (see Performance below).

## Workflow Overview

The FLoRe pipeline consists of four main stages:

### Fingerprint Generation
Each input read is converted into a sequence of k-mer fingerprints. Every k-mer (subsequence of length k) is mapped to a unique integer in the range 0 to 4^k–1 using a 2-bit encoding per nucleotide. (For example, A=00, C=01, G=10, T=11 in binary.) This bijective encoding guarantees that no two different k-mers produce the same fingerprint. The resulting fingerprint list preserves the order of k-mers along the read, allowing fast, direct comparisons between reads using integer operations. To further condense the data, FLoRe applies run-length encoding to the fingerprint list, collapsing any run of identical fingerprint IDs into a single instance. This step removes redundancies from low-complexity regions (like long runs of AAAA... or tandem repeats) which would otherwise produce many repeated k-mer entries.

### Index Construction & Candidate Retrieval
As fingerprints are generated, FLoRe builds a compact inverted index that records, for each unique fingerprint value, the list of read IDs (and orientation) in which it appears. Once all reads are indexed, FLoRe uses this structure to quickly find candidate overlapping pairs. Essentially, if two reads share a fingerprint ID (meaning they have an identical k-mer in some position), they are potential overlap candidates. The index allows FLoRe to retrieve all reads sharing a given fingerprint in O(1) time, and by intersecting the fingerprint sets of reads, it identifies pairs that share a sufficient number of fingerprints in common. This pre-filtering stage prunes the search space drastically: only read pairs with enough common k-mer signals pass to the next stage, while unrelated reads (no significant shared k-mers) are discarded early. Efficient set intersection algorithms (such as hybrid sorted intersection) are used to perform these operations quickly, even for large sets.

### Overlap Detection (Seed-and-Extend)
Each candidate pair from Stage 2 is then examined in detail to determine if an actual sequence overlap exists and to compute its exact coordinates. FLoRe applies a chain of seed-and-extend algorithms of increasing thoroughness:

**Fingerprint-Guided Overlap Extension (FGOE):** A rapid check that tries to extend an overlap assuming the reads align at their ends. For example, FGOE will compare the last few fingerprint IDs of read A with the first few of read B (or vice versa) to see if they match contiguously. If a short overlap (above the minimum length) is found immediately at the ends, FGOE reports it. This method is extremely fast and strict; if the required minimum overlap length is more than 3 fingerprints, FGOE might not find anything. It’s mainly useful for catching very short overlaps or obvious alignments.

**Adaptive Overlap Extension (AOE):** If FGOE fails, FLoRe uses AOE, which scans more of the fingerprint lists to find a longer matching segment that could indicate an overlap. AOE is essentially a more flexible search across the two reads’ fingerprints, able to handle overlaps starting at internal positions. It can find overlaps longer than the FGOE limit but is computationally heavier (worst-case O(m·n) for fingerprint list lengths m, n).

**Progressive Seed Hopping (PSH):** This strategy attempts to bridge gaps between matching regions by progressively lowering the seed resolution (using smaller k-mers as needed). In practice, PSH can detect overlaps when there are mismatches or errors splitting the alignment into pieces – it “hops” through the reads, trying shorter k-mer matches to connect partial alignments. (In the thesis, this is referred to as k-mer hopping.)

**Combined Overlap Refinement Engine (CORE):** Finally, CORE takes the best candidate alignments (seeds) identified by FGOE, AOE, and PSH and merges them to produce a refined overlap alignment. It might select up to three candidate segments (one from each method) and perform a light-weight dynamic programming merge (similar to stitching subsequences) to see if they form a consistent overlap. CORE essentially ensures that if there are partial overlaps found by earlier methods, they are integrated into one longer overlap, maximizing the overlap length.

After these methods, if an overlap is confirmed by any of them, FLoRe records it. If none of the above methods yield a valid overlap (which could happen if the overlap is highly fragmented by errors), FLoRe triggers a suffix automaton fallback. In this step, a suffix automaton is built for one read’s sequence and used to find the Longest Common Substring (LCS) between the two reads. The LCS corresponds to the longest possible overlap. This fallback runs in linear time relative to the read lengths and guarantees that no true overlap is missed, albeit at higher computational cost, so it is used only as a last resort. By combining these approaches from fastest/most specific to most general, Stage 3 achieves both speed and high recall.

### Output and Post-processing
Each verified overlap is then output in a structured format (see Usage below). FLoRe records the read identifiers, overlap positions, length, and other details for each overlap found. It also notes which algorithm confirmed the overlap (e.g., FGOE, AOE, PSH, CORE, or suffix automaton) as part of the result metadata. This is useful for debugging or analyzing which overlaps required the heavier fallback approach. The tool’s built-in profiling module also logs the execution time and peak memory used for the run, providing insight into performance. Optionally, the JSON results can be converted to a human-friendly Excel spreadsheet automatically if Python is available (this is handled by the provided wrapper script). After this stage, the set of overlaps can be used in downstream applications (such as assembly graph construction or error correction).

## Installation and Requirements

**Language and Standards:** FLoRe is written in C++ (C++17 standard) for high performance. You will need a C++17-compatible compiler such as g++ or clang++ to build it  
GitHub  
The code has been tested primarily on Linux (e.g., Ubuntu) and macOS. It uses the C++ Standard Library threads for parallelism and OS-specific calls for memory usage statistics. Windows support is partially included (there are code paths for Windows memory profiling  
GitHub  
), but building on Windows may require a POSIX compatibility layer or using WSL (Windows Subsystem for Linux) because the code uses POSIX APIs like <getopt.h> for command-line parsing.

**Prerequisites:** No external libraries are required for core functionality – the implementation relies only on the standard C++ library (threads, data structures, etc.). For convenience, a Python 3 installation is recommended if you wish to use the automated JSON-to-Excel conversion (pandas and openpyxl will be used by the script, but they are installed on-the-fly if missing). Python is not required to find overlaps, only for formatting results into Excel.

**Supported Platforms:** The tool is developed and tested on 64-bit Linux and should compile on macOS as well (it includes macOS-specific memory queries via Mach APIs  
GitHub  
). On Windows, compilation with MinGW or MSYS2 may work, but due to the use of Unix-specific functions, the simplest method is to compile and run FLoRe under WSL or a Linux container.

## Building from Source

**Clone the Repository:**

```bash
git clone https://github.com/manuelsica/FLoRe.git && cd FLoRe
```

**Compile: A Makefile is provided. Simply run:**

```bash
make
```

This will produce an executable binary named FLORE_bin in the project directory. The Makefile uses g++ with -std=c++17 -O3 optimizations enabled  
All source files (for fingerprinting, indexing, overlap algorithms, etc.) will be compiled and linked into FLORE_bin. If the build succeeds, the make process will also set the execute permission on the binary automatically  
**(Optional) Install:** You can copy the FLORE_bin binary and the FLORE runner script to a directory in your $PATH for system-wide use. However, this is not required; you can also run it from the build directory.

No special installation steps are needed beyond compilation. The repository includes all necessary source code. Ensure your system has enough memory for the dataset you plan to analyze (long-read overlap detection can be memory-intensive for large genomes). For example, overlapping a few thousand long reads (length 10kbp+) may require a few GB of RAM.

## Usage

FLoRe operates on genomic read sequences provided by the user and outputs the overlaps it finds. The basic usage is:

```bash
./FLORE -f <reads.fasta> [options]
```

You can also use the provided FLORE wrapper script (note the capital letters) which invokes FLORE_bin and additionally converts the JSON output to an Excel file if Python is available. For example, running:

```bash
./FLORE -f reads.fasta -k 15 -t 8 -j overlaps.json
```

will process the reads from reads.fasta using k=15, employ 8 threads, and save the detected overlaps to overlaps.json (and also produce overlaps.json.xlsx if possible).

**Command-line Options:** FLoRe supports several parameters to tune its behavior (shown in the usage help). Key options include:

- **-f, --fasta <file> (required):** Input file in FASTA format containing the reads to process. Each read in the FASTA should have a unique identifier. (Multi-FASTA with many reads is supported.)
- **-k, --kmer <int>:** k-mer length to use for fingerprinting. Default is 15  
GitHub  
Shorter k-mers increase sensitivity (more potential overlaps found) at the cost of more candidates and possibly more false hits; longer k-mers make the fingerprints more specific. Typically 15-19 is a reasonable range for PacBio/ONT reads.
- **-m, --min_overlap <int>:** Minimum overlap length (in base pairs) to report. Default is 13 bp  
GitHub  
Overlaps shorter than this threshold will be ignored. You might raise this for very long reads to avoid trivial overlaps.
- **-r, --max_repeat_threshold <int>:** Maximum allowed run-length for a single fingerprint ID. Default 10  
GitHub  
This is essentially a filter for reads with very long homopolymers or low complexity – if a fingerprint value repeats more than this many times in a row in a read, FLoRe might break the fingerprint or flag it (to avoid overly repetitive candidate matches). In most cases the default is fine.
- **-t, --threads <int>:** Number of parallel threads to use. By default, FLoRe uses all available hardware threads (as reported by the system)  
GitHub  
You can cap this to avoid using all CPU cores. Multi-threading greatly speeds up Stage 3 (overlap validation), as different read pairs can be processed concurrently.
- **-j, --json <file>:** Output filename for the JSON results. Default is results.json  
GitHub  
You can specify a path to direct the overlap output to a particular location. The JSON will contain an array of overlap records (described below). If you use the wrapper script, it will also attempt to create an Excel (.xlsx) version of this file automatically.
- **-v, --verbose:** Enable verbose logging. By default, FLoRe prints minimal info to the console. With -v, it will log progress and debug information (and also save a detailed log to a file). This can help in understanding how many candidates were filtered, how many overlaps each method found, etc., but may slow down execution slightly due to I/O.
- **--solid_fingerprint:** Enable “solid” fingerprint mode. In this mode, FLoRe will ignore extremely rare or extremely frequent k-mers when building fingerprints. Specifically, you can set --solid_min_freq <n> and --solid_max_freq <m> to only consider k-mers that occur in at least n reads and at most m reads (across the dataset). By default, solid_min_freq=1 and solid_max_freq=100000 (effectively no filtering)  
GitHub  
For large datasets, increasing solid_min_freq to 2 (for example) filters out k-mers that occur in only one read (which are often error k-mers) and can greatly reduce false candidate overlaps at a negligible risk of missing true overlaps. Similarly, lowering solid_max_freq can remove overly common k-mers (e.g., from repetitive regions) that would cause many spurious matches. Use this feature with caution, as it may drop overlaps that rely on those k-mers if set too strictly.

**Fingerprint Variants (CFL/ICFL):** By default, FLoRe uses a CFL (Canonical Lyndon factor) compression in the fingerprint stage, an advanced technique that can further compress fingerprint sequences by leveraging Lyndon factorization of reads (a string decomposition approach). The tool offers experimental flags to control this:

- **--no_cfl** – disable Lyndon factor compression, using the plain run-length fingerprint instead (this might increase runtime/memory but could be more sensitive in some cases).
- **--cfl_threshold <int>** – length threshold for “long” Lyndon factors (default 30). Lyndon factors longer than this are treated specially in compression.
- **--icfl** – enable inverse Lyndon factor grouping (ICFL) as an alternative fingerprinting strategy.
- **--cfl_icfl** – combine both CFL and ICFL strategies for fingerprinting.  
These options correspond to variations of Stage 1 discussed in the thesis (they are mainly for research and ablation studies). In practice, the default settings (CFL on) should work well; turning these off or on might be useful for evaluating performance on extremely repetitive sequences.

- **-h, --help:** Show the full usage information and exit.

### Input Format
FLoRe expects reads in FASTA format (text file with sequences). Each read should have a header line beginning with > followed by an identifier, and one or more lines of nucleotide sequence. Multi-line sequences are handled, as are multi-read files. (FASTQ format reads can be converted to FASTA or provided directly; if a FASTQ is given, FLoRe will treat lines starting with @ as headers similarly to >). All DNA bases A,C,G,T (and lowercase) are recognized; other IUPAC nucleotide codes are not explicitly handled and should be converted or removed for best results. Ensure the reads are long-read sequences (e.g., 1,000 bp or longer) for the overlap problem to be meaningful – FLoRe is designed for long reads (PacBio, Oxford Nanopore, etc.), not short reads.

### Output Format
FLoRe outputs a JSON file listing all detected overlaps, where each overlap is a JSON object with detailed information. Each object (representing one overlap between a pair of reads) includes the following fields:

- **read1, read2:** Identifiers or indices of the two reads that overlap, as given in the input FASTA (e.g., read names or numbers). Each read’s orientation is also recorded (e.g., whether the read was used in forward or reverse-complement form for the overlap). In the JSON, this may be represented by separate fields like "orientation1" and "orientation2" (for example, "forward" or "reverse").

- **start1, end1, start2, end2:** The start and end coordinates (0-indexed or 1-indexed depending on implementation; typically 1-indexed in bioinformatics tools) of the overlapping region on each read. These define the exact span of the overlap on read1 and read2. For instance, "start1": 500, "end1": 1000 and "start2": 1, "end2": 501 would indicate the overlap covers positions 500–1000 of read1 and 1–501 of read2.

- **overlap_length:** The length of the overlap in base pairs (should be end1 - start1 + 1 and equal to end2 - start2 + 1). This gives a quick measure of how long the overlapping segment is.

- **overlap_region_read1, overlap_region_read2:** The actual nucleotide sequence of the overlapping segment from each read. These two strings should be identical (if the overlap alignment has no differences) or nearly identical if there are sequencing errors (note: FLoRe is detecting overlaps, not performing full error-corrected alignment, so these are just the raw sequences from each read over the overlapped interval – they may not perfectly match if there are mismatches). These fields help with downstream analysis or verification of the overlaps.

- **fingerprint_read1, fingerprint_read2:** The sequence of fingerprint IDs corresponding to the overlap region on each read. Essentially, this is the encoded k-mer representation of the overlapping segment. These can be useful for debugging or for analyzing the k-mer matching pattern. They will typically start and end with the same values if the overlap is aligned on k-mer boundaries.

- **used_algorithm:** Which algorithm or strategy confirmed this overlap. Possible values are "FGOE", "AOE", "PSH", "CORE", or "SuffixAutomaton" (for the fallback). This indicates the most permissive method that was needed – e.g., if FGOE found it, the others were not used for that pair; if all failed and suffix automaton found it, this field will read "SuffixAutomaton". This information is mostly for analysis; all overlaps listed are valid overlaps regardless of method.

An example snippet of the JSON output (for illustration) might look like:

```json
[
  {
    "read1": "readA", "orientation1": "forward",
    "start1": 100, "end1": 850, "len_read1": 1000,
    "read2": "readB", "orientation2": "reverse",
    "start2": 50, "end2": 800, "len_read2": 900,
    "overlap_length": 751,
    "overlap_region_read1": "ACCTG...TGA",
    "overlap_region_read2": "ACCTG...TGA",
    "fingerprint_read1": [12345, 67890],
    "fingerprint_read2": [12345, 67890],
    "used_algorithm": "CORE"
  }
]
```

(Note: The exact field names in the actual JSON may differ slightly, for example len_read1 and len_read2 might appear if the program includes full read lengths, etc. The above is for conceptual clarity. Use the actual output format as reference when parsing programmatically.)

This JSON format is both human-readable and easy to parse with scripts. It can be readily converted to other formats if needed. In fact, the provided FLORE wrapper will automatically create an Excel .xlsx spreadsheet from the JSON using Python, listing each overlap as a row (with nicely formatted columns). This is convenient for manual inspection or reporting. If you prefer CSV, you can similarly convert the JSON via tools like jq or Python/pandas.

Example: After running FLoRe on a test FASTA, you might get an output JSON indicating overlaps like read1 overlaps read5 by 500bp, read2 overlaps read7 by 1200bp, etc. You can then use those overlaps in a genome assembly pipeline (overlaps are the input to string graph assemblers), or simply analyze how many overlaps each read has (coverage) and so on.

## Performance

FLoRe was designed to be efficient in both runtime and memory, taking advantage of algorithmic pruning and parallelism. It has been benchmarked against a baseline overlap detection method (referred to as LROD in the thesis) and shows significant improvements. For instance, on a test of 100 long reads, FLoRe completed in 0.58 seconds, about 2.6× faster than the baseline (which took 1.52 s), and it correctly detected 8 true overlaps whereas the baseline found only 7. This indicates FLoRe not only runs faster but also achieved higher sensitivity in that case.

When scaling up to larger datasets, FLoRe maintains performance benefits. On 1,000 reads, enabling the default CFL compression in the fingerprint stage greatly reduced memory usage – peak RAM was around 915 MB (with no loss in recall) – and shaved off an additional 1.2 seconds of runtime compared to the uncompressed mode. Multi-threading further accelerates the overlap computation stage linearly with the number of cores (up to the point of I/O or memory bandwidth limits). The integrated profiling will report timings for each major phase so you can see how the workload scales.

In terms of complexity, the initial fingerprinting and indexing steps scale roughly linearly with the total input size (number of reads × average read length). The candidate filtering uses efficient intersection, which is typically sub-linear in the number of total comparisons due to skipping unmatched pairs. The overlap extension algorithms are designed to quickly accept or reject a candidate; in the best cases (clear overlaps or clear non-overlaps) one of the fast methods (FGOE/AOE) resolves the pair, and only difficult cases invoke CORE or suffix automaton. This layered approach keeps the average runtime per candidate low. In the worst-case scenario (very repetitive reads where many methods are tried), FLoRe might spend more time per pair, but those are exactly the cases where an overlap is genuinely hard to determine and requires thorough analysis.

Experiments in the thesis also compared FLoRe’s resource usage to other state-of-the-art tools. FLoRe’s memory footprint was observed to be lower than older overlap tools (like MHAP) and on par with or slightly higher than Minimap2, trading a bit of memory for speed and exhaustive overlap finding. The 915 MB for 1000 reads mentioned above shows it can handle moderate-size sets on a standard machine. For very large datasets (e.g., all reads of a mammalian genome), you may need a server-class machine, but FLoRe can be combined with read partitioning or downsampling strategies if needed.

In summary, FLoRe provides fast overlap detection by narrowing down candidates aggressively and leveraging efficient algorithms for validation. Its multi-threaded design and optimizations make it suitable for large-scale projects, and its accuracy remains high thanks to the comprehensive validation pipeline. Users can refer to the profiling output and logs to get detailed performance metrics for their specific run.
