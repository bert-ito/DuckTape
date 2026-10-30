# DuckTape

DuckTape discovers **natural antisense transcripts (NATs)** in strand-specific RNA-seq data
and produces a high-quality assembly of them.

* **Stage 1** (`ducktape_stage1.sh`) prepares the annotation and counts sense and antisense
  reads per gene. It then runs [NASTIseq](#references) to find *orphan* genes. These are genes
  with significant antisense transcription but no annotated antisense partner. It also builds
  a benchmark set from known sense–antisense pairs.
* **Stage 2** (`ducktape_stage2.sh`) removes antisense signal that comes from library
  preparation, calls antisense peaks around the orphan genes, and assembles them into NAT
  transcript models. In **test** mode, it hides the antisense partner of known pairs and
  measures how well DuckTape recovers it.

The method is described in detail in the DuckTape paper (_TODO: add reference_). If you use
DuckTape, please cite the paper and the tools it is built on (see [References](#references)).

## Requirements

* Linux with `bash` ≥ 4 and GNU coreutils. macOS is not supported out of the box.
* [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/). We recommend
  [Miniforge](https://github.com/conda-forge/miniforge).

## Installation

```bash
git clone https://github.com/bert-ito/DuckTape.git
cd DuckTape

conda env create -f environment.yml      # or: mamba env create -f environment.yml
conda activate ducktape

bash ducktape_stage1.sh --help
```

If the environment fails to solve, set strict channel priority as Bioconda recommends
(`conda config --set channel_priority strict`) and try again.

Both scripts use the active environment. Alternatively, pass `--env ducktape` (a name or a
path) and the script activates it for you.

## Input data

| Input | Description |
|---|---|
| BAM files | One `*.bam` file per sample, all in one directory. Use strand-specific RNA-seq aligned with a splice-aware aligner, coordinate-sorted (indexing is recommended). StringTie expects the `XS` strand tag on spliced reads, e.g. HISAT2 `--dta` or STAR `--outSAMstrandField intronMotif`. |
| Annotation | GFF3 or GTF, optionally gzipped, for the same genome assembly as the BAM files. `-g` can be given several times to merge annotations. |
| Library type | `RF` (dUTP / fr-firststrand, e.g. Illumina TruSeq Stranded; the default) or `FR` (fr-secondstrand). Unstranded data are not supported. |
| Read type | `PE` (paired-end; the default) or `SE` (single-end). |

The sample name is the BAM file name without `.bam`. Give Stage 2 the same BAM files as
Stage 1, because it reuses the Stage 1 read counts.

> **Note:** NASTIseq cuts sample names at the first dot. Give every BAM a name that is unique
> before its first dot: `S1.rep1.bam` and `S1.rep2.bam` would be merged as sample `S1`.

## Quick start

```bash
conda activate ducktape

# Stage 1
bash ducktape_stage1.sh \
    -g genome_annotation.gff3 \
    -b bams/ \
    -o stage1_results \
    --threads 16 --mem-gb 64

# Stage 2, analysis mode: final NAT models in stage2_results/FINAL_NATS.gtf
bash ducktape_stage2.sh \
    -s1 stage1_results \
    -b bams/ \
    -o stage2_results \
    --mode analysis \
    --threads 16 --mem-gb 64

# Stage 2, test mode (optional): benchmark on known sense-antisense pairs
bash ducktape_stage2.sh \
    -s1 stage1_results \
    -b bams/ \
    -o stage2_test \
    --mode test \
    --threads 16 --mem-gb 64
```

> **Note:** each run deletes and recreates `<output>/tmp_work`. Use a separate output directory
> for every run you want to keep.

## Stage 1

### Options

| Option | Default | Description |
|---|---|---|
| `-g, --annotation FILE` | *required* | GTF/GFF3 annotation (optionally `.gz`). Repeat to merge several annotations. |
| `-b, --bam-dir DIR` | *required* | Directory with `*.bam` files. |
| `-o, --output DIR` | *required* | Output directory. |
| `--library-type RF\|FR` | `RF` | Library strandedness. |
| `--read-type PE\|SE` | `PE` | Paired-end or single-end reads. |
| `--test-expression FLOAT` | `20` | Minimum total sense count (summed over samples) for both genes of a known pair. |
| `--query-fraction FLOAT` | `0.05` | Maximum query/test expression ratio for a known pair. |
| `--threads INT` | all available CPUs | Number of CPUs. |
| `--mem-gb INT` | available memory | Memory budget in GB. |
| `--env NAME\|PATH` | current environment | Conda environment to activate. |
| `--pipeline-dir DIR` | `scripts/` | Location of the helper scripts. |
| `--no-benchmark` | | Do not record [resource usage](#resource-usage). |
| `--benchmark-interval SEC` | `1` | Memory sampling interval for the resource usage tables. |

### Outputs

| File | Description |
|---|---|
| `annotation.gff3` | Normalised (and merged) annotation in GFF3. |
| `reference.gtf` | The same annotation in GTF. |
| `featurecounts/<sample>.sense.txt` | featureCounts sense-strand counts. |
| `featurecounts/<sample>.antisense.txt` | featureCounts antisense-strand counts. |
| `orphans.bed` | Orphan genes reported by NASTIseq. |
| `TEST_KNOWNPAIRS_expression_filtered.bed` | Test genes of the expression-filtered known pairs. |
| `TEST_KNOWNPAIRS_QUERY_expression_filtered.bed` | Query genes of the same pairs, in the same order. |
| `reference_test.gff3` | `annotation.gff3` without the query genes, used by Stage 2 test mode. |
| `benchmark/` | [Resource usage](#resource-usage) per step and per command. |
| `tmp_work/` | Intermediate files and logs, including the full NASTIseq results. |

## Stage 2

### Options

| Option | Default | Description |
|---|---|---|
| `-s1, --stage1-dir DIR` | *required* | Output directory of Stage 1. |
| `-b, --bam-dir DIR` | *required* | Directory with the same `*.bam` files used in Stage 1. |
| `-o, --output DIR` | *required* | Output directory. |
| `--mode analysis\|test` | `analysis` | Discover NATs, or benchmark on known pairs. |
| `--library-type RF\|FR` | `RF` | Library strandedness. |
| `--read-type PE\|SE` | `PE` | Paired-end or single-end reads. |
| `--target-class all\|coding\|noncoding` | `all` | Which target genes to use. |
| `--orphan-buffer INT` | `10000` | Extension (bp) around each target before the strand is flipped. |
| `--alpha-factor FLOAT` | `4` | Multiplier applied to the antisense leakage (α) in the coverage correction. |
| `--intersect-frac FLOAT` | `0.2` | Minimum fraction of a read that must overlap a peak. |
| `--threads INT` | all available CPUs | Number of CPUs. |
| `--mem-gb INT` | available memory | Memory budget in GB. |
| `--env NAME\|PATH` | current environment | Conda environment to activate. |
| `--pipeline-dir DIR` | `scripts/` | Location of the helper scripts. |
| `--no-benchmark` | | Do not record [resource usage](#resource-usage). |
| `--benchmark-interval SEC` | `1` | Memory sampling interval for the resource usage tables. |

### Outputs

| File | Mode | Description |
|---|---|---|
| `FINAL_NATS.gtf` | analysis | Final NAT transcript models. |
| `orphans.bed` | analysis | Orphan genes used as targets. |
| `FINAL_QUERY.gtf` | test | NAT models assembled for the test genes. |
| `test_knownpairs.bed` | test | Test genes used as targets. |
| `final_queries.stats` | test | gffcompare accuracy of `FINAL_QUERY.gtf` against the removed query genes. |
| `NAT_tracking_x.tsv` | both | gffcompare tracking lines (class code `x`) linking each NAT to its sense gene. |
| `benchmark/` | both | [Resource usage](#resource-usage) per step and per command. |
| `tmp_work/` | both | Intermediate files and tool logs. This can be large and may be deleted once the run is checked. |

## Resource usage

Both stages record the elapsed time, CPU time, memory and I/O of every step, and of every
tool call within it: featureCounts, NASTIseq, bamCoverage, samtools, bedtools, StringTie,
gffcompare, AGAT, gffread and the helper scripts. The figures include every process a tool
starts, such as the forked R workers inside NASTIseq. A summary per step and per tool is
printed at the end of each stage, and the full tables are written to `<output>/benchmark/`:

| File | Contents |
|---|---|
| `steps.tsv` | One row per step, plus a `TOTAL` row covering the whole stage from start to finish. |
| `commands.tsv` | One row per tool call, with its step, tool name, exit status and full command line. |

| Column | Meaning |
|---|---|
| `wall_s` | Elapsed time, in seconds. |
| `cpu_user_s`, `cpu_sys_s`, `cpu_s` | CPU time in user mode, in kernel mode, and in total, in seconds. |
| `avg_cores` | `cpu_s / wall_s`: the average number of CPUs kept busy. |
| `peak_mem_mb` | Peak memory of all processes together, as the sum of their PSS. PSS shares out pages used by several processes, so this is the figure to base memory requests on. |
| `peak_rss_mb` | The same peak, summing RSS instead. Shared pages are counted once per process, so this is an upper bound. |
| `max_rss_mb` | Peak RSS of the largest single process. |
| `read_mb`, `write_mb` | Data read and written, including pipes and cached files. |
| `disk_read_mb`, `disk_write_mb` | Data read from and written to storage. May be 0 on network file systems. |
| `max_procs`, `max_threads` | Most processes and threads running at the same time. |

`commands.tsv` also reports major page faults and voluntary and involuntary context switches.
Memory is in MB (1 MB = 2<sup>20</sup> bytes).

**How it is measured.** CPU time, `max_rss_mb`, page faults and context switches come from the
Linux kernel's accounting for each command and all of its child processes, so they are exact.
I/O comes from `/proc/<pid>/io`. The memory totals are sampled from `/proc` every
`--benchmark-interval` seconds, so a peak shorter than the interval can be missed. For commands
shorter than the interval, `peak_mem_mb` falls back to `max_rss_mb`. Commands too small to list
individually (`awk`, `cut`, `grep`) are still included in the step and stage totals. Monitoring
requires Linux.

**Notes.**

* `max_rss_mb` is never below about 13 MB, the size of the monitoring process. Linux carries a
  process's memory peak over into the programs it starts, and GNU `time` shows the same effect.
* In our tests, sampling once per second added about 5% to the run time of NASTIseq, the
  step with the most processes to sample. For long runs, `--benchmark-interval 5` reduces this
  further. `--no-benchmark` turns monitoring off; the pipeline's results are identical either
  way.
* StringTie numbers its genes (`STRG.*`) in the order its threads finish, so these IDs can
  differ between runs when it uses more than one thread. The transcript models themselves are
  the same.

## References

Please cite the following software when you use DuckTape.

### Methods and tools

| Software | Used for | Reference |
|---|---|---|
| NASTIseq | Antisense scoring, orphan and known-pair detection | Li S, Liberman LM, Mukherjee N, Benfey PN, Ohler U. Integrated detection of natural antisense transcripts using strand-specific RNA sequencing data. *Genome Res.* 2013;23(10):1730–1739. [doi:10.1101/gr.149310.112](https://doi.org/10.1101/gr.149310.112) |
| AGAT | Annotation normalisation and merging | Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format. Zenodo. [doi:10.5281/zenodo.3552717](https://doi.org/10.5281/zenodo.3552717) (cite the version you used) |
| GffRead, GffCompare | Annotation conversion/filtering; classification of assembled transcripts | Pertea G, Pertea M. GFF Utilities: GffRead and GffCompare. *F1000Research.* 2020;9:304. [doi:10.12688/f1000research.23297.2](https://doi.org/10.12688/f1000research.23297.2) |
| featureCounts (Subread) | Sense and antisense read counting | Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics.* 2014;30(7):923–930. [doi:10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656) |
| SAMtools | BAM filtering, sorting, merging and indexing | Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. *GigaScience.* 2021;10(2):giab008. [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008) |
| BEDTools | Peak merging and intersection | Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics.* 2010;26(6):841–842. [doi:10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033) |
| deepTools (`bamCoverage`) | Strand-specific coverage tracks | Ramírez F, Ryan DP, Grüning B, Bhardwaj V, Kilpert F, Richter AS, Heyne S, Dündar F, Manke T. deepTools2: a next generation web server for deep-sequencing data analysis. *Nucleic Acids Res.* 2016;44(W1):W160–W165. [doi:10.1093/nar/gkw257](https://doi.org/10.1093/nar/gkw257) |
| UCSC `bigWigToBedGraph` | bigWig to bedGraph conversion | Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. BigWig and BigBed: enabling browsing of large distributed datasets. *Bioinformatics.* 2010;26(17):2204–2207. [doi:10.1093/bioinformatics/btq351](https://doi.org/10.1093/bioinformatics/btq351) |
| StringTie | Transcript assembly | Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT, Salzberg SL. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. *Nat Biotechnol.* 2015;33(3):290–295. [doi:10.1038/nbt.3122](https://doi.org/10.1038/nbt.3122) <br> Kovaka S, Zimin AV, Pertea GM, Razaghi R, Salzberg SL, Pertea M. Transcriptome assembly from long-read RNA-seq alignments with StringTie2. *Genome Biol.* 2019;20:278. [doi:10.1186/s13059-019-1910-1](https://doi.org/10.1186/s13059-019-1910-1) |

### R and Bioconductor packages

| Software | Used for | Reference |
|---|---|---|
| R | NASTIseq implementation | R Core Team. R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/> |
| DESeq2 | Size factors and dispersion estimates | Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol.* 2014;15(12):550. [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8) |
| GenomicRanges, IRanges | Gene ranges and overlaps | Lawrence M, Huber W, Pagès H, Aboyoun P, Carlson M, Gentleman R, Morgan MT, Carey VJ. Software for computing and annotating genomic ranges. *PLoS Comput Biol.* 2013;9(8):e1003118. [doi:10.1371/journal.pcbi.1003118](https://doi.org/10.1371/journal.pcbi.1003118) |
| rtracklayer | Annotation import | Lawrence M, Gentleman R, Carey V. rtracklayer: an R package for interfacing with genome browsers. *Bioinformatics.* 2009;25(14):1841–1842. [doi:10.1093/bioinformatics/btp328](https://doi.org/10.1093/bioinformatics/btp328) |
| S4Vectors | Bioconductor containers (overlap hits, metadata columns) | Pagès H, Lawrence M, Aboyoun P. S4Vectors: Foundation of vector-like and list-like containers in Bioconductor. R package. [doi:10.18129/B9.bioc.S4Vectors](https://doi.org/10.18129/B9.bioc.S4Vectors) |
| ROCR | FDR threshold for NASTIseq scores | Sing T, Sander O, Beerenwinkel N, Lengauer T. ROCR: visualizing classifier performance in R. *Bioinformatics.* 2005;21(20):3940–3941. [doi:10.1093/bioinformatics/bti623](https://doi.org/10.1093/bioinformatics/bti623) |
| data.table | Count table handling | Barrett T, Dowle M, Srinivasan A, Gorecki J, Chirico M, Hocking T, Schwendinger B, Krylov I. data.table: Extension of 'data.frame'. R package. <https://CRAN.R-project.org/package=data.table> |
| dplyr | Result filtering | Wickham H, François R, Henry L, Müller K, Vaughan D. dplyr: A Grammar of Data Manipulation. R package. <https://CRAN.R-project.org/package=dplyr> |
| optparse | Command-line parsing | Davis TL. optparse: Command Line Option Parser. R package. <https://CRAN.R-project.org/package=optparse> |
| pbmcapply | Parallel scoring with progress bar | Kuang K, Kong Q, Napolitano F. pbmcapply: Tracking the Progress of Mc*pply with Progress Bar. R package. <https://CRAN.R-project.org/package=pbmcapply> |
| openxlsx | Excel export of orphans and known pairs | Schauberger P, Walker A. openxlsx: Read, Write and Edit xlsx Files. R package. <https://CRAN.R-project.org/package=openxlsx> |

### Python packages

| Software | Used for | Reference |
|---|---|---|
| Python | Helper scripts | Python Software Foundation. Python Language Reference. <https://www.python.org/> |
| pandas | Table handling | McKinney W. Data structures for statistical computing in Python. In: *Proc. 9th Python in Science Conf.* 2010:56–61. [doi:10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a) <br> The pandas development team. pandas-dev/pandas: Pandas. Zenodo. [doi:10.5281/zenodo.3509134](https://doi.org/10.5281/zenodo.3509134) |
| NumPy | Numerical arrays | Harris CR, Millman KJ, van der Walt SJ, *et al.* Array programming with NumPy. *Nature.* 2020;585:357–362. [doi:10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2) |
| SciPy | Student's t distribution for the correlation p-value | Virtanen P, Gommers R, Oliphant TE, *et al.* SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat Methods.* 2020;17:261–272. [doi:10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2) |
| pyBigWig | Reading and writing bigWig files | Ryan DP, *et al.* pyBigWig: a Python extension for quick access to bigWig and bigBed files. Zenodo. [doi:10.5281/zenodo.45238](https://doi.org/10.5281/zenodo.45238) |

### Software distribution

| Software | Reference |
|---|---|
| Bioconda | Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; The Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. *Nat Methods.* 2018;15(7):475–476. [doi:10.1038/s41592-018-0046-7](https://doi.org/10.1038/s41592-018-0046-7) |
| conda-forge | conda-forge community. The conda-forge Project: Community-based Software Distribution Built on the conda Package Format and Ecosystem. Zenodo. 2015. [doi:10.5281/zenodo.4774216](https://doi.org/10.5281/zenodo.4774216) |

To get the exact versions installed in your environment for a methods section, run
`conda list -n ducktape`.

## License

_TODO: choose a license._ `scripts/core.R` and `scripts/helper.R` are adapted from NASTIseq.
Check that the NASTIseq license allows redistribution and is compatible with the license you
choose.
