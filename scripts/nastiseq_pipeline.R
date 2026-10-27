# =============================================================================
# NASTIseq full pipeline — featureCounts → NAT prediction → export
# Requirements: DuckTape conda environment (environment.yml): DESeq2, GenomicRanges, rtracklayer, etc.
#
# Usage:
#   Rscript nastiseq_pipeline.R \
#     --fc-dir  /path/to/featurecounts/ \
#     --gtf     /path/to/annotation.gtf \
#     --outdir  /path/to/output/
#
# Optional:
#   --cores    16       (default: all cores minus 1)
#   --fdr      0.05     (default: 0.05)
#   --overlap  0.40     (default: 0.40 minimum reciprocal overlap fraction)
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DESeq2)
  library(GenomicRanges)
  library(rtracklayer)
  library(parallel)
  library(pbmcapply)
  library(ROCR)
  library(dplyr)
  library(openxlsx)
})

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

option_list <- list(
  make_option(c("-i", "--fc-dir"),
              type    = "character",
              default = NULL,
              help    = "Directory containing featureCounts output files
                         (sorted.sense.txt / sorted.antisense.txt) [required]",
              metavar = "PATH"),

  make_option(c("-g", "--gtf"),
              type    = "character",
              default = NULL,
              help    = "Path to GTF or GFF annotation file [required]",
              metavar = "FILE"),

  make_option(c("-o", "--outdir"),
              type    = "character",
              default = NULL,
              help    = "Output directory (created if it does not exist) [required]",
              metavar = "PATH"),

  make_option(c("-c", "--cores"),
              type    = "integer",
              default = max(1L, parallel::detectCores() - 1L),
              help    = "Number of cores for parallel scoring [default: all-1]",
              metavar = "INT"),

  make_option(c("-f", "--fdr"),
              type    = "double",
              default = 0.05,
              help    = "FDR threshold for NAT calling [default: 0.05]",
              metavar = "FLOAT"),

  make_option(c("-v", "--overlap"),
              type    = "double",
              default = 0.40,
              help    = "Minimum reciprocal overlap fraction for NAT pairs [default: 0.40]",
              metavar = "FLOAT")
)

parser <- OptionParser(
  usage       = "%prog -i FC_DIR -g GTF -o OUTDIR [options]",
  option_list = option_list,
  description = "NASTIseq pipeline: featureCounts counts -> NAT prediction"
)

opt <- parse_args(parser)

missing_args <- c()
if (is.null(opt[["fc-dir"]])) missing_args <- c(missing_args, "--fc-dir")
if (is.null(opt[["gtf"]]))    missing_args <- c(missing_args, "--gtf")
if (is.null(opt[["outdir"]])) missing_args <- c(missing_args, "--outdir")
if (length(missing_args) > 0) {
  print_help(parser)
  stop(sprintf("Missing required argument(s): %s\n", paste(missing_args, collapse = ", ")),
       call. = FALSE)
}

FC_DIR           <- normalizePath(opt[["fc-dir"]], mustWork = TRUE)
GTF_FILE         <- normalizePath(opt[["gtf"]],    mustWork = TRUE)
OUT_DIR          <- opt[["outdir"]]
N_CORES          <- opt[["cores"]]
FDR_THR          <- opt[["fdr"]]
MIN_OVERLAP_FRAC <- opt[["overlap"]]

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
OUT_DIR    <- normalizePath(OUT_DIR)
OUT_PREFIX <- file.path(OUT_DIR,
                        sub("\\.(gtf|gff|gff3)(\\.gz)?$", "",
                            basename(GTF_FILE), ignore.case = TRUE))

SCRIPT_DIR <- tryCatch({
  args        <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args[grep("--file=", args)])
  if (length(script_path) == 1) dirname(normalizePath(script_path)) else getwd()
}, error = function(e) getwd())

cat("══════════════════════════════════════════════════════════════════════════\n")
cat("NASTIseq pipeline\n")
cat("══════════════════════════════════════════════════════════════════════════\n")
cat(sprintf("  featureCounts dir : %s\n", FC_DIR))
cat(sprintf("  annotation file   : %s\n", GTF_FILE))
cat(sprintf("  output directory  : %s\n", OUT_DIR))
cat(sprintf("  cores             : %d\n", N_CORES))
cat(sprintf("  FDR threshold     : %.3f\n", FDR_THR))
cat(sprintf("  min overlap frac  : %.2f\n", MIN_OVERLAP_FRAC))
cat("══════════════════════════════════════════════════════════════════════════\n\n")

suppressWarnings(source(file.path(SCRIPT_DIR, "core.R")))
suppressWarnings(source(file.path(SCRIPT_DIR, "helper.R")))

# =============================================================================
# STEP 1 — Read featureCounts output files
# =============================================================================

cat("── Step 1: Reading featureCounts files ──────────────────────────────────\n")

read_fc <- function(f) {
  dt <- fread(f, skip = 1, showProgress = FALSE)
  dt <- dt[, .(Geneid, count = .SD[[ncol(dt)]])]
  dt[, Geneid := sub("^(gene|transcript|exon):", "", Geneid)]
  dt[, sample := sub("\\..*", "", basename(f))]
  dt
}

sense_files <- list.files(FC_DIR, pattern = "\\.sense\\.txt$",     full.names = TRUE)
antis_files <- list.files(FC_DIR, pattern = "\\.antisense\\.txt$", full.names = TRUE)

if (length(sense_files) == 0) stop("No sense files found in: ", FC_DIR)
if (length(antis_files) == 0) stop("No antisense files found in: ", FC_DIR)

cat(sprintf("Found %d sense files and %d antisense files\n",
            length(sense_files), length(antis_files)))

sense <- rbindlist(lapply(sense_files, read_fc))
antis <- rbindlist(lapply(antis_files, read_fc))

smat_dt  <- dcast(sense, Geneid ~ sample, value.var = "count")
asmat_dt <- dcast(antis, Geneid ~ sample, value.var = "count")
rm(sense, antis)

genes_shared <- intersect(smat_dt$Geneid, asmat_dt$Geneid)
smat  <- as.matrix(as.data.frame(smat_dt[match(genes_shared,  smat_dt$Geneid),  -1],
                                  row.names = genes_shared))
asmat <- as.matrix(as.data.frame(asmat_dt[match(genes_shared, asmat_dt$Geneid), -1],
                                  row.names = genes_shared))
rm(smat_dt, asmat_dt)

cat(sprintf("Count matrices: %d genes x %d samples\n", nrow(smat), ncol(smat)))

# =============================================================================
# STEP 2 — Import annotation and build NAT pairs
# =============================================================================
# Gene-ID-centric annotation parsing: works with any GTF or GFF3 regardless
# transcripts, giving the most complete gene set (handles GTF which lacks gene
# features, and GFF3 where some genes only have transcript entries).
# Gene IDs are normalized here (strip gene:/transcript: prefix) to match the
# featureCounts output. The same normalized IDs are reused in step 3 for
# genepos, guaranteeing nat_pairs and genepos are always consistent.
# =============================================================================

cat("\n── Step 2: Importing annotation and building NAT pairs ──────────────────\n")

# Import annotation once — reused across steps 2 and 3
gtf_gr     <- suppressWarnings(suppressMessages(import(GTF_FILE)))
mcol_names <- names(mcols(gtf_gr))
all_types  <- unique(as.character(gtf_gr$type))

cat(sprintf("Annotation loaded: %d features across %d feature types\n",
            length(gtf_gr), length(all_types)))

# ── Gene-ID-centric approach ─────────────────────────────────────────────────
# Instead of filtering by feature type (fragile: 'gene', 'ncRNA_gene',
# 'pseudogene', 'lincRNA_gene' differ across sources), we use gene_id as the
# universal key. In GTF every feature carries gene_id; in GFF3 (Ensembl) all
# gene-level features carry gene_id regardless of their type. Grouping ALL
# features by gene_id and taking min(start)/max(end) gives correct gene spans
# without any dependency on feature type names.
if ("gene_id" %in% mcol_names) {
  raw_ids <- as.character(mcols(gtf_gr)$gene_id)
} else if ("Name" %in% mcol_names) {
  raw_ids <- as.character(mcols(gtf_gr)$Name)
} else {
  raw_ids <- as.character(mcols(gtf_gr)$ID)
}

# Normalize: strip type: prefix (gene:AT1G01010 → AT1G01010) to match
# featureCounts output (step 1 applies the same stripping)
norm_ids <- sub("^(gene|transcript|exon):", "", raw_ids)

# Keep features that have a valid gene_id
keep   <- !is.na(norm_ids) & nchar(norm_ids) > 0
gr_sub <- gtf_gr[keep]
ids    <- norm_ids[keep]

# Group by gene_id: min(start) / max(end) across all features of that gene
gene_dt <- data.table(
  gene_id = ids,
  seqname = as.character(seqnames(gr_sub)),
  start   = start(gr_sub),
  end     = end(gr_sub),
  strand  = as.character(strand(gr_sub))
)

gene_ranges <- gene_dt[
  , .(seqname = seqname[1],
      start   = min(start),
      end     = max(end),
      strand  = strand[1]),
  by = .(gene_id)
]

gene_id_col <- gene_ranges$gene_id
gene_gr     <- GRanges(
  seqnames = gene_ranges$seqname,
  ranges   = IRanges(start = gene_ranges$start, end = gene_ranges$end),
  strand   = gene_ranges$strand
)
names(gene_gr) <- gene_id_col

cat(sprintf("Annotation: %d unique gene loci built from %d features\n",
            length(gene_gr), sum(keep)))

# ── Biotype (best-effort) ─────────────────────────────────────────────────────
# Look for gene-level feature types to get biotype; falls back to 'unknown'.
# This is purely informational and does NOT gate which genes are included.
gene_type_hits <- all_types[grepl("^gene$|_gene$|^pseudogene$", all_types)]
gene_features  <- gtf_gr[gtf_gr$type %in% gene_type_hits]
feat_mcols     <- names(mcols(gene_features))

if (length(gene_features) > 0 &&
    any(c("gene_biotype", "biotype") %in% feat_mcols)) {
  feat_ids <- sub("^(gene|transcript|exon):", "", as.character(
    if ("gene_id" %in% feat_mcols) gene_features$gene_id
    else if ("Name" %in% feat_mcols) gene_features$Name
    else gene_features$ID
  ))
  raw_bio     <- as.character(
    if ("gene_biotype" %in% feat_mcols) gene_features$gene_biotype
    else gene_features$biotype
  )
  bio_map     <- setNames(raw_bio, feat_ids)
  biotype_col <- unname(bio_map[gene_id_col])
  biotype_col[is.na(biotype_col)] <- "unknown"
  cat(sprintf("Biotype info found for %d / %d genes\n",
              sum(biotype_col != "unknown"), length(gene_id_col)))
} else {
  biotype_col <- rep("unknown", length(gene_gr))
  cat("No biotype info found — all genes set to 'unknown'\n")
}

# Find opposite-strand overlapping gene pairs (NATs)
ov  <- findOverlaps(gene_gr, ignore.strand = TRUE, minoverlap = 1)
ov  <- ov[queryHits(ov) != subjectHits(ov)]
opp <- strand(gene_gr)[queryHits(ov)] != strand(gene_gr)[subjectHits(ov)]
ov  <- ov[opp]

q          <- gene_gr[queryHits(ov)]
s          <- gene_gr[subjectHits(ov)]
ovlp_width <- width(pintersect(q, s, ignore.strand = TRUE))
frac_q     <- ovlp_width / width(q)
frac_s     <- ovlp_width / width(s)
ov         <- ov[frac_q >= MIN_OVERLAP_FRAC & frac_s >= MIN_OVERLAP_FRAC]

pairs <- data.frame(
  g1 = gene_id_col[queryHits(ov)],
  g2 = gene_id_col[subjectHits(ov)],
  stringsAsFactors = FALSE
)

nat_pairs <- unique(data.frame(
  geneA = pmin(pairs$g1, pairs$g2),
  geneB = pmax(pairs$g1, pairs$g2)
))

cat(sprintf("Found %d NAT pairs (>= %.0f%% reciprocal overlap)\n",
            nrow(nat_pairs), MIN_OVERLAP_FRAC * 100))

# =============================================================================
# STEP 3 — Build genepos annotation table (uses same gene_gr and IDs as step 2)
# =============================================================================

cat("\n── Step 3: Building gene position table ─────────────────────────────────\n")

# Use 'annot_source' internally to avoid conflict with base R source() function,
# then rename to 'source' after aggregation so column order matches GFF format
# (seqname, source, feature, start, end, score, strand, frame, attributes)
# which is what FindNATs expects.
genepos <- data.table(
  seqname      = as.character(seqnames(gene_gr)),
  annot_source = biotype_col,
  feature      = "gene",
  start        = start(gene_gr),
  end          = end(gene_gr),
  score        = ".",
  strand       = as.character(strand(gene_gr)),
  frame        = ".",
  attributes   = gene_id_col
)

# For genes with multiple entries, keep the largest span per gene_id.
# Group by .(attributes) WITHOUT repeating it in .() to avoid duplicate columns;
# reorder so attributes is the 9th column as FindNATs expects.
genepos <- genepos[
  , .(seqname      = seqname[1],
      annot_source = annot_source[1],
      feature      = "gene",
      start        = min(start),
      end          = max(end),
      score        = ".",
      strand       = strand[1],
      frame        = "."),
  by = .(attributes)
]

setnames(genepos, "annot_source", "source")
setcolorder(genepos, c("seqname", "source", "feature",
                        "start",  "end",    "score",
                        "strand", "frame",  "attributes"))

genepos           <- as.data.frame(genepos)
rownames(genepos) <- genepos$attributes

cat(sprintf("Gene position table: %d genes\n", nrow(genepos)))

# =============================================================================
# STEP 4 — Align matrices, pairs and annotation to shared genes
# =============================================================================

cat("\n── Step 4: Aligning genes across matrices, pairs and annotation ─────────\n")

# Keep only genepos rows whose attributes are in smat, drop NAs and duplicates
# (genes in smat missing from the annotation get silently skipped — this is
#  expected and does not affect scoring)
genepos <- genepos[genepos$attributes %in% rownames(smat), ]
genepos <- genepos[!is.na(genepos$attributes) & !duplicated(genepos$attributes), ]
rownames(genepos) <- genepos$attributes

# valid_genes = genes present in BOTH smat and genepos; used to filter all pairs
valid_genes <- intersect(rownames(smat), genepos$attributes)

nat_pairs_filtered <- nat_pairs[
  nat_pairs[, 1] %in% valid_genes &
  nat_pairs[, 2] %in% valid_genes, ]

cat(sprintf("nat_pairs after filtering: %d → %d\n",
            nrow(nat_pairs), nrow(nat_pairs_filtered)))
cat(sprintf("Genes in smat: %d  |  Genes in annotation: %d  |  Intersection: %d\n",
            nrow(smat), nrow(genepos), length(valid_genes)))

colnames(smat)  <- NULL
colnames(asmat) <- NULL

# =============================================================================
# STEP 5 — Save intermediate RDS files
# =============================================================================

cat("\n── Step 5: Saving intermediate RDS files ────────────────────────────────\n")

saveRDS(smat,               file.path(OUT_DIR, "smat.rds"))
saveRDS(asmat,              file.path(OUT_DIR, "asmat.rds"))
saveRDS(nat_pairs_filtered, file.path(OUT_DIR, "nat_pairs.rds"))
saveRDS(genepos,            file.path(OUT_DIR, "genepos.rds"))
cat(sprintf("Saved intermediate RDS files to %s/\n", OUT_DIR))

# =============================================================================
# STEP 6 — Define fast parallel getNASTIscore (DESeq2-based)
# =============================================================================

getNASTIscore_fast <- function(smat, asmat,
                               llanti  = llanti_full,
                               llsense = llsense_full,
                               n_cores = max(1L, detectCores() - 1L)) {
  labels <- rep("wr", ncol(smat))
  allmat <- smat + asmat

  colnames(allmat) <- paste0("s", seq_len(ncol(allmat)))
  col_data <- data.frame(condition = factor(labels), row.names = colnames(allmat))

  dds <- suppressMessages(
    DESeqDataSetFromMatrix(countData = allmat, colData = col_data, design = ~ 1)
  )
  dds        <- suppressMessages(estimateSizeFactors(dds))
  dds        <- suppressMessages(estimateDispersions(dds, quiet = TRUE))
  SF         <- sizeFactors(dds)
  rawVarFunc <- dispersionFunction(dds)

  p      <- colSums(asmat) / (colSums(smat) + colSums(asmat))
  ngenes <- nrow(smat)
  nconds <- length(labels)

  low_count_genes <- rowSums((smat + asmat) < 5) > 1
  cat(sprintf("Pre-filtered %d / %d low-count genes\n", sum(low_count_genes), ngenes))

  score_one_gene <- function(i) {
    if (low_count_genes[i]) return(-Inf)
    N1     <- smat[i, ];  N2 <- asmat[i, ]
    lowb   <- rep(0, nconds * 2)
    highb  <- c(N2, N1)
    startx <- pmin(c(floor(N1 * p), floor(N2 * p)), highb)
    sol <- tryCatch(
      optim(par     = startx,
            fn      = function(x) {
              x_int <- pmax(lowb, pmin(highb, round(x)))
              val   <- llanti(x_int, c(N1, N2), p, SF, rawVarFunc)
              if (is.na(val) || is.infinite(val)) return(1e20)
              -val
            },
            method  = "Nelder-Mead",
            control = list(maxit = 500, reltol = 1e-6)),
      error = function(e) NULL)
    if (is.null(sol)) return(-Inf)
    sol_clipped <- pmax(lowb, pmin(highb, round(sol$par)))
    BICanti(sol_clipped, c(N1, N2), p, SF, rawVarFunc, llanti, llsense)
  }

  scores           <- pbmclapply(seq_len(ngenes), score_one_gene, mc.cores = n_cores)
  BICcormat        <- unlist(scores)
  names(BICcormat) <- rownames(smat)
  sort(BICcormat, decreasing = TRUE)
}

defineFDR_robust <- function(pred, xthreshold = 0.05) {
  totalp  <- pred@tp[[1]] + pred@fp[[1]]
  fdr     <- pred@fp[[1]] / totalp
  cutoffs <- pred@cutoffs[[1]]
  valid   <- is.finite(cutoffs) & is.finite(fdr)
  fdr     <- fdr[valid];  cutoffs <- cutoffs[valid]
  ord     <- order(fdr);  fdr <- fdr[ord];  cutoffs <- cutoffs[ord]
  if (xthreshold < min(fdr)) return(max(cutoffs))
  if (xthreshold > max(fdr)) return(min(cutoffs))
  suppressWarnings(approx(x = fdr, y = cutoffs, xout = xthreshold)$y)
}

# =============================================================================
# STEP 7 — Run NASTIseq
# =============================================================================

cat("\n── Step 7: Running NASTIseq ──────────────────────────────────────────────\n")

WRscore  <- getNASTIscore_fast(smat, asmat, n_cores = N_CORES)
negpairs <- getnegativepairs(genepos)
negpairs    <- negpairs[
  negpairs[, 1] %in% valid_genes &
  negpairs[, 2] %in% valid_genes, ]

bicout  <- EvalBIC(WRscore, smat, asmat, nat_pairs_filtered, negpairs)
WRpred  <- list(predictions = bicout[, 1], labels = bicout[, 2])

pred        <- sign(WRpred$predictions) * log1p(abs(WRpred$predictions))
WRpred_rocr <- prediction(pred, WRpred$labels)
thr         <- defineFDR_robust(WRpred_rocr, FDR_THR)

cat(sprintf("FDR threshold at %.0f%%: %.4f\n", FDR_THR * 100, thr))

WR_names <- FindNATs(WRscore, thr, nat_pairs_filtered, genepos)

cat("\n── Results ──────────────────────────────────────────────────────────────\n")
cat("Known pairs: ", nrow(WR_names$knownpairs), "\n")
cat("New pairs:   ", nrow(WR_names$newpairs),   "\n")
cat("Orphans:     ", length(WR_names$neworphan), "\n")

# =============================================================================
# STEP 8 — Save results
# =============================================================================

cat("\n── Step 8: Saving results ───────────────────────────────────────────────\n")

saveRDS(
  list(WRscore = WRscore, WRpred = WRpred, WR_names = WR_names, thr = thr),
  paste0(OUT_PREFIX, "_results.rds")
)
cat(sprintf("Saved %s_results.rds\n", OUT_PREFIX))

# =============================================================================
# STEP 9 — Export orphans to Excel and BED
# =============================================================================

cat("\n── Step 9: Exporting orphans ────────────────────────────────────────────\n")

orphans <- genepos %>%
  filter(attributes %in% WR_names$neworphan) %>%
  select(seqname, start, end, attributes, score, strand)

write.xlsx(orphans, paste0(OUT_PREFIX, "_orphans.xlsx"), rowNames = FALSE)

write.table(
  orphans,
  paste0(OUT_PREFIX, "_orphans.bed"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

cat(sprintf("Exported %d orphans to:\n  %s_orphans.xlsx\n  %s_orphans.bed\n",
            nrow(orphans), OUT_PREFIX, OUT_PREFIX))
            
# =============================================================================
# STEP 10 — Export new pairs to Excel and BED
# =============================================================================

cat("\n── Step 10: Exporting knownpairs ────────────────────────────────────────────\n")

knownpairs <- genepos %>%
  filter(attributes %in% WR_names$knownpairs) %>%
  select(seqname, start, end, attributes, score, strand)

write.xlsx(knownpairs, paste0(OUT_PREFIX, "_knownpairs.xlsx"), rowNames = FALSE)

write.table(
  knownpairs,
  paste0(OUT_PREFIX, "_knownpairs.bed"),
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
)

cat(sprintf("Exported %d knownpairs to:\n  %s_knownpairs.xlsx\n  %s_knownpairs.bed\n",
            nrow(knownpairs), OUT_PREFIX, OUT_PREFIX))


cat("\n══════════════════════════════════════════════════════════════════════════\n")
cat("Pipeline complete\n")
cat("══════════════════════════════════════════════════════════════════════════\n")
