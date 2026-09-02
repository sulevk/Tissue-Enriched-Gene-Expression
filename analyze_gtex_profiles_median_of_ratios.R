resolve_project_dir <- function() {
  candidates <- c(
    Sys.getenv("WORKSPACE_DIR"),
    Sys.getenv("PROJECT_ROOT"),
    Sys.getenv("PWD"),
    getwd(),
    "."
  )

  for (candidate in candidates) {
    if (nzchar(candidate) && dir.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = FALSE))
    }
  }

  for (root in c("/scratch", "/mnt", "/workspace", "/tmp")) {
    if (dir.exists(root)) {
      return(normalizePath(root, winslash = "/", mustWork = FALSE))
    }
  }

  normalizePath(".", winslash = "/", mustWork = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(limma)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

project_dir <- resolve_project_dir()
rdata_file <- file.path(project_dir, "GTExWithTissues.RData")
onehot_file <- file.path(project_dir, "matrix1_onehot_tissue.txt")
out_dir <- file.path(project_dir, "gtex_profile_analysis_median_of_ratios")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

checkpoint <- function(obj, name, dir_path = out_dir) {
  saveRDS(obj, file.path(dir_path, paste0(name, ".rds")))
  invisible(obj)
}

if (!file.exists(rdata_file)) {
  stop("Missing GTEx input file: ", rdata_file)
}

load(rdata_file)
if (!exists("gct_obj2", inherits = FALSE)) {
  stop("gct_obj2 not found in RData")
}

counts <- as.matrix(gct_obj2@mat)
gene_ids <- as.character(gct_obj2@rid)
sample_ids <- as.character(gct_obj2@cid)

if (length(sample_ids) == 0L || nrow(counts) == 0L) {
  stop("Expression matrix is empty.")
}

onehot <- fread(onehot_file, sep = "\t", header = TRUE, check.names = FALSE)
if (!"SAMPID" %in% names(onehot)) {
  stop("Metadata must contain an SAMPID column")
}

idx <- match(sample_ids, onehot$SAMPID)
if (anyNA(idx)) {
  stop("Some expression-matrix sample IDs are missing from metadata")
}

onehot_aligned <- onehot[idx, ]
tissue_cols <- setdiff(names(onehot_aligned), c("SAMPID", "tissue"))
if (length(tissue_cols) == 0L) {
  stop("No tissue assignment columns found in metadata")
}

X <- as.matrix(onehot_aligned[, ..tissue_cols])
storage.mode(X) <- "double"
if (!all(rowSums(X) == 1)) {
  stop("Each sample must have exactly one tissue assignment")
}

tissue <- factor(tissue_cols[max.col(X, ties.method = "first")], levels = tissue_cols)
col_data <- data.frame(
  tissue = tissue,
  row.names = sample_ids,
  stringsAsFactors = FALSE
)

keep <- rowSums(counts >= 10L) >= 2L
counts <- counts[keep, , drop = FALSE]
filtered_gene_ids <- gene_ids[keep]
storage.mode(counts) <- "integer"

if (nrow(counts) == 0L) {
  stop("No genes remain after low-count filtering.")
}

checkpoint(counts, "filtered_counts")

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = col_data,
  design = ~ 1
)
dds <- estimateSizeFactors(dds, type = "ratio")
normalized_counts <- counts(dds, normalized = TRUE)
checkpoint(dds, "dds_after_size_factors")
checkpoint(normalized_counts, "normalized_counts")

annot_dt <- data.table(
  gene_id = filtered_gene_ids,
  ensembl_id = sub("\\..*$", "", filtered_gene_ids)
)
annot_dt[, hgnc_symbol := unname(mapIds(
  org.Hs.eg.db, keys = ensembl_id, keytype = "ENSEMBL",
  column = "SYMBOL", multiVals = "first"
))]
annot_dt[, gene_name := unname(mapIds(
  org.Hs.eg.db, keys = ensembl_id, keytype = "ENSEMBL",
  column = "GENENAME", multiVals = "first"
))]

tissue_counts <- colSums(X)
mean_matrix <- sweep(normalized_counts %*% X, 2L, tissue_counts, "/")
mean_dt <- cbind(data.table(gene_id = filtered_gene_ids), as.data.table(mean_matrix))
setnames(mean_dt, c("gene_id", tissue_cols))
mean_dt <- merge(
  annot_dt[, .(gene_id, gene_name, hgnc_symbol)], mean_dt,
  by = "gene_id", all.y = TRUE, sort = FALSE
)
fwrite(mean_dt, file.path(out_dir, "gene_expression_tissue_means.tsv"), sep = "\t", quote = FALSE)

var_dt <- data.table(
  gene_id = filtered_gene_ids,
  sd_across_tissues = apply(mean_matrix, 1L, sd)
)
var_dt <- merge(
  annot_dt[, .(gene_id, gene_name, hgnc_symbol)], var_dt,
  by = "gene_id", all.y = TRUE, sort = FALSE
)
setorder(var_dt, -sd_across_tissues)
fwrite(var_dt, file.path(out_dir, "gene_variability_across_tissues.tsv"), sep = "\t", quote = FALSE)
fwrite(var_dt[1:min(1000, .N)], file.path(out_dir, "top1000_variable_genes_across_tissues.tsv"), sep = "\t", quote = FALSE)

log_expr <- log2(normalized_counts + 1)
checkpoint(log_expr, "log2_normalized_counts")

design <- X
colnames(design) <- tissue_cols
n_tissues <- length(tissue_cols)
if (n_tissues < 2L) {
  stop("At least two tissues are required for tissue-vs-rest contrasts.")
}

contrast_mat <- matrix(
  -1 / (n_tissues - 1), nrow = n_tissues, ncol = n_tissues,
  dimnames = list(tissue_cols, paste0("vs_rest_", make.names(tissue_cols)))
)
diag(contrast_mat) <- 1

fit <- lmFit(log_expr, design)
fit2 <- eBayes(contrasts.fit(fit, contrast_mat), trend = TRUE, robust = TRUE)
all_tissue_sig <- vector("list", n_tissues)
sig_counts <- integer(n_tissues)

summarise_tissue_results <- function(model_fit, tissue_name, annotation_dt, alpha = 0.05) {
  coef_name <- paste0("vs_rest_", make.names(tissue_name))
  tt <- as.data.table(topTable(
    model_fit, coef = coef_name, number = Inf,
    adjust.method = "BH", sort.by = "P"
  ), keep.rownames = "gene_id")

  tt <- merge(
    annotation_dt[, .(gene_id, gene_name, hgnc_symbol)], tt,
    by = "gene_id", all.y = TRUE, sort = FALSE
  )
  tt <- tt[!is.na(adj.P.Val) & adj.P.Val <= alpha]
  tt[, tissue := tissue_name]
  setcolorder(tt, c("tissue", "gene_id", "gene_name", "hgnc_symbol"))
  tt
}

for (j in seq_along(tissue_cols)) {
  tissue_name <- tissue_cols[[j]]
  tt <- summarise_tissue_results(fit2, tissue_name, annot_dt)
  sig_counts[[j]] <- nrow(tt)
  safe_name <- gsub("[^A-Za-z0-9]+", "_", tissue_name)
  fwrite(tt, file.path(out_dir, paste0("significant_FDR0.05_", safe_name, ".tsv")), sep = "\t", quote = FALSE)
  fwrite(tt[logFC > 0], file.path(out_dir, paste0("significant_FDR0.05_", safe_name, "_upregulated.csv")), sep = ",", quote = FALSE)
  all_tissue_sig[[j]] <- tt
}

all_sig_dt <- rbindlist(all_tissue_sig, fill = TRUE)
fwrite(all_sig_dt, file.path(out_dir, "significant_FDR0.05_all_tissues.tsv"), sep = "\t", quote = FALSE)
fwrite(all_sig_dt[logFC > 0], file.path(out_dir, "significant_FDR0.05_all_tissues_upregulated.csv"), sep = ",", quote = FALSE)
checkpoint(all_sig_dt, "all_significant_tissue_results")

fwrite(data.table(
  metric = c("n_genes_input", "n_genes_after_filter", "n_samples", "n_tissues", "normalization"),
  value = c(nrow(gct_obj2@mat), nrow(counts), ncol(counts), n_tissues, "DESeq2 median-of-ratios")
), file.path(out_dir, "analysis_summary.tsv"), sep = "\t", quote = FALSE)

fwrite(data.table(
  tissue = tissue_cols,
  n_samples = as.integer(tissue_counts),
  n_significant_FDR0.05 = as.integer(sig_counts)
), file.path(out_dir, "tissue_sample_counts.tsv"), sep = "\t", quote = FALSE)

fwrite(data.table(
  sample_id = sample_ids,
  size_factor = sizeFactors(dds)
), file.path(out_dir, "sample_size_factors.tsv"), sep = "\t", quote = FALSE)

cat("Analysis complete. Output directory:", out_dir, "\n")
