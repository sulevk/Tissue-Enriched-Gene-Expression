# ============================================================
# Calculate conventional Tau tissue-specificity scores
# E-MTAB-2836: 32 human tissues
#
# Input 1: gene-by-sample expression matrix
#           First column = gene ID
#           Remaining columns = sample IDs
#
# Input 2: sample metadata
#           Columns: sample_id, tissue
#
# Output:
#   1. Classical Tau score for each gene
#   2. Mean linear expression by tissue
#   3. Summary of Tau scores
#
# IMPORTANT:
# Expression values used for Tau must be non-negative and
# on a linear scale. Do not use log-transformed values or
# DESeq2 log2FC values.
# ============================================================

# Install DESeq2 first if needed:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("DESeq2")

suppressPackageStartupMessages({
    library(DESeq2)
})

# ---- 1. User settings ---------------------------------------

expression_file <- "EMTAB2836_expression_matrix.tsv"
metadata_file   <- "EMTAB2836_sample_metadata.tsv"

# Choose one:
# "normalized" = already-normalized, linear-scale expression
# "raw_counts" = raw, non-negative integer read counts
input_type <- "normalized"

gene_id_column <- 1
sample_id_column <- "sampleID"
tissue_column <- "tissue"

output_prefix <- "EMTAB2836_classical_Tau"

# ---- 2. Read expression matrix ------------------------------

expr_df <- read.delim(
    expression_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

gene_ids <- as.character(expr_df[[gene_id_column]])

expr <- as.matrix(expr_df[, -gene_id_column, drop = FALSE])
storage.mode(expr) <- "numeric"

rownames(expr) <- gene_ids

if (anyDuplicated(rownames(expr))) {
    stop("Duplicate gene IDs detected. Resolve them before analysis.")
}

if (anyDuplicated(colnames(expr))) {
    stop("Duplicate sample IDs detected in expression matrix.")
}

if (any(!is.finite(expr)) || any(expr < 0)) {
    stop("Expression matrix contains missing, non-finite, or negative values.")
}

# ---- 3. Read and check sample metadata -----------------------

meta <- read.delim(
    metadata_file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

required_columns <- c(sample_id_column, tissue_column)

if (!all(required_columns %in% colnames(meta))) {
    stop("Metadata must contain columns named: ",
         paste(required_columns, collapse = ", "))
}

meta[[sample_id_column]] <- as.character(meta[[sample_id_column]])
meta[[tissue_column]] <- as.character(meta[[tissue_column]])

if (anyDuplicated(meta[[sample_id_column]])) {
    stop("Duplicate sample IDs detected in metadata.")
}

# Keep samples present in both files and align their order.
common_samples <- colnames(expr)[
    colnames(expr) %in% meta[[sample_id_column]]
]

if (length(common_samples) == 0) {
    stop("No matching sample IDs between expression matrix and metadata.")
}

if (length(common_samples) != ncol(expr)) {
    warning("Some expression-matrix samples are absent from metadata.")
}

if (length(common_samples) != nrow(meta)) {
    warning("Some metadata samples are absent from expression matrix.")
}

expr <- expr[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta[[sample_id_column]]), ,
             drop = FALSE]

stopifnot(
    identical(colnames(expr), meta[[sample_id_column]])
)

meta[[tissue_column]] <- factor(meta[[tissue_column]])

cat("Genes:", nrow(expr), "\n")
cat("Samples matched:", ncol(expr), "\n")
cat("Tissues represented:", nlevels(meta[[tissue_column]]), "\n")
print(table(meta[[tissue_column]]))

if (nlevels(meta[[tissue_column]]) != 32) {
    warning("The matched data do not contain exactly 32 tissues.")
}

# ---- 4. Normalize if input is raw counts --------------------

if (input_type == "raw_counts") {

    if (any(expr != round(expr))) {
        stop("Raw counts must be integer-valued.")
    }

    dds <- DESeqDataSetFromMatrix(
        countData = round(expr),
        colData = data.frame(
            tissue = meta[[tissue_column]],
            row.names = common_samples
        ),
        design = ~ tissue
    )

    dds <- estimateSizeFactors(dds)

    expr <- counts(dds, normalized = TRUE)

} else if (input_type != "normalized") {
    stop("input_type must be 'normalized' or 'raw_counts'.")
}

# ---- 5. Calculate arithmetic mean expression per tissue -----

tissue_levels <- levels(meta[[tissue_column]])

mean_expr <- vapply(
    tissue_levels,
    function(t) {
        samples <- which(meta[[tissue_column]] == t)
        rowMeans(expr[, samples, drop = FALSE])
    },
    FUN.VALUE = numeric(nrow(expr))
)

rownames(mean_expr) <- rownames(expr)
colnames(mean_expr) <- tissue_levels

# Save the tissue-level expression matrix for inspection.
write.table(
    data.frame(
        gene_id = rownames(mean_expr),
        mean_expr,
        check.names = FALSE
    ),
    paste0(output_prefix, "_mean_expression_by_tissue.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ---- 6. Calculate conventional Tau --------------------------

max_expr <- apply(mean_expr, 1, max)

# Tau is undefined when a gene has zero expression in all
# tissues. Leave these genes as NA rather than assigning Tau = 1.
tau <- rep(NA_real_, nrow(mean_expr))
names(tau) <- rownames(mean_expr)

expressed <- is.finite(max_expr) & max_expr > 0

tau[expressed] <- rowSums(
    1 - sweep(
        mean_expr[expressed, , drop = FALSE],
        1,
        max_expr[expressed],
        FUN = "/"
    )
) / (ncol(mean_expr) - 1)

# Small numerical deviations can occur through floating-point math.
tau[!is.na(tau)] <- pmin(1, pmax(0, tau[!is.na(tau)]))

# Identify the tissue with the highest mean expression.
top_index <- max.col(mean_expr, ties.method = "first")

results <- data.frame(
    gene_id = rownames(mean_expr),
    Tau = unname(tau),
    top_tissue = tissue_levels[top_index],
    maximum_mean_expression = max_expr,
    all_tissues_zero = !expressed,
    stringsAsFactors = FALSE
)

results <- results[order(results$Tau, decreasing = TRUE,
                         na.last = TRUE), ]

write.table(
    results,
    paste0(output_prefix, "_gene_scores.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ---- 7. Summary and basic checks -----------------------------

cat("\nTau calculation complete.\n")
cat("Genes analysed:", nrow(results), "\n")
cat("Genes with defined Tau:", sum(!is.na(results$Tau)), "\n")
cat("Genes zero across all tissues:", sum(!expressed), "\n")
cat("\nTau summary:\n")
print(summary(results$Tau))

cat("\nNumber of genes by top-expression tissue:\n")
print(sort(table(results$top_tissue), decreasing = TRUE))

write.table(
    data.frame(
        statistic = c(
            "Genes_total",
            "Genes_with_defined_Tau",
            "Genes_zero_in_all_tissues",
            "Tissues"
        ),
        value = c(
            nrow(results),
            sum(!is.na(results$Tau)),
            sum(!expressed),
            ncol(mean_expr)
        )
    ),
    paste0(output_prefix, "_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# End of script
# ============================================================