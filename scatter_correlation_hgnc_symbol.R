#!/usr/bin/env Rscript
# Scatterplot and Pearson correlation of log2FoldChange between resSig<TISSUE>FDR005.xlsx
# and GTEXresSig<TISSUE>FDR005.xlsx, matched by hgnc_symbol.
#
# Usage:
#   Rscript scatter_correlation_hgnc_symbol.R <TISSUE> "<Label>"
# Example:
#   Rscript scatter_correlation_hgnc_symbol.R PROSTATE "Prostate"
#
# Duplicate hgnc_symbol entries within each file are collapsed by averaging
# their log2FoldChange before matching.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(writexl)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript scatter_correlation_hgnc_symbol.R <TISSUE> [\"Label\"]")
}

tissue <- args[1]
label <- if (length(args) >= 2) args[2] else tools::toTitleCase(tolower(tissue))

base <- getwd()
f1 <- file.path(base, paste0('resSig', tissue, 'FDR005.xlsx'))
f2 <- file.path(base, paste0('GTEXresSig', tissue, 'FDR005.xlsx'))

d1 <- read_excel(f1)
d2 <- read_excel(f2)

collapse_symbol <- function(d) {
  d %>%
    transmute(
      hgnc_symbol = trimws(as.character(hgnc_symbol)),
      log2FoldChange = suppressWarnings(as.numeric(log2FoldChange))
    ) %>%
    filter(!is.na(hgnc_symbol), hgnc_symbol != '', !is.na(log2FoldChange)) %>%
    group_by(hgnc_symbol) %>%
    summarize(log2FoldChange = mean(log2FoldChange), .groups = 'drop')
}

keep1 <- collapse_symbol(d1)
keep2 <- collapse_symbol(d2)

merged <- inner_join(keep1, keep2, by = 'hgnc_symbol', suffix = c('_resSig', '_GTEx'))

n <- nrow(merged)
ct <- cor.test(merged$log2FoldChange_resSig, merged$log2FoldChange_GTEx, method = 'pearson')
r <- unname(ct$estimate)
p <- ct$p.value

out_png <- file.path(base, paste0('resSig_vs_GTEX_', tissue, '_scatter_hgnc_symbol.png'))
plot <- ggplot(merged, aes(x = log2FoldChange_resSig, y = log2FoldChange_GTEx)) +
  geom_point(alpha = 0.45, size = 1.2, color = '#1f78b4') +
  geom_smooth(method = 'lm', se = FALSE, color = '#e31a1c', linewidth = 0.8) +
  theme_classic(base_size = 13) +
  labs(
    title = label,
    x = 'E-MTAB-2836 log2FC',
    y = 'GTEx log2FC',
    subtitle = paste0('n = ', n, ', Pearson r = ', sprintf('%.4f', r), ', p = ', format(p, scientific = TRUE, digits = 3))
  )
ggsave(out_png, plot, width = 7.2, height = 5.6, dpi = 300)

out_xlsx <- file.path(base, paste0('matched_log2FoldChange_', tissue, '_hgnc_symbol.xlsx'))
merged_out <- merged %>%
  rename(
    !!paste0('log2FoldChange_resSig', tissue) := log2FoldChange_resSig,
    !!paste0('log2FoldChange_GTEX', tissue) := log2FoldChange_GTEx
  ) %>%
  arrange(hgnc_symbol)
write_xlsx(merged_out, out_xlsx)

cat(sprintf('n=%d\n', n))
cat(sprintf('r=%.10f\n', r))
cat(sprintf('p=%.16g\n', p))
cat(sprintf('plot=%s\n', out_png))
cat(sprintf('matched_table=%s\n', out_xlsx))
