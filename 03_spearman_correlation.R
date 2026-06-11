library(edgeR)
library(pheatmap)


marks <- list(
  ac = list(
    rds       = "~/svd_analysis/counts/featurecounts_H3K27ac.rds",
    name_gsub = "unique_pbmc_ac_12_12_2025_|\\.sorted\\.bam",
    title     = "H3K27ac",
    out_pdf   = "~/svd_analysis/figures/svg figure versions/QC_spearman_heatmap_ac.svg"
  ),
  k4me3 = list(
    rds       = "~/svd_analysis/counts/featurecounts_H3K4me3.rds",
    name_gsub = "unique_pbmc_k4me3_12_12_2025_|\\.sorted\\.bam",
    title     = "H3K4me3",
    out_pdf   = "~/svd_analysis/figures/svg figure versions/QC_spearman_heatmap_k4me3.svg"
  ),
  k27me3 = list(
    rds       = "~/svd_analysis/counts/featurecounts_H3K27me3.rds",
    name_gsub = "unique_pbmc_k27me3_12_12_2025_|\\.sorted\\.bam",
    title     = "H3K27me3",
    out_pdf   = "~/svd_analysis/figures/svg figure versions/QC_spearman_heatmap_k27me3.svg"
  )
)


dir.create("~/svd_analysis/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("~/svd_analysis/figures/svg figure versions", showWarnings = FALSE, recursive = TRUE)
dir.create("~/svd_analysis/qc", showWarnings = FALSE, recursive = TRUE)




for (mark_name in names(marks)) {

  
  # Load counts
  fc <- readRDS(m$rds)
  counts <- fc$counts
  colnames(counts) <- gsub(m$name_gsub, "", colnames(counts))
  
  # Normalise to log-CPM
  logcpm <- cpm(counts, log = TRUE, prior.count = 1)
  
  cat("Matrix:", nrow(logcpm), "peaks x", ncol(logcpm), "samples\n")
  
  # Spearman correlation
  cor_matrix <- cor(logcpm, method = "spearman")
  
  # Plot heatmap
  svg(m$out_pdf, width = 10, height = 8)
  pheatmap(cor_matrix, display_numbers = TRUE, number_format = "%.2f",
           main = paste0(m$title, ": Spearman Correlation Between Samples"))
  dev.off()
  cat("Heatmap saved:", m$out_pdf, "\n")
  
  # Summary statistics
  mean_cor <- rowMeans(cor_matrix)
  min_cor <- apply(cor_matrix, 1, function(x) min(x[x < 1]))
  
  cat("\nMean Spearman correlation per sample:\n")
  print(round(sort(mean_cor), 4))
  
  cat("\nMinimum pairwise correlation per sample:\n")
  print(round(sort(min_cor), 4))
  
  flagged_mean <- names(mean_cor[mean_cor < 0.75])
  flagged_min  <- names(min_cor[min_cor < 0.65])
  
  cat("\nFlagged (mean r < 0.75):",
      ifelse(length(flagged_mean) == 0, "none", paste(flagged_mean, collapse = ", ")), "\n")
  cat("Flagged (any pairwise r < 0.65):",
      ifelse(length(flagged_min) == 0, "none", paste(flagged_min, collapse = ", ")), "\n")
  
  # Log summary
  log_lines <- c(log_lines, "",
                 paste0("--- ", mark_name, " ---"),
                 paste("Peaks:", nrow(logcpm)),
                 paste("Samples:", ncol(logcpm)),
                 "",
                 "Mean Spearman correlation per sample (sorted):",
                 capture.output(print(round(sort(mean_cor), 4))),
                 "",
                 "Minimum pairwise correlation per sample (sorted):",
                 capture.output(print(round(sort(min_cor), 4))),
                 "",
                 paste("Overall mean correlation:", round(mean(cor_matrix[lower.tri(cor_matrix)]), 4)),
                 paste("Overall min correlation:", round(min(cor_matrix[lower.tri(cor_matrix)]), 4)),
                 paste("Overall max correlation:", round(max(cor_matrix[lower.tri(cor_matrix)]), 4)),
                 paste("Flagged (mean r < 0.75):",
                       ifelse(length(flagged_mean) == 0, "none", paste(flagged_mean, collapse = ", "))),
                 paste("Flagged (any pairwise r < 0.65):",
                       ifelse(length(flagged_min) == 0, "none", paste(flagged_min, collapse = ", ")))
  )
}


writeLines(log_lines, "~/svd_analysis/qc/03_spearman_summary.txt")
cat("\n=== Summary saved to ~/svd_analysis/qc/03_spearman_summary.txt ===\n")