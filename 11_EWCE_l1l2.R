library(EWCE)
library(ggplot2)

set.seed(1234)

# Paths
ctd_rds_path <- "~/svd_analysis/reference/terekhova_ctd/terekhova_pbmc_ctd_L1L2.rds"

ctd <- readRDS("~/svd_analysis/reference/terekhova_ctd/terekhova_pbmc_ctd_L1L2.rds")

peak_genes <- readRDS("~/svd_analysis/qc/09_chromatin_state_gene_lists_intersected.rds")

cat("L2 cell types:", ncol(ctd[[2]]$specificity), "\n")
cat(paste(colnames(ctd[[2]]$specificity), collapse = "\n"), "\n")

reps <- 100000
results <- list()

for (mark_name in names(peak_genes)) {
  cat("\n--- Testing", mark_name, "at L2\n")
  hits <- peak_genes[[mark_name]]
  
  full_results <- EWCE::bootstrap_enrichment_test(
    sct_data = ctd,
    hits = hits,
    reps = reps,
    annotLevel = 2,
    sctSpecies = "human",
    genelistSpecies = "human"
  )
  
  res <- full_results$results
  res$list <- mark_name
  results[[mark_name]] <- res
  
  sig <- res[res$q < 0.05, ]
  cat("  Significant (q < 0.05):", nrow(sig), "\n")
  if (nrow(sig) > 0) {
    sig <- sig[order(sig$p), ]
    for (i in 1:min(10, nrow(sig))) {
      cat(sprintf("    %s: FC=%.3f, SD=%.2f, q=%.2e\n",
                  sig$CellType[i], sig$fold_change[i],
                  sig$sd_from_mean[i], sig$q[i]))
    }
    if (nrow(sig) > 10) cat("    ... and", nrow(sig) - 10, "more\n")
  }
}

combined_res <- do.call(rbind, results)
write.csv(combined_res,
          "~/svd_analysis/qc/08c_ewce_terekhova_L2_results.csv",
          row.names = FALSE)
dir.create("~/svd_analysis/figures", showWarnings = FALSE, recursive = TRUE)

for (mark_name in names(results)) {
  res <- results[[mark_name]]
  res$list <- NULL
  
  plot_list <- EWCE::ewce_plot(
    total_res = res,
    mtc_method = "BH",
    ctd = ctd
  )
  
  ggsave(paste0("~/svd_analysis/figures/QC_ewce_terekhova_L2_", mark_name, ".png"),
         plot_list$plain, width = 16, height = 6, dpi = 300)
}

combined_res_plot <- do.call(rbind, results)
plot_combined <- EWCE::ewce_plot(
  total_res = combined_res_plot,
  mtc_method = "BH",
  ctd = ctd
)

ggsave("~/svd_analysis/figures/QC_ewce_terekhova_L2_combined.png",
       plot_combined$plain, width = 20, height = 14, dpi = 300)

cat("Plots saved\n")
cat("\nResults saved\n")