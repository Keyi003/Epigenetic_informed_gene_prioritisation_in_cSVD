library(EWCE)
library(ggplot2)

set.seed(1234)

# Paths
ctd_rds_path <- "~/svd_analysis/reference/terekhova_ctd/terekhova_pbmc_ctd.rds"

ctd <- readRDS("~/svd_analysis/reference/terekhova_ctd/terekhova_pbmc_ctd.rds")

peak_genes <- readRDS("~/svd_analysis/qc/09_chromatin_state_gene_lists_intersected.rds")

reps <- 100000
results <- list()

for (mark_name in names(peak_genes)) {
  hits <- peak_genes[[mark_name]]
  
  full_results <- EWCE::bootstrap_enrichment_test(
    sct_data = ctd,
    hits = hits,
    reps = reps,
    annotLevel = 1,
    sctSpecies = "human",
    genelistSpecies = "human"
  )
  
  res <- full_results$results
  res$list <- mark_name
  results[[mark_name]] <- res
  
  cat("  Significant (q < 0.05):\n")
  sig <- res[res$q < 0.05, ]
  if (nrow(sig) > 0) {
    sig <- sig[order(sig$p), ]
    for (i in 1:nrow(sig)) {
      cat(sprintf("    %s: FC=%.3f, SD=%.2f, q=%.2e\n",
                  sig$CellType[i], sig$fold_change[i],
                  sig$sd_from_mean[i], sig$q[i]))
    }
  } else {
    cat("    None\n")
  }
}


dir.create("~/svd_analysis/figures", showWarnings = FALSE, recursive = TRUE)

for (mark_name in names(results)) {
  res <- results[[mark_name]]
  res$list <- NULL  # Remove list column for individual (non-faceted) plots
  
  plot_list <- EWCE::ewce_plot(
    total_res = res,
    mtc_method = "BH",
    ctd = ctd         # Needed for dendrogram computation
  )
  
  
  ggsave(paste0("~/svd_analysis/figures/QC_ewce_terekhova_combined_marks_intersected", mark_name, ".png"),
         plot_list$plain, width = 10, height = 5, dpi = 300)
}

combined_res <- do.call(rbind, results)

write.csv(combined_res,
          "~/svd_analysis/qc/08c_ewce_terekhova_results_combined_marks_intersected.csv",
          row.names = FALSE)
cat("  Results saved to: ~/svd_analysis/qc/08c_ewce_terekhova_results_combined_marks_intersected.csv\n")