
# 08_celltype_verification.R
# Cell-type enrichment of CUT&Tag consensus peak genes
# TEST RUN: using ewceData mouse brain CTD as negative control
# Following Nathan Skene's EWCE workflow


library(EWCE)
library(ewceData)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggplot2)
library(patchwork)

set.seed(1234)

dir.create("~/svd_analysis/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("~/svd_analysis/qc", showWarnings = FALSE, recursive = TRUE)


cat("=== PART A: Peak-to-gene annotation ===\n")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
all_genes_gr <- genes(txdb)

annotate_peaks <- function(bed_path, mark_name) {
  cat("--- Annotating", mark_name, "---\n")
  
  # Read consensus peaks
  peaks <- import(bed_path)
  cat("  Peaks loaded:", length(peaks), "\n")
  
  # Find nearest gene for each peak
  nearest_idx <- nearest(peaks, all_genes_gr)
  
  # Remove peaks with no nearest gene
  valid <- !is.na(nearest_idx)
  nearest_entrez <- all_genes_gr[nearest_idx[valid]]$gene_id
  
  # Convert Entrez IDs to gene symbols
  symbols <- mapIds(org.Hs.eg.db,
                    keys = nearest_entrez,
                    keytype = "ENTREZID",
                    column = "SYMBOL",
                    multiVals = "first")
  
  genes <- unique(na.omit(symbols))
  cat("  Unique genes mapped:", length(genes), "\n")
  return(genes)
}

peak_genes <- list(
  H3K27ac  = annotate_peaks("~/svd_analysis/counts/H3K27ac_consensus_peaks.bed", "H3K27ac"),
  H3K4me3  = annotate_peaks("~/svd_analysis/counts/H3K4me3_consensus_peaks.bed", "H3K4me3"),
  H3K27me3 = annotate_peaks("~/svd_analysis/counts/H3K27me3_consensus_peaks.bed", "H3K27me3")
)

# Save gene lists for later use
saveRDS(peak_genes, "~/svd_analysis/qc/08_peak_gene_lists.rds")


cat("\n=== PART B: Loading mouse brain CTD (negative control) ===\n")

# This is the Zeisel et al. 2015 mouse cortex/hippocampus dataset
# Cell types: astrocytes, interneurons, microglia, oligodendrocytes,
#             pyramidal CA1, pyramidal SS, endothelial/mural

ctd <- ewceData::ctd()

cat("CTD loaded successfully\n")
cat("Cell types in CTD:\n")
print(names(ctd[[1]]$specificity[1,]))



cat("\n=== PART C: Running EWCE bootstrap enrichment ===\n")

reps <- 100
annotLevel <- 1

results <- list()

for (mark_name in names(peak_genes)) {
  cat("\n--- Testing", mark_name, "---\n")
  
  hits <- peak_genes[[mark_name]]
  
  full_results <- EWCE::bootstrap_enrichment_test(
    sct_data = ctd,
    hits = hits,
    reps = reps,
    annotLevel = annotLevel,
    sctSpecies = "mouse",
    genelistSpecies = "human"
  )
  
  # Store results and tag with mark name
  res <- full_results$results
  res$list <- mark_name
  results[[mark_name]] <- res
  
  cat("  Significant cell types (q < 0.05):\n")
  sig <- res[res$q < 0.05, ]
  if (nrow(sig) > 0) {
    print(sig[order(sig$p), c("CellType", "p", "fold_change", "sd_from_mean", "q")])
  } else {
    cat("  None (expected for brain CTD vs PBMC peaks)\n")
  }
}


for (mark_name in names(results)) {
  res <- results[[mark_name]]
  res$list <- NULL  # remove list column for individual plots
  
  plot_list <- EWCE::ewce_plot(
    total_res = res,
    mtc_method = "BH"
  )
  
  ggsave(
    paste0("~/svd_analysis/figures/QC_ewce_brain_", mark_name, ".png"),
    plot_list$plain,
    width = 8, height = 5, dpi = 300
  )
  cat("  Saved plot for", mark_name, "\n")
}


combined_res <- do.call(rbind, results)

plot_combined <- EWCE::ewce_plot(
  total_res = combined_res,
  mtc_method = "BH"
)

ggsave(
  "~/svd_analysis/figures/QC_ewce_brain_combined.png",
  plot_combined$plain,
  width = 12, height = 8, dpi = 300
)
cat("  Saved combined plot\n")



write.csv(combined_res, "~/svd_analysis/qc/08_ewce_brain_results.csv", row.names = FALSE)

