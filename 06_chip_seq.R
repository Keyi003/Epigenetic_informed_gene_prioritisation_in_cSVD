library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


marks <- list(
  ac     = "~/svd_analysis/counts/H3K27ac_consensus_peaks.bed",
  k4me3  = "~/svd_analysis/counts/H3K4me3_consensus_peaks.bed",
  k27me3 = "~/svd_analysis/counts/H3K27me3_consensus_peaks.bed"
)

peak_genes <- list()

for (mark_name in names(marks)) {
  cat("=== Annotating", mark_name, "===\n")
  
  peaks <- readPeakFile(marks[[mark_name]])
  
  peakAnno <- annotatePeak(peaks, TxDb = txdb, level = "gene",
                           annoDb = "org.Hs.eg.db",
                           tssRegion = c(-2000, 500))
  
  # Save annotation distribution plot
  png(paste0("~/svd_analysis/figures/QC_chipseeker_", mark_name, ".png"),
      width = 8, height = 6, units = "in", res = 300)
  print(plotAnnoPie(peakAnno))
  dev.off()
  
  # Extract gene symbols
  anno_df <- as.data.frame(peakAnno)
  genes <- unique(na.omit(anno_df$SYMBOL))
  peak_genes[[mark_name]] <- genes
  
  cat("  Peaks annotated:", nrow(anno_df), "\n")
  cat("  Unique genes:", length(genes), "\n")
}