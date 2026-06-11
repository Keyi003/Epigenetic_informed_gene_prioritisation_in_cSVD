#09_chromatin_states.R
#script to define chromatin states: active promoter region is H3K4me3 and H3K27ac and close to TSS region 
#Enhancer is h3k27ac without k4me3 overlap
#Bivalent region is h3k4me3 and h3k27me3 overlap near TSS
#Logic: identify promoters and enhancers from consensus peaks (all 20 samples merged)

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

ac_gr  <- readPeakFile("~/svd_analysis/counts/H3K27ac_consensus_peaks.bed")
k4_gr  <- readPeakFile("~/svd_analysis/counts/H3K4me3_consensus_peaks.bed")
k27_gr <- readPeakFile("~/svd_analysis/counts/H3K27me3_consensus_peaks.bed")


tss_gr <- promoters(txdb, upstream = 2000, downstream = 500)

hits <- findOverlaps(k4_gr, ac_gr)
active_promoter_gr <- pintersect(k4_gr[queryHits(hits)], ac_gr[subjectHits(hits)])
mcols(active_promoter_gr) <- NULL
active_promoter_gr <- subsetByOverlaps(active_promoter_gr, tss_gr)

ac_no_k4 <- subsetByOverlaps(ac_gr, k4_gr, invert = TRUE)
enhancer_gr <- subsetByOverlaps(ac_no_k4, tss_gr, invert = TRUE)

hits_biv <- findOverlaps(k4_gr, k27_gr)
bivalent_gr <- pintersect(k4_gr[queryHits(hits_biv)], k27_gr[subjectHits(hits_biv)])
mcols(bivalent_gr) <- NULL
bivalent_gr <- subsetByOverlaps(bivalent_gr, tss_gr)

cat("  Active promoters:", length(active_promoter_gr), "\n")
cat("  Active enhancers:", length(enhancer_gr), "\n")
cat("  Bivalent promoters:", length(bivalent_gr), "\n\n")

states <- list(
  active_promoter = active_promoter_gr,
  enhancer        = enhancer_gr,
  bivalent        = bivalent_gr
)

promoter_window <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)

tagMatrixList <- list(
  active_promoter = getTagMatrix(active_promoter_gr, windows = promoter_window),
  enhancer        = getTagMatrix(enhancer_gr, windows = promoter_window),
  bivalent        = getTagMatrix(bivalent_gr, windows = promoter_window)
)

png("~/svd_analysis/figures/QC_covplot_active_promoter_intersected.png",
    width = 10, height = 8, units = "in", res = 300)
covplot(active_promoter_gr)
dev.off()

png("~/svd_analysis/figures/QC_covplot_enhancer.png",
    width = 10, height = 8, units = "in", res = 300)
covplot(enhancer_gr)
dev.off()

png("~/svd_analysis/figures/QC_covplot_bivalent_intersected.png",
    width = 10, height = 8, units = "in", res = 300)
covplot(bivalent_gr)
dev.off()

png("~/svd_analysis/figures/QC_avgprof_chromatin_states_intersected.png",
    width = 8, height = 5, units = "in", res = 300)
plotAvgProf(tagMatrixList, xlim = c(-3000, 3000),
            xlab = "Genomic Region (5'->3')",
            ylab = "Peak Frequency")
dev.off()

png("~/svd_analysis/figures/QC_tagheatmap_chromatin_states_intersected.png",
    width = 12, height = 8, units = "in", res = 300)
tagHeatmap(tagMatrixList)
dev.off()

peak_genes <- list()

for (state_name in names(states)) {
  cat("Annotating", state_name, "===\n")
  
  peaks <- states[[state_name]]
  
  peakAnno <- annotatePeak(peaks, TxDb = txdb, level = "gene",
                           annoDb = "org.Hs.eg.db",
                           tssRegion = c(-2000, 500))
  
  png(paste0("~/svd_analysis/figures/QC_chipseeker_intersected_", state_name, ".png"),
      width = 8, height = 6, units = "in", res = 300)
  print(plotAnnoPie(peakAnno))
  dev.off()
  
  anno_df <- as.data.frame(peakAnno)
  genes <- unique(na.omit(anno_df$SYMBOL))
  peak_genes[[state_name]] <- genes
  
  cat("  Peaks annotated:", nrow(anno_df), "\n")
  cat("  Unique genes:", length(genes), "\n")
}

saveRDS(peak_genes, "~/svd_analysis/qc/09_chromatin_state_gene_lists_intersected.rds")
cat("\nGene lists saved to: ~/svd_analysis/qc/09_chromatin_state_gene_lists_intersected.rds\n")
