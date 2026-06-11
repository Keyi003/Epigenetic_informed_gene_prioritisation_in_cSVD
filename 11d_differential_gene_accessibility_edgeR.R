library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GenomicRanges)
library(data.table)
library(ggplot2)

candidates <- fread("~/svd_analysis/results/edgeR_3rd_run/ac_candidates_p0.05_lfc1.csv")

cat("Number of candidate peaks:", nrow(candidates), "\n")
head(candidates)

bed <- fread("~/svd_analysis/counts/H3K27ac_consensus_peaks.bed",
             col.names = c("chr", "start", "end", "name"))

cat("BED rows:", nrow(bed), "\n")
head(bed)

candidate_bed <- bed[bed$name %in% candidates$genes, ]

cat("Candidate BED rows:", nrow(candidate_bed), "\n")

gr <- GRanges(seqnames = candidate_bed$chr,
              ranges   = IRanges(start = candidate_bed$start,
                                 end   = candidate_bed$end),
              name     = candidate_bed$name)

peak_anno <- annotatePeak(gr,
                          TxDb      = TxDb.Hsapiens.UCSC.hg38.knownGene,
                          level     = "gene",
                          annoDb    = "org.Hs.eg.db",
                          tssRegion = c(-2000, 500),
                          verbose   = FALSE)

anno_df <- as.data.frame(peak_anno)

candidate_genes <- unique(na.omit(anno_df$SYMBOL))

cat("Candidate genes (unique):", length(candidate_genes), "\n")
head(candidate_genes)

gr_all <- GRanges(seqnames = bed$chr,
                  ranges   = IRanges(start = bed$start,
                                     end   = bed$end),
                  name     = bed$name)

bg_anno <- annotatePeak(gr_all,
                        TxDb      = TxDb.Hsapiens.UCSC.hg38.knownGene,
                        level     = "gene",
                        annoDb    = "org.Hs.eg.db",
                        tssRegion = c(-2000, 500),
                        verbose   = FALSE)

bg_genes <- unique(na.omit(as.data.frame(bg_anno)$SYMBOL))

cat("Background genes:", length(bg_genes), "\n")

cat("Candidate genes also in background:",
    sum(candidate_genes %in% bg_genes),
    "/", length(candidate_genes), "\n")

ego <- enrichGO(gene          = candidate_genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                universe      = bg_genes,
                readable      = TRUE)

res <- as.data.frame(ego)

cat("Terms tested:", nrow(res), "\n")
cat("Terms at q<0.05:", sum(res$p.adjust < 0.05), "\n")

head(res[, c("Description", "GeneRatio", "BgRatio",
             "pvalue", "p.adjust", "Count")], 20)

fwrite(res, "~/svd_analysis/results/GO_enrichment/ac_GO_BP_full.csv",
       sep = "\t")

plot_df <- head(res[order(res$p.adjust), ], 15)

plot_df$Description <- factor(plot_df$Description,
                              levels = rev(plot_df$Description))

plot_df$GeneCount <- as.integer(sub("/.*", "", plot_df$GeneRatio))

p <- ggplot(plot_df,
            aes(x = -log10(p.adjust),
                y = Description,
                size = GeneCount,
                fill = -log10(p.adjust))) +
  geom_point(shape = 21, colour = "grey20", stroke = 0.3) +
  scale_size_continuous(name = "Gene count", range = c(2, 8)) +
  scale_fill_gradient(name = expression(-log[10](q)),
                      low = "#FDB863", high = "#7F0000") +
  labs(title = "H3K27ac candidate peaks - GO:BP enrichment",
       x = expression(-log[10]("BH-adjusted q")),
       y = NULL) +
  theme_classic(base_size = 11)

ggsave("~/svd_analysis/figures/GO_enrichment/ac_GO_BP.png",
       p, width = 9, height = 6, dpi = 300)


