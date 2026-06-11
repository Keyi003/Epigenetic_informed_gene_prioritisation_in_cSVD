# Combined PCA of all three marks to check for mislabelling


library(ggplot2)
library(ggrepel)
library(edgeR)

dir.create("~/svd_analysis/figures/svg figure versions", showWarnings = FALSE, recursive = TRUE)


fc_ac    <- readRDS("~/svd_analysis/counts/featurecounts_H3K27ac.rds")
fc_k4    <- readRDS("~/svd_analysis/counts/featurecounts_H3K4me3.rds")
fc_k27   <- readRDS("~/svd_analysis/counts/featurecounts_H3K27me3.rds")


counts_ac  <- fc_ac$counts
counts_k4  <- fc_k4$counts
counts_k27 <- fc_k27$counts

colnames(counts_ac)  <- paste0(gsub("unique_pbmc_ac_12_12_2025_|\\.sorted\\.bam", "", colnames(counts_ac)), "_ac")
colnames(counts_k4)  <- paste0(gsub("unique_pbmc_k4me3_12_12_2025_|\\.sorted\\.bam", "", colnames(counts_k4)), "_k4me3")
colnames(counts_k27) <- paste0(gsub("unique_pbmc_k27me3_12_12_2025_|\\.sorted\\.bam", "", colnames(counts_k27)), "_k27me3")


peaks_ac  <- rownames(counts_ac)
peaks_k4  <- rownames(counts_k4)
peaks_k27 <- rownames(counts_k27)
all_peaks <- unique(c(peaks_ac, peaks_k4, peaks_k27))

fill_matrix <- function(counts, all_peaks) {
  mat <- matrix(0, nrow = length(all_peaks), ncol = ncol(counts))
  rownames(mat) <- all_peaks
  colnames(mat) <- colnames(counts)
  mat[rownames(counts), ] <- counts
  mat
}

full_ac  <- fill_matrix(counts_ac, all_peaks)
full_k4  <- fill_matrix(counts_k4, all_peaks)
full_k27 <- fill_matrix(counts_k27, all_peaks)

# Combine: 60 columns (20 samples x 3 marks)
combined <- cbind(full_ac, full_k4, full_k27)


logcpm <- cpm(combined, log = TRUE, prior.count = 1)


pca <- prcomp(t(logcpm), scale. = TRUE)
var_explained <- summary(pca)$importance[2, 1:5] * 100


pca_df <- data.frame(
  Sample_Mark = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

pca_df$Sample <- gsub("_(ac|k4me3|k27me3)$", "", pca_df$Sample_Mark)
pca_df$Mark <- gsub("^R[0-9]+_", "", pca_df$Sample_Mark)
pca_df$Mark <- factor(pca_df$Mark, levels = c("ac", "k4me3", "k27me3"),
                      labels = c("H3K27ac", "H3K4me3", "H3K27me3"))


mark_colours <- c("H3K27ac"  = "#D4726A",
                  "H3K4me3"  = "#5B9E8F",
                  "H3K27me3" = "#6B85A3")

p <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Mark, label = Sample)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, show.legend = FALSE, seed = 42, max.overlaps = 30) +
  scale_colour_manual(values = mark_colours) +
  labs(title = "Combined PCA — all three histone marks",
       colour = "Histone\nmark",
       x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  theme_minimal(base_size = 13)

ggsave("~/svd_analysis/figures/svg figure versions/QC_pca_all_marks.svg", p,
       width = 9, height = 6)

cat("Saved to ~/svd_analysis/figures/svg figure versions/QC_pca_all_marks.svg\n")
cat("Variance explained — PC1:", round(var_explained[1], 1),
    "% PC2:", round(var_explained[2], 1), "%\n")