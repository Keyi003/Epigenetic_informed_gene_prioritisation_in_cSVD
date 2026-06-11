library(edgeR)
library(readxl)
library(data.table)
library(ggplot2)
library(tidyr)

svg_dir <- "~/svd_analysis/figures/svg figure versions"
png_dir <- "~/svd_analysis/figures"
dir.create(svg_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

# Suffix tags every output filename so original/relaxed scripts
# don't collide.
SUFFIX <- "original"

marks <- list(
  ac     = list(rds       = "~/svd_analysis/counts/featurecounts_H3K27ac.rds",
                name_gsub = "unique_pbmc_ac_12_12_2025_|\\.sorted\\.bam",
                qc_sheet  = "AC",
                file_gsub = "pbmc_ac_12_12_2025_",
                label     = "H3K27ac"),
  k4me3  = list(rds       = "~/svd_analysis/counts/featurecounts_H3K4me3.rds",
                name_gsub = "unique_pbmc_k4me3_12_12_2025_|\\.sorted\\.bam",
                qc_sheet  = "k4me3",
                file_gsub = "pbmc_k4me3_12_12_2025_",
                label     = "H3K4me3"),
  k27me3 = list(rds       = "~/svd_analysis/counts/featurecounts_H3K27me3.rds",
                name_gsub = "unique_pbmc_k27me3_12_12_2025_|\\.sorted\\.bam",
                qc_sheet  = "k27me3",
                file_gsub = "pbmc_ke27me3_12_12_2025_",
                label     = "H3K27me3")
)

dir.create("~/svd_analysis/results", showWarnings = FALSE, recursive = TRUE)

# Thresholds for this script (original / strict)
P_THRESH   <- 0.05
LFC_THRESH <- 1

wmh_sheet <- read_excel("~/svd_analysis/data/WMH_disease_groups.xlsx",
                        sheet = "WMH Disease Groups")

wmh_groups <- data.frame(
  MSSB_ID = wmh_sheet$`Sample ID`,
  Group   = factor(wmh_sheet$Group, levels = c("Low", "High"))
)

cat("WMH groups loaded:\n")
print(table(wmh_groups$Group))

summary_list <- list()

for (mark_name in names(marks)) {
  
  
  qc <- read_excel("~/svd_analysis/data/QC_decision_sheet.xlsx",
                   sheet = m$qc_sheet)
  qc$R_number <- gsub(m$file_gsub, "", qc$FILE)
  
  groups <- merge(wmh_groups, qc[, c("ID", "R_number")],
                  by.x = "MSSB_ID", by.y = "ID")
  
  cat("Samples with group assignment:", nrow(groups), "\n")
  
  fc     <- readRDS(m$rds)
  counts <- fc$counts
  colnames(counts) <- gsub(m$name_gsub, "", colnames(counts))
  
  cat("Peaks loaded:", nrow(counts), "; samples:", ncol(counts), "\n")
  
  group <- groups$Group[match(colnames(counts), groups$R_number)]
  stopifnot(!any(is.na(group)))
  
  cat("Group sizes: Low =", sum(group == "Low"),
      ", High =", sum(group == "High"), "\n")

  
  dgList <- DGEList(
    counts = counts,
    genes  = rownames(counts),
    group  = group
  )
  
  # Keep peaks with CPM > 1 in at least 2 samples
  countsPerMillion <- cpm(dgList)
  countCheck       <- countsPerMillion > 1
  keep             <- which(rowSums(countCheck) >= 2)
  dgList           <- dgList[keep, ]
  cat("Peaks after filtering:", nrow(dgList), "\n")
  
  dgList <- calcNormFactors(dgList, method = "TMM")
  cat("TMM norm factors:\n")
  print(round(dgList$samples$norm.factors, 3))
  
  design <- model.matrix(~ 0 + group, data = dgList$samples)
  colnames(design) <- gsub("group", "", colnames(design))
  cat("Design matrix columns:", paste(colnames(design), collapse = ", "), "\n")
  
  dgList <- estimateDisp(dgList, design, robust = TRUE)
  cat("Common dispersion:", round(dgList$common.dispersion, 4), "\n")
  
  fit <- glmQLFit(dgList, design, robust = TRUE)
  
  contrast <- makeContrasts(High - Low, levels = design)
  qlf      <- glmQLFTest(fit, contrast = contrast)
  
  results <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  cat("Results extracted:", nrow(results), "peaks\n")
  

  
  # Full results table (every peak with its stats)
  out_path <- paste0("~/svd_analysis/results/", mark_name,
                     "_edgeR_WMH_High_vs_Low_", SUFFIX, ".csv")
  fwrite(results, file = out_path, sep = "\t")
  cat("Results saved:", out_path, "\n")
  
  # FDR-based candidate list (the actual hits used in this script)
  candidates_fdr <- results[
    results$FDR < P_THRESH & abs(results$logFC) > LFC_THRESH,
  ]
  cand_path <- paste0("~/svd_analysis/results/", mark_name,
                      "_candidates_FDR", P_THRESH, "_lfc", LFC_THRESH,
                      "_", SUFFIX, ".csv")
  fwrite(candidates_fdr, file = cand_path, sep = "\t")
  cat("Candidates (FDR<", P_THRESH, ", |logFC|>", LFC_THRESH, "):",
      nrow(candidates_fdr), "saved to", cand_path, "\n")
  
  # Normalised log-CPM matrix (for downstream visualisation / GSEA)
  logcpm      <- cpm(dgList, normalized.lib.sizes = TRUE, log = TRUE)
  logcpm_path <- paste0("~/svd_analysis/results/", mark_name,
                        "_logCPM_normalised_", SUFFIX, ".csv")
  fwrite(as.data.frame(logcpm), file = logcpm_path, sep = "\t")
  cat("Log-CPM matrix saved:", logcpm_path, "\n")
  
 
  
  n_total   <- nrow(results)
  n_fdr05   <- sum(results$FDR < 0.05 & abs(results$logFC) > LFC_THRESH)
  n_fdr10   <- sum(results$FDR < 0.10 & abs(results$logFC) > LFC_THRESH)
  n_nom_p05 <- sum(results$PValue < 0.05 & abs(results$logFC) > LFC_THRESH)
  n_up      <- sum(results$FDR < 0.05 & results$logFC >  LFC_THRESH)
  n_down    <- sum(results$FDR < 0.05 & results$logFC < -LFC_THRESH)
  med_fdr   <- round(median(results$FDR), 4)
  
  cat("\n--- Summary:", m$label, "---\n")
  cat("Total peaks tested:                     ", n_total,   "\n")
  cat("Significant FDR < 0.05, |logFC| > 1:    ", n_fdr05,   "\n")
  cat("Significant FDR < 0.10, |logFC| > 1:    ", n_fdr10,   "\n")
  cat("Nominal p < 0.05, |logFC| > 1:          ", n_nom_p05, "\n")
  cat("Up (higher in High WMH, logFC > 1):     ", n_up,      "\n")
  cat("Down (lower in High WMH, logFC < -1):   ", n_down,    "\n")
  cat("Median FDR across all peaks:            ", med_fdr,   "\n")
  
  summary_list[[mark_name]] <- data.frame(
    Mark              = m$label,
    Total_peaks       = n_total,
    Sig_FDR05_lfc1    = n_fdr05,
    Sig_FDR10_lfc1    = n_fdr10,
    Nominal_p05_lfc1  = n_nom_p05,
    Up_High_WMH       = n_up,
    Down_High_WMH     = n_down,
    Median_FDR        = med_fdr
  )
  

  
  grp_col <- ifelse(group == "High", "#D6604D", "#2166AC")
  
  # 1. MDS plot (top 500 most variable peaks — edgeR default)
  svg(file.path(svg_dir, paste0("11_MDS_", mark_name, "_", SUFFIX, ".svg")),
      width = 7, height = 6)
  plotMDS(dgList, col = grp_col, pch = 16, cex = 1.4,
          main = paste(m$label, "- MDS (top 500 variable peaks)"),
          xlab = "Leading logFC dim 1", ylab = "Leading logFC dim 2")
  legend("topright", legend = c("Low WMH", "High WMH"),
         col = c("#2166AC", "#D6604D"), pch = 16, pt.cex = 1.2, bty = "n")
  dev.off()
  png(file.path(png_dir, paste0("11_MDS_", mark_name, "_", SUFFIX, ".png")),
      width = 7, height = 6, units = "in", res = 300)
  plotMDS(dgList, col = grp_col, pch = 16, cex = 1.4,
          main = paste(m$label, "- MDS (top 500 variable peaks)"),
          xlab = "Leading logFC dim 1", ylab = "Leading logFC dim 2")
  legend("topright", legend = c("Low WMH", "High WMH"),
         col = c("#2166AC", "#D6604D"), pch = 16, pt.cex = 1.2, bty = "n")
  dev.off()
  
  # 2. BCV plot — biological coefficient of variation across peaks
  svg(file.path(svg_dir, paste0("11_BCV_", mark_name, "_", SUFFIX, ".svg")),
      width = 7, height = 5)
  plotBCV(dgList,
          main = paste(m$label, "- Biological CV (dispersion estimation)"))
  dev.off()
  png(file.path(png_dir, paste0("11_BCV_", mark_name, "_", SUFFIX, ".png")),
      width = 7, height = 5, units = "in", res = 300)
  plotBCV(dgList,
          main = paste(m$label, "- Biological CV (dispersion estimation)"))
  dev.off()
  
  # 3. QL dispersion plot
  svg(file.path(svg_dir, paste0("11_QLDisp_", mark_name, "_", SUFFIX, ".svg")),
      width = 7, height = 5)
  plotQLDisp(fit,
             main = paste(m$label, "- QL dispersions (EB-squeezed)"))
  dev.off()
  png(file.path(png_dir, paste0("11_QLDisp_", mark_name, "_", SUFFIX, ".png")),
      width = 7, height = 5, units = "in", res = 300)
  plotQLDisp(fit,
             main = paste(m$label, "- QL dispersions (EB-squeezed)"))
  dev.off()
  
  # 4. P-value histogram
  p_hist <- ggplot(results, aes(x = PValue)) +
    geom_histogram(bins = 50, fill = "#4393C3",
                   colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0.05, linetype = "dashed", colour = "black") +
    annotate("text", x = 0.07, y = Inf, label = "p = 0.05",
             hjust = 0, vjust = 1.5, size = 3.5, colour = "black") +
    labs(title    = paste(m$label, "- P-value distribution"),
         subtitle = "Flat = no signal; spike near 0 = differential peaks present",
         x = "Nominal p-value (glmQLFTest)", y = "Number of peaks") +
    theme_classic(base_size = 12)
  ggsave(file.path(svg_dir, paste0("11_pval_hist_", mark_name, "_", SUFFIX, ".svg")),
         p_hist, width = 7, height = 5)
  ggsave(file.path(png_dir, paste0("11_pval_hist_", mark_name, "_", SUFFIX, ".png")),
         p_hist, width = 7, height = 5, dpi = 300)
  
  # 5. Volcano plot — categorical colours by FDR + logFC
  results$Category <- "Not significant"
  results$Category[results$FDR < 0.05 & results$logFC >  1] <- "Up in High WMH (FDR<0.05)"
  results$Category[results$FDR < 0.05 & results$logFC < -1] <- "Down in High WMH (FDR<0.05)"
  results$Category[results$FDR >= 0.05 & results$PValue < 0.05 &
                     abs(results$logFC) > 1] <- "Nominal p<0.05 only"
  results$Category <- factor(results$Category, levels = c(
    "Not significant", "Nominal p<0.05 only",
    "Up in High WMH (FDR<0.05)", "Down in High WMH (FDR<0.05)"
  ))
  
  cat_cols <- c(
    "Not significant"              = "grey75",
    "Nominal p<0.05 only"          = "#FDB863",
    "Up in High WMH (FDR<0.05)"    = "#D6604D",
    "Down in High WMH (FDR<0.05)"  = "#2166AC"
  )
  
  fdr_line <- -log10(0.05)
  
  volcano <- ggplot(results, aes(x = logFC, y = -log10(FDR),
                                 colour = Category)) +
    geom_point(size = 0.7, alpha = 0.55) +
    scale_colour_manual(values = cat_cols, name = NULL) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = fdr_line, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    annotate("text", x = max(results$logFC, na.rm = TRUE),
             y = fdr_line + 0.05, label = "FDR = 0.05",
             hjust = 1, vjust = 0, size = 3, colour = "grey30") +
    labs(title    = paste(m$label, "- Volcano Plot"),
         subtitle = paste0("High vs Low WMH | n=", n_fdr05,
                           " peaks FDR<0.05 & |logFC|>1"),
         x = "log2 fold change (High / Low WMH)",
         y = "-log10(FDR)") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9))
  ggsave(file.path(svg_dir, paste0("11_volcano_", mark_name, "_", SUFFIX, ".svg")),
         volcano, width = 8, height = 6)
  ggsave(file.path(png_dir, paste0("11_volcano_", mark_name, "_", SUFFIX, ".png")),
         volcano, width = 8, height = 6, dpi = 300)
  
  # 6. MA plot
  ma <- ggplot(results, aes(x = logCPM, y = logFC, colour = Category)) +
    geom_point(size = 0.6, alpha = 0.45) +
    scale_colour_manual(values = cat_cols, name = NULL) +
    geom_hline(yintercept = c(-1, 0, 1),
               linetype = c("dashed", "solid", "dashed"),
               colour = "grey40", linewidth = 0.4) +
    labs(title    = paste(m$label, "- MA Plot"),
         subtitle = "Average expression vs fold change; expect symmetry around 0",
         x = "Average log2 CPM (normalised)",
         y = "log2 fold change (High / Low WMH)") +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9))
  ggsave(file.path(svg_dir, paste0("11_MA_", mark_name, "_", SUFFIX, ".svg")),
         ma, width = 8, height = 6)
  ggsave(file.path(png_dir, paste0("11_MA_", mark_name, "_", SUFFIX, ".png")),
         ma, width = 8, height = 6, dpi = 300)
  
  cat("Plots saved for", m$label, "\n")
  
} # end mark loop


summary_df <- do.call(rbind, summary_list)
print(summary_df, row.names = FALSE)

fwrite(summary_df,
       file = paste0("~/svd_analysis/results/edgeR_summary_all_marks_",
                     SUFFIX, ".csv"),
       sep  = "\t")
cat("\nSummary saved to ~/svd_analysis/results/edgeR_summary_all_marks_",
    SUFFIX, ".csv\n", sep = "")

# Summary bar chart
summary_long <- pivot_longer(
  summary_df,
  cols      = c(Sig_FDR05_lfc1, Sig_FDR10_lfc1, Nominal_p05_lfc1),
  names_to  = "Threshold",
  values_to = "N_peaks"
)
summary_long$Threshold <- factor(
  summary_long$Threshold,
  levels = c("Nominal_p05_lfc1", "Sig_FDR10_lfc1", "Sig_FDR05_lfc1"),
  labels = c("Nominal p<0.05", "FDR<0.10", "FDR<0.05")
)

bar <- ggplot(summary_long, aes(x = Mark, y = N_peaks, fill = Threshold)) +
  geom_col(position = "dodge", width = 0.65) +
  geom_text(aes(label = N_peaks), position = position_dodge(width = 0.65),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("Nominal p<0.05" = "#FDB863",
                               "FDR<0.10"       = "#74ADD1",
                               "FDR<0.05"       = "#2166AC"),
                    name = NULL) +
  labs(title    = "Differential peak accessibility: High vs Low WMH volume",
       subtitle = "edgeR glmQLFTest | all thresholds also require |logFC| > 1",
       x = NULL, y = "Number of peaks") +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(svg_dir, paste0("11_summary_bar_all_marks_", SUFFIX, ".svg")),
       bar, width = 7, height = 5)
ggsave(file.path(png_dir, paste0("11_summary_bar_all_marks_", SUFFIX, ".png")),
       bar, width = 7, height = 5, dpi = 300)
cat("Summary bar chart saved.\n")