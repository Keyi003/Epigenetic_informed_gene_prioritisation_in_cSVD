library(edgeR)
library(readxl)
library(data.table)
library(ggplot2)
library(tidyr)
library(ggthemes)



results_dir <- "~/svd_analysis/results/edgeR_with_covariates"
fig_dir     <- "~/svd_analysis/figures/edgeR_with_covariates"
svg_dir     <- file.path(fig_dir, "svg_versions")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(svg_dir,     showWarnings = FALSE, recursive = TRUE)



P_THRESH   <- 0.05
LFC_THRESH <- 1

mark_cols <- c(
  H3K27ac  = "#D6604D",
  H3K4me3  = "#2166AC",
  H3K27me3 = "#4DAC26"
)


base_theme <- function(size = 12) {
  theme_base(base_size = size, base_family = "Helvetica") +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = "black"),
      legend.background = element_rect(fill = "white", colour = NA),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(colour = "grey25",
                                      face = "italic", size = 10)
    )
}



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



wmh_sheet <- read_excel("~/svd_analysis/data/WMH_disease_groups.xlsx",
                        sheet = "WMH Disease Groups")
wmh_groups <- data.frame(
  MSSB_ID = wmh_sheet$`Sample ID`,
  Group   = factor(wmh_sheet$Group, levels = c("Low", "High"))
)
cat("WMH groups loaded:\n")
print(table(wmh_groups$Group))

# Age + Sex from the same metadata file the limma script reads
meta_all <- read_excel("~/svd_analysis/data/sample_metadata.xlsx",
                       sheet = "metadata")
cat("Clinical metadata rows:", nrow(meta_all), "\n")
cat("Age range:", round(min(meta_all$Age), 1),
    "-", round(max(meta_all$Age), 1), "years\n")
cat("Sex distribution:",
    paste(names(table(meta_all$Sex)), table(meta_all$Sex),
          sep = "=", collapse = ", "), "\n")


summary_list <- list()

for (mark_name in names(marks)) {
  
 
  # Sample-to-group mapping via the QC decision sheet
  qc <- read_excel("~/svd_analysis/data/QC_decision_sheet.xlsx",
                   sheet = m$qc_sheet)
  qc$R_number <- gsub(m$file_gsub, "", qc$FILE)
  
  groups <- merge(wmh_groups, qc[, c("ID", "R_number")],
                  by.x = "MSSB_ID", by.y = "ID")
  cat("Samples with group assignment:", nrow(groups), "\n")
  
  # Load counts and clean column names to R-numbers
  fc     <- readRDS(m$rds)
  counts <- fc$counts
  colnames(counts) <- gsub(m$name_gsub, "", colnames(counts))
  cat("Peaks loaded:", nrow(counts), "; samples:", ncol(counts), "\n")
  
  # Group factor (High/Low) aligned to count matrix column order
  group <- groups$Group[match(colnames(counts), groups$R_number)]
  stopifnot(!any(is.na(group)))
  
  # Age + Sex aligned to same column order
  age_sex <- meta_all[match(colnames(counts), meta_all$R_number),
                      c("Age", "Sex")]
  stopifnot(!any(is.na(age_sex$Age)),
            !any(is.na(age_sex$Sex)))
  
  cat("Group sizes: Low =", sum(group == "Low"),
      ", High =", sum(group == "High"), "\n")
  cat("Mean age (Low):  ", round(mean(age_sex$Age[group == "Low"]), 1),
      "  | Mean age (High):", round(mean(age_sex$Age[group == "High"]), 1),
      "\n")
  cat("Sex (Low):       ",
      paste(names(table(age_sex$Sex[group == "Low"])),
            table(age_sex$Sex[group == "Low"]),
            sep = "=", collapse = ","),
      " | Sex (High):",
      paste(names(table(age_sex$Sex[group == "High"])),
            table(age_sex$Sex[group == "High"]),
            sep = "=", collapse = ","),
      "\n")
  

  
  dgList <- DGEList(
    counts = counts,
    genes  = rownames(counts),
    group  = group
  )
  
  # CPM > 1 in at least 2 samples
  countsPerMillion <- cpm(dgList)
  countCheck       <- countsPerMillion > 1
  keep             <- which(rowSums(countCheck) >= 2)
  dgList           <- dgList[keep, , keep.lib.sizes = FALSE]
  cat("Peaks after filtering:", nrow(dgList), "\n")
  
  # TMM normalisation
  dgList <- calcNormFactors(dgList, method = "TMM")
  cat("TMM norm factors range:",
      round(min(dgList$samples$norm.factors), 3), "-",
      round(max(dgList$samples$norm.factors), 3), "\n")
  
  # Design matrix WITH covariates.
  # ~ 0 + group gives one column per group (High, Low) which we
  # contrast directly. Adding + Age + Sex adds those as additional
  # columns whose coefficients absorb the variance they explain
  # without changing the High-vs-Low contrast definition.
  # Numeric Age, factor Sex (model.matrix handles dummy coding).
  Age <- as.numeric(age_sex$Age)
  Sex <- factor(age_sex$Sex)
  
  design <- model.matrix(~ 0 + group + Age + Sex,
                         data = dgList$samples)
  colnames(design) <- gsub("group", "", colnames(design))
  cat("Design columns:", paste(colnames(design), collapse = ", "), "\n")
  cat("Design rank:", qr(design)$rank, "/", ncol(design),
      "(equal = no collinearity)\n")
  
  # Dispersion estimation + QL fit
  dgList <- estimateDisp(dgList, design, robust = TRUE)
  cat("Common dispersion:", round(dgList$common.dispersion, 4), "\n")
  
  fit      <- glmQLFit(dgList, design, robust = TRUE)
  contrast <- makeContrasts(High - Low, levels = design)
  qlf      <- glmQLFTest(fit, contrast = contrast)
  
  results <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  cat("Results extracted:", nrow(results), "peaks\n")
  

  
  # Full results table
  out_path <- file.path(results_dir,
                        paste0(mark_name, "_edgeR_covar_full_results.csv"))
  fwrite(results, file = out_path, sep = "\t")
  cat("Full results saved:", out_path, "\n")
  
  # Candidates: nominal p < 0.05 AND logFC > 1
  candidates <- results[
    results$PValue < P_THRESH & abs(results$logFC) > LFC_THRESH,
  ]
  cand_path <- file.path(results_dir,
                         paste0(mark_name, "_candidates_covar_p",
                                P_THRESH, "_lfc", LFC_THRESH, ".csv"))
  fwrite(candidates, file = cand_path, sep = "\t")
  cat("Candidates (p<", P_THRESH, ", |logFC|>", LFC_THRESH, "): ",
      nrow(candidates), " saved\n", sep = "")
  cat("  Up in High WMH:  ", sum(candidates$logFC > 0), "\n")
  cat("  Down in High WMH:", sum(candidates$logFC < 0), "\n")
  
  # Normalised log-CPM matrix
  logcpm      <- cpm(dgList, normalized.lib.sizes = TRUE, log = TRUE)
  logcpm_path <- file.path(results_dir,
                           paste0(mark_name, "_logCPM_covar.csv"))
  fwrite(as.data.frame(logcpm), file = logcpm_path, sep = "\t")

  
  n_total      <- nrow(results)
  n_candidates <- nrow(candidates)
  n_up         <- sum(candidates$logFC > 0)
  n_down       <- sum(candidates$logFC < 0)
  n_fdr05_lfc1 <- sum(results$FDR < 0.05 & abs(results$logFC) > LFC_THRESH)
  med_fdr      <- round(median(results$FDR), 4)
  min_fdr      <- signif(min(results$FDR), 3)
  min_p        <- signif(min(results$PValue), 3)
  med_p        <- round(median(results$PValue), 4)
  
 
  
  summary_list[[mark_name]] <- data.frame(
    Mark                = m$label,
    Total_peaks         = n_total,
    Candidates_p05_lfc1 = n_candidates,
    Up_High_WMH         = n_up,
    Down_High_WMH       = n_down,
    FDR05_lfc1_FYI      = n_fdr05_lfc1,
    Min_pvalue          = min_p,
    Median_pvalue       = med_p,
    Min_FDR             = min_fdr,
    Median_FDR          = med_fdr
  )
 
  
  mark_colour <- mark_cols[m$label]
  grp_col <- ifelse(group == "High", "#D6604D", "#2166AC")
  

  svg(file.path(svg_dir, paste0("11_MDS_", mark_name, "_covar.svg")),
      width = 7, height = 6, family = "Helvetica")
  par(family = "Helvetica")
  plotMDS(dgList, col = grp_col, pch = 16, cex = 1.4,
          main = paste(m$label, "- MDS (top 500 variable peaks, +Age+Sex)"),
          xlab = "Leading logFC dim 1", ylab = "Leading logFC dim 2")
  legend("topright", legend = c("Low WMH", "High WMH"),
         col = c("#2166AC", "#D6604D"), pch = 16, pt.cex = 1.2, bty = "n")
  dev.off()
  png(file.path(fig_dir, paste0("11_MDS_", mark_name, "_covar.png")),
      width = 7, height = 6, units = "in", res = 300, family = "Helvetica")
  par(family = "Helvetica")
  plotMDS(dgList, col = grp_col, pch = 16, cex = 1.4,
          main = paste(m$label, "- MDS (top 500 variable peaks, +Age+Sex)"),
          xlab = "Leading logFC dim 1", ylab = "Leading logFC dim 2")
  legend("topright", legend = c("Low WMH", "High WMH"),
         col = c("#2166AC", "#D6604D"), pch = 16, pt.cex = 1.2, bty = "n")
  dev.off()
  

  if (n_candidates >= 10) {
    cand_ids    <- candidates$genes
    dgList_cand <- dgList[rownames(dgList) %in% cand_ids, ]
    
    svg(file.path(svg_dir, paste0("11_MDS_", mark_name,
                                  "_candidates_covar.svg")),
        width = 7, height = 6, family = "Helvetica")
    par(family = "Helvetica")
    plotMDS(dgList_cand, top = nrow(dgList_cand),
            col = grp_col, pch = 16, cex = 1.4,
            main = paste0(m$label,
                          " - MDS (candidate peaks only, n=",
                          nrow(dgList_cand), ", +Age+Sex)"),
            xlab = "Leading logFC dim 1",
            ylab = "Leading logFC dim 2")
    legend("topright", legend = c("Low WMH", "High WMH"),
           col = c("#2166AC", "#D6604D"), pch = 16,
           pt.cex = 1.2, bty = "n")
    dev.off()
    
    png(file.path(fig_dir, paste0("11_MDS_", mark_name,
                                  "_candidates_covar.png")),
        width = 7, height = 6, units = "in", res = 300,
        family = "Helvetica")
    par(family = "Helvetica")
    plotMDS(dgList_cand, top = nrow(dgList_cand),
            col = grp_col, pch = 16, cex = 1.4,
            main = paste0(m$label,
                          " - MDS (candidate peaks only, n=",
                          nrow(dgList_cand), ", +Age+Sex)"),
            xlab = "Leading logFC dim 1",
            ylab = "Leading logFC dim 2")
    legend("topright", legend = c("Low WMH", "High WMH"),
           col = c("#2166AC", "#D6604D"), pch = 16,
           pt.cex = 1.2, bty = "n")
    dev.off()
    cat("Candidate-peak MDS saved.\n")
  }
  

  svg(file.path(svg_dir, paste0("11_BCV_", mark_name, "_covar.svg")),
      width = 7, height = 5, family = "Helvetica")
  par(family = "Helvetica")
  plotBCV(dgList,
          main = paste(m$label, "- Biological CV (+Age+Sex)"))
  dev.off()
  png(file.path(fig_dir, paste0("11_BCV_", mark_name, "_covar.png")),
      width = 7, height = 5, units = "in", res = 300,
      family = "Helvetica")
  par(family = "Helvetica")
  plotBCV(dgList,
          main = paste(m$label, "- Biological CV (+Age+Sex)"))
  dev.off()
  

  svg(file.path(svg_dir, paste0("11_QLDisp_", mark_name, "_covar.svg")),
      width = 7, height = 5, family = "Helvetica")
  par(family = "Helvetica")
  plotQLDisp(fit,
             main = paste(m$label, "- QL dispersions (+Age+Sex)"))
  dev.off()
  png(file.path(fig_dir, paste0("11_QLDisp_", mark_name, "_covar.png")),
      width = 7, height = 5, units = "in", res = 300,
      family = "Helvetica")
  par(family = "Helvetica")
  plotQLDisp(fit,
             main = paste(m$label, "- QL dispersions (+Age+Sex)"))
  dev.off()

  p_hist <- ggplot(results, aes(x = PValue)) +
    geom_histogram(bins = 50, fill = mark_colour,
                   colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = P_THRESH, linetype = "dashed", colour = "black") +
    annotate("text", x = P_THRESH + 0.02, y = Inf,
             label = paste0("p = ", P_THRESH),
             hjust = 0, vjust = 1.5, size = 3.5, colour = "black") +
    labs(title    = paste(m$label, "(+Age+Sex) - P-value distribution"),
         subtitle = "Flat = no signal; spike near 0 = differential peaks present",
         x = "Nominal p-value (glmQLFTest)", y = "Number of peaks") +
    base_theme()
  ggsave(file.path(svg_dir, paste0("11_pval_hist_", mark_name, "_covar.svg")),
         p_hist, width = 7, height = 5)
  ggsave(file.path(fig_dir, paste0("11_pval_hist_", mark_name, "_covar.png")),
         p_hist, width = 7, height = 5, dpi = 300)
  
  # Y-axis is -log10(PValue) since nominal p defines candidates.
  # Three categories: Not significant, Up in High WMH, Down in High WMH.
  results$Category <- "Not significant"
  results$Category[results$PValue < P_THRESH &
                     results$logFC >  LFC_THRESH] <- "Up in High WMH"
  results$Category[results$PValue < P_THRESH &
                     results$logFC < -LFC_THRESH] <- "Down in High WMH"
  results$Category <- factor(results$Category,
                             levels = c("Not significant",
                                        "Up in High WMH",
                                        "Down in High WMH"))
  
  cat_cols <- c("Not significant"   = "grey80",
                "Up in High WMH"    = "#D6604D",
                "Down in High WMH"  = "#2166AC")
  
  pv_line <- -log10(P_THRESH)
  
  volcano <- ggplot(results, aes(x = logFC, y = -log10(PValue),
                                 colour = Category)) +
    geom_point(size = 0.7, alpha = 0.55) +
    scale_colour_manual(values = cat_cols, name = NULL) +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = pv_line, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    annotate("text", x = max(results$logFC, na.rm = TRUE),
             y = pv_line + 0.05,
             label = paste0("p = ", P_THRESH),
             hjust = 1, vjust = 0, size = 3, colour = "grey30") +
    labs(title    = paste(m$label, "(+Age+Sex) - Volcano Plot"),
         subtitle = paste0("High vs Low WMH | Candidates: p<", P_THRESH,
                           " & |logFC|>", LFC_THRESH,
                           " (n=", n_candidates, " peaks)"),
         x = "log2 fold change (High / Low WMH)",
         y = "-log10(nominal p)") +
    base_theme() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9))
  ggsave(file.path(svg_dir, paste0("11_volcano_", mark_name, "_covar.svg")),
         volcano, width = 8, height = 6)
  ggsave(file.path(fig_dir, paste0("11_volcano_", mark_name, "_covar.png")),
         volcano, width = 8, height = 6, dpi = 300)
  

  ma <- ggplot(results, aes(x = logCPM, y = logFC, colour = Category)) +
    geom_point(size = 0.6, alpha = 0.45) +
    scale_colour_manual(values = cat_cols, name = NULL) +
    geom_hline(yintercept = c(-LFC_THRESH, 0, LFC_THRESH),
               linetype = c("dashed", "solid", "dashed"),
               colour = "grey40", linewidth = 0.4) +
    labs(title    = paste(m$label, "(+Age+Sex) - MA Plot"),
         subtitle = paste0("Candidate bands at |logFC| = ", LFC_THRESH),
         x = "Average log2 CPM (normalised)",
         y = "log2 fold change (High / Low WMH)") +
    base_theme() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9))
  ggsave(file.path(svg_dir, paste0("11_MA_", mark_name, "_covar.svg")),
         ma, width = 8, height = 6)
  ggsave(file.path(fig_dir, paste0("11_MA_", mark_name, "_covar.png")),
         ma, width = 8, height = 6, dpi = 300)
  
  cat("Plots saved for", m$label, "\n")
  
} # end mark loop




summary_df <- do.call(rbind, summary_list)
print(summary_df, row.names = FALSE)

fwrite(summary_df,
       file = file.path(results_dir, "edgeR_covar_summary_all_marks.csv"),
       sep  = "\t")
cat("\nSummary saved.\n")



summary_long <- pivot_longer(
  summary_df,
  cols      = c(Up_High_WMH, Down_High_WMH),
  names_to  = "Direction",
  values_to = "N_peaks"
)
summary_long$Direction <- factor(
  summary_long$Direction,
  levels = c("Up_High_WMH", "Down_High_WMH"),
  labels = c("Up in High WMH", "Down in High WMH")
)

bar <- ggplot(summary_long, aes(x = Mark, y = N_peaks, fill = Direction)) +
  geom_col(position = "dodge", width = 0.65) +
  geom_text(aes(label = N_peaks), position = position_dodge(width = 0.65),
            vjust = -0.4, size = 3.5, family = "Helvetica") +
  scale_fill_manual(values = c("Up in High WMH"   = "#D6604D",
                               "Down in High WMH" = "#2166AC"),
                    name = NULL) +
  labs(title    = "Candidate peaks (+Age+Sex): High vs Low WMH",
       subtitle = paste0("edgeR glmQLFTest | nominal p<", P_THRESH,
                         " & |log2FC|>", LFC_THRESH),
       x = NULL, y = "Number of candidate peaks") +
  base_theme(size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(svg_dir, "11_summary_bar_covar.svg"),
       bar, width = 7, height = 5)
ggsave(file.path(fig_dir, "11_summary_bar_covar.png"),
       bar, width = 7, height = 5, dpi = 300)
cat("Summary bar chart saved.\n")