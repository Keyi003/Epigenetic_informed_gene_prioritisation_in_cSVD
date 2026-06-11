library(ggplot2)
library(ggrepel)
library(edgeR)
library(readxl)

#Configuration
marks <- list(
  ac = list(
    rds        = "~/svd_analysis/counts/featurecounts_H3K27ac.rds",
    name_gsub  = "unique_pbmc_ac_12_12_2025_|\\.sorted\\.bam",
    qc_sheet   = "AC",
    file_gsub  = "pbmc_ac_12_12_2025_",
    mark_label = "H3K27ac",
    fig_prefix = "QC_pca_ac"
  ),
  k4me3 = list(
    rds        = "~/svd_analysis/counts/featurecounts_H3K4me3.rds",
    name_gsub  = "unique_pbmc_k4me3_12_12_2025_|\\.sorted\\.bam",
    qc_sheet   = "k4me3",
    file_gsub  = "pbmc_k4me3_12_12_2025_",
    mark_label = "H3K4me3",
    fig_prefix = "QC_pca_k4me3"
  ),
  k27me3 = list(
    rds        = "~/svd_analysis/counts/featurecounts_H3K27me3.rds",
    name_gsub  = "unique_pbmc_k27me3_12_12_2025_|\\.sorted\\.bam",
    qc_sheet   = "k27me3",
    file_gsub  = "pbmc_ke27me3_12_12_2025_",
    mark_label = "H3K27me3",
    fig_prefix = "QC_pca_k27me3"
  )
)

#Output directories
dir.create("~/svd_analysis/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("~/svd_analysis/figures/svg figure versions", showWarnings = FALSE, recursive = TRUE)
dir.create("~/svd_analysis/qc", showWarnings = FALSE, recursive = TRUE)
dir.create("~/svd_analysis/data", showWarnings = FALSE, recursive = TRUE)

#Load clinical metadata (once) 
clinical <- read_excel("~/svd_analysis/data/YWu_20250925.xlsx")

#Determine MoCA column to use
# f07 = 1 month (primary), f09 = 1 year (fallback)
f07_available <- sum(!is.na(clinical$f07_MoCA) & clinical$f07_MoCA != "")
f09_available <- sum(!is.na(clinical$f09_MoCA) & clinical$f09_MoCA != "")

cat("MoCA coverage check:\n")
cat("  f07_MoCA (1 month):", f07_available, "out of", nrow(clinical), "\n")
cat("  f09_MoCA (1 year):", f09_available, "out of", nrow(clinical), "\n")

if (f07_available >= f09_available) {
  moca_col <- "f07_MoCA"
  moca_label <- "MoCA (1 month)"
  cat("Using f07_MoCA (1 month) — better or equal coverage\n")
} else {
  moca_col <- "f09_MoCA"
  moca_label <- "MoCA (1 year)"
  cat("Using f09_MoCA (1 year) — better coverage\n")
}


#Process each mark
for (mark_name in names(marks)) {
  
  m <- marks[[mark_name]]
  cat("\n=== Processing", mark_name, "===\n")
  
  #Load counts
  fc <- readRDS(m$rds)
  counts <- fc$counts
  colnames(counts) <- gsub(m$name_gsub, "", colnames(counts))
  
  #Normalise to log-CPM
  logcpm <- cpm(counts, log = TRUE, prior.count = 1)
  
  #Build sample mapping
  qc_sheet <- read_excel("~/svd_analysis/data/QC_decision_sheet.xlsx", sheet = m$qc_sheet)
  qc_sheet$R_number <- gsub(m$file_gsub, "", qc_sheet$FILE)
  
  sample_map <- data.frame(
    R_number = qc_sheet$R_number,
    MSSB_ID  = qc_sheet$ID
  )
  
  #Merge with clinical data
  metadata <- merge(sample_map, clinical, by.x = "MSSB_ID", by.y = "ID")
  
  #Run PCA
  pca <- prcomp(t(logcpm), scale. = TRUE)
  var_explained <- summary(pca)$importance[2, 1:5] * 100
  
  cat("Variance explained (PC1-5):", round(var_explained, 1), "\n")
  
  #Build PCA data frame
  pca_df <- data.frame(
    R_number = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3]
  )
  
  pca_df <- merge(pca_df, metadata, by = "R_number")
  
  #Prepare clinical variables
  # f01 baseline clinical
  pca_df$Sex <- factor(pca_df$f01_Sex)
  pca_df$Age <- as.numeric(pca_df$f01_Age)
  
  #f03 baseline imaging
  pca_df$Fazekas <- factor(pca_df$f03_Fazekas)
  pca_df$WMH_raw <- as.numeric(pca_df$f03_WMH_ml)
  pca_df$ICV     <- as.numeric(pca_df$f03_ICV_ml)
  
  #Normalise WMH to ICV (Joanna's preference: 100 * WMH / ICV)
  pca_df$WMH_norm <- 100 * pca_df$WMH_raw / pca_df$ICV
  
  cat("WMH_norm range:", round(min(pca_df$WMH_norm, na.rm = TRUE), 2), "-",
      round(max(pca_df$WMH_norm, na.rm = TRUE), 2), "%\n")
  
  #MoCA (chosen column)
  pca_df$MoCA <- suppressWarnings(as.numeric(pca_df[[moca_col]]))
  n_moca_missing <- sum(is.na(pca_df$MoCA))
  n_moca_present <- sum(!is.na(pca_df$MoCA))
  cat("MoCA available:", n_moca_present, "/ missing:", n_moca_missing, "\n")
  
  #Plot 1: Basic
  svg(paste0("~/svd_analysis/figures/svg figure versions/", m$fig_prefix, "_basic.svg"), width = 8, height = 6)
  print(
    ggplot(pca_df, aes(x = PC1, y = PC2, label = R_number)) +
      geom_point(size = 3) +
      geom_text_repel(size = 3) +
      labs(title = paste(m$mark_label, "PCA"),
           x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
           y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
      theme_minimal()
  )
  dev.off()
  
  #Plot 2: Fazekas
  svg(paste0("~/svd_analysis/figures/svg figure versions/", m$fig_prefix, "_fazekas.svg"), width = 9, height = 6)
  print(
    ggplot(pca_df, aes(x = PC1, y = PC2, colour = Fazekas, label = R_number)) +
      geom_point(size = 4) +
      geom_text_repel(size = 3, show.legend = FALSE) +
      scale_colour_brewer(palette = "RdYlBu", direction = -1) +
      labs(title = paste(m$mark_label, "PCA - coloured by Fazekas score (f03)"),
           colour = "Fazekas\nscore",
           x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
           y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
      theme_minimal(base_size = 13)
  )
  dev.off()
  
  #Plot 3: Sex
  svg(paste0("~/svd_analysis/figures/svg figure versions/", m$fig_prefix, "_sex.svg"), width = 9, height = 6)
  print(
    ggplot(pca_df, aes(x = PC1, y = PC2, colour = Sex, label = R_number)) +
      geom_point(size = 4) +
      geom_text_repel(size = 3, show.legend = FALSE) +
      scale_colour_manual(values = c("F" = "#E41A1C", "M" = "#377EB8")) +
      labs(title = paste(m$mark_label, "PCA - coloured by sex (f01)"),
           colour = "Sex",
           x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
           y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
      theme_minimal(base_size = 13)
  )
  dev.off()
  
  #Plot 4: Age
  svg(paste0("~/svd_analysis/figures/svg figure versions/", m$fig_prefix, "_age.svg"), width = 9, height = 6)
  print(
    ggplot(pca_df, aes(x = PC1, y = PC2, colour = Age, label = R_number)) +
      geom_point(size = 4) +
      geom_text_repel(size = 3, show.legend = FALSE) +
      scale_colour_viridis_c(option = "plasma") +
      labs(title = paste(m$mark_label, "PCA - coloured by age (f01)"),
           colour = "Age",
           x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
           y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
      theme_minimal(base_size = 13)
  )
  dev.off()
  
  #Plot 5: WMH normalised
  svg(paste0("~/svd_analysis/figures/svg figure versions/", m$fig_prefix, "_wmh.svg"), width = 9, height = 6)
  print(
    ggplot(pca_df, aes(x = PC1, y = PC2, colour = WMH_norm, label = R_number)) +
      geom_point(size = 4) +
      geom_text_repel(size = 3, show.legend = FALSE) +
      scale_colour_viridis_c(option = "inferno") +
      labs(title = paste(m$mark_label, "PCA - coloured by WMH volume (% ICV, f03)"),
           colour = "WMH\n(% ICV)",
           x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
           y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
      theme_minimal(base_size = 13)
  )
  dev.off()
  
  #Plot 6: MoCA (only if >=10 samples have data)
  if (n_moca_present >= 10) {
    svg(paste0("~/svd_analysis/figures/svg figure versions/", m$fig_prefix, "_moca.svg"), width = 9, height = 6)
    print(
      ggplot(pca_df[!is.na(pca_df$MoCA), ],
             aes(x = PC1, y = PC2, colour = MoCA, label = R_number)) +
        geom_point(size = 4) +
        geom_text_repel(size = 3, show.legend = FALSE) +
        scale_colour_viridis_c(option = "viridis", direction = -1) +
        labs(title = paste(m$mark_label, "PCA - coloured by", moca_label),
             colour = "MoCA\nscore",
             x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
             y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
        theme_minimal(base_size = 13)
    )
    dev.off()
    cat("MoCA PCA plot saved\n")
  } else {
    cat("Skipping MoCA plot - only", n_moca_present, "samples with data (need >=10)\n")
  }
  
  #Flag outliers
  pc1_z <- abs(scale(pca_df$PC1))
  pc2_z <- abs(scale(pca_df$PC2))
  outliers_pc1 <- pca_df$R_number[pc1_z > 3]
  outliers_pc2 <- pca_df$R_number[pc2_z > 3]
  
  cat("PC1 outliers (>3 SD):", ifelse(length(outliers_pc1) == 0, "none",
                                      paste(outliers_pc1, collapse = ", ")), "\n")
  cat("PC2 outliers (>3 SD):", ifelse(length(outliers_pc2) == 0, "none",
                                      paste(outliers_pc2, collapse = ", ")), "\n")
  
  #Save PCA data
  write.csv(pca_df, paste0("~/svd_analysis/qc/pca_with_clinical_", mark_name, ".csv"),
            row.names = FALSE)

  log_lines <- c(log_lines, "",
                 paste0("--- ", mark_name, " ---"),
                 paste("Peaks:", nrow(logcpm)),
                 paste("Samples:", ncol(logcpm)),
                 "",
                 "Variance explained:",
                 paste(paste0("  PC", 1:5, ": ", round(var_explained, 1), "%"), collapse = "\n"),
                 "",
                 paste("PC1 outliers (>3 SD):", ifelse(length(outliers_pc1) == 0, "none",
                                                       paste(outliers_pc1, collapse = ", "))),
                 paste("PC2 outliers (>3 SD):", ifelse(length(outliers_pc2) == 0, "none",
                                                       paste(outliers_pc2, collapse = ", "))),
                 "",
                 "Clinical variable ranges:",
                 paste("  Fazekas levels:", paste(sort(unique(pca_df$Fazekas)), collapse = ", ")),
                 paste("  Sex: F =", sum(pca_df$Sex == "F"), ", M =", sum(pca_df$Sex == "M")),
                 paste("  Age:", round(min(pca_df$Age), 1), "-", round(max(pca_df$Age), 1)),
                 paste("  WMH raw:", round(min(pca_df$WMH_raw, na.rm = TRUE), 1), "-",
                       round(max(pca_df$WMH_raw, na.rm = TRUE), 1), "ml"),
                 paste("  WMH normalised:", round(min(pca_df$WMH_norm, na.rm = TRUE), 2), "-",
                       round(max(pca_df$WMH_norm, na.rm = TRUE), 2), "% ICV"),
                 paste("  MoCA:", n_moca_present, "available,", n_moca_missing, "missing")
  )
  
  if (n_moca_present > 0) {
    moca_vals <- pca_df$MoCA[!is.na(pca_df$MoCA)]
    log_lines <- c(log_lines,
                   paste("  MoCA range:", min(moca_vals), "-", max(moca_vals)),
                   paste("  MoCA median:", median(moca_vals))
    )
  }
}


write.csv(metadata, "~/svd_analysis/data/sample_metadata_merged.csv", row.names = FALSE)


writeLines(log_lines, "~/svd_analysis/qc/04_pca_summary.txt")
cat("\n=== Summary saved to ~/svd_analysis/qc/04_pca_summary.txt ===\n")