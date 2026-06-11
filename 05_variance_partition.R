library(reformulas)
library(variancePartition)
library(edgeR)
library(readxl)
library(ggplot2)

#Configuration
marks <- list(
  ac = list(
    rds        = "~/svd_analysis/counts/featurecounts_H3K27ac.rds",
    name_gsub  = "unique_pbmc_ac_12_12_2025_|\\.sorted\\.bam",
    qc_sheet   = "AC",
    file_gsub  = "pbmc_ac_12_12_2025_",
    mark_label = "H3K27ac",
    out_pdf    = "~/svd_analysis/figures/svg figure versions/QC_variance_partition_ac.svg",
    out_csv    = "~/svd_analysis/qc/variance_partition_H3K27ac.csv",
    out_rds    = "~/svd_analysis/qc/variance_partition_H3K27ac.rds"
  ),
  k4me3 = list(
    rds        = "~/svd_analysis/counts/featurecounts_H3K4me3.rds",
    name_gsub  = "unique_pbmc_k4me3_12_12_2025_|\\.sorted\\.bam",
    qc_sheet   = "k4me3",
    file_gsub  = "pbmc_k4me3_12_12_2025_",
    mark_label = "H3K4me3",
    out_pdf    = "~/svd_analysis/figures/svg figure versions/QC_variance_partition_k4me3.svg",
    out_csv    = "~/svd_analysis/qc/variance_partition_H3K4me3.csv",
    out_rds    = "~/svd_analysis/qc/variance_partition_H3K4me3.rds"
  ),
  k27me3 = list(
    rds        = "~/svd_analysis/counts/featurecounts_H3K27me3.rds",
    name_gsub  = "unique_pbmc_k27me3_12_12_2025_|\\.sorted\\.bam",
    qc_sheet   = "k27me3",
    file_gsub  = "pbmc_ke27me3_12_12_2025_",
    mark_label = "H3K27me3",
    out_pdf    = "~/svd_analysis/figures/svg figure versions/QC_variance_partition_k27me3.svg",
    out_csv    = "~/svd_analysis/qc/variance_partition_H3K27me3.csv",
    out_rds    = "~/svd_analysis/qc/variance_partition_H3K27me3.rds"
  )
)

#Output directories
dir.create("~/svd_analysis/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("~/svd_analysis/figures/svg figure versions", showWarnings = FALSE, recursive = TRUE)
dir.create("~/svd_analysis/qc", showWarnings = FALSE, recursive = TRUE)

#Load clinical metadata 
clinical <- read_excel("~/svd_analysis/data/YWu_20250925.xlsx")

#Determine MoCA column 
f07_available <- sum(!is.na(clinical$f07_MoCA) & clinical$f07_MoCA != "")
f09_available <- sum(!is.na(clinical$f09_MoCA) & clinical$f09_MoCA != "")

if (f07_available >= f09_available) {
  moca_col <- "f07_MoCA"
  moca_label <- "MoCA (1 month)"
} else {
  moca_col <- "f09_MoCA"
  moca_label <- "MoCA (1 year)"
}
cat("Using", moca_col, "for MoCA\n")



#Process each mark
for (mark_name in names(marks)) {
  
  m <- marks[[mark_name]]
  cat("\n=== Processing", mark_name, "===\n")
  
  #Load counts
  fc <- readRDS(m$rds)
  counts <- fc$counts
  colnames(counts) <- gsub(m$name_gsub, "", colnames(counts))
  
  #Build sample mapping
  qc_sheet <- read_excel("~/svd_analysis/data/QC_decision_sheet.xlsx", sheet = m$qc_sheet)
  qc_sheet$R_number <- gsub(m$file_gsub, "", qc_sheet$FILE)
  
  sample_map <- data.frame(
    R_number = qc_sheet$R_number,
    MSSB_ID  = qc_sheet$ID
  )
  
  #Merge with clinical data
  metadata <- merge(sample_map, clinical, by.x = "MSSB_ID", by.y = "ID")
  
  #Align metadata to count matrix
  metadata <- metadata[match(colnames(counts), metadata$R_number), ]
  rownames(metadata) <- metadata$R_number
  cat("Samples match:", all(colnames(counts) == rownames(metadata)), "\n")
  
  #Prepare variables
  metadata$Sex     <- factor(metadata$f01_Sex)
  metadata$Age     <- as.numeric(metadata$f01_Age)
  metadata$Fazekas <- factor(metadata$f03_Fazekas)
  
  #WMH normalised to ICV
  metadata$WMH_raw  <- as.numeric(metadata$f03_WMH_ml)
  metadata$ICV      <- as.numeric(metadata$f03_ICV_ml)
  metadata$WMH_norm <- 100 * metadata$WMH_raw / metadata$ICV
  
  #MoCA
  metadata$MoCA <- suppressWarnings(as.numeric(metadata[[moca_col]]))
  
  #Check missing values
  cat("Missing values:\n")
  cat("  Sex:", sum(is.na(metadata$Sex)), "\n")
  cat("  Age:", sum(is.na(metadata$Age)), "\n")
  cat("  Fazekas:", sum(is.na(metadata$Fazekas)), "\n")
  cat("  WMH_norm:", sum(is.na(metadata$WMH_norm)), "\n")
  cat("  MoCA:", sum(is.na(metadata$MoCA)), "\n")
  
  #Decide formula based on MoCA availability
  n_moca_complete <- sum(complete.cases(metadata[, c("Sex", "Age", "Fazekas", "WMH_norm", "MoCA")]))
  
  if (n_moca_complete >= 15) {
    cat("Including MoCA (", n_moca_complete, "complete samples)\n")
    vars_to_check <- c("Sex", "Age", "Fazekas", "WMH_norm", "MoCA")
    form <- ~ (1|Sex) + (1|Fazekas) + Age + WMH_norm + MoCA
    use_moca <- TRUE
  } else {
    cat("Excluding MoCA (only", n_moca_complete, "complete samples, need >=15)\n")
    vars_to_check <- c("Sex", "Age", "Fazekas", "WMH_norm")
    form <- ~ (1|Sex) + (1|Fazekas) + Age + WMH_norm
    use_moca <- FALSE
  }
  
  #Filter to complete cases
  complete <- complete.cases(metadata[, vars_to_check])
  n_dropped <- sum(!complete)
  
  if (!all(complete)) {
    dropped <- metadata$R_number[!complete]
    cat("Dropping", n_dropped, "samples with missing data:", paste(dropped, collapse = ", "), "\n")
    metadata <- metadata[complete, ]
    counts <- counts[, metadata$R_number]
  }
  
  cat("Samples used:", nrow(metadata), "\n")
  
  #Normalise with voom 
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~ 1, data = metadata)
  voom_data <- voom(dge, design)
  
  #Fit variance partition
  cat("Fitting variance partition model (formula:", deparse(form), ")...\n")
  vp <- fitExtractVarPartModel(voom_data, form, metadata)
  
  #Plot
  svg(m$out_pdf, width = 10, height = 6)
  plotVarPart(sortCols(vp))
  dev.off()
  cat("Plot saved:", m$out_pdf, "\n")
  
  #Summary statistics
  vp_means   <- colMeans(vp)
  vp_medians <- apply(vp, 2, median)
  
  cat("\nMean variance explained per variable:\n")
  print(round(sort(vp_means, decreasing = TRUE), 4))
  
  cat("\nMedian variance explained per variable:\n")
  print(round(sort(vp_medians, decreasing = TRUE), 4))
  
  #Top peaks per variable of interest
  top_fazekas <- head(sort(vp$Fazekas, decreasing = TRUE), 20)
  top_wmh     <- head(sort(vp$WMH_norm, decreasing = TRUE), 20)
  
  cat("\nTop 20 peaks by Fazekas variance:\n")
  print(round(top_fazekas, 3))
  
  cat("\nTop 20 peaks by WMH_norm variance:\n")
  print(round(top_wmh, 3))
  
  #Save
  write.csv(as.data.frame(vp), m$out_csv, row.names = TRUE)
  saveRDS(vp, m$out_rds)
  
  #Log summary
  log_lines <- c(log_lines, "",
                 paste0("--- ", mark_name, " ---"),
                 paste("Formula:", deparse(form)),
                 paste("MoCA included:", use_moca),
                 paste("Samples used:", nrow(metadata)),
                 paste("Samples dropped:", n_dropped),
                 paste("Peaks:", nrow(vp)),
                 "",
                 "Mean variance explained:",
                 capture.output(print(round(sort(vp_means, decreasing = TRUE), 4))),
                 "",
                 "Median variance explained:",
                 capture.output(print(round(sort(vp_medians, decreasing = TRUE), 4))),
                 "",
                 "Top 20 peaks by Fazekas variance:",
                 capture.output(print(round(top_fazekas, 3))),
                 "",
                 "Top 20 peaks by WMH_norm variance:",
                 capture.output(print(round(top_wmh, 3)))
  )
  
  if (use_moca) {
    top_moca <- head(sort(vp$MoCA, decreasing = TRUE), 20)
    cat("\nTop 20 peaks by MoCA variance:\n")
    print(round(top_moca, 3))
    
    log_lines <- c(log_lines, "",
                   "Top 20 peaks by MoCA variance:",
                   capture.output(print(round(top_moca, 3)))
    )
  }
}

#Save summary logs
writeLines(log_lines, "~/svd_analysis/qc/05_variance_partition_summary.txt")
cat("\n=== Summary saved to ~/svd_analysis/qc/05_variance_partition_summary.txt ===\n")