library(data.table)
library(Rsubread)

#Configuration
marks <- list(
  ac = list(
    bam_dir     = "~/pbmc/ac/",
    bam_pattern = "\\.sorted\\.bam$",
    saf_file    = "~/svd_analysis/counts/H3K27ac_consensus_peaks.saf",
    out_rds     = "~/svd_analysis/counts/featurecounts_H3K27ac.rds",
    out_csv     = "~/svd_analysis/counts/counts_H3K27ac.csv",
    name_gsub   = "unique_pbmc_ac_12_12_2025_|\\.sorted\\.bam"
  ),
  k4me3 = list(
    bam_dir     = "~/pbmc/k4me3/",
    bam_pattern = "\\.sorted\\.bam$",
    saf_file    = "~/svd_analysis/counts/H3K4me3_consensus_peaks.saf",
    out_rds     = "~/svd_analysis/counts/featurecounts_H3K4me3.rds",
    out_csv     = "~/svd_analysis/counts/counts_H3K4me3.csv",
    name_gsub   = "unique_pbmc_k4me3_12_12_2025_|\\.sorted\\.bam"
  ),
  k27me3 = list(
    bam_dir     = "~/pbmc/k27me3/",
    bam_pattern = "\\.sorted\\.bam$",
    saf_file    = "~/svd_analysis/counts/H3K27me3_consensus_peaks.saf",
    out_rds     = "~/svd_analysis/counts/featurecounts_H3K27me3.rds",
    out_csv     = "~/svd_analysis/counts/counts_H3K27me3.csv",
    name_gsub   = "unique_pbmc_k27me3_12_12_2025_|\\.sorted\\.bam"
  )
)

#Output directories
dir.create("~/svd_analysis/counts", showWarnings = FALSE, recursive = TRUE)
dir.create("~/svd_analysis/qc", showWarnings = FALSE, recursive = TRUE)

#Summary log
log_lines <- c(
  "======================================",
  "02_featurecounts.R - Summary",
  paste("Date:", Sys.time()),
  "======================================",
  "",
)

#Process each mark
for (mark_name in names(marks)) {
  
  m <- marks[[mark_name]]
  cat("\n=== Processing", mark_name, "===\n")
  
  #Find BAM files
  bam_files <- list.files(m$bam_dir, pattern = m$bam_pattern, full.names = TRUE)
  cat("Found", length(bam_files), "BAM files\n")
  
  #Load SAF as data frame (lab approach: fread)
  annot.ext <- fread(input = m$saf_file, sep = "\t")
  cat("Consensus peaks in SAF:", nrow(annot.ext), "\n")
  
  #Run featureCounts
  cat("Running featureCounts (this may take 15-30 min)...\n")
  fc <- featureCounts(
    files                  = bam_files,
    annot.ext              = annot.ext,
    isGTFAnnotationFile    = FALSE,
    isPairedEnd            = TRUE,
    countMultiMappingReads = FALSE,
    nthreads               = 4
  )
  
  #Save full featureCounts object
  saveRDS(fc, m$out_rds)
  cat("Saved RDS:", m$out_rds, "\n")
  
  #Clean column names and save count matrix
  counts <- fc$counts
  colnames(counts) <- gsub(m$name_gsub, "", colnames(counts))
  write.csv(counts, m$out_csv)
  cat("Saved CSV:", m$out_csv, "\n")
  
  #Assignment summary
  stat <- fc$stat
  assigned <- as.numeric(stat[stat$Status == "Assigned", -1])
  total <- colSums(stat[, -1])
  pct_assigned <- round(assigned / total * 100, 1)
  sample_names <- gsub(m$name_gsub, "", colnames(stat)[-1])
  
  cat("\nAssignment rates:\n")
  for (i in seq_along(sample_names)) {
    cat(sprintf("  %s: %s assigned / %s total (%.1f%%)\n",
                sample_names[i], format(assigned[i], big.mark = ","),
                format(total[i], big.mark = ","), pct_assigned[i]))
  }
  
  #Log summary
  log_lines <- c(log_lines, "",
                 paste0("--- ", mark_name, " ---"),
                 paste("BAM files:", length(bam_files)),
                 paste("Consensus peaks:", nrow(annot.ext)),
                 paste("Assignment rate range:", paste0(min(pct_assigned), "% - ", max(pct_assigned), "%")),
                 paste("Median assignment rate:", paste0(median(pct_assigned), "%")),
                 paste("Mean assignment rate:", paste0(round(mean(pct_assigned), 1), "%")),
                 "",
                 "Per-sample assignment rates:"
  )
  for (i in seq_along(sample_names)) {
    log_lines <- c(log_lines,
                   sprintf("  %s: %s / %s (%.1f%%)", sample_names[i],
                           format(assigned[i], big.mark = ","),
                           format(total[i], big.mark = ","), pct_assigned[i]))
  }
  
  log_lines <- c(log_lines, "", "Full assignment summary:",
                 capture.output(print(stat)))
}

#Save summary log
writeLines(log_lines, "~/svd_analysis/qc/02_featurecounts_summary.txt")
cat("\n=== Summary saved to ~/svd_analysis/qc/02_featurecounts_summary.txt ===\n")