library(data.table)
library(GenomicRanges)
library(IRanges)


marks <- list(
  ac = list(
    peak_dir = "~/pbmc/ac/unique_peaks/",
    pattern  = "narrowPeak$",
    out_saf  = "~/svd_analysis/counts/H3K27ac_consensus_peaks.saf",
    out_bed  = "~/svd_analysis/counts/H3K27ac_consensus_peaks.bed"
  ),
  k4me3 = list(
    peak_dir = "~/pbmc/k4me3/unique_peaks/",
    pattern  = "broadPeak$",
    out_saf  = "~/svd_analysis/counts/H3K4me3_consensus_peaks.saf",
    out_bed  = "~/svd_analysis/counts/H3K4me3_consensus_peaks.bed"
  ),
  k27me3 = list(
    peak_dir = "~/pbmc/k27me3/unique_peaks/",
    pattern  = "narrowPeak$",
    out_saf  = "~/svd_analysis/counts/H3K27me3_consensus_peaks.saf",
    out_bed  = "~/svd_analysis/counts/H3K27me3_consensus_peaks.bed"
  )
)


dir.create("~/svd_analysis/counts", showWarnings = FALSE, recursive = TRUE)
dir.create("~/svd_analysis/qc", showWarnings = FALSE, recursive = TRUE)

#Summary log
log_lines <- c(,
  "01_consensus_peaks.R - Summary",
  paste("Date:", Sys.time()),s,
  "",
)

#Process each mark
for (mark_name in names(marks)) {
  
  m <- marks[[mark_name]]
  cat("\n=== Processing", mark_name, "===\n")
  
  #Find all Peak files
  peak_paths <- list.files(path = m$peak_dir, pattern = m$pattern,
                           all.files = TRUE, full.names = TRUE)
  cat("Found", length(peak_paths), "peak files\n")
  
  #Read all peak files into one table (lab approach: rbindlist + fread)
  consensus_peaks <- rbindlist(l = lapply(peak_paths, fread))
  total_raw <- nrow(consensus_peaks)
  consensus_peaks <- unique(consensus_peaks)
  total_unique <- nrow(consensus_peaks)
  
  cat("Total peaks (raw):", total_raw, "\n")
  cat("Total peaks (deduplicated):", total_unique, "\n")
  
  #Filter: keep chr1-22 only (remove chrX, chrY, and contigs)
  removed_sex <- sum(consensus_peaks$V1 %in% c("chrX", "chrY"))
  removed_contigs <- sum(!(consensus_peaks$V1 %in% c(paste0("chr", 1:22), "chrX", "chrY")))
  consensus_peaks <- subset(consensus_peaks, V1 %in% paste0("chr", 1:22))
  
  cat("Removed chrX/chrY peaks:", removed_sex, "\n")
  cat("Removed contig peaks:", removed_contigs, "\n")
  cat("Peaks after filtering:", nrow(consensus_peaks), "\n")
  
  #Create GRanges object
  consensus_peaks_gr <- GRanges(
    seqnames = Rle(consensus_peaks$V1),
    ranges   = IRanges(start = consensus_peaks$V2, end = consensus_peaks$V3)
  )
  
  #Combine overlapping regions with at least 30bp overlap (lab standard)
  consensus_combined_gr <- GenomicRanges::reduce(
    x = consensus_peaks_gr, min.gapwidth = 31
  )
  
  consensus_combined <- as.data.frame(consensus_combined_gr)
  
  cat("Consensus peaks after reduce:", nrow(consensus_combined), "\n")
  
  #Build SAF format (GeneID, Chr, Start, End, Strand)
  saf <- data.frame(
    GeneID = paste0("PeakCount_", seq_len(nrow(consensus_combined))),
    Chr    = as.character(consensus_combined$seqnames),
    Start  = consensus_combined$start,
    End    = consensus_combined$end,
    Strand = "-"
  )
  
  #Save SAF
  fwrite(x = saf, file = m$out_saf, sep = "\t")
  cat("SAF saved:", m$out_saf, "\n")
  
  #Save BED
  bed <- data.frame(
    chr   = saf$Chr,
    start = saf$Start,
    end   = saf$End,
    name  = saf$GeneID
  )
  fwrite(x = bed, file = m$out_bed, sep = "\t", col.names = FALSE)
  cat("BED saved:", m$out_bed, "\n")
  
  #Chromosome distribution
  chr_table <- sort(table(saf$Chr))
  cat("\nChromosome distribution:\n")
  print(chr_table)
  
  #Peak width statistics
  widths <- saf$End - saf$Start
  cat("\nPeak width statistics:\n")
  cat("  Min:", min(widths), "bp\n")
  cat("  Median:", median(widths), "bp\n")
  cat("  Mean:", round(mean(widths)), "bp\n")
  cat("  Max:", max(widths), "bp\n")
  
  #Log summary
  log_lines <- c(log_lines, "",
                 paste0("--- ", mark_name, " ---"),
                 paste("Peak files found:", length(peak_paths)),
                 paste("Total raw peaks:", total_raw),
                 paste("After deduplication:", total_unique),
                 paste("chrX/chrY peaks removed:", removed_sex),
                 paste("Contig peaks removed:", removed_contigs),
                 paste("Peaks after filtering (chr1-22 only):", nrow(consensus_peaks)),
                 paste("Consensus peaks (after reduce):", nrow(saf)),
                 paste("Peak width - min:", min(widths), "bp"),
                 paste("Peak width - median:", median(widths), "bp"),
                 paste("Peak width - mean:", round(mean(widths)), "bp"),
                 paste("Peak width - max:", max(widths), "bp"),
                 paste("SAF file:", m$out_saf),
                 paste("BED file:", m$out_bed),
                 "",
                 "Chromosome distribution:",
                 capture.output(print(chr_table))
  )
}

#Save summary log
writeLines(log_lines, "~/svd_analysis/qc/01_consensus_peaks_summary.txt")
cat("\n=== Summary saved to ~/svd_analysis/qc/01_consensus_peaks_summary.txt ===\n")