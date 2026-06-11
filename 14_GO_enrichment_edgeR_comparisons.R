library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GenomicRanges)
library(data.table)
library(ggplot2)
library(ggthemes)



results_dir_out <- "~/svd_analysis/results/GO_enrichment_comparison"
fig_dir         <- "~/svd_analysis/figures/GO_enrichment_comparison"
svg_dir         <- file.path(fig_dir, "svg_versions")

dir.create(results_dir_out, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,         showWarnings = FALSE, recursive = TRUE)
dir.create(svg_dir,         showWarnings = FALSE, recursive = TRUE)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Per-mark colour for plot accents
mark_cols <- c(
  H3K27ac  = "#D6604D",
  H3K4me3  = "#2166AC",
  H3K27me3 = "#4DAC26"
)


base_theme <- function(size = 11) {
  theme_base(base_size = size, base_family = "Helvetica") +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = "black"),
      legend.background = element_rect(fill = "white", colour = NA),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(colour = "grey25",
                                      face = "italic", size = 9)
    )
}



marks <- list(
  ac = list(
    label  = "H3K27ac",
    bed    = "~/svd_analysis/counts/H3K27ac_consensus_peaks.bed",
    csv_nocovar = "~/svd_analysis/results/edgeR_3rd_run/ac_candidates_p0.05_lfc1.csv",
    csv_covar   = "~/svd_analysis/results/edgeR_with_covariates/ac_candidates_covar_p0.05_lfc1.csv"
  ),
  k4me3 = list(
    label  = "H3K4me3",
    bed    = "~/svd_analysis/counts/H3K4me3_consensus_peaks.bed",
    csv_nocovar = "~/svd_analysis/results/edgeR_3rd_run/k4me3_candidates_p0.05_lfc1.csv",
    csv_covar   = "~/svd_analysis/results/edgeR_with_covariates/k4me3_candidates_covar_p0.05_lfc1.csv"
  ),
  k27me3 = list(
    label  = "H3K27me3",
    bed    = "~/svd_analysis/counts/H3K27me3_consensus_peaks.bed",
    csv_nocovar = "~/svd_analysis/results/edgeR_3rd_run/k27me3_candidates_p0.05_lfc1.csv",
    csv_covar   = "~/svd_analysis/results/edgeR_with_covariates/k27me3_candidates_covar_p0.05_lfc1.csv"
  )
)



annotate_peaks_to_genes <- function(bed_path, peak_ids = NULL) {
  bed <- fread(bed_path,
               col.names = c("chr", "start", "end", "name"))
  
  if (!is.null(peak_ids)) {
    bed <- bed[bed$name %in% peak_ids, ]
    if (nrow(bed) < length(peak_ids)) {
      warning("Mismatch: ", length(peak_ids) - nrow(bed),
              " peak IDs not in BED. Has consensus been regenerated?")
    }
  }
  
  if (nrow(bed) == 0) return(list(genes = character(0), anno = data.frame()))
  
  gr <- GRanges(seqnames = bed$chr,
                ranges   = IRanges(start = bed$start, end = bed$end),
                name     = bed$name)
  
  peak_anno <- annotatePeak(gr,
                            TxDb      = txdb,
                            level     = "gene",
                            annoDb    = "org.Hs.eg.db",
                            tssRegion = c(-2000, 500),
                            verbose   = FALSE)
  
  anno_df <- as.data.frame(peak_anno)
  genes   <- unique(na.omit(anno_df$SYMBOL))
  return(list(genes = genes, anno = anno_df))
}


run_enrichGO <- function(fg_genes, bg_genes, label) {
  if (length(fg_genes) < 5) {
    cat("  Too few foreground genes (", length(fg_genes),
        ") for", label, "-- skipping.\n")
    return(NULL)
  }
  ego <- enrichGO(
    gene          = fg_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    universe      = bg_genes,
    readable      = TRUE
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat("  No enrichment results for", label, "\n")
    return(NULL)
  }
  return(ego)
}


dedup_terms <- function(ego_df) {
  ego_dt <- as.data.table(ego_df)
  ego_dt[, gene_key := vapply(strsplit(geneID, "/"),
                              function(x) paste(sort(x), collapse = ","),
                              character(1))]
  keep_idx <- ego_dt[, .I[which.max(p.adjust)], by = gene_key]$V1
  dedup    <- ego_dt[keep_idx][order(pvalue)]
  dedup[, gene_key := NULL]
  return(dedup)
}



make_dotplot <- function(dedup_df, mark_label, model_label,
                         fg_n, bg_n, peak_n, mark_colour,
                         out_png, out_svg, top_n = 15) {
  
  plot_n  <- min(top_n, nrow(dedup_df))
  plot_df <- as.data.frame(dedup_df[seq_len(plot_n)])
  
  # Truncate long term names
  plot_df$term_short <- ifelse(
    nchar(plot_df$Description) > 55,
    paste0(substr(plot_df$Description, 1, 52), "..."),
    plot_df$Description
  )
  plot_df$term_short <- factor(plot_df$term_short,
                               levels = rev(plot_df$term_short))
  
  # GeneRatio is k/n -- parse numerator for sizing
  plot_df$GeneCount <- as.integer(sub("/.*", "", plot_df$GeneRatio))
  
  # Significance reference lines
  p_005 <- -log10(0.05)
  p_001 <- -log10(0.01)
  
  # Q status for subtitle
  n_q05 <- sum(dedup_df$p.adjust < 0.05)
  status_line <- if (n_q05 == 0) {
    paste0(peak_n, " candidate peaks \u2192 ", fg_n, " genes (",
           bg_n, " bg). No q<0.05; ranked by raw p.")
  } else {
    paste0(n_q05, " terms at q<0.05 (", peak_n, " peaks, ",
           fg_n, " genes, ", bg_n, " bg).")
  }
  
  p <- ggplot(plot_df,
              aes(x = -log10(pvalue), y = term_short,
                  size = GeneCount, fill = -log10(pvalue))) +
    geom_point(shape = 21, colour = "grey20", stroke = 0.3) +
    geom_vline(xintercept = p_005, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = p_001, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    annotate("text", x = p_005, y = 0.45,
             label = "p=0.05", hjust = -0.1, vjust = 0,
             size = 2.8, colour = "grey30", family = "Helvetica") +
    annotate("text", x = p_001, y = 0.45,
             label = "p=0.01", hjust = -0.1, vjust = 0,
             size = 2.8, colour = "grey30", family = "Helvetica") +
    scale_size_continuous(name = "Gene count", range = c(2, 8)) +
    scale_fill_gradient(name = expression(-log[10](p)),
                        low = "#FDB863", high = mark_colour) +
    labs(title    = paste0(mark_label, " \u2014 GO:BP enrichment (", model_label, ")"),
         subtitle = status_line,
         x = expression(-log[10]("nominal p-value")),
         y = NULL) +
    base_theme() +
    theme(axis.text.y     = element_text(size = 9),
          legend.position = "right")
  
  ggsave(out_png, p, width = 9, height = 6, dpi = 300)
  ggsave(out_svg, p, width = 9, height = 6)
}


summary_list <- list()

for (mark_name in names(marks)) {
  

  
  # Background only needs annotating once per mark (same BED for
  # both models). Slow step (~minutes for H3K27ac's 68k peaks).
  cat("Annotating background (full consensus peak set)...\n")
  bg <- annotate_peaks_to_genes(m$bed, peak_ids = NULL)
  cat("  Background genes:", length(bg$genes), "\n")
  
  for (model in c("nocovar", "covar")) {
    
    model_label <- if (model == "nocovar") "no covariates" else "+ Age + Sex"
    csv_path    <- if (model == "nocovar") m$csv_nocovar else m$csv_covar
    
    cat("\n--- Model:", model_label, "---\n")
    
    candidates <- fread(csv_path)
    peak_n     <- nrow(candidates)
    cat("Candidates:", peak_n, "\n")
    
    if (peak_n < 10) {
      cat("  Too few candidates -- skipping.\n")
      next
    }
    
    # Foreground: annotate candidate peaks to unique gene symbols.
    fg <- annotate_peaks_to_genes(m$bed, peak_ids = candidates$genes)
    cat("  Foreground genes:", length(fg$genes), "\n")
    
    # Enrichment
    ego <- run_enrichGO(fg_genes = fg$genes,
                        bg_genes = bg$genes,
                        label    = paste(m$label, model_label))
    if (is.null(ego)) next
    
    res <- as.data.frame(ego)
    cat("  Terms tested:", nrow(res), "\n")
    cat("  Terms at q<0.05:", sum(res$p.adjust < 0.05), "\n")
    
    # Save full results
    full_path <- file.path(results_dir_out,
                           paste0(mark_name, "_", model, "_GO_BP_full.csv"))
    fwrite(res, file = full_path, sep = "\t")
    
    # Deduplicate
    dedup <- dedup_terms(res)
    dedup_path <- file.path(results_dir_out,
                            paste0(mark_name, "_", model, "_GO_BP_dedup.csv"))
    fwrite(dedup, file = dedup_path, sep = "\t")
    cat("  Unique gene-set signatures:", nrow(dedup),
        "(removed", nrow(res) - nrow(dedup), "redundant)\n")
    
    # Dotplot
    make_dotplot(
      dedup_df    = dedup,
      mark_label  = m$label,
      model_label = model_label,
      fg_n        = length(fg$genes),
      bg_n        = length(bg$genes),
      peak_n      = peak_n,
      mark_colour = mark_cols[m$label],
      out_png = file.path(fig_dir,
                          paste0("12_GO_", mark_name, "_", model, ".png")),
      out_svg = file.path(svg_dir,
                          paste0("12_GO_", mark_name, "_", model, ".svg"))
    )
    cat("  Dotplot saved.\n")
    
    # Summary row -- top term by raw p after dedup
    summary_list[[paste0(mark_name, "_", model)]] <- data.frame(
      Mark       = m$label,
      Model      = model_label,
      Peaks      = peak_n,
      FG_genes   = length(fg$genes),
      BG_genes   = length(bg$genes),
      Terms      = nrow(res),
      Terms_dedup = nrow(dedup),
      Sig_q05    = sum(res$p.adjust < 0.05),
      Top_term   = if (nrow(dedup) > 0) dedup$Description[1] else NA_character_,
      Top_p      = if (nrow(dedup) > 0) signif(dedup$pvalue[1],   3) else NA_real_,
      Top_q      = if (nrow(dedup) > 0) signif(dedup$p.adjust[1], 3) else NA_real_
    )
  }
}




summary_df <- do.call(rbind, summary_list)
print(summary_df, row.names = FALSE)

fwrite(summary_df,
       file = file.path(results_dir_out, "GO_BP_summary_comparison.csv"),
       sep  = "\t")
cat("\nSummary saved.\n")