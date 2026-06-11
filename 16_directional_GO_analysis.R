library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(GenomicRanges)
library(data.table)
library(ggplot2)
library(ggthemes)




results_dir_out <- "~/svd_analysis/results/GO_enrichment_directional"
fig_dir         <- "~/svd_analysis/figures/GO_enrichment_directional"
svg_dir         <- file.path(fig_dir, "svg_versions")

dir.create(results_dir_out, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,         showWarnings = FALSE, recursive = TRUE)
dir.create(svg_dir,         showWarnings = FALSE, recursive = TRUE)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

mark_cols <- c(
  H3K27ac  = "#D6604D",
  H3K4me3  = "#2166AC",
  H3K27me3 = "#4DAC26"
)

# Direction colours used in the merged dotplot.
# Up-in-High-WMH = red (chromatin/gene activity gain with disease)
# Down-in-High-WMH = blue (chromatin/gene activity loss with disease)
dir_cols <- c(UP = "#B2182B", DOWN = "#2166AC")

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
  if (nrow(bed) == 0) return(list(bed = bed, genes = character(0)))
  
  gr <- GRanges(seqnames = bed$chr,
                ranges   = IRanges(start = bed$start, end = bed$end),
                name     = bed$name)
  
  peak_anno <- annotatePeak(gr,
                            tssRegion  = c(-2000, 500),
                            TxDb       = txdb,
                            level      = "gene",
                            annoDb     = "org.Hs.eg.db",
                            verbose    = FALSE)
  anno_df <- as.data.frame(peak_anno@anno)
  # Use SYMBOL; drop NA/blank
  genes <- unique(anno_df$SYMBOL[!is.na(anno_df$SYMBOL) & anno_df$SYMBOL != ""])
  list(bed = bed, genes = genes)
}



run_enrichGO <- function(fg_genes, bg_genes, label = "") {
  if (length(fg_genes) < 5) {
    cat("  ", label, ": <5 foreground genes -- skipping enrichGO\n")
    return(NULL)
  }
  tryCatch(
    enrichGO(gene          = fg_genes,
             universe      = bg_genes,
             OrgDb         = org.Hs.eg.db,
             keyType       = "SYMBOL",
             ont           = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff  = 1,    # keep all terms; threshold at plot
             qvalueCutoff  = 1,
             readable      = FALSE),
    error = function(e) {
      cat("  ", label, ": enrichGO error --", conditionMessage(e), "\n")
      NULL
    }
  )
}




dedup_terms <- function(ego_df) {
  ego_dt <- as.data.table(ego_df)
  if (nrow(ego_dt) == 0) return(ego_dt)
  ego_dt[, gene_key := vapply(strsplit(geneID, "/"),
                              function(x) paste(sort(x), collapse = ","),
                              character(1))]
  keep_idx <- ego_dt[, .I[which.max(p.adjust)], by = gene_key]$V1
  dedup    <- ego_dt[keep_idx][order(pvalue)]
  dedup[, gene_key := NULL]
  return(dedup)
}



make_dotplot <- function(dedup_df, mark_label, model_label, direction,
                         fg_n, bg_n, peak_n, mark_colour,
                         out_png, out_svg, top_n = 15) {
  if (nrow(dedup_df) == 0) {
    cat("  No terms to plot for", direction, "\n")
    return(invisible(NULL))
  }
  plot_n  <- min(top_n, nrow(dedup_df))
  plot_df <- as.data.frame(dedup_df[seq_len(plot_n)])
  
  plot_df$term_short <- ifelse(
    nchar(plot_df$Description) > 55,
    paste0(substr(plot_df$Description, 1, 52), "..."),
    plot_df$Description
  )
  plot_df$term_short <- factor(plot_df$term_short,
                               levels = rev(plot_df$term_short))
  plot_df$GeneCount <- as.integer(sub("/.*", "", plot_df$GeneRatio))
  
  p_005 <- -log10(0.05)
  p_001 <- -log10(0.01)
  
  n_q05 <- sum(dedup_df$p.adjust < 0.05)
  status_line <- if (n_q05 == 0) {
    paste0(peak_n, " ", direction, " peaks \u2192 ", fg_n, " genes (",
           bg_n, " bg). No q<0.05; ranked by raw p.")
  } else {
    paste0(n_q05, " terms at q<0.05 (", peak_n, " ", direction,
           " peaks, ", fg_n, " genes, ", bg_n, " bg).")
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
    labs(title    = paste0(mark_label, " \u2014 GO:BP enrichment (",
                           model_label, ", ", direction, ")"),
         subtitle = status_line,
         x = expression(-log[10]("nominal p-value")),
         y = NULL) +
    base_theme() +
    theme(axis.text.y     = element_text(size = 9),
          legend.position = "right")
  
  ggsave(out_png, p, width = 9, height = 6, dpi = 300)
  ggsave(out_svg, p, width = 9, height = 6)
}




make_merged_dotplot <- function(up_df, down_df,
                                mark_label, model_label,
                                up_peaks, down_peaks,
                                up_genes, down_genes, bg_n,
                                out_png, out_svg, top_n = 10) {
  build_subset <- function(df, direction) {
    if (nrow(df) == 0) return(NULL)
    n  <- min(top_n, nrow(df))
    sub <- as.data.frame(df[seq_len(n)])
    sub$Direction  <- direction
    sub$term_short <- ifelse(
      nchar(sub$Description) > 50,
      paste0(substr(sub$Description, 1, 47), "..."),
      sub$Description
    )
    sub$GeneCount <- as.integer(sub("/.*", "", sub$GeneRatio))
    sub
  }
  up_sub   <- build_subset(up_df,   "UP")
  down_sub <- build_subset(down_df, "DOWN")
  plot_df  <- rbind(up_sub, down_sub)
  if (is.null(plot_df) || nrow(plot_df) == 0) return(invisible(NULL))
  
  # Order terms by -log10(p) within direction so the dotplot reads
  # cleanly top-to-bottom for each colour block.
  plot_df <- plot_df[order(plot_df$Direction,
                           -(-log10(plot_df$pvalue))), ]
  plot_df$term_short <- factor(plot_df$term_short,
                               levels = rev(unique(plot_df$term_short)))
  
  p_005 <- -log10(0.05)
  p_001 <- -log10(0.01)
  
  status_line <- paste0("UP: ", up_peaks, " peaks \u2192 ", up_genes,
                        " genes  |  DOWN: ", down_peaks, " peaks \u2192 ",
                        down_genes, " genes  |  bg ", bg_n)
  
  p <- ggplot(plot_df,
              aes(x = -log10(pvalue), y = term_short,
                  size = GeneCount, fill = Direction)) +
    geom_point(shape = 21, colour = "grey20", stroke = 0.3) +
    geom_vline(xintercept = p_005, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = p_001, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    scale_size_continuous(name = "Gene count", range = c(2, 8)) +
    scale_fill_manual(name   = "Direction",
                      values = dir_cols,
                      labels = c(UP   = "UP in High WMH",
                                 DOWN = "DOWN in High WMH")) +
    labs(title    = paste0(mark_label,
                           " \u2014 GO:BP enrichment by direction (",
                           model_label, ")"),
         subtitle = status_line,
         x = expression(-log[10]("nominal p-value")),
         y = NULL) +
    base_theme() +
    theme(axis.text.y     = element_text(size = 9),
          legend.position = "right")
  
  ggsave(out_png, p, width = 10, height = 7, dpi = 300)
  ggsave(out_svg, p, width = 10, height = 7)
}



summary_list <- list()

for (mark_name in names(marks)) {
  

  bg <- annotate_peaks_to_genes(m$bed, peak_ids = NULL)
  cat("  Background genes:", length(bg$genes), "\n")
  
  for (model in c("nocovar", "covar")) {
    
    model_label <- if (model == "nocovar") "no covariates" else "+ Age + Sex"
    csv_path    <- if (model == "nocovar") m$csv_nocovar else m$csv_covar
    
    cat("\n--- Model:", model_label, "---\n")
    
    candidates <- fread(csv_path)
    cat("Total candidates:", nrow(candidates), "\n")
    

    up_ids   <- candidates$genes[candidates$logFC > 0]
    down_ids <- candidates$genes[candidates$logFC < 0]
    all_ids  <- candidates$genes
    cat("  UP in High WMH  :", length(up_ids), "peaks\n")
    cat("  DOWN in High WMH:", length(down_ids), "peaks\n")
    

    fg_all  <- annotate_peaks_to_genes(m$bed, peak_ids = all_ids)
    fg_up   <- annotate_peaks_to_genes(m$bed, peak_ids = up_ids)
    fg_down <- annotate_peaks_to_genes(m$bed, peak_ids = down_ids)
    cat("  Foreground genes -- ALL:", length(fg_all$genes),
        "| UP:", length(fg_up$genes),
        "| DOWN:", length(fg_down$genes), "\n")
    

      all  = list(genes = fg_all$genes,  peaks = length(all_ids)),
      up   = list(genes = fg_up$genes,   peaks = length(up_ids)),
      down = list(genes = fg_down$genes, peaks = length(down_ids))
    )
    
    dedup_store <- list()  # holds dedup tables per direction for the merged plot
    
    for (direction in names(direction_sets)) {
      d <- direction_sets[[direction]]
      cat("\n  Direction:", direction, "(", d$peaks, "peaks,",
          length(d$genes), "genes )\n")
      
      if (length(d$genes) < 5) {
        cat("    <5 genes -- skipping.\n")
        dedup_store[[direction]] <- data.table()
        next
      }
      
      ego <- run_enrichGO(fg_genes = d$genes,
                          bg_genes = bg$genes,
                          label    = paste(m$label, model_label, direction))
      if (is.null(ego)) {
        dedup_store[[direction]] <- data.table()
        next
      }
      
      res <- as.data.frame(ego)
      cat("    Terms tested:", nrow(res), "\n")
      cat("    Terms at q<0.05:", sum(res$p.adjust < 0.05), "\n")
      
      full_path  <- file.path(results_dir_out,
                              paste0(mark_name, "_", model, "_", direction,
                                     "_GO_BP_full.csv"))
      fwrite(res, file = full_path, sep = "\t")
      
      dedup <- dedup_terms(res)
      dedup_path <- file.path(results_dir_out,
                              paste0(mark_name, "_", model, "_", direction,
                                     "_GO_BP_dedup.csv"))
      fwrite(dedup, file = dedup_path, sep = "\t")
      cat("    Unique gene-set signatures:", nrow(dedup),
          "(removed", nrow(res) - nrow(dedup), ")\n")
      dedup_store[[direction]] <- dedup
      
      # Single-direction dotplot.
      make_dotplot(
        dedup_df    = dedup,
        mark_label  = m$label,
        model_label = model_label,
        direction   = direction,
        fg_n        = length(d$genes),
        bg_n        = length(bg$genes),
        peak_n      = d$peaks,
        mark_colour = mark_cols[m$label],
        out_png = file.path(fig_dir,
                            paste0("15_GO_", mark_name, "_", model,
                                   "_", direction, ".png")),
        out_svg = file.path(svg_dir,
                            paste0("15_GO_", mark_name, "_", model,
                                   "_", direction, ".svg"))
      )
      
      # Summary row -- top term by raw p
      summary_list[[paste0(mark_name, "_", model, "_", direction)]] <-
        data.frame(
          Mark       = m$label,
          Model      = model_label,
          Direction  = direction,
          Peaks      = d$peaks,
          FG_genes   = length(d$genes),
          BG_genes   = length(bg$genes),
          Terms      = nrow(res),
          Terms_dedup = nrow(dedup),
          Sig_q05    = sum(res$p.adjust < 0.05),
          Top_term   = if (nrow(dedup) > 0) dedup$Description[1] else NA_character_,
          Top_p      = if (nrow(dedup) > 0) signif(dedup$pvalue[1],   3) else NA_real_,
          Top_q      = if (nrow(dedup) > 0) signif(dedup$p.adjust[1], 3) else NA_real_
        )
    }  # end direction loop
    
    # Merged UP+DOWN dotplot -- the slide-ready figure.
    make_merged_dotplot(
      up_df       = dedup_store$up,
      down_df     = dedup_store$down,
      mark_label  = m$label,
      model_label = model_label,
      up_peaks    = length(up_ids),
      down_peaks  = length(down_ids),
      up_genes    = length(fg_up$genes),
      down_genes  = length(fg_down$genes),
      bg_n        = length(bg$genes),
      out_png     = file.path(fig_dir,
                              paste0("15_GO_", mark_name, "_", model,
                                     "_MERGED.png")),
      out_svg     = file.path(svg_dir,
                              paste0("15_GO_", mark_name, "_", model,
                                     "_MERGED.svg"))
    )
  }  # end model loop
}    # end mark loop




summary_df <- do.call(rbind, summary_list)
print(summary_df, row.names = FALSE)

fwrite(summary_df,
       file = file.path(results_dir_out, "GO_BP_directional_summary.csv"),
       sep  = "\t")
cat("\nDirectional summary saved.\n")