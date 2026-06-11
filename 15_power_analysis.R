suppressPackageStartupMessages({
  library(edgeR)
  library(readxl)
  library(RNASeqPower)   # Hart 2013 score-test formula
  library(ggplot2)
  library(data.table)
})


marks <- list(
  ac     = list(rds       = "~/svd_analysis/counts/featurecounts_H3K27ac.rds",
                name_gsub = "unique_pbmc_ac_12_12_2025_|\\.sorted\\.bam",
                qc_sheet  = "AC",
                file_gsub = "pbmc_ac_12_12_2025_",
                label     = "H3K27ac",
                colour    = "#f33c84"),
  k4me3  = list(rds       = "~/svd_analysis/counts/featurecounts_H3K4me3.rds",
                name_gsub = "unique_pbmc_k4me3_12_12_2025_|\\.sorted\\.bam",
                qc_sheet  = "k4me3",
                file_gsub = "pbmc_k4me3_12_12_2025_",
                label     = "H3K4me3",
                colour    = "#049093"),
  k27me3 = list(rds       = "~/svd_analysis/counts/featurecounts_H3K27me3.rds",
                name_gsub = "unique_pbmc_k27me3_12_12_2025_|\\.sorted\\.bam",
                qc_sheet  = "k27me3",
                file_gsub = "pbmc_ke27me3_12_12_2025_",
                label     = "H3K27me3",
                colour    = "#531c92")
)

results_dir <- "~/svd_analysis/results/power_analysis"
fig_dir     <- "~/svd_analysis/figures/power_analysis"
gp_dir      <- file.path(results_dir, "graphpad")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(gp_dir,      showWarnings = FALSE, recursive = TRUE)




# Discrete grid for the headline table 
fold_changes <- c(1.25, 1.5, 1.75, 2.0, 2.5, 3.0)
powers       <- c(0.80, 0.90)

# Dense grid for smooth GraphPad curves (1.10 to 3.00 in 0.05 steps)
fold_dense   <- seq(1.10, 3.00, by = 0.05)

# Reverse calculation: which fold change can we detect at fixed n?
n_per_group_grid <- c(10, 20, 30, 50, 75, 100, 150, 200)

# How many true differential peaks are we realistically aiming for?
# For peripheral epigenomic studies, 100 is a generous target.
target_positives <- 100



mark_inputs <- list()

for (mark_name in names(marks)) {
  
  m <- marks[[mark_name]]

  
  # Load counts (same as 11_*.R)
  fc <- readRDS(m$rds)
  counts <- fc$counts
  colnames(counts) <- gsub(m$name_gsub, "", colnames(counts))
  
  # Identical filtering to differential analysis
  dgList <- DGEList(counts = counts, genes = rownames(counts))
  cpm_mat <- cpm(dgList)
  keep    <- which(rowSums(cpm_mat > 1) >= 2)
  dgList  <- dgList[keep, , keep.lib.sizes = FALSE]
  
  # TMM normalisation
  dgList <- calcNormFactors(dgList, method = "TMM")
  
  
  design_null <- model.matrix(~ 1, data = dgList$samples)
  dgList <- estimateDisp(dgList, design_null, robust = TRUE)
  
  bcv_common  <- sqrt(dgList$common.dispersion)
  bcv_trended <- sqrt(median(dgList$trended.dispersion, na.rm = TRUE))
  bcv_tagwise <- sqrt(median(dgList$tagwise.dispersion, na.rm = TRUE))
  
  # Lambda: median expected count per peak. RNASeqPower wants the
  # baseline (lower group) mean count, so we use the median across peaks
  # of the median count across samples.
  median_lib <- median(dgList$samples$lib.size)
  median_cpm <- 2 ^ median(aveLogCPM(dgList))
  lambda     <- median_cpm * median_lib / 1e6
  

  mark_inputs[[mark_name]] <- list(
    label   = m$label,
    colour  = m$colour,
    npeaks  = nrow(dgList),
    bcv     = bcv_common,
    lambda  = lambda
  )
  
  # Save BCV diagnostic plot
  pdf(file.path(fig_dir, paste0("15_BCV_", mark_name, ".pdf")),
      width = 7, height = 5)
  plotBCV(dgList,
          main = paste(m$label, "- BCV used for power analysis"))
  abline(h = bcv_common, col = "red", lty = 2)
  legend("topright",
         legend = c(paste("Common BCV =", round(bcv_common, 3))),
         col = "red", lty = 2, bty = "n")
  dev.off()
}



inputs_df <- do.call(rbind, lapply(mark_inputs, function(x) {
  data.frame(Mark    = x$label,
             N_peaks = x$npeaks,
             BCV     = round(x$bcv, 4),
             Lambda  = round(x$lambda, 1))
}))
fwrite(inputs_df,
       file.path(results_dir, "15_power_inputs.csv"))
cat("\nInputs to power analysis:\n")
print(inputs_df, row.names = FALSE)




power_grid <- list()

for (mark_name in names(mark_inputs)) {
  mi  <- mark_inputs[[mark_name]]
  # Bonferroni-equivalent alpha for FDR<0.05 with ~100 expected positives
  alpha_fdr <- 0.05 * target_positives / mi$npeaks
  
  for (rho in fold_changes) {
    for (pw in powers) {
      # rnapower returns samples per group
      n_nom <- ceiling(rnapower(depth = mi$lambda,
                                cv    = mi$bcv,
                                effect = rho,
                                alpha  = 0.05,
                                power  = pw))
      n_fdr <- ceiling(rnapower(depth = mi$lambda,
                                cv    = mi$bcv,
                                effect = rho,
                                alpha  = alpha_fdr,
                                power  = pw))
      power_grid[[length(power_grid) + 1]] <- data.frame(
        Mark                    = mi$label,
        log2FC                  = round(log2(rho), 2),
        FoldChange              = rho,
        Power                   = pw,
        n_per_group_nominal_p05 = n_nom,
        n_per_group_FDR05       = n_fdr,
        n_total_nominal         = n_nom * 2,
        n_total_FDR05           = n_fdr * 2,
        alpha_FDR_eq            = signif(alpha_fdr, 3)
      )
    }
  }
}
power_df <- do.call(rbind, power_grid)
fwrite(power_df, file.path(results_dir, "15_power_curves_long.csv"))




dense_grid <- list()
for (mark_name in names(mark_inputs)) {
  mi  <- mark_inputs[[mark_name]]
  alpha_fdr <- 0.05 * target_positives / mi$npeaks
  for (rho in fold_dense) {
    for (pw in powers) {
      n_nom <- ceiling(rnapower(depth=mi$lambda, cv=mi$bcv,
                                effect=rho, alpha=0.05, power=pw))
      n_fdr <- ceiling(rnapower(depth=mi$lambda, cv=mi$bcv,
                                effect=rho, alpha=alpha_fdr, power=pw))
      dense_grid[[length(dense_grid) + 1]] <- data.frame(
        Mark = mi$label, FoldChange = rho, Power = pw,
        n_per_group_nominal_p05 = n_nom,
        n_per_group_FDR05       = n_fdr
      )
    }
  }
}
dense_df <- do.call(rbind, dense_grid)

# GraphPad XY wide: one X (FoldChange), three Y columns (marks)
for (pw in powers) {
  for (alpha_tier in c("nominal_p05", "FDR05")) {
    sub <- dense_df[dense_df$Power == pw, ]
    val_col <- paste0("n_per_group_", alpha_tier)
    wide <- reshape(sub[, c("FoldChange", "Mark", val_col)],
                    timevar = "Mark", idvar = "FoldChange",
                    direction = "wide")
    names(wide) <- gsub(paste0(val_col, "\\."), "", names(wide))
    fn <- file.path(gp_dir,
                    sprintf("graphpad_xy_%s_power%02d.csv",
                            alpha_tier, pw * 100))
    fwrite(wide, fn)
  }
}



mdfc_rows <- list()
for (mark_name in names(mark_inputs)) {
  mi <- mark_inputs[[mark_name]]
  alpha_fdr <- 0.05 * target_positives / mi$npeaks
  for (n in n_per_group_grid) {
    for (pw in powers) {
      # rnapower can solve for effect when n is supplied
      mdfc_nom <- tryCatch(
        rnapower(depth = mi$lambda, n = n, cv = mi$bcv,
                 alpha = 0.05, power = pw),
        error = function(e) NA
      )
      mdfc_fdr <- tryCatch(
        rnapower(depth = mi$lambda, n = n, cv = mi$bcv,
                 alpha = alpha_fdr, power = pw),
        error = function(e) NA
      )
      mdfc_rows[[length(mdfc_rows) + 1]] <- data.frame(
        Mark        = mi$label,
        n_per_group = n,
        Power       = pw,
        MDFC_nominal_p05  = round(mdfc_nom, 3),
        MDFC_FDR05_genome = round(mdfc_fdr, 3),
        log2_MDFC_nominal = round(log2(mdfc_nom), 3),
        log2_MDFC_FDR     = round(log2(mdfc_fdr), 3)
      )
    }
  }
}
mdfc_df <- do.call(rbind, mdfc_rows)
fwrite(mdfc_df, file.path(results_dir, "15_min_detectable_FC.csv"))

# GraphPad: X = n_per_group, Y = MDFC per mark, at 80% power
for (alpha_tier in c("nominal_p05", "FDR05")) {
  val_col <- paste0("MDFC_", ifelse(alpha_tier == "nominal_p05",
                                    "nominal_p05", "FDR05_genome"))
  sub <- mdfc_df[mdfc_df$Power == 0.80, c("n_per_group", "Mark", val_col)]
  wide <- reshape(sub, timevar = "Mark", idvar = "n_per_group",
                  direction = "wide")
  names(wide) <- gsub(paste0(val_col, "\\."), "", names(wide))
  fn <- file.path(gp_dir,
                  sprintf("graphpad_xy_MDFC_%s_power80.csv", alpha_tier))
  fwrite(wide, fn)
}



mark_cols <- setNames(
  vapply(mark_inputs, function(x) x$colour, character(1)),
  vapply(mark_inputs, function(x) x$label,  character(1))
)

p_dense <- ggplot(dense_df,
                  aes(x = FoldChange, y = n_per_group_nominal_p05,
                      colour = Mark, linetype = factor(Power))) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 10, linetype = "dotted", colour = "grey50") +
  annotate("text", x = 2.8, y = 11.5, label = "current pilot (n=10/group)",
           size = 3, colour = "grey50") +
  scale_y_log10(breaks = c(2, 5, 10, 20, 50, 100, 200, 500)) +
  scale_colour_manual(values = mark_cols) +
  labs(title    = "Sample size required to detect a fold change",
       subtitle = "Hart 2013 score test, nominal p < 0.05 per peak",
       x = "Fold change (High / Low WMH)",
       y = "Samples per group (log scale)",
       colour = "Mark", linetype = "Power") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")
ggsave(file.path(fig_dir, "15_power_curves_nominal.pdf"),
       p_dense, width = 8, height = 5)
ggsave(file.path(fig_dir, "15_power_curves_nominal.png"),
       p_dense, width = 8, height = 5, dpi = 300)

p_fdr <- p_dense %+% dense_df +
  aes(y = n_per_group_FDR05) +
  labs(subtitle = "FDR < 0.05 genome-wide (Bonferroni-equivalent)")
ggsave(file.path(fig_dir, "15_power_curves_FDR.pdf"),
       p_fdr, width = 8, height = 5)
ggsave(file.path(fig_dir, "15_power_curves_FDR.png"),
       p_fdr, width = 8, height = 5, dpi = 300)


