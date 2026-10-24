suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
out_dir <- Sys.getenv("OUT_DIR", "multistudy_validation_ibd")
s <- read.csv(file.path(out_dir, "multistudy_validation_summary.csv"))
col <- "#2a78d6"
th <- theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(),
  panel.grid.major = element_line(colour = "#ecebe8"), axis.title = element_text(colour = "#52514e"),
  axis.text = element_text(colour = "#52514e"), plot.background = element_rect(fill = "white", colour = NA))
pan <- function(exp, met, xlab, ylab, scale = 1, ref = NULL) {
  d <- s[s$experiment == exp & s$metric == met, ]; d$mean <- d$mean * scale; d$se <- d$se * scale
  g <- ggplot(d, aes(level, mean))
  if (!is.null(ref)) g <- g + ref
  g + geom_ribbon(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), fill = col, alpha = 0.18) +
    geom_line(colour = col, linewidth = 0.9) + geom_point(colour = col, size = 2.2) +
    labs(x = xlab, y = ylab) + th
}
f1 <- pan("A. effect.cor", "coef_cor", "Effect correlation setting", "Realized coefficient correlation",
          ref = geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "#9a9893")) +
      pan("A. effect.cor", "transfer_R2", "Effect correlation setting",
          expression("Prediction " * R^2 * " in held-out cohorts"))
ggsave(file.path(out_dir, "ibd_sweep_fig1_effect_correlation.png"), f1, width = 9, height = 4.2, dpi = 200, bg = "white")
noise <- mean(s$mean[s$experiment == "D. batch.loc.sd" & s$metric == "within_cor_noise"])
d2 <- s[s$experiment == "D. batch.loc.sd" & s$metric == "within_cor_change", ]
s <- rbind(s, transform(d2, metric = "within_cor_rel", mean = mean / noise, se = se / noise))
f2 <- pan("D. batch.loc.sd", "study_R2", "Batch effect strength", "Feature variance explained by cohort", scale = 100) +
        scale_y_continuous(labels = function(x) paste0(x, "%")) +
      pan("D. batch.loc.sd", "within_cor_rel", "Batch effect strength",
          "Change in within-cohort correlations\n(relative to replicate-to-replicate noise)",
          ref = geom_hline(yintercept = 1, linetype = "dotted", colour = "#9a9893")) +
        scale_y_continuous(limits = c(0, 1.05))
ggsave(file.path(out_dir, "ibd_sweep_fig2_batch_effects.png"), f2, width = 9, height = 4.2, dpi = 200, bg = "white")
