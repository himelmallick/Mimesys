suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
r <- readRDS(file.path(Sys.getenv("OUT_DIR", "ibd_validation"), "ibd_validation_results.rds"))
coh <- c("PRISM", "Validation", "LloydPrice")
grp <- c(KS_taxa = "Marginal distributions", KS_metab = "Marginal distributions",
  zero_taxa = "Zero fractions (taxa)", cor_within_taxa = "Correlation structure",
  cor_within_metab = "Correlation structure", cor_cross_view = "Correlation structure")
agg <- function(d, kind) { d$group <- grp[d$metric]; a <- aggregate(value ~ rep + real + group, d, mean); a$kind <- kind; a }
f <- r$fid
d <- rbind(agg(f[f$real == f$sim, ], "Simulated from its own template"),
           agg(f[f$real != f$sim, ], "Simulated from other cohorts' templates"),
           agg(r$floor, "Sampling noise floor"))
s <- aggregate(value ~ real + group + kind, d, function(v) c(m = mean(v), se = sd(v) / sqrt(length(v))))
s <- do.call(data.frame, s); names(s)[4:5] <- c("mean", "se")
s$real <- factor(s$real, coh); s$group <- factor(s$group, unique(grp))
s$kind <- factor(s$kind, c("Simulated from its own template", "Simulated from other cohorts' templates", "Sampling noise floor"))
cols <- c("Simulated from its own template" = "#2a78d6", "Simulated from other cohorts' templates" = "#eb6834", "Sampling noise floor" = "#8a8984")
th <- theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(colour = "#ecebe8"), axis.title = element_text(colour = "#52514e"),
  axis.text = element_text(colour = "#52514e"), strip.text = element_text(face = "bold", size = 11),
  legend.position = "top", legend.title = element_blank(), plot.background = element_rect(fill = "white", colour = NA))
pd <- position_dodge(width = 0.55)
g1 <- ggplot(s, aes(real, mean, colour = kind, shape = kind)) +
  geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), width = 0, position = pd, linewidth = 0.8) +
  geom_point(size = 3, position = pd) +
  scale_colour_manual(values = cols) + scale_shape_manual(values = c(16, 16, 4)) +
  facet_wrap(~group, scales = "free_y") + expand_limits(y = 0) +
  labs(x = NULL, y = "Distance from the real cohort") + th + guides(colour = guide_legend(nrow = 1))
ggsave(file.path(Sys.getenv("OUT_DIR", "ibd_validation"), "ibd_fig1_fidelity.png"), g1, width = 11, height = 4.3, dpi = 200, bg = "white")

## Figure 2: between-cohort differences and cohort identity
ex <- r$smd_example; ex$pair <- factor(ex$pair, unique(ex$pair))
ex$view <- factor(ex$view, c("Microbiome", "Metabolome"))
g2a <- ggplot(ex, aes(real, sim, colour = view)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "#9a9893") +
  geom_point(size = 1.5, alpha = 0.75) +
  scale_colour_manual(values = c(Microbiome = "#1baf7a", Metabolome = "#6250d6")) +
  facet_wrap(~pair, scales = "free", nrow = 1) +
  labs(x = "Real difference between cohorts (standardized mean difference per feature)",
       y = "Simulated difference") + th +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
cl <- r$cls
cl2 <- rbind(data.frame(kind = cl$kind, dir = "Trained on real,\ntested on simulated", acc = cl$train_real_test_sim),
             data.frame(kind = cl$kind, dir = "Trained on simulated,\ntested on real", acc = cl$train_sim_test_real))
cs <- aggregate(acc ~ kind + dir, cl2, function(v) c(m = mean(v), se = sd(v) / sqrt(length(v)))); cs <- do.call(data.frame, cs); names(cs)[3:4] <- c("mean", "se")
cs$kind <- factor(cs$kind, c("Per-study templates", "One template"))
g2b <- ggplot(cs, aes(dir, mean, fill = kind)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), position = position_dodge(width = 0.7), width = 0, colour = "#52514e") +
  geom_hline(yintercept = r$chance, linetype = "dashed", colour = "#8a8984") +
  scale_fill_manual(values = c("Per-study templates" = "#2a78d6", "One template" = "#c3c2b7")) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Samples assigned to the correct cohort") + th + theme(panel.grid.major.y = element_line(colour = "#ecebe8"))
ggsave(file.path(Sys.getenv("OUT_DIR", "ibd_validation"), "ibd_fig2_cohort_differences.png"), g2a, width = 11, height = 4.2, dpi = 200, bg = "white")
ggsave(file.path(Sys.getenv("OUT_DIR", "ibd_validation"), "ibd_fig3_cohort_classifier.png"), g2b, width = 6.5, height = 4.6, dpi = 200, bg = "white")
cat("r per pair (example):\n"); print(sapply(split(ex, ex$pair), function(x) round(cor(x$real, x$sim), 3)))
print(cs, digits = 3); cat("chance", round(r$chance, 3), "\n")
