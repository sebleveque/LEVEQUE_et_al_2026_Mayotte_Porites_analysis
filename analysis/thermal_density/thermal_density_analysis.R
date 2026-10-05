# MAYOTTE THERMAL STRESS + SYMBIONT DENSITY ANALYSIS
#   1) quantify thermal exposure by sampling site and sampling timepoint
#   2) test whether symbiont density was higher in April 2024 than April 2022
#   3) retain May/June 2022 as seasonal context
#   4) export tables and figures


# 0. PACKAGES
required_pkgs <- c(
  "tidyverse", "lme4", "lmerTest", "emmeans", "scales"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      FUN.VALUE = logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them first with:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(scales)
})


old_contrasts <- options("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))

# 1. data and outrput folders
input_file <- file.path("data", "core", "metadata_mayotte_all_v2.txt")
out_dir <- file.path("outputs", "thermal_density")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
diag_dir <- file.path(out_dir, "diagnostics")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, showWarnings = FALSE)
dir.create(tab_dir, showWarnings = FALSE)
dir.create(diag_dir, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file,
       "\nRun the script from the repository root and ensure data/core/metadata_mayotte_all_v2.txt exists.")
}

# 2. READ DATA
meta <- read.delim(
  input_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

# Convert decimal commas to decimal points safely.
to_numeric_decimal_comma <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x), fixed = TRUE)))
}

numeric_cols <- c("density_count", "calibrated_SST", "DHW", "BAA", "Hotspot")
for (nm in numeric_cols) {
  if (!nm %in% names(meta)) stop("Required column missing: ", nm)
  meta[[nm]] <- to_numeric_decimal_comma(meta[[nm]])
}

required_cols <- c(
  "sample_type", "collection_site", "year", "month", "collection_date",
  "colonyID", "density_count", "calibrated_SST", "DHW", "Hotspot"
)
missing_cols <- setdiff(required_cols, names(meta))
if (length(missing_cols) > 0) {
  stop("Required columns missing from metadata: ", paste(missing_cols, collapse = ", "))
}

# coral tissue samples only.
coral <- meta %>%
  filter(sample_type == "coral_field") %>%
  mutate(
    year = as.integer(year),
    month = as.integer(month),
    site = factor(collection_site,
                  levels = c("Boueni", "Mtsamboro", "Musical", "Sakouli", "Tanaraki")),
    sampling = case_when(
      year == 2022 & month == 4 ~ "Apr 2022",
      year == 2022 & month == 5 ~ "May 2022",
      year == 2022 & month == 6 ~ "Jun 2022",
      year == 2024 & month == 4 ~ "Apr 2024",
      TRUE ~ NA_character_
    ),
    sampling = factor(sampling,
                      levels = c("Apr 2022", "May 2022", "Jun 2022", "Apr 2024")),
    log10_density = if_else(density_count > 0, log10(density_count), NA_real_)
  ) %>%
  filter(!is.na(sampling))

cat("\n============================================================\n")
cat("MAYOTTE THERMAL + DENSITY ANALYSIS\n")
cat("============================================================\n")
cat("Coral rows retained:", nrow(coral), "\n")
cat("Density measurements available:", sum(!is.na(coral$density_count)), "\n")
cat("Unique colonies:", n_distinct(coral$colonyID), "\n\n")

# 3. THERMAL SUMMARY TABLE
# Thermal values unique for each Site x Sampling occasion.
thermal_check <- coral %>%
  group_by(site, sampling) %>%
  summarise(
    n_SST = n_distinct(calibrated_SST[!is.na(calibrated_SST)]),
    n_Hotspot = n_distinct(Hotspot[!is.na(Hotspot)]),
    n_DHW = n_distinct(DHW[!is.na(DHW)]),
    .groups = "drop"
  )

if (any(thermal_check$n_SST > 1 | thermal_check$n_Hotspot > 1 | thermal_check$n_DHW > 1)) {
  warning("More than one thermal value was found for at least one Site x Sampling combination. Inspect thermal_uniqueness_check.csv")
}
write.csv(thermal_check, file.path(tab_dir, "thermal_uniqueness_check.csv"), row.names = FALSE)

thermal_table <- coral %>%
  group_by(site, sampling) %>%
  summarise(
    calibrated_SST_C = first(na.omit(calibrated_SST)),
    HotSpot_C = first(na.omit(Hotspot)),
    DHW_C_weeks = first(na.omit(DHW)),
    .groups = "drop"
  ) %>%
  arrange(site, sampling)

write.csv(thermal_table,
          file.path(tab_dir, "Table_thermal_by_site_timepoint.csv"),
          row.names = FALSE)

# Wide version
thermal_wide <- thermal_table %>%
  pivot_wider(
    names_from = sampling,
    values_from = c(calibrated_SST_C, HotSpot_C, DHW_C_weeks),
    names_glue = "{.value}_{sampling}"
  )
write.csv(thermal_wide,
          file.path(tab_dir, "Table_thermal_by_site_timepoint_WIDE.csv"),
          row.names = FALSE)

cat("THERMAL VALUES BY SITE AND SAMPLING OCCASION\n")
print(thermal_table, n = Inf)

#4. THERMAL FIGURE

thermal_long <- thermal_table %>%
  pivot_longer(
    cols = c(calibrated_SST_C, HotSpot_C, DHW_C_weeks),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("calibrated_SST_C", "HotSpot_C", "DHW_C_weeks"),
      labels = c(
        "(a)   Calibrated SST (°C)",
        "(b)   HotSpot (°C)",
        "(c)   DHW (°C-weeks)"
      )
    )
  )

site_cols <- c(
  "Boueni"     = "#0072B2",
  "Mtsamboro"  = "#E69F00",
  "Musical"    = "#009E73",
  "Sakouli"    = "#CC79A7",
  "Tanaraki"   = "#D55E00"
)

theme_pub <- function(base_size = 15) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", hjust = 0),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      plot.title = element_text(size = base_size + 2, face = "bold"),
      plot.subtitle = element_text(size = base_size - 1),
      plot.caption = element_text(size = base_size - 2, hjust = 0),
      panel.spacing = grid::unit(1.0, "lines"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

# Mean and parametric 95% CI
mean_ci95 <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- mean(x)
  if (n < 2) return(data.frame(y = m, ymin = NA_real_, ymax = NA_real_))
  se <- stats::sd(x) / sqrt(n)
  ci <- stats::qt(0.975, df = n - 1) * se
  data.frame(y = m, ymin = m - ci, ymax = m + ci)
}

p_thermal <- ggplot(
  thermal_long,
  aes(x = sampling, y = value, colour = site, group = site)
) +
  geom_line(linewidth = 1.05, alpha = 0.95) +
  geom_point(size = 3.7) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = site_cols, name = "Site", drop = FALSE) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_pub(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))

# High-resolution raster + vector output.
ggsave(file.path(fig_dir, "Figure_thermal_site_timepoint_PUBLICATION.png"),
       p_thermal, width = 8.2, height = 9.2, dpi = 600, bg = "white")
ggsave(file.path(fig_dir, "Figure_thermal_site_timepoint_PUBLICATION.pdf"),
       p_thermal, width = 8.2, height = 9.2, device = "pdf")

#5. DENSITY DESCRIPTIVE TABLE
density_summary <- coral %>%
  filter(!is.na(density_count), density_count > 0) %>%
  group_by(site, sampling) %>%
  summarise(
    n = n(),
    mean_density = mean(density_count),
    sd_density = sd(density_count),
    median_density = median(density_count),
    se_density = sd_density / sqrt(n),
    ci95_low = mean_density - qt(0.975, df = n - 1) * se_density,
    ci95_high = mean_density + qt(0.975, df = n - 1) * se_density,
    .groups = "drop"
  )

write.csv(density_summary,
          file.path(tab_dir, "Table_density_summary_by_site_timepoint.csv"),
          row.names = FALSE)

cat("\nDENSITY SUMMARY BY SITE AND SAMPLING OCCASION\n")
print(density_summary, n = Inf)

#6. FULL REPEATED-MEASURES MODEL
density_dat <- coral %>%
  filter(!is.na(log10_density), !is.na(site), !is.na(colonyID)) %>%
  droplevels()

# Full model retains the 2022 seasonal trajectory while testing Site x Sampling.
model_full <- lmer(
  log10_density ~ sampling * site + (1 | colonyID),
  data = density_dat,
  REML = TRUE
)

anova_full <- anova(model_full, type = 3, ddf = "Satterthwaite")
anova_full_df <- as.data.frame(anova_full) %>%
  rownames_to_column("term")
write.csv(anova_full_df,
          file.path(tab_dir, "Model_full_density_TypeIII_ANOVA.csv"),
          row.names = FALSE)

capture.output(summary(model_full),
               file = file.path(tab_dir, "Model_full_density_summary.txt"))

# Estimated marginal means by sampling occasion, averaging equally across sites.
emm_sampling <- emmeans(model_full, ~ sampling, weights = "equal")
write.csv(as.data.frame(emm_sampling),
          file.path(tab_dir, "Model_full_density_EMMeans_sampling.csv"),
          row.names = FALSE)

# All pairwise temporal contrasts
pair_sampling <- pairs(emm_sampling, adjust = "tukey")
write.csv(as.data.frame(pair_sampling),
          file.path(tab_dir, "Model_full_density_all_pairwise_sampling.csv"),
          row.names = FALSE)

# Targeted biologically relevant contrasts on log10 scale.
# Positive values below mean the later/2024 time point has higher density.
target_contrasts <- contrast(
  emm_sampling,
  method = list(
    "Apr2024_minus_Apr2022" = c(-1, 0, 0, 1),
    "May2022_minus_Apr2022" = c(-1, 1, 0, 0),
    "Jun2022_minus_Apr2022" = c(-1, 0, 1, 0),
    "Apr2024_minus_Jun2022" = c(0, 0, -1, 1)
  ),
  adjust = "holm"
)

target_df <- as.data.frame(confint(target_contrasts)) %>%
  mutate(
    fold_change = 10^estimate,
    fold_change_CI_low = 10^lower.CL,
    fold_change_CI_high = 10^upper.CL,
    interpretation = case_when(
      fold_change > 1 ~ "numerator timepoint has higher estimated density",
      fold_change < 1 ~ "numerator timepoint has lower estimated density",
      TRUE ~ "no difference"
    )
  )
write.csv(target_df,
          file.path(tab_dir, "Model_full_density_TARGET_contrasts_fold_change.csv"),
          row.names = FALSE)

# Site-specific April 2024 minus April 2022 effects.
emm_site_sampling <- emmeans(model_full, ~ sampling | site)
site_april_contrast <- contrast(
  emm_site_sampling,
  method = list("Apr2024_minus_Apr2022" = c(-1, 0, 0, 1)),
  adjust = "none"
)
site_april_df <- as.data.frame(confint(site_april_contrast)) %>%
  mutate(
    fold_change = 10^estimate,
    fold_change_CI_low = 10^lower.CL,
    fold_change_CI_high = 10^upper.CL
  )
write.csv(site_april_df,
          file.path(tab_dir, "Model_full_density_Apr2024_vs_Apr2022_BY_SITE.csv"),
          row.names = FALSE)

cat("\nFULL DENSITY MODEL — TYPE III TESTS\n")
print(anova_full)
cat("\nTARGET DENSITY CONTRASTS (fold-change included in CSV)\n")
print(target_df)

#7. PRIMARY MATCHED APRIL 2022 vs 2024
april_dat <- density_dat %>%
  filter(sampling %in% c("Apr 2022", "Apr 2024")) %>%
  mutate(april_year = factor(as.character(sampling),
                             levels = c("Apr 2022", "Apr 2024")))

# Restrict to colonies with density observations in both Apr 2022 and Apr 2024.
paired_ids <- april_dat %>%
  count(colonyID, april_year) %>%
  complete(colonyID, april_year, fill = list(n = 0)) %>%
  group_by(colonyID) %>%
  summarise(has_both = all(n > 0), .groups = "drop") %>%
  filter(has_both) %>%
  pull(colonyID)

april_paired <- april_dat %>%
  filter(colonyID %in% paired_ids) %>%
  droplevels()

cat("\nMatched colonies with both Apr 2022 and Apr 2024 density:",
    n_distinct(april_paired$colonyID), "\n")

paired_counts <- april_paired %>% count(site, april_year)
write.csv(paired_counts,
          file.path(tab_dir, "Matched_April_sample_sizes.csv"),
          row.names = FALSE)

model_april <- lmer(
  log10_density ~ april_year * site + (1 | colonyID),
  data = april_paired,
  REML = TRUE
)

anova_april <- anova(model_april, type = 3, ddf = "Satterthwaite")
anova_april_df <- as.data.frame(anova_april) %>%
  rownames_to_column("term")
write.csv(anova_april_df,
          file.path(tab_dir, "PRIMARY_April2022_vs_April2024_TypeIII_ANOVA.csv"),
          row.names = FALSE)

capture.output(summary(model_april),
               file = file.path(tab_dir, "PRIMARY_April2022_vs_April2024_model_summary.txt"))

# Overall April effect averaged equally across sites.
emm_april <- emmeans(model_april, ~ april_year, weights = "equal")
april_effect <- contrast(
  emm_april,
  method = list("Apr2024_minus_Apr2022" = c(-1, 1)),
  adjust = "none"
)

april_effect_df <- as.data.frame(confint(april_effect)) %>%
  mutate(
    fold_change = 10^estimate,
    fold_change_CI_low = 10^lower.CL,
    fold_change_CI_high = 10^upper.CL
  )
write.csv(april_effect_df,
          file.path(tab_dir, "PRIMARY_April2024_vs_April2022_fold_change.csv"),
          row.names = FALSE)

# Site-specific matched April effects.
emm_april_site <- emmeans(model_april, ~ april_year | site)
april_site_effect <- contrast(
  emm_april_site,
  method = list("Apr2024_minus_Apr2022" = c(-1, 1)),
  adjust = "none"
)
april_site_effect_df <- as.data.frame(confint(april_site_effect)) %>%
  mutate(
    fold_change = 10^estimate,
    fold_change_CI_low = 10^lower.CL,
    fold_change_CI_high = 10^upper.CL
  )
write.csv(april_site_effect_df,
          file.path(tab_dir, "PRIMARY_April2024_vs_April2022_fold_change_BY_SITE.csv"),
          row.names = FALSE)

cat("\nPRIMARY MATCHED APRIL MODEL — TYPE III TESTS\n")
print(anova_april)
cat("\nPRIMARY APRIL 2024 / APRIL 2022 EFFECT\n")
print(april_effect_df)
cat("\nSITE-SPECIFIC APRIL 2024 / APRIL 2022 EFFECTS\n")
print(april_site_effect_df)

#8. DIAGNOSTICS
png(file.path(diag_dir, "diagnostics_full_model.png"),
    width = 1800, height = 900, res = 180)
par(mfrow = c(1, 2))
plot(fitted(model_full), resid(model_full),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Full model: residuals vs fitted")
abline(h = 0, lty = 2)
qqnorm(resid(model_full), main = "Full model: residual Q-Q")
qqline(resid(model_full))
dev.off()

png(file.path(diag_dir, "diagnostics_primary_april_model.png"),
    width = 1800, height = 900, res = 180)
par(mfrow = c(1, 2))
plot(fitted(model_april), resid(model_april),
     xlab = "Fitted values", ylab = "Residuals",
     main = "April model: residuals vs fitted")
abline(h = 0, lty = 2)
qqnorm(resid(model_april), main = "April model: residual Q-Q")
qqline(resid(model_april))
dev.off()
par(mfrow = c(1, 1))

#9. DENSITY FIGURES

site_labeller <- c(
  "Boueni" = "Boueni",
  "Mtsamboro" = "Mtsamboro",
  "Musical" = "Musical",
  "Sakouli" = "Sakouli",
  "Tanaraki" = "Tanaraki"
)

# Figure A: all sampling occasions, per site.
# Raw observations are retained. Large diamonds + 95% CI show the site/timepoint mean.
# Colour encodes site consistently with the thermal figure.
p_density_all <- ggplot(
  density_dat,
  aes(x = sampling, y = density_count, colour = site)
) +
  geom_jitter(width = 0.11, height = 0, alpha = 0.42, size = 2.0,
              show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", size = 4.2, shape = 18,
               show.legend = FALSE) +
  stat_summary(fun.data = mean_ci95, geom = "errorbar", width = 0.11,
               linewidth = 0.95, show.legend = FALSE) +
  facet_wrap(
    ~ site, nrow = 1,
    labeller = as_labeller(site_labeller)
  ) +
  scale_colour_manual(values = site_cols, drop = FALSE) +
  scale_y_log10(
    labels = label_number(scale = 1e-6, suffix = " M", accuracy = 0.1)
  ) +
  labs(
    x = NULL,
    y = expression("Symbiodiniaceae density (cells; log"[10]~"scale)")
  ) +
  theme_pub(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(file.path(fig_dir, "Figure_density_all_timepoints_by_site_PUBLICATION.png"),
       p_density_all, width = 14.2, height = 5.7, dpi = 600, bg = "white")
ggsave(file.path(fig_dir, "Figure_density_all_timepoints_by_site_PUBLICATION.pdf"),
       p_density_all, width = 14.2, height = 5.7, device = "pdf")

# Figure B: primary matched April comparison.
# Grey lines connect the same tagged colony; colored points identify the site.
# Large diamonds are observed site means. The raw paired response remains visible.
p_density_paired <- ggplot(
  april_paired,
  aes(x = april_year, y = density_count, group = colonyID)
) +
  geom_line(colour = "grey55", alpha = 0.50, linewidth = 0.65) +
  geom_point(aes(colour = site), alpha = 0.78, size = 2.5,
             show.legend = FALSE) +
  stat_summary(aes(group = 1, colour = site), fun = mean, geom = "line",
               linewidth = 1.35, show.legend = FALSE) +
  stat_summary(aes(group = 1, colour = site), fun = mean, geom = "point",
               size = 4.6, shape = 18, show.legend = FALSE) +
  facet_wrap(
    ~ site, nrow = 1,
    labeller = as_labeller(site_labeller)
  ) +
  scale_colour_manual(values = site_cols, drop = FALSE) +
  scale_y_log10(
    labels = label_number(scale = 1e-6, suffix = " M", accuracy = 0.1)
  ) +
  labs(
    x = NULL,
    y = expression("Symbiodiniaceae density (cells; log"[10]~"scale)")
  ) +
  theme_pub(base_size = 15) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave(file.path(fig_dir, "Figure_density_MATCHED_Apr2022_vs_Apr2024_by_site_PUBLICATION.png"),
       p_density_paired, width = 14.2, height = 5.7, dpi = 600, bg = "white")
ggsave(file.path(fig_dir, "Figure_density_MATCHED_Apr2022_vs_Apr2024_by_site_PUBLICATION.pdf"),
       p_density_paired, width = 14.2, height = 5.7, device = "pdf")

# Figure C: model derived April 2024 / April 2022 fold change by site.

fold_plot_dat <- april_site_effect_df %>%
  mutate(
    site = factor(site, levels = levels(density_dat$site)),
    site_label = factor(
      as.character(site),
      levels = levels(density_dat$site),
      labels = unname(site_labeller[levels(density_dat$site)])
    )
  )

p_fold <- ggplot(
  fold_plot_dat,
  aes(x = site_label, y = fold_change, colour = site)
) +
  geom_errorbar(
    aes(ymin = fold_change_CI_low, ymax = fold_change_CI_high),
    width = 0.12, linewidth = 1.0
  ) +
  geom_point(size = 4.3) +
  scale_colour_manual(values = site_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "April 2024 / April 2022 density ratio"
  ) +
  theme_pub(base_size = 15) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(fig_dir, "Figure_density_foldchange_by_site_PUBLICATION.png"),
       p_fold, width = 8.4, height = 6.0, dpi = 600, bg = "white")
ggsave(file.path(fig_dir, "Figure_density_foldchange_by_site_PUBLICATION.pdf"),
       p_fold, width = 8.4, height = 6.0, device = "pdf")

#10. MANUSCRIPT tablers 
results_key <- list(
  thermal = thermal_table,
  density_primary_effect = april_effect_df,
  density_primary_typeIII = anova_april_df,
  density_primary_by_site = april_site_effect_df,
  density_full_target_contrasts = target_df
)
saveRDS(results_key, file.path(out_dir, "key_results_for_manuscript.rds"))

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("Outputs written to:", normalizePath(out_dir), "\n")
cat("\nMost important files:\n")
cat("  tables/Table_thermal_by_site_timepoint.csv\n")
cat("  tables/PRIMARY_April2022_vs_April2024_TypeIII_ANOVA.csv\n")
cat("  tables/PRIMARY_April2024_vs_April2022_fold_change.csv\n")
cat("  tables/PRIMARY_April2024_vs_April2022_fold_change_BY_SITE.csv\n")
cat("  figures/Figure_thermal_site_timepoint_PUBLICATION.png\n")
cat("  figures/Figure_density_all_timepoints_by_site_PUBLICATION.png\n")
cat("  figures/Figure_density_MATCHED_Apr2022_vs_Apr2024_by_site_PUBLICATION.png\n")
cat("  figures/Figure_density_foldchange_by_site_PUBLICATION.png\n")
cat("============================================================\n")
options(contrasts = old_contrasts$contrasts)
