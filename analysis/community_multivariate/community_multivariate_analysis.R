# Mayotte Symbiodiniaceae community and multivariate analyses
#
# Main analysis choices preserved from the submitted analysis script
#   - P5 read-depth filtering
#   - 1% rare-taxon filter
#   - Hellinger transformation followed by Bray-Curtis dissimilarity
#   - NMDS, PERMANOVA, PERMDISP, envfit, SIMPER and indicator-species analyses

# 0. Packages, repository paths and output folders

library(vegan)
library(indicspecies)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(patchwork)

options(stringsAsFactors = FALSE)

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_core_dir <- file.path(repo_root, "data", "core")
outputs_dir <- file.path(repo_root, "outputs")
figures_dir <- file.path(outputs_dir, "figures")
tables_dir <- file.path(outputs_dir, "tables")
rds_dir <- file.path(outputs_dir, "rds")
main_tables_dir <- file.path(outputs_dir, "manuscript_tables")
supp_tables_dir <- file.path(outputs_dir, "supplementary_tables")
diagnostics_dir <- file.path(outputs_dir, "diagnostics")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(main_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)

required_input_files <- file.path(
  data_core_dir,
  c("postMED_rel_abund_only.txt",
    "postMED_abs_abund_only.txt",
    "metadata_mayotte_all_v2.txt")
)
missing_input_files <- required_input_files[!file.exists(required_input_files)]
if (length(missing_input_files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(paste0(" - ", missing_input_files), collapse = "\n"),
    "\nRun the script from the repository root."
  )
}

# 1. Load raw data

postMED <- read.table(file.path(data_core_dir, "postMED_rel_abund_only.txt"), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
postMED_abs <- read.table(file.path(data_core_dir, "postMED_abs_abund_only.txt"), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
meta <- read.table(file.path(data_core_dir, "metadata_mayotte_all_v2.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

required_cols <- c(
  "sample_uid", "sample_type", "collection_site", "timepoint",
  "year", "month", "calibrated_SST", "DHW", "BAA", "Hotspot",
  "NDTI", "density_count"
)
stopifnot(all(required_cols %in% names(meta)))

rownames(meta) <- meta$sample_uid
shared_ids <- Reduce(intersect, list(rownames(meta), rownames(postMED), rownames(postMED_abs)))
meta <- meta[shared_ids, , drop = FALSE]
postMED <- postMED[shared_ids, , drop = FALSE]
postMED_abs <- postMED_abs[shared_ids, , drop = FALSE]

# 2. P5 filtering and 1% rare-taxon filtering

read_depth <- rowSums(postMED_abs, na.rm = TRUE)
th_p5 <- as.numeric(quantile(read_depth, 0.05, na.rm = TRUE))
keep_p5 <- read_depth >= th_p5

postMED_p5 <- postMED[keep_p5, , drop = FALSE]
postMED_abs_p5 <- postMED_abs[keep_p5, , drop = FALSE]
meta_p5 <- meta[rownames(postMED_p5), , drop = FALSE]

lost_p5_table <- meta[!keep_p5, c("sample_uid", "sample_type", "collection_site", "timepoint", "year", "month"), drop = FALSE]
lost_p5_table$read_depth <- read_depth[!keep_p5]
write.csv(lost_p5_table, file.path(tables_dir, "lost_samples_P5.csv"), row.names = FALSE)

rare_threshold <- 0.01
max_abund <- apply(postMED_p5, 2, max, na.rm = TRUE)
keep_taxa <- max_abund >= rare_threshold

postMED_filt <- postMED_p5[, keep_taxa, drop = FALSE]
postMED_abs_filt <- postMED_abs_p5[, colnames(postMED_filt), drop = FALSE]
meta_filt <- meta_p5[rownames(postMED_filt), , drop = FALSE]

cat("P5 threshold:", th_p5, "reads\n")
cat("Samples before P5:", nrow(postMED), "\n")
cat("Samples after P5:", nrow(postMED_filt), "\n")
cat("Taxa before 1% filter:", ncol(postMED_p5), "\n")
cat("Taxa after 1% filter:", ncol(postMED_filt), "\n")

# 3. Clean metadata and environmental variables by compartment

meta_filt$sample_type <- factor(meta_filt$sample_type)
meta_filt$collection_site <- factor(meta_filt$collection_site)
meta_filt$timepoint <- factor(meta_filt$timepoint)
meta_filt$year <- factor(as.character(meta_filt$year))
meta_filt$month <- factor(as.character(meta_filt$month), levels = c("4", "5", "6"))

meta_filt$calibrated_SST <- as.numeric(gsub(",", ".", trimws(as.character(meta_filt$calibrated_SST))))
meta_filt$DHW <- as.numeric(gsub(",", ".", trimws(as.character(meta_filt$DHW))))
meta_filt$BAA <- as.numeric(gsub(",", ".", trimws(as.character(meta_filt$BAA))))
meta_filt$Hotspot <- as.numeric(gsub(",", ".", trimws(as.character(meta_filt$Hotspot))))
meta_filt$NDTI <- as.numeric(gsub(",", ".", trimws(as.character(meta_filt$NDTI))))
meta_filt$density_count <- as.numeric(gsub(",", ".", trimws(as.character(meta_filt$density_count))))

# Seawater NDTI_for_analysis from mean coral NDTI at same site x timepoint.
coral_ndti <- aggregate(
  NDTI ~ collection_site + timepoint,
  data = meta_filt[meta_filt$sample_type == "coral_field", ],
  FUN = mean,
  na.rm = TRUE
)
coral_ndti$key <- paste(coral_ndti$collection_site, coral_ndti$timepoint, sep = "__")
meta_filt$key <- paste(meta_filt$collection_site, meta_filt$timepoint, sep = "__")
ndti_lookup <- setNames(coral_ndti$NDTI, coral_ndti$key)
meta_filt$NDTI_for_analysis <- meta_filt$NDTI
meta_filt$NDTI_for_analysis[meta_filt$sample_type == "seawater"] <- ndti_lookup[meta_filt$key[meta_filt$sample_type == "seawater"]]

# Standardized thermal metrics for comparable model scales.
meta_filt$thermal_cal_SST <- as.numeric(scale(meta_filt$calibrated_SST))
meta_filt$thermal_DHW <- as.numeric(scale(meta_filt$DHW))
meta_filt$thermal_BAA <- as.numeric(scale(meta_filt$BAA))
meta_filt$thermal_Hotspot <- as.numeric(scale(meta_filt$Hotspot))

# Thermal PC1, oriented so high values = high thermal stress.
thermal_mat <- meta_filt[, c("calibrated_SST", "DHW", "BAA", "Hotspot")]
thermal_ok <- complete.cases(thermal_mat)
thermal_pc <- prcomp(thermal_mat[thermal_ok, ], center = TRUE, scale. = TRUE)
meta_filt$thermal_PC1 <- NA_real_
meta_filt$thermal_PC1[thermal_ok] <- thermal_pc$x[, 1]
pc_cor <- cor(meta_filt$thermal_PC1, meta_filt$DHW, use = "pairwise.complete.obs")
meta_filt$thermal_PC1 <- meta_filt$thermal_PC1 * sign(pc_cor)
write.csv(thermal_pc$rotation, file.path(tables_dir, "thermal_PCA_loadings.csv"))
write.csv(summary(thermal_pc)$importance, file.path(tables_dir, "thermal_PCA_importance.csv"))

# Extreme sediment samples removed from plots only.
extreme_sediment_samples <- c("264610")

# 4. Manual functions

prepare_dataset <- function(comm, meta_df) {
  ids <- rownames(meta_df)
  comm <- comm[ids, , drop = FALSE]
  keep_rows <- rowSums(comm, na.rm = TRUE) > 0
  comm <- comm[keep_rows, , drop = FALSE]
  meta_df <- meta_df[rownames(comm), , drop = FALSE]
  keep_cols <- colSums(comm, na.rm = TRUE) > 0
  comm <- comm[, keep_cols, drop = FALSE]
  list(comm = comm, meta = droplevels(meta_df))
}

cb_site_cols <- c(
  "Boueni" = "#0072B2",
  "Mtsamboro" = "#D55E00",
  "Musical" = "#009E73",
  "Sakouli" = "#CC79A7",
  "Tanaraki" = "#E69F00"
)

time_cols <- c(
  "2022" = "#0072B2",
  "2024" = "#D55E00",
  "April 2022" = "#0072B2",
  "May 2022" = "#009E73",
  "June 2022" = "#CC79A7",
  "2022_04" = "#0072B2",
  "2022_05" = "#009E73",
  "2022_06" = "#CC79A7"
)

time_pch <- c(
  "2022" = 16,
  "2024" = 15,
  "April 2022" = 16,
  "May 2022" = 15,
  "June 2022" = 18,
  "2022_04" = 16,
  "2022_05" = 15,
  "2022_06" = 18
)

simple_env_names <- c(
  "thermal_cal_SST" = "cal SST",
  "thermal_DHW" = "DHW",
  "thermal_BAA" = "BAA",
  "thermal_Hotspot" = "Hotspot",
  "NDTI_for_analysis" = "NDTI",
  "density_count" = "density"
)

run_permanova <- function(dat, rhs, label) {
  vars <- all.vars(as.formula(paste("~", rhs)))
  keep <- complete.cases(dat$meta[, vars, drop = FALSE])
  comm <- dat$comm[keep, , drop = FALSE]
  meta <- droplevels(dat$meta[keep, , drop = FALSE])
  f <- as.formula(paste("comm ~", rhs))
  out <- adonis2(f, data = meta, permutations = 999, method = "bray", by = "margin")
  tab <- as.data.frame(out)
  tab$term <- rownames(tab)
  tab$analysis <- label
  write.csv(tab, file.path(tables_dir, paste0("PERMANOVA_", gsub("[^A-Za-z0-9]+", "_", label), ".csv")), row.names = FALSE)
  print(out)
  out
}

run_betadisper_test <- function(dat, group_var, label) {
  keep <- !is.na(dat$meta[[group_var]])
  comm <- dat$comm[keep, , drop = FALSE]
  meta <- droplevels(dat$meta[keep, , drop = FALSE])
  bd <- betadisper(vegdist(comm, method = "bray"), group = meta[[group_var]])
  a <- anova(bd)
  p <- permutest(bd, permutations = 999)
  write.csv(as.data.frame(a), file.path(tables_dir, paste0("BETADISPER_anova_", gsub("[^A-Za-z0-9]+", "_", label), ".csv")))
  capture.output(p, file = file.path(tables_dir, paste0("BETADISPER_permutest_", gsub("[^A-Za-z0-9]+", "_", label), ".txt")))
  print(a)
  print(p)
  bd
}

run_envfit_table <- function(dat, env_vars, label) {
  vars <- env_vars
  keep <- complete.cases(dat$meta[, vars, drop = FALSE])
  comm <- dat$comm[keep, , drop = FALSE]
  meta <- dat$meta[keep, , drop = FALSE]
  comm_hell <- decostand(comm, method = "hellinger")
  set.seed(123)
  nmds <- metaMDS(comm_hell, distance = "bray", k = 2, trymax = 300, autotransform = FALSE, trace = FALSE)
  fit <- envfit(nmds, meta[, vars, drop = FALSE], permutations = 999, na.rm = TRUE)
  vectors <- as.data.frame(scores(fit, display = "vectors"))
  vectors$variable <- rownames(vectors)
  vectors$r2 <- fit$vectors$r
  vectors$p_value <- fit$vectors$pvals
  vectors$label <- simple_env_names[vectors$variable]
  write.csv(vectors, file.path(tables_dir, paste0("ENVFIT_", gsub("[^A-Za-z0-9]+", "_", label), ".csv")), row.names = FALSE)
  print(fit)
  fit
}

run_simper_summary <- function(dat, group_var, label, top_n = 20) {
  keep <- !is.na(dat$meta[[group_var]])
  comm <- dat$comm[keep, , drop = FALSE]
  groups <- droplevels(dat$meta[[group_var]][keep])
  sim <- simper(comm, groups, permutations = 999)
  sim_names <- names(sim)
  sim_tables <- lapply(sim_names, function(nm) {
    x <- as.data.frame(sim[[nm]])
    x$taxon <- rownames(x)
    x$comparison <- nm
    x <- x[order(x$average, decreasing = TRUE), , drop = FALSE]
    head(x, top_n)
  })
  out <- do.call(rbind, sim_tables)
  out$analysis <- label
  write.csv(out, file.path(tables_dir, paste0("SIMPER_", gsub("[^A-Za-z0-9]+", "_", label), ".csv")), row.names = FALSE)
  out
}

run_indicator_summary <- function(dat, group_var, label) {
  keep <- !is.na(dat$meta[[group_var]])
  comm <- dat$comm[keep, , drop = FALSE]
  groups <- droplevels(dat$meta[[group_var]][keep])
  ind <- multipatt(comm, groups, func = "IndVal.g", control = how(nperm = 999))
  out <- as.data.frame(ind$sign)
  out$taxon <- rownames(out)
  out$analysis <- label
  out <- out[order(out$p.value), , drop = FALSE]
  write.csv(out, file.path(tables_dir, paste0("INDICATOR_", gsub("[^A-Za-z0-9]+", "_", label), ".csv")), row.names = FALSE)
  print(summary(ind))
  ind
}

# 5. Manual NMDS objects layer by layer plotting

make_manual_nmds_object <- function(dat, color_var = "collection_site", symbol_var = "time_group",
                                    ellipse_var = "time_group", env_vars = character(0),
                                    remove_plot_samples = character(0),
                                    trajectory_var = NULL, trajectory_both_groups = NULL,
                                    seed = 123, trymax = 500) {

  stopifnot(color_var %in% names(dat$meta))
  stopifnot(symbol_var %in% names(dat$meta))
  stopifnot(ellipse_var %in% names(dat$meta))

  plot_meta <- dat$meta[!rownames(dat$meta) %in% remove_plot_samples, , drop = FALSE]
  plot_comm <- dat$comm[rownames(plot_meta), , drop = FALSE]
  plot_dat <- prepare_dataset(plot_comm, plot_meta)
  plot_comm <- plot_dat$comm
  plot_meta <- plot_dat$meta

  comm_hell <- vegan::decostand(plot_comm, method = "hellinger")

  set.seed(seed)
  nmds <- vegan::metaMDS(
    comm_hell,
    distance = "bray",
    k = 2,
    trymax = trymax,
    autotransform = FALSE,
    trace = FALSE
  )

  scr <- as.data.frame(vegan::scores(nmds, display = "sites"))
  scr$sample_id <- rownames(scr)
  scr <- cbind(scr, plot_meta[rownames(scr), , drop = FALSE])

  color_groups <- factor(scr[[color_var]])
  symbol_groups <- factor(scr[[symbol_var]])
  ellipse_groups <- factor(scr[[ellipse_var]])

  # Point colors
  color_levels <- levels(color_groups)
  color_palette <- cb_site_cols[color_levels]
  missing_cols <- is.na(color_palette)
  fallback_cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")
  color_palette[missing_cols] <- rep(fallback_cols, length.out = sum(missing_cols))
  names(color_palette) <- color_levels
  point_cols <- color_palette[as.character(color_groups)]

  # Point symbols
  symbol_levels <- levels(symbol_groups)
  symbol_palette <- time_pch[symbol_levels]
  missing_pch <- is.na(symbol_palette)
  fallback_pch <- c(16, 15, 18, 17, 8, 3, 4, 7)
  symbol_palette[missing_pch] <- rep(fallback_pch, length.out = sum(missing_pch))
  names(symbol_palette) <- symbol_levels
  point_pchs <- symbol_palette[as.character(symbol_groups)]

  # Ellipse colors
  ellipse_levels <- levels(ellipse_groups)
  ellipse_palette <- time_cols[ellipse_levels]
  missing_ell <- is.na(ellipse_palette)
  fallback_ell <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000")
  ellipse_palette[missing_ell] <- rep(fallback_ell, length.out = sum(missing_ell))
  names(ellipse_palette) <- ellipse_levels

  # Manual ellipse coordinates
  ellipse_list <- list()
  for (g in ellipse_levels) {
    idx <- ellipse_groups == g
    if (sum(idx, na.rm = TRUE) >= 3) {
      pts <- as.matrix(scr[idx, c("NMDS1", "NMDS2")])
      cen <- colMeans(pts, na.rm = TRUE)
      S <- stats::cov(pts, use = "complete.obs")
      if (all(is.finite(S))) {
        theta <- seq(0, 2 * pi, length.out = 160)
        circle <- cbind(cos(theta), sin(theta))
        eig <- eigen(S)
        radius <- sqrt(stats::qchisq(0.95, df = 2))
        ell <- sweep(circle %*% diag(sqrt(pmax(eig$values, 0))) %*% t(eig$vectors) * radius,
                     2, cen, "+")
        ellipse_list[[g]] <- ell
      }
    }
  }

  # Envfit arrows
  envfit_obj <- NULL
  env_arrows <- NULL
  env_vars_use <- env_vars[env_vars %in% names(plot_meta)]
  if (length(env_vars_use) > 0 && nrow(plot_meta) >= 6) {
    envfit_obj <- vegan::envfit(nmds, plot_meta[, env_vars_use, drop = FALSE], permutations = 999, na.rm = TRUE)
    ef <- as.data.frame(vegan::scores(envfit_obj, display = "vectors"))
    if (nrow(ef) > 0) {
      ef$variable <- rownames(ef)
      ef$label <- simple_env_names[ef$variable]
      ef$label[is.na(ef$label)] <- ef$variable[is.na(ef$label)]
      mult <- vegan:::ordiArrowMul(as.matrix(ef[, c("NMDS1", "NMDS2")]))
      ef$NMDS1s <- ef$NMDS1 * mult
      ef$NMDS2s <- ef$NMDS2 * mult
      env_arrows <- ef
    }
  }

  # Trajectory coordinates for paired coral colonies
  traj <- NULL
  if (!is.null(trajectory_var) && trajectory_var %in% names(scr) && !is.null(trajectory_both_groups)) {
    traj_list <- list()
    ids <- unique(na.omit(scr[[trajectory_var]]))
    for (id in ids) {
      sub <- scr[scr[[trajectory_var]] == id & as.character(scr[[symbol_var]]) %in% trajectory_both_groups, , drop = FALSE]
      if (length(unique(as.character(sub[[symbol_var]]))) == 2) {
        xy <- aggregate(cbind(NMDS1, NMDS2) ~ time_tmp,
                        data = transform(sub, time_tmp = as.character(sub[[symbol_var]])),
                        FUN = mean)
        xy <- xy[match(trajectory_both_groups, xy$time_tmp), , drop = FALSE]
        if (nrow(xy) == 2 && all(!is.na(xy$NMDS1))) {
          traj_list[[id]] <- data.frame(
            trajectory_id = id,
            x0 = xy$NMDS1[1], y0 = xy$NMDS2[1],
            x1 = xy$NMDS1[2], y1 = xy$NMDS2[2],
            stringsAsFactors = FALSE
          )
        }
      }
    }
    if (length(traj_list) > 0) traj <- do.call(rbind, traj_list)
  }

  out <- list(
    nmds = nmds,
    scores = scr,
    color_var = color_var,
    symbol_var = symbol_var,
    ellipse_var = ellipse_var,
    point_cols = point_cols,
    point_pchs = point_pchs,
    color_palette = color_palette,
    symbol_palette = symbol_palette,
    ellipse_palette = ellipse_palette,
    ellipse_list = ellipse_list,
    envfit = envfit_obj,
    env_arrows = env_arrows,
    trajectories = traj
  )
  class(out) <- "manual_nmds_layers"
  out
}

manual_nmds_limits <- function(obj, expand = 0.12, include_ellipses = TRUE) {
  x <- obj$scores$NMDS1
  y <- obj$scores$NMDS2
  if (include_ellipses && length(obj$ellipse_list) > 0) {
    ell <- do.call(rbind, obj$ellipse_list)
    x <- c(x, ell[, 1])
    y <- c(y, ell[, 2])
  }
  xr <- range(x, na.rm = TRUE)
  yr <- range(y, na.rm = TRUE)
  xr <- xr + c(-1, 1) * diff(xr) * expand
  yr <- yr + c(-1, 1) * diff(yr) * expand
  list(xlim = xr, ylim = yr)
}

open_manual_nmds_plot <- function(obj, title, xlim = NULL, ylim = NULL, expand = 0.12,
                                  main_cex = 1, mar = c(5, 5, 4, 2)) {
  lim <- manual_nmds_limits(obj, expand = expand, include_ellipses = TRUE)
  if (is.null(xlim)) xlim <- lim$xlim
  if (is.null(ylim)) ylim <- lim$ylim
  par(mar = mar, xpd = NA)
  plot(obj$scores$NMDS1, obj$scores$NMDS2,
       type = "n",
       xlab = "NMDS1",
       ylab = "NMDS2",
       xlim = xlim,
       ylim = ylim,
       main = paste0(title, "\nStress = ", round(obj$nmds$stress, 3)),
       cex.main = main_cex)
  invisible(list(xlim = xlim, ylim = ylim))
}

add_manual_ellipses <- function(obj, alpha = 0.18, lwd = 2.2, draw_border = TRUE) {
  for (g in names(obj$ellipse_list)) {
    ell <- obj$ellipse_list[[g]]
    border_col <- if (draw_border) obj$ellipse_palette[g] else NA
    polygon(ell[, 1], ell[, 2],
            col = adjustcolor(obj$ellipse_palette[g], alpha.f = alpha),
            border = border_col,
            lwd = lwd)
  }
}

add_manual_trajectories <- function(obj, col = "grey35", alpha = 0.35,
                                    lwd = 0.65, arrow_length = 0.045) {
  if (is.null(obj$trajectories)) return(invisible(NULL))
  with(obj$trajectories, arrows(x0, y0, x1, y1,
                                col = adjustcolor(col, alpha.f = alpha),
                                lwd = lwd,
                                length = arrow_length))
}

add_manual_points <- function(obj, point_cex = 1.35, point_lwd = 2, bg = NULL) {
  points(obj$scores$NMDS1, obj$scores$NMDS2,
         col = obj$point_cols,
         pch = obj$point_pchs,
         cex = point_cex,
         lwd = point_lwd,
         bg = bg)
}

add_manual_envfit <- function(obj, arrow_lwd = 1.4, arrow_length = 0.08,
                              arrow_col = "black", label_cex = 0.82,
                              label_font = 2, label_offset = 1.28,
                              label_jitter_y = 0.035) {
  ef <- obj$env_arrows
  if (is.null(ef) || nrow(ef) == 0) return(invisible(NULL))
  arrows(0, 0, ef$NMDS1s, ef$NMDS2s,
         length = arrow_length,
         lwd = arrow_lwd,
         col = arrow_col)
  jitter_y <- seq(-label_jitter_y, label_jitter_y, length.out = nrow(ef))
  text(ef$NMDS1s * label_offset,
       ef$NMDS2s * label_offset + jitter_y,
       labels = ef$label,
       cex = label_cex,
       font = label_font)
}

add_manual_legends <- function(obj, site_pos = "topright", time_pos = "bottomright",
                               ellipse_pos = NULL, legend_cex = 0.8,
                               site_title = "Site", time_title = "Time",
                               ellipse_title = "Ellipse") {
  legend(site_pos,
         legend = names(obj$color_palette),
         col = obj$color_palette,
         pch = 16,
         pt.cex = 1.1,
         title = site_title,
         bty = "n",
         cex = legend_cex)

  legend(time_pos,
         legend = names(obj$symbol_palette),
         pch = obj$symbol_palette,
         pt.cex = 1.1,
         title = time_title,
         bty = "n",
         cex = legend_cex)

  if (!is.null(ellipse_pos)) {
    legend(ellipse_pos,
           legend = names(obj$ellipse_palette),
           fill = adjustcolor(obj$ellipse_palette, alpha.f = 0.30),
           border = obj$ellipse_palette,
           title = ellipse_title,
           bty = "n",
           cex = legend_cex)
  }
}

save_manual_nmds_png <- function(filename, plot_code, width = 2200, height = 1900, res = 220) {
  png(filename, width = width, height = height, res = res)
  force(plot_code)
  dev.off()
}

# 6. NMDS PLOTTING FUNCTION

plot_nmds_base <- function(obj,
                           main_title = "NMDS",
                           xlim = NULL,
                           ylim = NULL,
                           expand = 0.15,
                           mar = c(5, 5, 4, 2),
                           point_cex = 1.35,
                           point_lwd = 2,
                           ellipse_alpha = 0.20,
                           ellipse_lwd = 2.2,
                           trajectory_lwd = 0.65,
                           trajectory_alpha = 0.35,
                           trajectory_arrow_length = 0.045,
                           trajectory_col = "grey35",
                           env_arrow_lwd = 1.4,
                           env_arrow_length = 0.08,
                           env_arrow_col = "black",
                           env_label_cex = 0.82,
                           env_label_font = 2,
                           env_label_offset = 1.28,
                           env_label_jitter_y = 0.035,
                           legend_cex = 0.80,
                           site_legend_pos = "topright",
                           time_legend_pos = "bottomright",
                           ellipse_legend_pos = NULL,
                           site_legend_inset = 0,
                           time_legend_inset = 0,
                           show_ellipses = TRUE,
                           show_trajectories = FALSE,
                           show_points = TRUE,
                           show_envfit = TRUE,
                           show_legend = TRUE) {

  scr <- obj$scores

  # Plot limits include ellipses and envfit arrows/labels by default
  # so polygons/arrows/text are not clipped.
  x <- scr$NMDS1
  y <- scr$NMDS2

  # Include ellipse coordinates.
  if (show_ellipses && length(obj$ellipse_list) > 0) {
    ell_all <- do.call(rbind, obj$ellipse_list)
    x <- c(x, ell_all[, 1])
    y <- c(y, ell_all[, 2])
  }

  # Include envfit arrow tips and label positions.
  if (show_envfit && !is.null(obj$env_arrows) && nrow(obj$env_arrows) > 0) {
    ef_lim <- obj$env_arrows
    jitter_y_lim <- seq(-env_label_jitter_y, env_label_jitter_y, length.out = nrow(ef_lim))

    x <- c(x, 0, ef_lim$NMDS1s, ef_lim$NMDS1s * env_label_offset)
    y <- c(y, 0, ef_lim$NMDS2s, ef_lim$NMDS2s * env_label_offset + jitter_y_lim)
  }

  if (is.null(xlim)) {
    xr <- range(x, na.rm = TRUE)
    xpad <- diff(xr) * expand
    if (xpad == 0) xpad <- 0.5
    xlim <- c(xr[1] - xpad, xr[2] + xpad)
  }

  if (is.null(ylim)) {
    yr <- range(y, na.rm = TRUE)
    ypad <- diff(yr) * expand
    if (ypad == 0) ypad <- 0.5
    ylim <- c(yr[1] - ypad, yr[2] + ypad)
  }

  par(mar = mar, xpd = NA)

  plot(
    scr$NMDS1,
    scr$NMDS2,
    type = "n",
    xlim = xlim,
    ylim = ylim,
    xlab = "NMDS1",
    ylab = "NMDS2",
    main = main_title
    #main = paste0(main_title, "\nStress = ", round(obj$nmds$stress, 3))
  )

  # Colored ellipses
  if (show_ellipses && length(obj$ellipse_list) > 0) {
    for (g in names(obj$ellipse_list)) {
      ell <- obj$ellipse_list[[g]]
      polygon(
        ell[, 1],
        ell[, 2],
        col = adjustcolor(obj$ellipse_palette[g], alpha.f = ellipse_alpha),
        border = obj$ellipse_palette[g],
        lwd = ellipse_lwd
      )
    }
  }

  # Optional paired-colony trajectories
  if (show_trajectories && !is.null(obj$trajectories)) {
    with(
      obj$trajectories,
      arrows(
        x0, y0, x1, y1,
        col = adjustcolor(trajectory_col, alpha.f = trajectory_alpha),
        lwd = trajectory_lwd,
        length = trajectory_arrow_length
      )
    )
  }

  # Points
  if (show_points) {
    points(
      scr$NMDS1,
      scr$NMDS2,
      col = obj$point_cols,
      pch = obj$point_pchs,
      cex = point_cex,
      lwd = point_lwd
    )
  }

  # Envfit arrows and labels
  if (show_envfit && !is.null(obj$env_arrows) && nrow(obj$env_arrows) > 0) {
    ef <- obj$env_arrows

    arrows(
      0, 0,
      ef$NMDS1s,
      ef$NMDS2s,
      length = env_arrow_length,
      lwd = env_arrow_lwd,
      col = env_arrow_col
    )

    jitter_y <- seq(-env_label_jitter_y, env_label_jitter_y, length.out = nrow(ef))

    text(
      ef$NMDS1s * env_label_offset,
      ef$NMDS2s * env_label_offset + jitter_y,
      labels = ef$label,
      cex = env_label_cex,
      font = env_label_font
    )
  }

  # Legends
  if (show_legend) {
    legend(
      site_legend_pos,
      legend = names(obj$color_palette),
      col = obj$color_palette,
      pch = 16,
      pt.cex = 1.1,
      title = "Site",
      bty = "n",
      cex = legend_cex,
      inset = site_legend_inset
    )

    legend(
      time_legend_pos,
      legend = names(obj$symbol_palette),
      pch = obj$symbol_palette,
      pt.cex = 1.1,
      title = "Time",
      bty = "n",
      cex = legend_cex,
      inset = time_legend_inset
    )

    if (!is.null(ellipse_legend_pos)) {
      
      ellipse_labels <- names(obj$ellipse_palette)
      ellipse_labels[ellipse_labels == "coral_field"] <- "coral"
      
      legend(
        ellipse_legend_pos,
        legend = ellipse_labels,
        fill = adjustcolor(obj$ellipse_palette, alpha.f = 0.30),
        border = obj$ellipse_palette,
        title = "Ellipse",
        bty = "n",
        cex = legend_cex
      )
    }
  }

  invisible(list(xlim = xlim, ylim = ylim))
}

# 7. datasets to spearate year comparison and seasonal comparison 

# All compartments: April 2022 vs April 2024
meta_all_april_22_24 <- meta_filt[meta_filt$month == "4" & meta_filt$year %in% c("2022", "2024"), , drop = FALSE]
meta_all_april_22_24$time_group <- factor(ifelse(meta_all_april_22_24$year == "2022", "2022", "2024"), levels = c("2022", "2024"))
dat_all_april_22_24 <- prepare_dataset(postMED_filt, meta_all_april_22_24)

# All compartments: April vs May vs June 2022
meta_all_season_2022 <- meta_filt[meta_filt$year == "2022" & meta_filt$month %in% c("4", "5", "6"), , drop = FALSE]
meta_all_season_2022$time_group <- factor(
  c("4" = "April 2022", "5" = "May 2022", "6" = "June 2022")[as.character(meta_all_season_2022$month)],
  levels = c("April 2022", "May 2022", "June 2022")
)
dat_all_season_2022 <- prepare_dataset(postMED_filt, meta_all_season_2022)

# Coral: April 2022 vs April 2024
meta_coral_april_22_24 <- meta_filt[meta_filt$sample_type == "coral_field" & meta_filt$month == "4" & meta_filt$year %in% c("2022", "2024"), , drop = FALSE]
meta_coral_april_22_24$time_group <- factor(ifelse(meta_coral_april_22_24$year == "2022", "2022", "2024"), levels = c("2022", "2024"))
dat_coral_april_22_24 <- prepare_dataset(postMED_filt, meta_coral_april_22_24)

# Sediment: April 2022 vs April 2024
meta_sed_april_22_24 <- meta_filt[meta_filt$sample_type == "sediment" & meta_filt$month == "4" & meta_filt$year %in% c("2022", "2024"), , drop = FALSE]
meta_sed_april_22_24$time_group <- factor(ifelse(meta_sed_april_22_24$year == "2022", "2022", "2024"), levels = c("2022", "2024"))
dat_sed_april_22_24 <- prepare_dataset(postMED_filt, meta_sed_april_22_24)

# Coral: April vs May vs June 2022
meta_coral_season_2022 <- meta_filt[meta_filt$sample_type == "coral_field" & meta_filt$year == "2022" & meta_filt$month %in% c("4", "5", "6"), , drop = FALSE]
meta_coral_season_2022$time_group <- factor(
  c("4" = "April 2022", "5" = "May 2022", "6" = "June 2022")[as.character(meta_coral_season_2022$month)],
  levels = c("April 2022", "May 2022", "June 2022")
)
dat_coral_season_2022 <- prepare_dataset(postMED_filt, meta_coral_season_2022)

# Sediment: April vs May vs June 2022
meta_sed_season_2022 <- meta_filt[meta_filt$sample_type == "sediment" & meta_filt$year == "2022" & meta_filt$month %in% c("4", "5", "6"), , drop = FALSE]
meta_sed_season_2022$time_group <- factor(
  c("4" = "April 2022", "5" = "May 2022", "6" = "June 2022")[as.character(meta_sed_season_2022$month)],
  levels = c("April 2022", "May 2022", "June 2022")
)
dat_sed_season_2022 <- prepare_dataset(postMED_filt, meta_sed_season_2022)

# Seawater: April vs May vs June 2022
meta_sw_season_2022 <- meta_filt[meta_filt$sample_type == "seawater" & meta_filt$year == "2022" & meta_filt$month %in% c("4", "5", "6"), , drop = FALSE]
meta_sw_season_2022$time_group <- factor(
  c("4" = "April 2022", "5" = "May 2022", "6" = "June 2022")[as.character(meta_sw_season_2022$month)],
  levels = c("April 2022", "May 2022", "June 2022")
)
dat_sw_season_2022 <- prepare_dataset(postMED_filt, meta_sw_season_2022)

# 8. Environmental variable sets

env_all <- c("thermal_cal_SST", "thermal_DHW", "thermal_BAA", "thermal_Hotspot", "NDTI_for_analysis")
env_coral <- c("thermal_cal_SST", "thermal_DHW", "thermal_BAA", "thermal_Hotspot", "NDTI_for_analysis", "density_count")
env_sed <- c("thermal_cal_SST", "thermal_DHW", "thermal_BAA", "thermal_Hotspot", "NDTI_for_analysis")
env_sw <- c("thermal_cal_SST", "thermal_DHW", "thermal_BAA", "thermal_Hotspot", "NDTI_for_analysis")

write.csv(meta_filt, file.path(tables_dir, "metadata_after_P5_1pct_environment.csv"), row.names = TRUE)
write.csv(postMED_filt, file.path(tables_dir, "postMED_after_P5_1pct.csv"), row.names = TRUE)

# 9. ANALYSIS SECTIONS

# 9A. ALL COMPARTMENTS: April 2022 vs April 2024

run_permanova(dat_all_april_22_24, "sample_type + collection_site + time_group", "all_compartments_April_2022_vs_2024")
run_betadisper_test(dat_all_april_22_24, "time_group", "all_compartments_April_2022_vs_2024_by_time")
run_betadisper_test(dat_all_april_22_24, "sample_type", "all_compartments_April_2022_vs_2024_by_compartment")
run_envfit_table(dat_all_april_22_24, env_all, "all_compartments_April_2022_vs_2024")
run_simper_summary(dat_all_april_22_24, "time_group", "all_compartments_April_2022_vs_2024")
run_indicator_summary(dat_all_april_22_24, "time_group", "all_compartments_April_2022_vs_2024")

# Manual NMDS object: All compartments April 2022 vs April 2024
nmds_all_april_22_24 <- make_manual_nmds_object(
  dat_all_april_22_24,
  color_var = "collection_site",
  symbol_var = "time_group",
  ellipse_var = "sample_type",
  env_vars = env_all,
  remove_plot_samples = extreme_sediment_samples,
  trajectory_var = NULL,
  trajectory_both_groups = NULL
)

# Plot preview
plot_nmds_base(
  nmds_all_april_22_24,
  main_title = "All compartments April 2022 vs April 2024",
  expand = 0.15,
  point_cex = 1.35,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = TRUE,
  env_label_offset = 1.25,
  env_label_jitter_y = 0.23,
  env_label_cex = 0.75,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)

# Save PNG
png(file.path(figures_dir, "NMDS_all_compartments_April_2022_vs_2024.png"), width = 3000, height =2200, res = 300)
plot_nmds_base(
  nmds_all_april_22_24,
  main_title = "All compartments April 2022 vs April 2024",
  expand = 0.15,
  point_cex = 1.35,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = TRUE,
  env_label_offset = 1.25,
  env_label_jitter_y = 0.23,
  env_label_cex = 0.75,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)
dev.off()

write.csv(nmds_all_april_22_24$scores, file.path(tables_dir, "NMDS_scores_all_compartments_April_2022_vs_2024.csv"), row.names = FALSE)

# 9B. ALL COMPARTMENTS: April vs May vs June 2022

run_permanova(dat_all_season_2022, "sample_type + collection_site + time_group", "all_compartments_seasonal_2022")
run_betadisper_test(dat_all_season_2022, "time_group", "all_compartments_seasonal_2022_by_time")
run_betadisper_test(dat_all_season_2022, "sample_type", "all_compartments_seasonal_2022_by_compartment")
run_envfit_table(dat_all_season_2022, env_all, "all_compartments_seasonal_2022")
run_simper_summary(dat_all_season_2022, "time_group", "all_compartments_seasonal_2022")
run_indicator_summary(dat_all_season_2022, "time_group", "all_compartments_seasonal_2022")

# Manual NMDS object: All compartments April May June 2022
nmds_all_season_2022 <- make_manual_nmds_object(
  dat_all_season_2022,
  color_var = "collection_site",
  symbol_var = "time_group",
  ellipse_var = "sample_type",
  env_vars = env_all,
  remove_plot_samples = extreme_sediment_samples,
  trajectory_var = NULL,
  trajectory_both_groups = NULL
)

# Plot preview
plot_nmds_base(
  nmds_all_season_2022,
  main_title = "All compartments April May June 2022",
  expand = 0.15,
  point_cex = 1.35,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = TRUE,
  env_label_offset = 1.55,
  env_label_jitter_y = 0.23,
  env_label_cex = 0.75,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)

# Save PNG
png(file.path(figures_dir, "NMDS_all_compartments_April_May_June_2022.png"), width = 3000, height = 2200, res = 300)
plot_nmds_base(
  nmds_all_season_2022,
  main_title = "All compartments April May June 2022",
  expand = 0.15,
  point_cex = 1.35,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = TRUE,
  env_label_offset = 1.55,
  env_label_jitter_y = 0.23,
  env_label_cex = 0.75,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)
dev.off()

write.csv(nmds_all_season_2022$scores, file.path(tables_dir, "NMDS_scores_all_compartments_April_May_June_2022.csv"), row.names = FALSE)

# 9C. CORAL: April 2022 vs April 2024

run_permanova(dat_coral_april_22_24, "collection_site + time_group", "coral_April_2022_vs_2024")
run_permanova(dat_coral_april_22_24, "collection_site + time_group + thermal_cal_SST + NDTI_for_analysis + density_count", "coral_April_2022_vs_2024_calSST")
run_permanova(dat_coral_april_22_24, "collection_site + time_group + thermal_DHW + NDTI_for_analysis + density_count", "coral_April_2022_vs_2024_DHW")
run_permanova(dat_coral_april_22_24, "collection_site + time_group + thermal_BAA + NDTI_for_analysis + density_count", "coral_April_2022_vs_2024_BAA")
run_permanova(dat_coral_april_22_24, "collection_site + time_group + thermal_Hotspot + NDTI_for_analysis + density_count", "coral_April_2022_vs_2024_Hotspot")
run_betadisper_test(dat_coral_april_22_24, "time_group", "coral_April_2022_vs_2024_by_time")
run_envfit_table(dat_coral_april_22_24, env_coral, "coral_April_2022_vs_2024")
run_simper_summary(dat_coral_april_22_24, "time_group", "coral_April_2022_vs_2024")
run_indicator_summary(dat_coral_april_22_24, "time_group", "coral_April_2022_vs_2024")

# Manual NMDS object: Coral April 2022 vs April 2024
nmds_coral_april_22_24 <- make_manual_nmds_object(
  dat_coral_april_22_24,
  color_var = "collection_site",
  symbol_var = "time_group",
  ellipse_var = "time_group",
  env_vars = env_coral,
  remove_plot_samples = character(0),
  trajectory_var = NULL,
  trajectory_both_groups = NULL
)

# Plot preview
plot_nmds_base(
  nmds_coral_april_22_24,
  main_title = "Coral April 2022 vs April 2024",
  expand = 0.15,
  point_cex = 1.35,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = TRUE,
  env_label_offset = 1.15,
  env_label_jitter_y = 0.25,
  env_label_cex = 0.75,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)

# Save PNG
png(file.path(figures_dir, "NMDS_coral_April_2022_vs_2024.png"), width = 3000, height = 2200, res = 300)
plot_nmds_base(
  nmds_coral_april_22_24,
  main_title = "Coral April 2022 vs April 2024",
  expand = 0.15,
  point_cex = 1.35,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = TRUE,
  env_label_offset = 1.15,
  env_label_jitter_y = 0.25,
  env_label_cex = 0.75,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)
dev.off()

write.csv(nmds_coral_april_22_24$scores, file.path(tables_dir, "NMDS_scores_coral_April_2022_vs_2024.csv"), row.names = FALSE)


# Coral paired-colony trajectory plot
colony_var <- "colony_number"
dat_coral_april_22_24$meta$paired_colony_key <- paste(
  dat_coral_april_22_24$meta$collection_site,
  dat_coral_april_22_24$meta[[colony_var]],
  sep = "__"
)
paired_tab <- table(dat_coral_april_22_24$meta$paired_colony_key,
                    dat_coral_april_22_24$meta$time_group)
paired_ids <- rownames(paired_tab)[paired_tab[, "2022"] > 0 & paired_tab[, "2024"] > 0]
meta_coral_paired <- dat_coral_april_22_24$meta[
  dat_coral_april_22_24$meta$paired_colony_key %in% paired_ids,
  , drop = FALSE
]
dat_coral_paired_april <- prepare_dataset(dat_coral_april_22_24$comm, meta_coral_paired)

# Manual NMDS object: Coral paired colonies April 2022 to April 2024
nmds_coral_paired_april_22_24 <- make_manual_nmds_object(
  dat_coral_paired_april,
  color_var = "collection_site",
  symbol_var = "time_group",
  ellipse_var = "time_group",
  env_vars = env_coral,
  remove_plot_samples = character(0),
  trajectory_var = "paired_colony_key",
  trajectory_both_groups = c("2022", "2024")
)

# Plot preview
plot_nmds_base(
  nmds_coral_paired_april_22_24,
  main_title = "Coral paired colonies April 2022 to April 2024",
  expand = 0.15,
  point_cex = 1.35,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = FALSE,
  show_trajectories = TRUE,
  trajectory_lwd = 1.3,
  trajectory_alpha = 0.35,
  trajectory_arrow_length = 0.045,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)

# Save PNG
png(file.path(figures_dir, "NMDS_coral_paired_colonies_April_2022_to_2024.png"), width = 3000, height = 2200, res = 300)
plot_nmds_base(
  nmds_coral_paired_april_22_24,
  main_title = "Coral paired colonies April 2022 to April 2024",
  expand = 0.15,
  point_cex = 1.35,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = FALSE,
  show_trajectories = TRUE,
  trajectory_lwd = 1.3,
  trajectory_alpha = 0.35,
  trajectory_arrow_length = 0.045,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)
dev.off()

write.csv(nmds_coral_paired_april_22_24$scores, file.path(tables_dir, "NMDS_scores_coral_paired_colonies_April_2022_to_2024.csv"), row.names = FALSE)

# 9D. SEDIMENT: April 2022 vs April 2024


run_permanova(dat_sed_april_22_24, "collection_site + time_group", "sediment_April_2022_vs_2024")
run_permanova(dat_sed_april_22_24, "collection_site + time_group + thermal_cal_SST + NDTI_for_analysis", "sediment_April_2022_vs_2024_calSST")
run_permanova(dat_sed_april_22_24, "collection_site + time_group + thermal_DHW + NDTI_for_analysis", "sediment_April_2022_vs_2024_DHW")
run_permanova(dat_sed_april_22_24, "collection_site + time_group + thermal_BAA + NDTI_for_analysis", "sediment_April_2022_vs_2024_BAA")
run_permanova(dat_sed_april_22_24, "collection_site + time_group + thermal_Hotspot + NDTI_for_analysis", "sediment_April_2022_vs_2024_Hotspot")
run_betadisper_test(dat_sed_april_22_24, "time_group", "sediment_April_2022_vs_2024_by_time")
run_envfit_table(dat_sed_april_22_24, env_sed, "sediment_April_2022_vs_2024")
run_simper_summary(dat_sed_april_22_24, "time_group", "sediment_April_2022_vs_2024")
run_indicator_summary(dat_sed_april_22_24, "time_group", "sediment_April_2022_vs_2024")

# Manual NMDS object: Sediment April 2022 vs April 2024
nmds_sed_april_22_24 <- make_manual_nmds_object(
  dat_sed_april_22_24,
  color_var = "collection_site",
  symbol_var = "time_group",
  ellipse_var = "time_group",
  env_vars = env_sed,
  remove_plot_samples = extreme_sediment_samples,
  trajectory_var = NULL,
  trajectory_both_groups = NULL
)

# Plot preview
plot_nmds_base(
  nmds_sed_april_22_24,
  main_title = "Sediment April 2022 vs April 2024",
  expand = 0.15,
  point_cex = 1.35,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = TRUE,
  env_label_offset = 1.6,
  env_label_jitter_y = 0.33,
  env_label_cex = 0.75,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)

# Save PNG
png(file.path(figures_dir, "NMDS_sediment_April_2022_vs_2024.png"), width = 3000, height = 2200, res = 300)
plot_nmds_base(
  nmds_sed_april_22_24,
  main_title = "Sediment April 2022 vs April 2024",
  expand = 0.15,
  point_cex = 1.35,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = TRUE,
  env_label_offset = 1.6,
  env_label_jitter_y = 0.33,
  env_label_cex = 0.75,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)
dev.off()

write.csv(nmds_sed_april_22_24$scores, file.path(tables_dir, "NMDS_scores_sediment_April_2022_vs_2024.csv"), row.names = FALSE)

# 9E. CORAL: April vs May vs June 2022

run_permanova(dat_coral_season_2022, "collection_site + time_group", "coral_seasonal_2022")
run_betadisper_test(dat_coral_season_2022, "time_group", "coral_seasonal_2022_by_time")
run_envfit_table(dat_coral_season_2022, env_coral, "coral_seasonal_2022")
run_simper_summary(dat_coral_season_2022, "time_group", "coral_seasonal_2022")
run_indicator_summary(dat_coral_season_2022, "time_group", "coral_seasonal_2022")

# Manual NMDS object: Coral April May June 2022
nmds_coral_season_2022 <- make_manual_nmds_object(
  dat_coral_season_2022,
  color_var = "collection_site",
  symbol_var = "time_group",
  ellipse_var = "time_group",
  env_vars = env_coral,
  remove_plot_samples = character(0),
  trajectory_var = NULL,
  trajectory_both_groups = NULL
)

# Plot preview
plot_nmds_base(
  nmds_coral_season_2022,
  main_title = "Coral April May June 2022",
  expand = 0.15,
  point_cex = 1.35,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = TRUE,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)

# Save PNG
png(file.path(figures_dir, "NMDS_coral_April_May_June_2022.png"), width = 3000, height = 2200, res = 300)
plot_nmds_base(
  nmds_coral_season_2022,
  main_title = "Coral April May June 2022",
  expand = 0.15,
  point_cex = 1.35,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = TRUE,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)
dev.off()

write.csv(nmds_coral_season_2022$scores, file.path(tables_dir, "NMDS_scores_coral_April_May_June_2022.csv"), row.names = FALSE)

# 9F. SEDIMENT: April vs May vs June 2022

run_permanova(dat_sed_season_2022, "collection_site + time_group", "sediment_seasonal_2022")
run_betadisper_test(dat_sed_season_2022, "time_group", "sediment_seasonal_2022_by_time")
run_envfit_table(dat_sed_season_2022, env_sed, "sediment_seasonal_2022")
run_simper_summary(dat_sed_season_2022, "time_group", "sediment_seasonal_2022")
run_indicator_summary(dat_sed_season_2022, "time_group", "sediment_seasonal_2022")

# Manual NMDS object: Sediment April May June 2022
nmds_sed_season_2022 <- make_manual_nmds_object(
  dat_sed_season_2022,
  color_var = "collection_site",
  symbol_var = "time_group",
  ellipse_var = "time_group",
  env_vars = env_sed,
  remove_plot_samples = extreme_sediment_samples,
  trajectory_var = NULL,
  trajectory_both_groups = NULL
)

# Plot preview
plot_nmds_base(
  nmds_sed_season_2022,
  main_title = "Sediment April May June 2022",
  expand = 0.15,
  point_cex = 1.35,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = TRUE,
  env_label_offset = 1.3,
  env_label_jitter_y = 0.43,
  env_label_cex = 0.75,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)

# Save PNG
png(file.path(figures_dir, "NMDS_sediment_April_May_June_2022.png"), width = 3000, height = 2200, res = 300)
plot_nmds_base(
  nmds_sed_season_2022,
  main_title = "Sediment April May June 2022",
  expand = 0.15,
  point_cex = 1.35,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = TRUE,
  env_label_offset = 1.3,
  env_label_jitter_y = 0.43,
  env_label_cex = 0.75,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)
dev.off()

write.csv(nmds_sed_season_2022$scores, file.path(tables_dir, "NMDS_scores_sediment_April_May_June_2022.csv"), row.names = FALSE)

# 9G. SEAWATER: April vs May vs June 2022

run_permanova(dat_sw_season_2022, "collection_site + time_group", "seawater_seasonal_2022")
run_betadisper_test(dat_sw_season_2022, "time_group", "seawater_seasonal_2022_by_time")
run_envfit_table(dat_sw_season_2022, env_sw, "seawater_seasonal_2022")
run_simper_summary(dat_sw_season_2022, "time_group", "seawater_seasonal_2022")
run_indicator_summary(dat_sw_season_2022, "time_group", "seawater_seasonal_2022")

# Manual NMDS object: Seawater April May June 2022
nmds_sw_season_2022 <- make_manual_nmds_object(
  dat_sw_season_2022,
  color_var = "collection_site",
  symbol_var = "time_group",
  ellipse_var = "time_group",
  env_vars = env_sw,
  remove_plot_samples = character(0),
  trajectory_var = NULL,
  trajectory_both_groups = NULL
)

# Plot preview
plot_nmds_base(
  nmds_sw_season_2022,
  main_title = "Seawater April May June 2022",
  expand = 0.15,
  point_cex = 1.35,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = TRUE,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)

# Save PNG
png(file.path(figures_dir, "NMDS_seawater_April_May_June_2022.png"), width = 3000, height = 2200, res = 300)
plot_nmds_base(
  nmds_sw_season_2022,
  main_title = "Seawater April May June 2022",
  expand = 0.15,
  point_cex = 1.35,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = TRUE,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)
dev.off()

write.csv(nmds_sw_season_2022$scores, file.path(tables_dir, "NMDS_scores_seawater_April_May_June_2022.csv"), row.names = FALSE)

# 10. Pairwise PERMANOVA among sites

run_pairwise_site_permanova <- function(dat,
                                        label,
                                        group_var = "collection_site",
                                        permutations = 999,
                                        p_adjust_method = "BH") {
  
  meta <- dat$meta
  comm <- dat$comm
  
  keep <- !is.na(meta[[group_var]])
  meta <- droplevels(meta[keep, , drop = FALSE])
  comm <- comm[keep, , drop = FALSE]
  
  groups <- droplevels(factor(meta[[group_var]]))
  
  if (length(levels(groups)) < 2) {
    message("Skipping pairwise site PERMANOVA for ", label, ": fewer than two groups.")
    return(NULL)
  }
  
  pair_list <- combn(levels(groups), 2, simplify = FALSE)
  
  out <- lapply(pair_list, function(pair) {
    
    idx <- groups %in% pair
    
    comm_sub <- comm[idx, , drop = FALSE]
    meta_sub <- droplevels(meta[idx, , drop = FALSE])
    
    f <- as.formula(paste("comm_sub ~", group_var))
    
    ad <- vegan::adonis2(
      f,
      data = meta_sub,
      method = "bray",
      permutations = permutations,
      by = "margin"
    )
    
    data.frame(
      analysis = label,
      group_1 = pair[1],
      group_2 = pair[2],
      Df = ad$Df[1],
      SumOfSqs = ad$SumOfSqs[1],
      R2 = ad$R2[1],
      F_value = ad$F[1],
      p_value = ad$`Pr(>F)`[1],
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, out)
  out$p_adjusted <- p.adjust(out$p_value, method = p_adjust_method)
  
  write.csv(
    out,
    paste0(
      file.path(tables_dir, "PAIRWISE_SITE_PERMANOVA_"),
      gsub("[^A-Za-z0-9]+", "_", label),
      ".csv"
    ),
    row.names = FALSE
  )
  
  print(out)
  out
}

pairwise_site_all_april_22_24 <- run_pairwise_site_permanova(
  dat_all_april_22_24,
  "all_compartments_April_2022_vs_2024"
)

pairwise_site_all_season_2022 <- run_pairwise_site_permanova(
  dat_all_season_2022,
  "all_compartments_seasonal_2022"
)

pairwise_site_coral_april_22_24 <- run_pairwise_site_permanova(
  dat_coral_april_22_24,
  "coral_April_2022_vs_2024"
)

pairwise_site_coral_season_2022 <- run_pairwise_site_permanova(
  dat_coral_season_2022,
  "coral_seasonal_2022"
)

pairwise_site_sed_april_22_24 <- run_pairwise_site_permanova(
  dat_sed_april_22_24,
  "sediment_April_2022_vs_2024"
)

pairwise_site_sed_season_2022 <- run_pairwise_site_permanova(
  dat_sed_season_2022,
  "sediment_seasonal_2022"
)

pairwise_site_sw_season_2022 <- run_pairwise_site_permanova(
  dat_sw_season_2022,
  "seawater_seasonal_2022"
)


# 10A. Pairwise site PERMANOVA heatmap

library(ggplot2)

pairwise_site <- read.csv(
  file.path(tables_dir, "PAIRWISE_SITE_PERMANOVA_coral_April_2022_vs_2024.csv"),
  stringsAsFactors = FALSE
)

pairwise_site$significance <- cut(
  pairwise_site$p_adjusted,
  breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
  labels = c("***", "**", "*", ".", "ns")
)

# Preview
ggplot(
  pairwise_site,
  aes(
    x = group_1,
    y = group_2,
    fill = R2
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  
  geom_text(
    aes(
      label = paste0(
        "p=", round(p_adjusted, 3),
        "\n",
        significance
      )
    ),
    size = 4.6,
    fontface = "bold"
  ) +
  
  scale_fill_gradientn(
    colours = c(
      "#fff5eb",
      "#fdae6b",
      "#e6550d",
      "#a63603"
    ),
    name = expression(R^2)
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    title = "Pairwise PERMANOVA among coral sites"
  ) +
  
  coord_fixed() +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 12
    ),
    
    axis.text.y = element_text(
      face = "bold",
      size = 12
    ),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 13
    ),
    
    legend.text = element_text(
      size = 11
    )
  )

# Save PNG
png(
  file.path(figures_dir, "PAIRWISE_PERMANOVA_heatmap_coral_sites.png"),
  width = 3000,
  height = 2200,
  res = 300
)

ggplot(
  pairwise_site,
  aes(
    x = group_1,
    y = group_2,
    fill = R2
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  
  geom_text(
    aes(
      label = paste0(
        "p=", round(p_adjusted, 3),
        "\n",
        significance
      )
    ),
    size = 5.2,
    fontface = "bold"
  ) +
  
  scale_fill_gradientn(
    colours = c(
      "#fff5eb",
      "#fdae6b",
      "#e6550d",
      "#a63603"
    ),
    name = expression(R^2)
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    title = ""
  ) +
  
  coord_fixed() +
  
  theme_bw(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 13
    ),
    
    axis.text.y = element_text(
      face = "bold",
      size = 13
    ),
    
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 17
    ),
    
    legend.title = element_text(
      face = "bold",
      size = 14
    ),
    
    legend.text = element_text(
      size = 12
    )
  )

dev.off()

# 10B. Within-site temporal PERMANOVA: coral April 2022 vs April 2024

run_within_site_temporal_permanova <- function(dat,
                                               label,
                                               site_var = "collection_site",
                                               time_var = "time_group",
                                               permutations = 999,
                                               p_adjust_method = "BH") {
  
  meta <- dat$meta
  comm <- dat$comm
  
  keep <- complete.cases(meta[, c(site_var, time_var), drop = FALSE])
  meta <- droplevels(meta[keep, , drop = FALSE])
  comm <- comm[keep, , drop = FALSE]
  
  sites <- levels(droplevels(factor(meta[[site_var]])))
  
  out <- lapply(sites, function(site_i) {
    
    idx <- meta[[site_var]] == site_i
    
    meta_site <- droplevels(meta[idx, , drop = FALSE])
    comm_site <- comm[idx, , drop = FALSE]
    
    if (nrow(meta_site) < 6 ||
        length(unique(meta_site[[time_var]])) < 2) {
      return(data.frame(
        analysis = label,
        site = site_i,
        n_samples = nrow(meta_site),
        Df = NA,
        SumOfSqs = NA,
        R2 = NA,
        F_value = NA,
        p_value = NA,
        stringsAsFactors = FALSE
      ))
    }
    
    f <- as.formula(paste("comm_site ~", time_var))
    
    ad <- vegan::adonis2(
      f,
      data = meta_site,
      method = "bray",
      permutations = permutations,
      by = "margin"
    )
    
    data.frame(
      analysis = label,
      site = site_i,
      n_samples = nrow(meta_site),
      Df = ad$Df[1],
      SumOfSqs = ad$SumOfSqs[1],
      R2 = ad$R2[1],
      F_value = ad$F[1],
      p_value = ad$`Pr(>F)`[1],
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, out)
  out$p_adjusted <- p.adjust(out$p_value, method = p_adjust_method)
  
  write.csv(
    out,
    paste0(
      file.path(tables_dir, "WITHIN_SITE_TEMPORAL_PERMANOVA_"),
      gsub("[^A-Za-z0-9]+", "_", label),
      ".csv"
    ),
    row.names = FALSE
  )
  
  print(out)
  out
}

within_site_coral_april_22_24 <- run_within_site_temporal_permanova(
  dat_coral_april_22_24,
  "coral_April_2022_vs_2024"
)

# 10C. Plot within-site temporal PERMANOVA results

within_site_coral_plot <- within_site_coral_april_22_24

within_site_coral_plot$significance <- cut(
  within_site_coral_plot$p_adjusted,
  breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
  labels = c("***", "**", "*", ".", "ns")
)

within_site_coral_plot$site <- factor(
  within_site_coral_plot$site,
  levels = within_site_coral_plot$site[order(within_site_coral_plot$R2)]
)

p_within_site_coral <- ggplot(
  within_site_coral_plot,
  aes(
    x = R2,
    y = site,
    fill = p_adjusted <= 0.05
  )
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(
      label = paste0(
        "p=", round(p_adjusted, 3),
        " ",
        significance
      )
    ),
    hjust = -0.05,
    size = 4.2,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = "#D55E00",
      "FALSE" = "grey70"
    ),
    labels = c(
      "TRUE" = "Adjusted p ≤ 0.05",
      "FALSE" = "Not significant"
    ),
    name = NULL
  ) +
  labs(
    x = expression(PERMANOVA~R^2),
    y = NULL,
    title = "Within-site coral community change: April 2022 vs April 2024"
  ) +
  coord_cartesian(
    xlim = c(0, max(within_site_coral_plot$R2, na.rm = TRUE) * 1.35)
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )

# Plot preview
p_within_site_coral

# Save PNG
png(
  file.path(figures_dir, "WITHIN_SITE_TEMPORAL_PERMANOVA_coral_April_2022_vs_2024.png"),
  width = 2400,
  height = 1600,
  res = 300
)

print(p_within_site_coral)

dev.off()

# 10D. Combined PERMANOVA figure
# (a) Pairwise site heatmap
# (b) Within-site temporal coral change

library(patchwork)

# (a) Pairwise site PERMANOVA heatmap

p_pairwise_heatmap <- ggplot(
  pairwise_site,
  aes(
    x = group_1,
    y = group_2,
    fill = R2
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  
  geom_text(
    aes(
      label = paste0(
        "p=", round(p_adjusted, 3),
        "\n",
        significance
      )
    ),
    size = 4.5,
    fontface = "bold"
  ) +
  
  scale_fill_gradientn(
    colours = c(
      "#fff5eb",
      "#fdae6b",
      "#e6550d",
      "#a63603"
    ),
    name = expression(R^2)
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  coord_fixed(clip = "off") +
  
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "(a)",
    hjust = 1.2,
    vjust = 1.5,
    fontface = "bold",
    size = 6
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    aspect.ratio = 1,
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 20
    ),
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 12
    ),
    
    axis.text.y = element_text(
      face = "bold",
      size = 12
    ),
    
    plot.title = element_blank(),
    
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    
    legend.text = element_text(size = 11)
  )

# (b) Within-site temporal PERMANOVA barplot

p_within_site_combined <- ggplot(
  within_site_coral_plot,
  aes(
    x = R2,
    y = site,
    fill = p_adjusted <= 0.05
  )
) +
  
  geom_col(width = 0.7) +
  
  geom_text(
    aes(
      label = paste0(
        "p=", round(p_adjusted, 3),
        " ",
        significance
      )
    ),
    hjust = -0.05,
    size = 4.2,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c(
      "TRUE" = "#D55E00",
      "FALSE" = "grey70"
    ),
    guide = "none"
  ) +
  
  labs(
    x = expression(PERMANOVA~R^2),
    y = NULL
  ) +
  
  coord_cartesian(
    xlim = c(
      0,
      max(within_site_coral_plot$R2, na.rm = TRUE) * 1.35
    ),
    clip = "off"
  ) +
  
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "(b)",
    hjust = 1.2,
    vjust = 1.5,
    fontface = "bold",
    size = 6
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    aspect.ratio = 1,
    panel.grid.major.y = element_blank(),
    
    axis.text.y = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.text.x = element_text(
      size = 11
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 12
    ),
    
    plot.title = element_blank()
  )

# 10E. Combine panels with equal square plots

combined_permanova_square <- p_pairwise_heatmap +
  p_within_site_combined +
  plot_layout(
    ncol = 2,
    widths = c(1, 1)
  )

# Plot preview
combined_permanova_square

# Save PNG

png(
  file.path(figures_dir, "Combined_pairwise_and_within_site_PERMANOVA_square.png"),
  width = 3600,
  height = 1800,
  res = 300
)

print(combined_permanova_square)

dev.off()

# 11. SIMPER diverging barplot: coral April 2022 vs April 2024

simper_coral_file <- file.path(tables_dir, "SIMPER_coral_April_2022_vs_2024.csv")

simper_coral <- read.csv(simper_coral_file, stringsAsFactors = FALSE)

prepare_simper_diverging <- function(simper_df, top_n = 12) {
  
  simper_df <- simper_df[order(simper_df$average, decreasing = TRUE), ]
  simper_df <- head(simper_df, top_n)
  
  simper_df$higher_in <- ifelse(
    simper_df$ava > simper_df$avb,
    "Higher in 2024",
    "Higher in 2022"
  )
  
  simper_df$effect <- ifelse(
    simper_df$higher_in == "Higher in 2024",
    simper_df$average,
    -simper_df$average
  )
  
  simper_df <- simper_df[order(simper_df$effect), ]
  simper_df
}

simper_coral_plot <- prepare_simper_diverging(simper_coral, top_n = 12)

plot_simper_diverging_base <- function(simper_df,
                                       main_title = "",
                                       xlab = "SIMPER contribution to Bray-Curtis dissimilarity",
                                       cex_names = 0.85,
                                       cex_axis = 0.85,
                                       cex_lab = 0.95,
                                       legend_cex = 0.80) {
  
  cols <- ifelse(
    simper_df$effect > 0,
    "#D55E00",
    "#0072B2"
  )
  
  max_abs <- max(abs(simper_df$effect), na.rm = TRUE)
  
  par(mar = c(5, 7, 3, 2))
  
  barplot(
    simper_df$effect,
    horiz = TRUE,
    names.arg = simper_df$taxon,
    col = cols,
    border = NA,
    las = 1,
    xlim = c(-max_abs * 1.25, max_abs * 1.25),
    xlab = xlab,
    main = main_title,
    cex.names = cex_names,
    cex.axis = cex_axis,
    cex.lab = cex_lab
  )
  
  #abline(v = 0, lty = 2, lwd = 1.2)
  
  legend(
    "bottomright",
    legend = c("Higher in 2022", "Higher in 2024"),
    fill = c("#0072B2", "#D55E00"),
    bty = "n",
    cex = legend_cex
  )
}

# Plot preview 
plot_simper_diverging_base(simper_coral_plot)

# Save SIMPER plot
png(
  file.path(figures_dir, "SIMPER_diverging_coral_April_2022_vs_2024.png"),
  width = 2200,
  height = 1800,
  res = 300
)
plot_simper_diverging_base(simper_coral_plot)
dev.off()

# 12. Coral NMDS with site ellipses + combined coral figure

# Coral NMDS with site-based ellipses
nmds_coral_april_22_24_site_ellipses <- make_manual_nmds_object(
  dat_coral_april_22_24,
  color_var = "collection_site",
  symbol_var = "time_group",
  ellipse_var = "collection_site",
  env_vars = env_coral,
  remove_plot_samples = character(0),
  trajectory_var = NULL,
  trajectory_both_groups = NULL
)

nmds_coral_april_22_24_site_ellipses$ellipse_palette <-
  nmds_coral_april_22_24_site_ellipses$color_palette[
    names(nmds_coral_april_22_24_site_ellipses$ellipse_palette)
  ]


# Save individual site-ellipse NMDS
png(
  file.path(figures_dir, "NMDS_coral_April_2022_vs_2024_site_ellipses.png"),
  width = 3000,
  height = 2200,
  res = 300
)

plot_nmds_base(
  nmds_coral_april_22_24_site_ellipses,
  main_title = "",
  expand = 0.15,
  point_cex = 1.15,
  point_lwd = 2,
  ellipse_alpha = 0.15,
  ellipse_lwd = 2.2,
  show_envfit = TRUE,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = NULL,
  legend_cex = 0.75
)

dev.off()

# 13. Combined coral NMDS + SIMPER figure

add_panel_label <- function(label,
                            x_frac = 0.03,
                            y_frac = 0.97,
                            cex = 1.4,
                            font = 2) {
  usr <- par("usr")
  x <- usr[1] + x_frac * diff(usr[1:2])
  y <- usr[3] + y_frac * diff(usr[3:4])
  text(x, y, labels = label, cex = cex, font = font, xpd = TRUE)
}


png(
  file.path(figures_dir, "NMDS_coral_combined_envfit_year_site_SIMPER.png"),
  width = 3400,
  height = 3000,
  res = 300
)

par(mfrow = c(2, 2))

# (a) Coral NMDS with envfit only
plot_nmds_base(
  nmds_coral_april_22_24,
  main_title = "",
  expand = 0.15,
  point_cex = 1.05,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_ellipses = FALSE,
  show_envfit = TRUE,
  env_label_offset = 1.1,
  env_label_jitter_y = 0.43,
  env_label_cex = 0.75,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = NULL,
  legend_cex = 1.0
)
add_panel_label("(a)", cex = 1.5)

# (b) Coral NMDS with April 2022 vs April 2024 ellipses
plot_nmds_base(
  nmds_coral_paired_april_22_24,
  main_title = "",
  expand = 0.15,
  point_cex = 1.05,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_ellipses = TRUE,
  show_envfit = FALSE,
  show_trajectories = TRUE,
  trajectory_lwd = 1.3,
  trajectory_alpha = 0.35,
  trajectory_arrow_length = 0.045,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)
add_panel_label("(b)", cex = 1.5)

# (c) Coral NMDS with site ellipses
plot_nmds_base(
  nmds_coral_april_22_24_site_ellipses,
  main_title = "",
  expand = 0.15,
  point_cex = 1.05,
  point_lwd = 2,
  ellipse_alpha = 0.15,
  ellipse_lwd = 2.2,
  show_ellipses = TRUE,
  show_envfit = FALSE,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)
add_panel_label("(c)", cex = 1.5)

# (d) SIMPER diverging barplot
plot_simper_diverging_base(
  simper_coral_plot,
  main_title = "",
  cex_names = 1.0,
  legend_cex = 1.0
)
add_panel_label("(d)", cex = 1.5)

par(mfrow = c(1, 1))

dev.off()

# 14. Combined compartment NMDS figure
# seasonality and April bleaching comparison by compartment

png(
  file.path(figures_dir, "NMDS_combined_compartments_seasonality_bleaching_3x2.png"),
  width = 4800,
  height = 3200,
  res = 300
)

par(mfrow = c(2, 3), mar = c(5, 5, 3, 2))

# (a) Coral seasonality
plot_nmds_base(
  nmds_coral_season_2022,
  main_title = "",
  expand = 0.15,
  point_cex = 1.05,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_ellipses = TRUE,
  show_envfit = FALSE,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.2
)
add_panel_label("(a)", cex = 1.5)

# (b) Sediment seasonality
plot_nmds_base(
  nmds_sed_season_2022,
  main_title = "",
  expand = 0.15,
  point_cex = 1.05,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_ellipses = TRUE,
  show_envfit = FALSE,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.2
)
add_panel_label("(b)", cex = 1.5)

# (c) Seawater seasonality
plot_nmds_base(
  nmds_sw_season_2022,
  main_title = "",
  expand = 0.15,
  point_cex = 1.05,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_ellipses = TRUE,
  show_envfit = FALSE,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.2
)
add_panel_label("(c)", cex = 1.5)

# (d) Coral April 2022 vs April 2024
plot_nmds_base(
  nmds_coral_april_22_24,
  main_title = "",
  expand = 0.15,
  point_cex = 1.05,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_ellipses = TRUE,
  show_envfit = FALSE,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.2
)
add_panel_label("(d)", cex = 1.5)

# (e) Sediment April 2022 vs April 2024
plot_nmds_base(
  nmds_sed_april_22_24,
  main_title = "",
  expand = 0.15,
  point_cex = 1.05,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_ellipses = TRUE,
  show_envfit = FALSE,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.2
)
add_panel_label("(e)", cex = 1.5)

# (f) Seawater April 2022 vs April 2024: no data
plot.new()
text(
  x = 0.5,
  y = 0.5,
  labels = "No Data",
  cex = 2,
  font = 2
)
text(
  x = 0.05,
  y = 0.95,
  labels = "(f)",
  cex = 1.5,
  font = 2,
  adj = c(0, 1),
  xpd = TRUE
)

par(mfrow = c(1, 1))

dev.off()

# 15. Summary tables

read_all_csv <- function(pattern) {
  
  files <- list.files(
    tables_dir,
    pattern = pattern,
    full.names = TRUE
  )
  
  if (length(files) == 0) return(NULL)
  
  out <- lapply(files, function(f) {
    
    x <- read.csv(
      f,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    x$source_file <- basename(f)
    
    x
  })
  
  dplyr::bind_rows(out)
}

clean_analysis_name <- function(x) {
  x <- gsub("^PERMANOVA_|^ENVFIT_|^SIMPER_|^INDICATOR_", "", x)
  x <- gsub("\\.csv$", "", x)
  x <- gsub("_", " ", x)
  x
}

format_p <- function(p) {
  ifelse(
    is.na(p), NA,
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

# 15A. TABLE 1. PERMANOVA summary table

permanova_all <- read_all_csv("^PERMANOVA_.*\\.csv$")

if (!is.null(permanova_all)) {
  
  table1_permanova <- permanova_all
  
  table1_permanova$analysis_clean <- clean_analysis_name(table1_permanova$source_file)
  
  # Keep only real model terms, remove residual/total rows
  table1_permanova <- table1_permanova[
    !tolower(table1_permanova$term) %in% c("residual", "residuals", "total"),
    ,
    drop = FALSE
  ]
  
  # Standardize column names depending on adonis2 output
  if ("R2" %in% names(table1_permanova)) {
    table1_permanova$R2_value <- table1_permanova$R2
  } else {
    table1_permanova$R2_value <- NA
  }
  
  if ("F" %in% names(table1_permanova)) {
    table1_permanova$F_value <- table1_permanova$F
  } else {
    table1_permanova$F_value <- NA
  }
  
  # Detect PERMANOVA p-value column
  p_col <- grep("Pr", names(table1_permanova), value = TRUE)
  
  if (length(p_col) > 0) {
    table1_permanova$p_value <- table1_permanova[[p_col[1]]]
  } else {
    table1_permanova$p_value <- NA
  }
  
  table1_permanova_out <- data.frame(
    Analysis = table1_permanova$analysis_clean,
    Term = table1_permanova$term,
    Df = table1_permanova$Df,
    SumOfSqs = round(table1_permanova$SumOfSqs, 4),
    R2 = round(table1_permanova$R2_value, 4),
    F_value = round(table1_permanova$F_value, 3),
    p_value = format_p(table1_permanova$p_value),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    table1_permanova_out,
    file.path(main_tables_dir, "Table_1_PERMANOVA_summary.csv"),
    row.names = FALSE
  )
}

# 15B. TABLE 2. Envfit summary table, all variables

envfit_all <- read_all_csv("^ENVFIT_.*\\.csv$")

if (!is.null(envfit_all)) {
  
  envfit_all$analysis_clean <- clean_analysis_name(envfit_all$source_file)
  
  table2_envfit <- data.frame(
    Analysis = envfit_all$analysis_clean,
    Variable = ifelse(!is.na(envfit_all$label) & envfit_all$label != "",
                      envfit_all$label,
                      envfit_all$variable),
    NMDS1 = round(envfit_all$NMDS1, 4),
    NMDS2 = round(envfit_all$NMDS2, 4),
    r2 = round(envfit_all$r2, 4),
    p_value = format_p(envfit_all$p_value),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    table2_envfit,
    file.path(main_tables_dir, "Table_2_ENVFIT_all_variables.csv"),
    row.names = FALSE
  )
}

# 15C. TABLE 3. SIMPER top 20 contributors

simper_all <- read_all_csv("^SIMPER_.*\\.csv$")

if (!is.null(simper_all)) {
  
  simper_all$analysis_clean <- clean_analysis_name(simper_all$source_file)
  
  # Add direction when possible
  simper_all$higher_in <- ifelse(
    simper_all$ava > simper_all$avb,
    sub("_.*", "", simper_all$comparison),
    sub(".*_", "", simper_all$comparison)
  )
  
  simper_split <- split(simper_all, simper_all$analysis_clean)
  
  table3_simper <- do.call(rbind, lapply(simper_split, function(x) {
    x <- x[order(x$average, decreasing = TRUE), , drop = FALSE]
    head(x, 20)
  }))
  
  table3_simper_out <- data.frame(
    Analysis = table3_simper$analysis_clean,
    Comparison = table3_simper$comparison,
    Taxon = table3_simper$taxon,
    Average_contribution = round(table3_simper$average, 5),
    Cumulative_contribution = round(table3_simper$cusum, 4),
    Mean_group_A = round(table3_simper$ava, 5),
    Mean_group_B = round(table3_simper$avb, 5),
    Higher_in = table3_simper$higher_in,
    p_value = format_p(table3_simper$p),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    table3_simper_out,
    file.path(main_tables_dir, "Table_3_SIMPER_top20.csv"),
    row.names = FALSE
  )
}

# 15D. TABLE 4. Indicator species summary

indicator_all <- read_all_csv("^INDICATOR_.*\\.csv$")

if (!is.null(indicator_all)) {
  
  indicator_all$analysis_clean <- clean_analysis_name(indicator_all$source_file)
  
  group_cols <- grep("^s\\.", names(indicator_all), value = TRUE)
  
  indicator_all$associated_group <- apply(
    indicator_all[, group_cols, drop = FALSE],
    1,
    function(z) {
      g <- group_cols[which(z == 1)]
      g <- gsub("^s\\.", "", g)
      paste(g, collapse = " + ")
    }
  )
  
  table4_indicator <- indicator_all[
    order(indicator_all$analysis_clean, indicator_all$p.value),
    ,
    drop = FALSE
  ]
  
  table4_indicator_out <- data.frame(
    Analysis = table4_indicator$analysis_clean,
    Taxon = table4_indicator$taxon,
    Associated_group = table4_indicator$associated_group,
    Indicator_value = round(table4_indicator$stat, 4),
    p_value = format_p(table4_indicator$p.value),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    table4_indicator_out,
    file.path(main_tables_dir, "Table_4_Indicator_species_summary.csv"),
    row.names = FALSE
  )
  
  # Significant indicators only
  table4_indicator_sig <- table4_indicator_out[
    as.numeric(gsub("<0.001", "0.0009", table4_indicator_out$p_value)) <= 0.05,
    ,
    drop = FALSE
  ]
  
  write.csv(
    table4_indicator_sig,
    file.path(main_tables_dir, "Table_4_Indicator_species_significant_only.csv"),
    row.names = FALSE
  )
}

# 15E. Supplementary tables

# Supplementary Table S3: PERMANOVA outputs used for the manuscript
if (!is.null(permanova_all)) {
  write.csv(
    permanova_all,
    file.path(supp_tables_dir, "Table_S3_PERMANOVA_results.csv"),
    row.names = FALSE
  )
}

# Supplementary Table S4: envfit results used for the manuscript
if (!is.null(envfit_all)) {
  write.csv(
    envfit_all,
    file.path(supp_tables_dir, "Table_S4_ENVFIT_results.csv"),
    row.names = FALSE
  )
}

# Supplementary Table S5: SIMPER results used for the manuscript
if (!is.null(simper_all)) {
  write.csv(
    simper_all,
    file.path(supp_tables_dir, "Table_S5_SIMPER_results.csv"),
    row.names = FALSE
  )
}

# Supplementary Table S6: indicator-species results used for the manuscript
if (!is.null(indicator_all)) {
  write.csv(
    indicator_all,
    file.path(supp_tables_dir, "Table_S6_indicator_species_results.csv"),
    row.names = FALSE
  )
}

# combined NMDS scores (not a manuscript supplementary table)
nmds_score_files <- list.files(
  tables_dir,
  pattern = "^NMDS_scores_.*\\.csv$",
  full.names = TRUE
)

if (length(nmds_score_files) > 0) {
  
  nmds_scores_all <- lapply(nmds_score_files, function(f) {
    
    x <- read.csv(
      f,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    x$source_file <- basename(f)
    
    # Convert all columns to character
    x[] <- lapply(x, as.character)
    
    x
  })
  
  nmds_scores_all <- dplyr::bind_rows(nmds_scores_all)
  
  write.csv(
    nmds_scores_all,
    file.path(diagnostics_dir, "NMDS_scores_all_analyses.csv"),
    row.names = FALSE
  )
}

# metadata after filtering (not a manuscript supplementary table)
metadata_file <- file.path(tables_dir, "metadata_after_P5_1pct_environment.csv")

if (file.exists(metadata_file)) {
  metadata_filtered <- read.csv(metadata_file, stringsAsFactors = FALSE)
  
  write.csv(
    metadata_filtered,
    file.path(diagnostics_dir, "metadata_after_filtering.csv"),
    row.names = FALSE
  )
}

cat("\nManuscript-ready tables written to:\n")
cat(" - ", main_tables_dir, "\n")

cat("\nSupplementary tables written to:\n")
cat(" - ", supp_tables_dir, "\n")

# 16. Combined all-compartment NMDS figure
#     April 2022 vs April 2024 + April, May, June 2022

add_panel_label <- function(label,
                            x_frac = 0.03,
                            y_frac = 0.97,
                            cex = 1.4,
                            font = 2) {
  usr <- par("usr")
  x <- usr[1] + x_frac * diff(usr[1:2])
  y <- usr[3] + y_frac * diff(usr[3:4])
  text(x, y, labels = label, cex = cex, font = font, xpd = TRUE)
}

# Plot preview
par(mfrow = c(1, 2), mar = c(5, 5, 3, 2))

plot_nmds_base(
  nmds_all_april_22_24,
  main_title = "",
  expand = 0.15,
  point_cex = 1.15,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = FALSE,
  env_label_offset = 1.3,
  env_label_jitter_y = 0.53,
  env_label_cex = 0.75,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)
add_panel_label("(a)")

plot_nmds_base(
  nmds_all_season_2022,
  main_title = "",
  expand = 0.15,
  point_cex = 1.15,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = FALSE,
  env_label_offset = 1.5,
  env_label_jitter_y = 0.38,
  env_label_cex = 0.75,
  show_trajectories = FALSE,
  site_legend_pos = "topright",
  time_legend_pos = "bottomright",
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.0
)
add_panel_label("(b)")

par(mfrow = c(1, 1))




# Save as PNG
png(
  file.path(figures_dir, "NMDS_combined_all_compartments_bleaching_and_seasonality.png"),
  width = 3600,
  height = 1700,
  res = 300
)

#par(mfrow = c(1, 2))

# Two plots with a narrow empty column between them
layout(
  matrix(c(1, 0, 2), nrow = 1),
  widths = c(1, 0.08, 1)
)

plot_nmds_base(
  nmds_all_april_22_24,
  main_title = "",
  expand = 0.15,
  mar = c(5, 5, 3, 1),
  point_cex = 1.05,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = FALSE,
  env_label_offset = 1.1,
  env_label_jitter_y = 0.53,
  env_label_cex = 0.75,
  show_trajectories = FALSE,
  
  # Both legends on right
  site_legend_pos = "bottomright",
  site_legend_inset = c(0.01, 0.13),
  time_legend_pos = "bottomright",
  time_legend_inset = c(0.045, 0),
  
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.2
)

add_panel_label("(a)", cex = 1.5)



plot_nmds_base(
  nmds_all_season_2022,
  main_title = "",
  expand = 0.15,
  mar = c(5, 1, 3, 2),
  point_cex = 1.05,
  point_lwd = 2,
  ellipse_alpha = 0.20,
  ellipse_lwd = 2.2,
  show_envfit = FALSE,
  env_label_offset = 1.5,
  env_label_jitter_y = 0.38,
  env_label_cex = 0.75,
  show_trajectories = FALSE,
  
  # Keep original Site position
  site_legend_pos = "topright",
  site_legend_inset = 0,
  
  time_legend_pos = "bottomright",
  time_legend_inset = 0,
  
  ellipse_legend_pos = "bottomleft",
  legend_cex = 1.2
)

add_panel_label("(b)", cex = 1.5)

#par(mfrow = c(1, 1))

dev.off()

# Save workspace object

saveRDS(
  list(
    meta = meta_filt,
    comm = postMED_filt,
    p5_threshold = th_p5,
    thermal_pc = thermal_pc,
    extreme_sediment_samples = extreme_sediment_samples
  ),
  file = file.path(rds_dir, "seasonality_bleaching_results_workspace.rds")
)

cat("\nCommunity/multivariate analysis complete. See outputs/ for figures, tables and diagnostics.\n")


