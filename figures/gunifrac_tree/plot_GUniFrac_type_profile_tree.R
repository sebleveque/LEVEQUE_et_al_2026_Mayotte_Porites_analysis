# GUniFrac ITS2 post-MED sample trees

input_dir <- file.path("data", "gunifrac")
figure_dir <- file.path("outputs", "figures")
diag_dir <- file.path("outputs", "gunifrac", "diagnostics")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

# 1. LOAD LIBRARIES

library(Biostrings)
library(ape)
library(kmer)
library(GUniFrac)
library(colorspace)
library(vegan)

# 2. INPUT FILES

its2_seqs <- readDNAStringSet(file.path(input_dir, "ITS2_postMED_seq.fasta"))

its2_tab <- read.table(
  file.path(input_dir, "ITS2_postMED_relative_abundances.txt"),
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

type_tab <- read.table(
  file.path(input_dir, "ITS2_type_profiles_relative_abundances.txt"),
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Use untransformed post-MED relative abundances; GUniFrac normalizes rows internally.
# Identify dominant profiles directly from their relative abundances.

# 3. CLEAN ITS2 TABLE

its2_tab <- its2_tab[, colSums(its2_tab) > 0, drop = FALSE]
its2_seqs <- its2_seqs[names(its2_seqs) %in% colnames(its2_tab)]
its2_tab <- its2_tab[, colnames(its2_tab) %in% names(its2_seqs), drop = FALSE]

cat("Number of samples in ITS2 table:", nrow(its2_tab), "\n")
cat("Number of ITS2 variants retained:", ncol(its2_tab), "\n")
cat("Number of ITS2 sequences retained:", length(its2_seqs), "\n\n")

# 4. k-mer ITS2 SEQUENCE TREE

build_kmer_tree <- function(dna_set, k = 7) {
  
  dna_char <- as.character(dna_set)
  names(dna_char) <- names(dna_set)
  
  max_len <- max(nchar(dna_char))
  
  dna_pad <- sapply(dna_char, function(x) {
    paste0(x, paste(rep("-", max_len - nchar(x)), collapse = ""))
  })
  
  dna_split <- do.call(rbind, strsplit(dna_pad, split = ""))
  rownames(dna_split) <- names(dna_char)
  
  dna_dnabin <- as.DNAbin(dna_split)
  
  dist_mat <- kdistance(
    dna_dnabin,
    k = k,
    method = "edgar",
    residues = "DNA"
  )
  
  dist_mat <- as.matrix(dist_mat)
  dist_mat[!is.finite(dist_mat)] <- 0
  diag(dist_mat) <- 0
  
  hc <- hclust(as.dist(dist_mat), method = "average")
  phy <- as.phylo(hc)
  
  phy$edge.length <- as.numeric(phy$edge.length)
  
  bad_edges <- is.na(phy$edge.length) |
    !is.finite(phy$edge.length) |
    phy$edge.length <= 0
  
  if (any(bad_edges)) {
    phy$edge.length[bad_edges] <- 1e-8
  }
  
  phy
}

its2_tree <- build_kmer_tree(its2_seqs, k = 7)

cat("ITS2 sequence tree built.\n")
cat("Tips:", Ntip(its2_tree), "\n")
cat("Nodes:", Nnode(its2_tree), "\n\n")

# 5. MATCH ITS2 TABLE AND TREE

common_otus <- intersect(colnames(its2_tab), its2_tree$tip.label)

its2_tab_guni <- its2_tab[, common_otus, drop = FALSE]

its2_tree_guni <- drop.tip(
  its2_tree,
  setdiff(its2_tree$tip.label, common_otus)
)

its2_tab_guni <- its2_tab_guni[, its2_tree_guni$tip.label, drop = FALSE]
its2_tab_guni <- data.matrix(its2_tab_guni)
its2_tab_guni <- its2_tab_guni[rowSums(its2_tab_guni) > 0, , drop = FALSE]

bad_edges <- is.na(its2_tree_guni$edge.length) |
  !is.finite(its2_tree_guni$edge.length) |
  its2_tree_guni$edge.length <= 0

if (any(bad_edges)) {
  its2_tree_guni$edge.length[bad_edges] <- 1e-8
}

cat("ITS2 OTU match check:",
    all(colnames(its2_tab_guni) == its2_tree_guni$tip.label), "\n")
cat("ITS2 OTUs retained for GUniFrac:", ncol(its2_tab_guni), "\n")
cat("Samples retained for GUniFrac:", nrow(its2_tab_guni), "\n\n")

# 6. GUniFrac CALCULATION

its2_guni <- GUniFrac(
  its2_tab_guni,
  its2_tree_guni,
  alpha = c(0, 0.5, 1)
)

its2_dist <- its2_guni$unifracs[, , "d_0.5"]

cat("GUniFrac distance matrix dimensions:\n")
print(dim(its2_dist))
cat("\n")

# 7. SAMPLE TYPE GROUPS

sample_names <- rownames(its2_tab_guni)

sample_group <- ifelse(grepl("^WS", sample_names), "WS",
                       ifelse(grepl("^WC", sample_names), "WC",
                              ifelse(grepl("^S", sample_names), "S", NA)))

keep_samples <- !is.na(sample_group)

its2_tab_guni <- its2_tab_guni[keep_samples, , drop = FALSE]
sample_group <- sample_group[keep_samples]
sample_names <- rownames(its2_tab_guni)

its2_dist <- its2_dist[sample_names, sample_names]

cat("Sample counts per group:\n")
print(table(sample_group))
cat("\n")

# 8. NEIGHBOR JOINING SAMPLE TREES


build_sample_tree <- function(dist_mat) {
  
  phy <- nj(as.dist(dist_mat))
  
  phy$edge.length <- as.numeric(phy$edge.length)
  
  bad_edges <- is.na(phy$edge.length) | !is.finite(phy$edge.length)
  
  if (any(bad_edges)) {
    phy$edge.length[bad_edges] <- 1e-8
  }
  
  phy
}

dist_S  <- its2_dist[sample_group == "S",  sample_group == "S",  drop = FALSE]
dist_WS <- its2_dist[sample_group == "WS", sample_group == "WS", drop = FALSE]
dist_WC <- its2_dist[sample_group == "WC", sample_group == "WC", drop = FALSE]

tree_S_nj  <- if (nrow(dist_S)  >= 2) build_sample_tree(dist_S)  else NULL
tree_WS_nj <- if (nrow(dist_WS) >= 2) build_sample_tree(dist_WS) else NULL
tree_WC_nj <- if (nrow(dist_WC) >= 2) build_sample_tree(dist_WC) else NULL

# 9. DOMINANT ITS2 TYPE PROFILE PER SAMPLE

common_type_samples <- intersect(
  rownames(type_tab),
  rownames(its2_tab_guni)
)

type_tab <- type_tab[common_type_samples, , drop = FALSE]
type_tab <- type_tab[, colSums(type_tab) > 0, drop = FALSE]
type_tab <- type_tab[rowSums(type_tab) > 0, , drop = FALSE]

dominant_type <- colnames(type_tab)[
  max.col(type_tab, ties.method = "first")
]

names(dominant_type) <- rownames(type_tab)

type_summary <- data.frame(
  Sample = names(dominant_type),
  Dominant_Type = dominant_type,
  stringsAsFactors = FALSE
)

write.table(
  type_summary,
  file = file.path(diag_dir, "Dominant_ITS2_type_profile_per_sample.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 10. TYPE PROFILE COLOR PALETTE

blue_grad   <- sequential_hcl(26, "Blues 3")
green_grad  <- sequential_hcl(10, "Greens 3")
purple_grad <- sequential_hcl(4, "Purples 3")
red_grad    <- sequential_hcl(4, "Reds 3")
grey_grad   <- sequential_hcl(11, "Light Grays")

gradient_list <- list(
  blue   = blue_grad,
  green  = green_grad,
  purple = purple_grad,
  red    = red_grad,
  grey   = grey_grad
)

unique_types <- sort(unique(dominant_type))

type_group <- rep("black", length(unique_types))

type_group[grepl("^C15adl-", unique_types)] <- "green"
type_group[grepl("^C15adl/", unique_types)] <- "purple"
type_group[grepl("^C15add", unique_types)]  <- "#BAAD50"
type_group[grepl("^C15ado", unique_types)]  <- "#A2B353"
type_group[grepl("^C15-", unique_types)]    <- "blue"
type_group[grepl("^C15/", unique_types)]    <- "red"

type_group[grepl("^C21", unique_types)]     <- "#00C0C1"
type_group[grepl("^C66", unique_types)]     <- "#CEA55C"
type_group[grepl("^C3aei", unique_types)]   <- "#7F3B08"
type_group[grepl("^C3aej", unique_types)]   <- "#B35806"
type_group[grepl("^C3u", unique_types)]     <- "#E08214"
type_group[grepl("^D1", unique_types)]      <- "grey"
type_group[grepl("^F4", unique_types)]      <- "#33C192"
type_group[grepl("^A2p", unique_types)]     <- "#87ACEB"
type_group[grepl("^C1", unique_types) & !grepl("^C15", unique_types)] <- "#13BCD4"

type_group_df <- data.frame(
  TypeProfile = unique_types,
  GradientGroup = type_group,
  stringsAsFactors = FALSE
)

profiles_by_group <- split(
  type_group_df$TypeProfile,
  type_group_df$GradientGroup
)

profile_colors <- c()

for (grp in names(profiles_by_group)) {
  
  these_profiles <- sort(profiles_by_group[[grp]])
  n_profiles <- length(these_profiles)
  
  if (grp %in% names(gradient_list)) {
    
    pal <- gradient_list[[grp]]
    
    if (n_profiles > length(pal)) {
      pal <- colorRampPalette(pal)(n_profiles)
    } else {
      pal <- pal[seq_len(n_profiles)]
    }
    
  } else {
    
    pal <- rep(grp, n_profiles)
  }
  
  names(pal) <- these_profiles
  profile_colors <- c(profile_colors, pal)
}

profile_color_table <- data.frame(
  TypeProfile = names(profile_colors),
  Color = unname(profile_colors),
  stringsAsFactors = FALSE
)

write.table(
  profile_color_table,
  file = file.path(diag_dir, "ITS2_type_profile_colors.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 11. PLOTTING FUNCTIONS

plot_fan_tree_equal_labels <- function(tree_obj,
                                       dominant_type,
                                       profile_colors,
                                       main_title = "",
                                       cex_tip = 0.8,
                                       tip_pch = 19,
                                       tip_cex = 1.2,
                                       label_cex = 0.5,
                                       radial_padding = 0.05,
                                       rotate_labels = TRUE,
                                       label_colored = TRUE,
                                       extra_space = 0.6,
                                       mar = c(0.5, 0.5, 2, 0.5)) {
  
  if (is.null(tree_obj)) {
    message("Tree is NULL")
    return(invisible(NULL))
  }
  
  tip_types <- dominant_type[
    match(tree_obj$tip.label, names(dominant_type))
  ]
  
  tip_cols <- profile_colors[tip_types]
  tip_cols[is.na(tip_cols)] <- "grey60"
  
  tree_obj <- ladderize(tree_obj)
  
  rt <- node.depth.edgelength(tree_obj)[1:Ntip(tree_obj)]
  max_depth <- max(rt, na.rm = TRUE)
  plot_limit <- max_depth + radial_padding + extra_space
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  par(xpd = NA, mar = mar)
  
  plot(
    tree_obj,
    type = "fan",
    use.edge.length = TRUE,
    tip.color = tip_cols,
    show.tip.label = FALSE,
    cex = cex_tip,
    no.margin = FALSE,
    x.lim = c(-plot_limit, plot_limit),
    y.lim = c(-plot_limit, plot_limit),
    main = main_title
  )
  
  tiplabels(
    pch = tip_pch,
    col = tip_cols,
    cex = tip_cex
  )
  
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
  tip_x <- lastPP$xx[1:Ntip(tree_obj)]
  tip_y <- lastPP$yy[1:Ntip(tree_obj)]
  
  tip_r <- sqrt(tip_x^2 + tip_y^2)
  tip_r[tip_r == 0] <- 1e-8
  
  max_r <- max(tip_r, na.rm = TRUE)
  label_r <- max_r + radial_padding
  
  x_lab <- tip_x / tip_r * label_r
  y_lab <- tip_y / tip_r * label_r
  
  angles <- atan2(y_lab, x_lab) * 180 / pi
  label_cols <- if (label_colored) tip_cols else "black"
  
  for (i in seq_along(tree_obj$tip.label)) {
    
    if (rotate_labels) {
      
      ang <- angles[i]
      
      if (ang < -90 || ang > 90) {
        srt_i <- ang + 180
        adj_i <- 1
      } else {
        srt_i <- ang
        adj_i <- 0
      }
      
      text(
        x_lab[i],
        y_lab[i],
        labels = tree_obj$tip.label[i],
        cex = label_cex,
        font = 2,
        srt = srt_i,
        adj = adj_i,
        col = label_cols[i]
      )
      
    } else {
      
      text(
        x_lab[i],
        y_lab[i],
        labels = tree_obj$tip.label[i],
        cex = label_cex,
        font = 2,
        adj = ifelse(x_lab[i] < 0, 1, 0),
        col = label_cols[i]
      )
    }
  }
  
  invisible(NULL)
}

shorten_type_label <- function(x, n = 3) {
  
  sapply(x, function(label) {
    
    parts <- unlist(
      strsplit(label, "(?<=[-/])|(?=[-/])", perl = TRUE)
    )
    
    non_sep <- parts[!parts %in% c("-", "/")]
    
    if (length(non_sep) <= n) {
      return(label)
    }
    
    keep_count <- 0
    out <- character()
    
    for (p in parts) {
      
      out <- c(out, p)
      
      if (!p %in% c("-", "/")) {
        keep_count <- keep_count + 1
      }
      
      if (keep_count == n) break
    }
    
    paste0(paste(out, collapse = ""), "...")
    
  }, USE.NAMES = FALSE)
}

draw_short_group_legend <- function(tree_obj,
                                    dominant_type,
                                    profile_colors,
                                    ncol = 6,
                                    cex = 0.8,
                                    pt.cex = 1.2) {
  
  if (is.null(tree_obj)) {
    message("Tree is NULL")
    return(invisible(NULL))
  }
  
  tip_types <- dominant_type[
    match(tree_obj$tip.label, names(dominant_type))
  ]
  
  legend_types <- sort(unique(tip_types))
  legend_types <- legend_types[!is.na(legend_types)]
  
  legend_cols <- profile_colors[legend_types]
  legend_labels <- shorten_type_label(legend_types)
  
  legend(
    "center",
    legend = legend_labels,
    col = legend_cols,
    pch = 19,
    pt.cex = pt.cex,
    cex = cex,
    ncol = ncol,
    bty = "n"
  )
  
  invisible(NULL)
}

plot_tree_with_short_legend <- function(tree_obj,
                                        dominant_type,
                                        profile_colors,
                                        legend_ncol = 6,
                                        legend_cex = 0.8,
                                        legend_pt_cex = 1.2,
                                        tree_panel = c(0, 1, 0.22, 1),
                                        legend_panel = c(0, 1, 0, 0.25),
                                        tree_mar = c(1, 1, 2, 1),
                                        legend_mar = c(0, 1, 0, 1),
                                        ...) {
  
  if (is.null(tree_obj)) {
    message("Tree is NULL")
    return(invisible(NULL))
  }
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  par(fig = tree_panel, mar = tree_mar, xpd = NA)
  
  plot_fan_tree_equal_labels(
    tree_obj,
    dominant_type,
    profile_colors,
    main_title = "",
    ...
  )
  
  par(fig = legend_panel, mar = legend_mar, xpd = NA, new = TRUE)
  plot.new()
  
  draw_short_group_legend(
    tree_obj,
    dominant_type,
    profile_colors,
    ncol = legend_ncol,
    cex = legend_cex,
    pt.cex = legend_pt_cex
  )
  
  invisible(NULL)
}

# 12. EXPORT FINAL TIFF FIGURES

# S tree (coral)
tiff(
  file.path(figure_dir, "Figure_S2a_coral_GUniFrac_tree.tiff"),
  width = 9000,
  height = 9500,
  res = 300,
  compression = "lzw"
)

plot_tree_with_short_legend(
  tree_S_nj,
  dominant_type,
  profile_colors,
  legend_ncol = 6,
  legend_cex = 1.75,
  tree_panel = c(0, 1, 0.10, 1),
  legend_panel = c(0, 1, 0, 0.12),
  tree_mar = c(1, 1, 2, 1),
  legend_mar = c(0, 0.5, 0.3, 0.5),
  cex_tip = 0.8,
  tip_cex = 1.5,
  label_cex = 1.5,
  radial_padding = 0.005,
  extra_space = 0.008,
  rotate_labels = TRUE
)

dev.off()

# WS tree (sediment)
tiff(
  file.path(figure_dir, "Figure_S2b_sediment_GUniFrac_tree.tiff"),
  width = 8500,
  height = 9000,
  res = 300,
  compression = "lzw"
)

plot_tree_with_short_legend(
  tree_WS_nj,
  dominant_type,
  profile_colors,
  legend_ncol = 7,
  legend_cex = 1.75,
  tree_panel = c(0, 1, 0.14, 1),
  legend_panel = c(0, 1, 0, 0.16),
  tree_mar = c(1, 1, 2, 1),
  legend_mar = c(0, 0.5, 0.3, 0.5),
  cex_tip = 0.8,
  tip_cex = 1.5,
  label_cex = 1.5,
  radial_padding = 0.04,
  extra_space = 0.10,
  rotate_labels = TRUE
)

dev.off()

# WC tree (seawater)
tiff(
  file.path(figure_dir, "Figure_S2c_seawater_GUniFrac_tree.tiff"),
  width = 6000,
  height = 6500,
  res = 300,
  compression = "lzw"
)

plot_tree_with_short_legend(
  tree_WC_nj,
  dominant_type,
  profile_colors,
  legend_ncol = 4,
  legend_cex = 1.75,
  tree_panel = c(0, 1, 0.18, 1),
  legend_panel = c(0, 1, 0, 0.20),
  tree_mar = c(1, 1, 2, 1),
  legend_mar = c(0, 0.5, 0.3, 0.5),
  cex_tip = 0.8,
  tip_cex = 1.5,
  label_cex = 1.5,
  radial_padding = 0.025,
  extra_space = 0.35,
  rotate_labels = TRUE
)

dev.off()

# 13. TREES .nwk AND DISTANCE MATRIX

if (!is.null(tree_S_nj)) {
  write.tree(tree_S_nj, file = file.path(diag_dir, "ITS2_sample_tree_S_NJ.nwk"))
}

if (!is.null(tree_WS_nj)) {
  write.tree(tree_WS_nj, file = file.path(diag_dir, "ITS2_sample_tree_WS_NJ.nwk"))
}

if (!is.null(tree_WC_nj)) {
  write.tree(tree_WC_nj, file = file.path(diag_dir, "ITS2_sample_tree_WC_NJ.nwk"))
}

write.tree(
  its2_tree_guni,
  file = file.path(diag_dir, "ITS2_sequence_tree_for_GUniFrac.nwk")
)

write.table(
  its2_dist,
  file = file.path(diag_dir, "ITS2_GUniFrac_d0.5_distance_matrix.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

cat("Done.\n")
