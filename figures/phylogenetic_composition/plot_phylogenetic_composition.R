# Figure S5 script
# Coral-field GUniFrac tree with metadata and matched Post-MED rings fan plot

# Rings:
# 1. Sites
# 2. Year
# 3. Coral_field Post-MED
# 4. Matched sediment Post-MED by colonyID + timepoint
# 5. Matched seawater Post-MED by collection_site + timepoint

input_core <- file.path("data", "core")
input_phylogeny <- file.path("data", "phylogeny")
figure_dir <- file.path("outputs", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# 1. packages

library(ape)
library(GUniFrac)
library(phangorn)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(cowplot)
library(grid)
library(svglite)

# 2. data

postMED <- read.table(
  file.path(input_core, "postMED_rel_abund_only.txt"),
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

meta <- read.table(
  file.path(input_core, "metadata_mayotte_all_v2.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

rownames(meta) <- meta$sample_uid

postMED_tree <- read.tree(file.path(input_phylogeny, "DIV_tree.nwk"))

# 3. Parameters

TOP_N_POSTMED <- 56

missing_match_label <- "No matched sample"
missing_match_col <- "grey90"
other_col <- "#444444"

# plots output format 
output_svg <- file.path(figure_dir, "Figure_S5_phylogenetic_composition.svg")
output_png <- file.path(figure_dir, "Figure_S5_phylogenetic_composition.png")

# 4. Filtering and Analysis functions
# Functions for post-MED taxon grouping, colour mapping, compartment filtering,
# sample matching, tree construction, and plot-data preparation.

get_prefix <- function(x) {
  sub("^([A-Za-z]+).*", "\\1", x)
}

get_core_family <- function(x) {
  ifelse(grepl("^C15", x), "C15",
         ifelse(grepl("^C116", x), "C116",
                ifelse(grepl("^C3", x), "C3",
                       ifelse(grepl("^[A-Za-z]+[0-9]+", x),
                              regmatches(x, regexpr("^[A-Za-z]+[0-9]+", x)),
                              get_prefix(x)))))
}

make_shades <- function(base_col, n) {
  if (n == 1) return(base_col)
  colorRampPalette(c(base_col, "#EDEDED"))(n)
}

build_postMED_color_map <- function(taxa_names) {
  
  taxa_names <- unique(taxa_names)
  prefixes <- get_prefix(taxa_names)
  core_families <- get_core_family(taxa_names)
  
  prefix_base_cols <- c(
    "A" = "#CCBB44",
    "B" = "#56B4E9",
    "C" = "#0072B2",
    "D" = "#E69F00",
    "F" = "#CC79A7",
    "G" = "#D55E00",
    "H" = "#999999"
  )
  
  missing_prefixes <- setdiff(unique(prefixes), names(prefix_base_cols))
  if (length(missing_prefixes) > 0) {
    extra_cols <- colorRampPalette(c("#444444", "#BBBBBB"))(length(missing_prefixes))
    names(extra_cols) <- missing_prefixes
    prefix_base_cols <- c(prefix_base_cols, extra_cols)
  }
  
  color_map <- character(length(taxa_names))
  names(color_map) <- taxa_names
  
  for (p in unique(prefixes)) {
    
    group_taxa <- sort(taxa_names[prefixes == p])
    
    if (p == "C") {
      
      c_core <- core_families[prefixes == "C"]
      
      for (fam in unique(c_core)) {
        
        fam_taxa <- sort(taxa_names[prefixes == "C" & core_families == fam])
        
        if (fam == "C15") {
          fam_cols <- colorRampPalette(c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF"))(length(fam_taxa))
        } else if (fam == "C3") {
          fam_cols <- colorRampPalette(c("#7F0000", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7"))(length(fam_taxa))
        } else if (fam == "C116") {
          fam_cols <- colorRampPalette(c("#762A83", "#9970AB", "#C2A5CF", "#E7D4E8"))(length(fam_taxa))
        } else {
          fam_cols <- make_shades("#88CCEE", length(fam_taxa))
        }
        
        names(fam_cols) <- fam_taxa
        color_map[fam_taxa] <- fam_cols
      }
      
    } else if (p == "A") {
      
      group_cols <- colorRampPalette(c("#6B4E00", "#A6761D", "#CCBB44", "#E6D76A", "#F2E6A2"))(length(group_taxa))
      names(group_cols) <- group_taxa
      color_map[group_taxa] <- group_cols
      
    } else {
      
      group_cols <- make_shades(prefix_base_cols[p], length(group_taxa))
      names(group_cols) <- group_taxa
      color_map[group_taxa] <- group_cols
    }
  }
  
  color_map
}

make_discrete_colors <- function(x, palette = "Dark2") {
  
  vals <- sort(unique(as.character(x)))
  vals <- vals[!is.na(vals)]
  
  n <- length(vals)
  
  if (n <= 8) {
    cols <- RColorBrewer::brewer.pal(max(3, n), palette)[seq_len(n)]
  } else {
    cols <- colorRampPalette(RColorBrewer::brewer.pal(8, palette))(n)
  }
  
  names(cols) <- vals
  cols
}

prepare_postMED_compartment <- function(compartment_name) {
  
  sample_ids <- rownames(meta)[meta$sample_type == compartment_name]
  postMED_comp <- postMED[sample_ids, , drop = FALSE]
  
  common_taxa <- intersect(colnames(postMED_comp), postMED_tree$tip.label)
  postMED_comp <- postMED_comp[, common_taxa, drop = FALSE]
  
  postMED_comp <- postMED_comp[, colSums(postMED_comp) > 0, drop = FALSE]
  postMED_comp <- postMED_comp[rowSums(postMED_comp) > 0, , drop = FALSE]
  
  tree_comp <- drop.tip(postMED_tree, setdiff(postMED_tree$tip.label, colnames(postMED_comp)))
  tree_comp <- root(tree_comp, outgroup = tree_comp$tip.label[1], resolve.root = TRUE)
  
  top_taxa <- names(sort(colMeans(postMED_comp), decreasing = TRUE))
  top_taxa <- top_taxa[seq_len(min(TOP_N_POSTMED, length(top_taxa)))]
  
  postMED_top <- postMED_comp[, top_taxa, drop = FALSE]
  postMED_other <- pmax(0, 1 - rowSums(postMED_top))
  postMED_plot <- cbind(postMED_top, Other = postMED_other)
  
  gu <- GUniFrac(as.matrix(postMED_comp), tree_comp, alpha = 0.5)
  # Generalized UniFrac with alpha = 0.5 (d_UW would select unweighted UniFrac).
  dist_guni <- gu$unifracs[, , "d_0.5"]
  
  tree_samples <- as.phylo(hclust(as.dist(dist_guni), method = "average"))
  
  list(
    compartment = compartment_name,
    tree = tree_samples,
    order = tree_samples$tip.label,
    postMED_plot = postMED_plot
  )
}

make_matched_ring <- function(target_order,
                              source_object,
                              match_keys,
                              missing_label = missing_match_label) {
  
  target_meta <- meta[target_order, , drop = FALSE]
  source_ids <- rownames(source_object$postMED_plot)
  source_meta <- meta[source_ids, , drop = FALSE]
  
  target_key <- apply(target_meta[, match_keys, drop = FALSE], 1, paste, collapse = "___")
  source_key <- apply(source_meta[, match_keys, drop = FALSE], 1, paste, collapse = "___")
  
  source_lookup <- source_ids[!duplicated(source_key)]
  names(source_lookup) <- source_key[!duplicated(source_key)]
  
  matched_source <- source_lookup[target_key]
  
  output_cols <- c(colnames(source_object$postMED_plot), missing_label)
  
  out <- matrix(
    0,
    nrow = length(target_order),
    ncol = length(output_cols),
    dimnames = list(target_order, output_cols)
  )
  
  for (i in seq_along(target_order)) {
    
    src <- matched_source[i]
    
    if (!is.na(src) && src %in% rownames(source_object$postMED_plot)) {
      out[target_order[i], colnames(source_object$postMED_plot)] <-
        as.numeric(source_object$postMED_plot[src, ])
    } else {
      out[target_order[i], missing_label] <- 1
    }
  }
  
  out
}

matrix_to_long <- function(mat, sample_order, taxon_levels) {
  
  as.data.frame(mat, check.names = FALSE) |>
    rownames_to_column("SampleID") |>
    pivot_longer(
      cols = -SampleID,
      names_to = "Taxon",
      values_to = "Abundance"
    ) |>
    mutate(
      SampleID = factor(SampleID, levels = sample_order),
      Taxon = factor(Taxon, levels = taxon_levels)
    )
}

# 5. Prepare compartments (coral/sediment/seawater)
# Process each compartment and extract the UniFrac-based coral sample order

coral_obj <- prepare_postMED_compartment("coral_field")
sediment_obj <- prepare_postMED_compartment("sediment")
seawater_obj <- prepare_postMED_compartment("seawater")

coral_order <- coral_obj$order

# 6. Color palettes maps

all_taxa <- names(sort(colMeans(postMED), decreasing = TRUE))

# Top 10 most abundant Post-MED sequences across the full dataset
top10_postMED <- head(all_taxa, 10)

postMED_cols_all <- build_postMED_color_map(all_taxa)

used_taxa <- unique(c(
  colnames(coral_obj$postMED_plot),
  colnames(sediment_obj$postMED_plot),
  colnames(seawater_obj$postMED_plot)
))

postMED_cols <- postMED_cols_all[used_taxa]
postMED_cols["Other"] <- other_col
postMED_cols[missing_match_label] <- missing_match_col

postMED_cols <- postMED_cols[!is.na(names(postMED_cols))]
taxon_levels <- names(postMED_cols)

# Legend labels: bold top 10 Post-MED sequences
postMED_legend_labels <- lapply(names(postMED_cols), function(x) {
  if (x %in% top10_postMED) {
    bquote(bold(.(x)))
  } else {
    x
  }
})

names(postMED_legend_labels) <- names(postMED_cols)

year_cols <- make_discrete_colors(meta$year, palette = "Set2")

site_levels <- sort(unique(as.character(meta$collection_site)))
site_palette <- c(
  "#0072B2",
  "#D55E00",
  "#009E73",
  "#CC79A7",
  "#E69F00",
  "#56B4E9",
  "#F0E442",
  "#000000"
)

site_cols <- setNames(
  rep(site_palette, length.out = length(site_levels)),
  site_levels
)

# 7. Build ring data
# Prepare metadata and abundance data for tree rings, matching sediment by
# colony and timepoint and seawater by site and timepoint to the coral order

tip_meta <- meta[coral_order, , drop = FALSE] |>
  rownames_to_column("label") |>
  mutate(
    collection_site = factor(collection_site, levels = names(site_cols)),
    year = factor(year, levels = names(year_cols))
  )

year_ring <- data.frame(
  SampleID = factor(coral_order, levels = coral_order),
  Year = factor(meta[coral_order, "year"], levels = names(year_cols)),
  Abundance = 1
)

coral_ring_long <- matrix_to_long(
  coral_obj$postMED_plot,
  sample_order = coral_order,
  taxon_levels = taxon_levels
)

sediment_matched <- make_matched_ring(
  target_order = coral_order,
  source_object = sediment_obj,
  match_keys = c("colonyID", "timepoint")
)

sediment_ring_long <- matrix_to_long(
  sediment_matched,
  sample_order = coral_order,
  taxon_levels = taxon_levels
)

seawater_matched <- make_matched_ring(
  target_order = coral_order,
  source_object = seawater_obj,
  match_keys = c("collection_site", "timepoint")
)

seawater_ring_long <- matrix_to_long(
  seawater_matched,
  sample_order = coral_order,
  taxon_levels = taxon_levels
)

# 8. Plot
# Build the circular UniFrac-based coral tree with collection site and year
# annotations and concentric post-MED composition rings for coral, sediment,
# and matched seawater samples

matched_plot <- ggtree(
  coral_obj$tree,
  layout = "fan",
  open.angle = 0
) %<+% tip_meta +
  
  xlim(0, 1.1) +
  
  geom_point2(
    aes(
      subset = isTip,
      x = x + 0.02,
      color = collection_site
    ),
    size = 2.3,
    alpha = 0.95
  ) +
  
  scale_color_manual(
    values = site_cols,
    breaks = names(site_cols),
    drop = FALSE,
    name = "Collection site"
  ) +
  
  geom_fruit(
    data = year_ring,
    geom = geom_col,
    mapping = aes(
      y = SampleID,
      x = Abundance,
      fill = Year
    ),
    orientation = "y",
    pwidth = 0.04,
    offset = 0.10,
    width = 0.9
  ) +
  
  scale_fill_manual(
    values = year_cols,
    name = "Year"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  geom_fruit(
    data = coral_ring_long,
    geom = geom_col,
    mapping = aes(
      y = SampleID,
      x = Abundance,
      fill = Taxon
    ),
    orientation = "y",
    pwidth = 0.26,
    offset = 0.04
  ) +
  
  geom_fruit(
    data = sediment_ring_long,
    geom = geom_col,
    mapping = aes(
      y = SampleID,
      x = Abundance,
      fill = Taxon
    ),
    orientation = "y",
    pwidth = 0.26,
    offset = 0.04
  ) +
  
  geom_fruit(
    data = seawater_ring_long,
    geom = geom_col,
    mapping = aes(
      y = SampleID,
      x = Abundance,
      fill = Taxon
    ),
    orientation = "y",
    pwidth = 0.26,
    offset = 0.04
  ) +
  
  scale_fill_manual(
    values = postMED_cols,
    breaks = names(postMED_cols),
    labels = postMED_legend_labels[names(postMED_cols)],
    name = "Post-MED sequences",
    guide = guide_legend(
      ncol = 10,
      byrow = TRUE,
      keyheight = unit(0.25, "cm"),
      keywidth = unit(0.25, "cm")
    )
  ) +
  
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE),
    fill = guide_legend(
      ncol = 10,
      byrow = TRUE,
      keyheight = unit(0.25, "cm"),
      keywidth = unit(0.25, "cm")
    )
  ) +
  
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.spacing.y = unit(0.4, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

# 9. Compact layout and export
# Separate and reposition the plot and legend for a compact final layout,
# export figure as SVG and PNG files

matched_legend <- cowplot::get_legend(
  matched_plot +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.spacing.y = unit(0.1, "cm"),
      legend.spacing.x = unit(0.1, "cm")
    )
)

matched_plot_noleg <- matched_plot +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )

matched_plot_big <- cowplot::ggdraw() +
  cowplot::draw_plot(
    matched_plot_noleg,
    x = 0,
    y = -0.08,
    width = 1,
    height = 1.15,
    scale = 1.15
  )

matched_legend_bg <- cowplot::ggdraw(matched_legend) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

matched_plot_final <- cowplot::plot_grid(
  matched_plot_big,
  matched_legend_bg,
  ncol = 1,
  rel_heights = c(1, 0.32)
)

ggsave(
  output_svg,
  matched_plot_final,
  width = 18,
  height = 17,
  device = svglite::svglite,
  bg = "white"
)

ggsave(
  output_png,
  matched_plot_final,
  width = 18,
  height = 17,
  dpi = 300,
  bg = "white"
)

