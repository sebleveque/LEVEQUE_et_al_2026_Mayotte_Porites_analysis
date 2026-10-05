# Figure 2 script
# Gridcell post-MED relative abundance barplots for Coral, sediment, and seawater 
# stacked vertically with shared postMED seqeunce legend

input_dir <- file.path("data", "core")
figure_dir <- file.path("outputs", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# 0. packages

library(RColorBrewer)

# 1. data

postMED <- read.table(
  file.path(input_dir, "postMED_rel_abund_only.txt"),
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

meta <- read.table(
  file.path(input_dir, "metadata_mayotte_all_v2.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

rownames(meta) <- meta$sample_uid
postMED <- postMED[rownames(meta), , drop = FALSE]

# 2. Parameters
# post-MED sequences displayed per grid cell, site order,
# and colours used for aggregated or missing data
TOP_N_POSTMED_CELL <- 15

site_order <- c("SAKOULI", "MUSICAL", "TANARAKI", "BOUENI", "MTSAMBORO")

postMED_other_col <- "#444444"
postMED_nodata_col <- "grey80"

# 3. Color palette selection functions
# Group post-MED sequences by lineage, assigned related colour gradients,
# and organised colours for plotting and legend display

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
  uniq_prefixes <- unique(prefixes)
  
  prefix_base_cols <- c(
    "A" = "#CCBB44",
    "B" = "#56B4E9",
    "C" = "#0072B2",
    "D" = "#E69F00",
    "F" = "#CC79A7",
    "G" = "#D55E00",
    "H" = "#999999",
    "I" = "#882255",
    "J" = "#44AA99",
    "K" = "#117733",
    "L" = "#88CCEE",
    "M" = "#DDCC77",
    "N" = "#AA4499",
    "O" = "#661100",
    "P" = "#6699CC",
    "Q" = "#AA4466",
    "R" = "#4477AA",
    "S" = "#228833",
    "T" = "#CCBB44",
    "U" = "#EE6677",
    "V" = "#BBBBBB"
  )
  
  missing_prefixes <- setdiff(uniq_prefixes, names(prefix_base_cols))
  
  if (length(missing_prefixes) > 0) {
    fallback_cols <- colorRampPalette(c("#444444", "#BBBBBB"))(length(missing_prefixes))
    names(fallback_cols) <- missing_prefixes
    prefix_base_cols <- c(prefix_base_cols, fallback_cols)
  }
  
  color_map <- character(length(taxa_names))
  names(color_map) <- taxa_names
  
  c_taxa <- taxa_names[prefixes == "C"]
  
  if (length(c_taxa) > 0) {
    
    c_core <- core_families[prefixes == "C"]
    
    for (fam in unique(c_core)) {
      
      fam_taxa <- sort(c_taxa[c_core == fam])
      
      if (fam == "C15") {
        fam_cols <- colorRampPalette(c(
          "#08306B", "#2171B5", "#6BAED6", "#C6DBEF"
        ))(length(fam_taxa))
      } else if (fam == "C3") {
        fam_cols <- colorRampPalette(c(
          "#7F0000", "#B2182B", "#D6604D", "#F4A582"
        ))(length(fam_taxa))
      } else {
        fam_cols <- make_shades("#CC79A7", length(fam_taxa))
      }
      
      names(fam_cols) <- fam_taxa
      color_map[fam_taxa] <- fam_cols[fam_taxa]
    }
  }
  
  non_c_prefixes <- setdiff(uniq_prefixes, "C")
  
  for (p in non_c_prefixes) {
    
    group_taxa <- sort(taxa_names[prefixes == p])
    
    if (p == "A") {
      group_cols <- colorRampPalette(c(
        "#6B4E00", "#A6761D", "#CCBB44", "#E6D76A", "#F2E6A2"
      ))(length(group_taxa))
    } else {
      group_cols <- make_shades(prefix_base_cols[p], length(group_taxa))
    }
    
    names(group_cols) <- group_taxa
    color_map[group_taxa] <- group_cols[group_taxa]
  }
  
  color_map
}

get_legend_group <- function(x) {
  ifelse(grepl("^C15", x), "C15",
         ifelse(grepl("^C116", x), "C116",
                ifelse(grepl("^C3", x), "C3",
                       ifelse(x %in% c("Other", "NoData"),
                              x,
                              get_prefix(x)))))
}

sort_colors_by_gradient <- function(cols_named) {
  
  taxa <- names(cols_named)
  cols <- unname(cols_named)
  grp <- get_legend_group(taxa)
  
  preferred_order <- c(
    "A", "B",
    "C3", "C15", "C116", "C",
    "D", "F", "G", "H",
    "I", "J", "K", "L", "M", "N", "O", "P",
    "Q", "R", "S", "T", "U", "V",
    "Other", "NoData"
  )
  
  grp <- factor(grp, levels = preferred_order)
  
  rgb <- grDevices::col2rgb(cols) / 255
  lab <- grDevices::convertColor(t(rgb), from = "sRGB", to = "Lab")
  brightness <- lab[, 1]
  
  ord <- order(grp, brightness, taxa, na.last = TRUE)
  cols_named[ord]
}

move_last <- function(cols_named, name_last) {
  if (name_last %in% names(cols_named)) {
    cols_named <- c(
      cols_named[names(cols_named) != name_last],
      cols_named[name_last]
    )
  }
  cols_named
}

# 4. Metadata preparation and global colors
# Standardized metadata fields, defined the global post-MED colour map,
# and established site and sampling-date structure

meta$collection_site <- toupper(trimws(meta$collection_site))
meta$collection_date <- as.integer(meta$collection_date)
meta$colony_number <- suppressWarnings(as.integer(as.character(meta$colony_number)))

all_postMED_taxa <- names(sort(colMeans(postMED, na.rm = TRUE), decreasing = TRUE))
postMED_color_map <- build_postMED_color_map(all_postMED_taxa)

meta_coral_ref <- meta[
  meta$sample_type == "coral_field" &
    meta$collection_site %in% site_order,
  , drop = FALSE
]

date_order <- sort(unique(meta_coral_ref$collection_date))
date_labels <- format(
  as.Date(as.character(date_order), format = "%Y%m%d"),
  "%Y-%m-%d"
)

N_SITE <- length(site_order)
N_DATE <- length(date_order)

# 5. Build compartment grids
#  coral, sediment, and seawater post-MED compositions organisation by site and date,
# retained most abundant sequences per grid cell and identification of missing samples
build_compartment_grid <- function(compartment_name) {
  
  meta_comp <- meta[
    meta$sample_type == compartment_name &
      meta$collection_site %in% site_order,
    , drop = FALSE
  ]
  
  postMED_comp <- postMED[rownames(meta_comp), , drop = FALSE]
  
  grid_cells <- matrix(
    vector("list", N_DATE * N_SITE),
    nrow = N_DATE,
    ncol = N_SITE
  )
  
  rownames(grid_cells) <- date_labels
  colnames(grid_cells) <- site_order
  
  for (r in seq_along(date_order)) {
    for (c in seq_along(site_order)) {
      
      date_k <- date_order[r]
      site_k <- site_order[c]
      
      coral_cell <- meta_coral_ref[
        meta_coral_ref$collection_site == site_k &
          meta_coral_ref$collection_date == date_k,
        , drop = FALSE
      ]
      
      coral_cell <- coral_cell[!is.na(coral_cell$colony_number), , drop = FALSE]
      if (nrow(coral_cell) == 0) next
      
      if (compartment_name %in% c("coral_field", "sediment")) {
        
        if (site_k == "BOUENI") {
          all_colonies <- seq(min(coral_cell$colony_number), 15)
        } else {
          all_colonies <- seq(
            min(coral_cell$colony_number),
            max(coral_cell$colony_number)
          )
        }
        
      } else {
        
        all_colonies <- "Seawater"
      }
      
      meta_cell <- meta_comp[
        meta_comp$collection_site == site_k &
          meta_comp$collection_date == date_k,
        , drop = FALSE
      ]
      
      if (compartment_name %in% c("coral_field", "sediment")) {
        meta_cell <- meta_cell[!is.na(meta_cell$colony_number), , drop = FALSE]
        meta_cell <- meta_cell[order(meta_cell$colony_number), , drop = FALSE]
      } else {
        meta_cell <- meta_cell[order(rownames(meta_cell)), , drop = FALSE]
      }
      
      if (nrow(meta_cell) > 0) {
        
        postMED_present_all <- postMED_comp[rownames(meta_cell), , drop = FALSE]
        postMED_present_all[is.na(postMED_present_all)] <- 0
        
        pm_means <- colMeans(postMED_present_all, na.rm = TRUE)
        pm_means <- pm_means[!names(pm_means) %in% c("Other", "NoData")]
        
        topN <- names(sort(pm_means, decreasing = TRUE))[
          seq_len(min(TOP_N_POSTMED_CELL, length(pm_means)))
        ]
        
      } else {
        topN <- character(0)
      }
      
      pm_plot <- matrix(
        0,
        nrow = length(all_colonies),
        ncol = length(topN) + 2
      )
      
      rownames(pm_plot) <- as.character(all_colonies)
      colnames(pm_plot) <- c(topN, "Other", "NoData")
      
      if (nrow(meta_cell) > 0 && length(topN) > 0) {
        
        if (compartment_name %in% c("coral_field", "sediment")) {
          
          plot_positions <- meta_cell$colony_number
          keep <- plot_positions %in% all_colonies
          
          plot_positions <- plot_positions[keep]
          meta_cell_kept <- meta_cell[keep, , drop = FALSE]
          
          postMED_present <- postMED_comp[
            rownames(meta_cell_kept),
            topN,
            drop = FALSE
          ]
          
          postMED_present[is.na(postMED_present)] <- 0
          other_present <- pmax(0, 1 - rowSums(postMED_present))
          
          pm_plot[as.character(plot_positions), topN] <- as.matrix(postMED_present)
          pm_plot[as.character(plot_positions), "Other"] <- other_present
          pm_plot[as.character(plot_positions), "NoData"] <- 0
          
        } else {
          
          postMED_present <- postMED_comp[
            rownames(meta_cell),
            topN,
            drop = FALSE
          ]
          
          postMED_present[is.na(postMED_present)] <- 0
          
          seawater_mean <- colMeans(postMED_present, na.rm = TRUE)
          other_mean <- pmax(0, 1 - sum(seawater_mean))
          
          pm_plot["Seawater", topN] <- seawater_mean
          pm_plot["Seawater", "Other"] <- other_mean
          pm_plot["Seawater", "NoData"] <- 0
        }
      }
      
      missing <- rowSums(
        pm_plot[, !colnames(pm_plot) %in% "NoData", drop = FALSE]
      ) == 0
      
      pm_plot[missing, "NoData"] <- 1
      
      pm_cols <- c(
        postMED_color_map[
          colnames(pm_plot)[!colnames(pm_plot) %in% c("Other", "NoData")]
        ],
        "Other" = postMED_other_col,
        "NoData" = postMED_nodata_col
      )
      
      grid_cells[[r, c]] <- list(
        postMED_plot = pm_plot,
        postMED_cols = pm_cols
      )
    }
  }
  
  grid_cells
}

coral_grid <- build_compartment_grid("coral_field")
sediment_grid <- build_compartment_grid("sediment")
seawater_grid <- build_compartment_grid("seawater")

# 6. Shared legend colors across compartment (coral/sediment/seawater)
# Combined post-MED sequences present across all compartments into a single,
# consistently ordered colour legend, including Other and No Data categories

make_common_legend_cols <- function(...) {
  
  grid_list <- list(...)
  
  legend_labels <- unique(unlist(lapply(grid_list, function(grid_cells) {
    unlist(lapply(as.vector(grid_cells), function(x) {
      if (is.null(x)) return(NULL)
      colnames(x$postMED_plot)
    }))
  })))
  
  legend_cols <- c(
    postMED_color_map[legend_labels[!legend_labels %in% c("Other", "NoData")]],
    "Other" = postMED_other_col,
    "NoData" = postMED_nodata_col
  )
  
  legend_cols <- legend_cols[names(legend_cols) %in% legend_labels]
  legend_cols <- legend_cols[!is.na(legend_cols)]
  
  legend_cols <- sort_colors_by_gradient(legend_cols)
  legend_cols <- move_last(legend_cols, "Other")
  legend_cols <- move_last(legend_cols, "NoData")
  
  legend_cols
}

# 7. Combined figure
# Arrange coral, sediment, and seawater post-MED barplots by site and date
# in vertical stacked panels with a shared sequence legend

plot_combined_three_compartments <- function(coral_grid,
                                             sediment_grid,
                                             seawater_grid,
                                             figure_title = "") {
  
  # Display dates as "Month year" in English
  dates_display <- as.Date(as.character(date_order), format = "%Y%m%d")
  
  date_labels_display <- paste(
    month.name[as.integer(format(dates_display, "%m"))],
    format(dates_display, "%Y")
  )
  
  legend_cols <- make_common_legend_cols(
    coral_grid,
    sediment_grid,
    seawater_grid
  )
  
  total_plot_rows <- N_DATE * 3
  total_plot_cols <- N_SITE
  
  plot_ids <- matrix(
    seq_len(total_plot_rows * total_plot_cols),
    nrow = total_plot_rows,
    ncol = total_plot_cols,
    byrow = TRUE
  )
  
  legend_id <- total_plot_rows * total_plot_cols + 1
  
  layout_mat <- cbind(
    plot_ids,
    rep(legend_id, total_plot_rows)
  )
  
  row_heights <- rep(1, total_plot_rows)
  
  sediment_last_row <- N_DATE * 2
  row_heights[sediment_last_row] <- 1.45
  
  layout(
    layout_mat,
    widths = c(rep(1.75, N_SITE), 2.0),
    heights = row_heights
  )
  
  par(
    oma = c(2.2, 4.4, 2.0, 0.5),
    xaxs = "i",
    yaxs = "i"
  )
  
  # (a) Coral tissue
  
  for (r in seq_len(N_DATE)) {
    for (c in seq_len(N_SITE)) {
      
      cell <- coral_grid[[r, c]]
      
      par(mar = c(0.6, 0.6, ifelse(r == 1, 1.0, 0.2), 0.2))
      
      if (is.null(cell)) {
        plot.new()
        box(lwd = 0.6)
      } else {
        barplot(
          t(cell$postMED_plot),
          col = cell$postMED_cols[colnames(cell$postMED_plot)],
          border = NA,
          space = 0,
          axes = FALSE,
          names.arg = rep("", nrow(cell$postMED_plot))
        )
        box(lwd = 0.6)
      }
      
      if (r == 1) {
        mtext(site_order[c], side = 3, line = 0.25, cex = 1.2, font = 2)
      }
      
      if (c == 1) {
        mtext(date_labels_display[r], side = 2, line = 0.25, cex = 1.2, font = 2)
      }
    }
  }
  
  mtext(
    expression(bold("(a) ") * bolditalic("Porites") * bold(" spp. tissue")),
    side = 2,
    outer = TRUE,
    line = 2,
    at = 0.84,
    cex = 1.4,
    font = 2
  )
  
  # (b) Interstitial sediment water
  
  for (r in seq_len(N_DATE)) {
    for (c in seq_len(N_SITE)) {
      
      cell <- sediment_grid[[r, c]]
      
      par(
        mar = c(
          ifelse(r == N_DATE, 4.0, 0.6),
          0.6,
          ifelse(r == 1, 1.0, 0.2),
          0.2
        )
      )
      
      if (is.null(cell)) {
        
        plot.new()
        box(lwd = 0.6)
        
      } else {
        
        bar_centers <- barplot(
          t(cell$postMED_plot),
          col = cell$postMED_cols[colnames(cell$postMED_plot)],
          border = NA,
          space = 0,
          axes = FALSE,
          names.arg = rep("", nrow(cell$postMED_plot))
        )
        
        box(lwd = 0.6)
        
        if (r == N_DATE) {
          axis(
            1,
            at = bar_centers,
            labels = rownames(cell$postMED_plot),
            las = 2,
            cex.axis = 1,
            tck = -0.02
          )
        }
      }
      
#      if (r == 1) {
#        mtext(site_order[c], side = 3, line = 0.25, cex = 1.0, font = 2)
#      }
      
      if (c == 1) {
        mtext(date_labels_display[r], side = 2, line = 0.25, cex = 1.2, font = 2)
      }
    }
  }
  
  mtext(
    "(b) Interstitial sediment water",
    side = 2,
    outer = TRUE,
    line = 2,
    at = 0.50,
    cex = 1.4,
    font = 2
  )
  
  # (c) Water column
  
  for (r in seq_len(N_DATE)) {
    for (c in seq_len(N_SITE)) {
      
      cell <- seawater_grid[[r, c]]
      
      par(mar = c(0.6, 0.6, ifelse(r == 1, 1.0, 0.2), 0.2))
      
      if (is.null(cell)) {
        
        plot.new()
        box(lwd = 0.6)
        
      } else {
        
        barplot(
          t(cell$postMED_plot),
          col = cell$postMED_cols[colnames(cell$postMED_plot)],
          border = NA,
          space = 0,
          axes = FALSE,
          names.arg = rep("", nrow(cell$postMED_plot))
        )
        
        box(lwd = 0.6)
      }
      
 #     if (r == 1) {
 #       mtext(site_order[c], side = 3, line = 0.25, cex = 1.0, font = 2)
 #     }
      
      if (c == 1) {
        mtext(date_labels_display[r], side = 2, line = 0.25, cex = 1.2, font = 2)
      }
    }
  }
  
  mtext(
    "(c) Water column",
    side = 2,
    outer = TRUE,
    line = 2,
    at = 0.17,
    cex = 1.4,
    font = 2
  )
  
  # Shared legend
  
  par(mar = c(0, 0.4, 0, 0))
  plot.new()
  
  legend(
    "center",
    legend = names(legend_cols),
    fill = unname(legend_cols),
    border = NA,
    bty = "n",
    cex = 1.6,
    ncol = 2,
    title = "Post-MED"
  )
  
  mtext(
    figure_title,
    outer = TRUE,
    side = 3,
    line = 0.4,
    cex = 1.2,
    font = 2
  )
}

# 8. Save combined figure
# Export as PNG and SVG files

png(
  filename = file.path(figure_dir, "Figure_2_postMED_composition.png"),
  width = 9500,
  height = 11000,
  res = 600
)

plot_combined_three_compartments(
  coral_grid,
  sediment_grid,
  seawater_grid,
  figure_title = ""
)

dev.off()

svg(
  filename = file.path(figure_dir, "Figure_2_postMED_composition.svg"),
  width = 15.8,
  height = 18.3
)

plot_combined_three_compartments(
  coral_grid,
  sediment_grid,
  seawater_grid,
  figure_title = ""
)

dev.off()
