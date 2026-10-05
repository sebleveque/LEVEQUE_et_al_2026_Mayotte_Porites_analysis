# General SymPortal sequence statistics

input_dir <- file.path("data", "core")
out_dir <- file.path("outputs", "general_statistics")
table_dir <- file.path(out_dir, "diagnostics")
manuscript_table_dir <- file.path("outputs", "supplementary_tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_table_dir, recursive = TRUE, showWarnings = FALSE)

# 1. data 
postMED_counts <- read.table(
  file.path(input_dir, "postMED_abs_abund_only.txt"),
  header = TRUE, sep = "\t",
  row.names = 1, check.names = FALSE
)

profiles_counts <- read.table(
  file.path(input_dir, "ITS2_profiles_absolute_abund_only.txt"),
  header = TRUE, sep = "\t",
  row.names = 1, check.names = FALSE
)

meta <- read.table(
  file.path(input_dir, "metadata_mayotte_all_v2.txt"),
  header = TRUE, sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 2. Sample IDs filtering
rownames(postMED_counts) <- as.character(rownames(postMED_counts))
rownames(profiles_counts) <- as.character(rownames(profiles_counts))
meta$sample_uid <- as.character(meta$sample_uid)

# 3. Keep only samples present in all files
common_samples <- Reduce(
  intersect,
  list(
    rownames(postMED_counts),
    rownames(profiles_counts),
    meta$sample_uid
  )
)

postMED_counts <- postMED_counts[common_samples, , drop = FALSE]
profiles_counts <- profiles_counts[common_samples, , drop = FALSE]
meta_sub <- meta[match(common_samples, meta$sample_uid), ]

# 4. Sample level calculations
sample_depth <- rowSums(postMED_counts)

profiles_per_sample <- rowSums(profiles_counts > 0)

majority_sequence <- apply(postMED_counts, 1, function(x) {
  colnames(postMED_counts)[which.max(x)]
})

majority_sequence_count <- apply(postMED_counts, 1, max)

profile_sequence_count <- rowSums(profiles_counts)

non_profile_sequence_count <- sample_depth - profile_sequence_count

prop_profile_sequences <- profile_sequence_count / sample_depth
prop_non_profile_sequences <- non_profile_sequence_count / sample_depth

# 5. Sample summary
sample_summary <- data.frame(
  sample_uid = common_samples,
  sample_name = meta_sub$sample_name,
  compartment = meta_sub$sample_type,
  host_genus = meta_sub$host_genus,
  total_postMED_sequences = sample_depth,
  n_ITS2_profiles = profiles_per_sample,
  majority_sequence = majority_sequence,
  majority_sequence_count = majority_sequence_count,
  profile_sequence_count = profile_sequence_count,
  non_profile_sequence_count = non_profile_sequence_count,
  prop_profile_sequences = prop_profile_sequences,
  prop_non_profile_sequences = prop_non_profile_sequences,
  stringsAsFactors = FALSE
)

# 6. Overall summary
overall_summary <- data.frame(
  total_samples = nrow(sample_summary),
  total_postMED_sequences = sum(sample_summary$total_postMED_sequences),
  mean_postMED_depth_per_sample = mean(sample_summary$total_postMED_sequences),
  median_postMED_depth_per_sample = median(sample_summary$total_postMED_sequences),
  min_postMED_depth_per_sample = min(sample_summary$total_postMED_sequences),
  max_postMED_depth_per_sample = max(sample_summary$total_postMED_sequences),
  total_ITS2_type_profiles = sum(colSums(profiles_counts > 0) > 0),
  stringsAsFactors = FALSE
)

# 7. Compartment summary
compartment_summary <- aggregate(
  cbind(
    total_postMED_sequences,
    n_ITS2_profiles,
    prop_profile_sequences,
    prop_non_profile_sequences
  ) ~ compartment,
  data = sample_summary,
  FUN = function(x) c(
    total = sum(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
)

compartment_summary <- do.call(data.frame, compartment_summary)

n_samples_compartment <- as.data.frame(table(sample_summary$compartment))
colnames(n_samples_compartment) <- c("compartment", "n_samples")

compartment_summary <- merge(
  n_samples_compartment,
  compartment_summary,
  by = "compartment",
  all.x = TRUE
)

# 8. Samples per host genus per compartment
genus_by_compartment <- as.data.frame(
  table(sample_summary$host_genus, sample_summary$compartment)
)

colnames(genus_by_compartment) <- c(
  "host_genus",
  "compartment",
  "n_samples"
)

genus_by_compartment <- genus_by_compartment[
  genus_by_compartment$n_samples > 0,
]

# 9. Majority sequence summary by compartment
majority_by_compartment <- as.data.frame(
  table(sample_summary$compartment, sample_summary$majority_sequence)
)

colnames(majority_by_compartment) <- c(
  "compartment",
  "majority_sequence",
  "n_samples"
)

majority_by_compartment <- majority_by_compartment[
  majority_by_compartment$n_samples > 0,
]

# 10. ITS2 type profile richness by compartment
profiles_present_by_compartment <- lapply(
  unique(sample_summary$compartment),
  function(comp) {
    samples_comp <- sample_summary$sample_uid[
      sample_summary$compartment == comp
    ]
    
    data.frame(
      compartment = comp,
      total_unique_ITS2_type_profiles = sum(
        colSums(profiles_counts[samples_comp, , drop = FALSE] > 0) > 0
      )
    )
  }
)

profile_richness_by_compartment <- do.call(
  rbind,
  profiles_present_by_compartment
)

# 11. Summary
overall_summary
compartment_summary
genus_by_compartment
majority_by_compartment
profile_richness_by_compartment

# 12. export as .txt files
write.table(
  sample_summary,
  file = file.path(table_dir, "symportal_sample_summary.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  overall_summary,
  file = file.path(table_dir, "symportal_overall_summary.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  compartment_summary,
  file = file.path(table_dir, "symportal_compartment_summary.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  genus_by_compartment,
  file = file.path(table_dir, "symportal_genus_by_compartment.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  majority_by_compartment,
  file = file.path(table_dir, "symportal_majority_sequence_by_compartment.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  profile_richness_by_compartment,
  file = file.path(table_dir, "symportal_profile_richness_by_compartment.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

list.files(table_dir, pattern = "symportal_.*\\.txt$")

# 13. Final formatted compartment table

# Relabel compartments
sample_summary$compartment_label <- ifelse(
  sample_summary$compartment == "coral_field", "Porites spp",
  ifelse(
    sample_summary$compartment == "sediment", "Interstitial sediment water",
    ifelse(
      sample_summary$compartment == "seawater", "Water column",
      sample_summary$compartment
    )
  )
)

# Desired column order
comp_order <- c(
  "Porites spp",
  "Interstitial sediment water",
  "Water column"
)

# Function for majority sequences per compartment
get_majority_sequences <- function(df) {
  tab <- sort(table(df$majority_sequence), decreasing = TRUE)
  max_n <- max(tab)
  paste(names(tab)[tab == max_n], collapse = "; ")
}

# Create summary table
final_summary <- data.frame(
  metric = c(
    "Total number of samples",
    "Total number of sequences",
    "Average sequencing depth per sample (# reads)",
    "Proportion of profile sequences",
    "Proportion of non-profile sequences",
    "Total number of ITS2 type profiles",
    "Majority sequence(s)"
  ),
  stringsAsFactors = FALSE
)

for (comp in comp_order) {
  
  samples_comp <- sample_summary$sample_uid[
    sample_summary$compartment_label == comp
  ]
  
  df_comp <- sample_summary[
    sample_summary$compartment_label == comp,
  ]
  
  values <- c(
    nrow(df_comp),
    sum(df_comp$total_postMED_sequences, na.rm = TRUE),
    mean(df_comp$total_postMED_sequences, na.rm = TRUE),
    mean(df_comp$prop_profile_sequences, na.rm = TRUE),
    mean(df_comp$prop_non_profile_sequences, na.rm = TRUE),
    sum(colSums(profiles_counts[samples_comp, , drop = FALSE] > 0) > 0),
    get_majority_sequences(df_comp)
  )
  
  final_summary[[comp]] <- values
}

# View table
final_summary

# Write as txt file
write.table(
  final_summary,
  file = file.path(table_dir, "symportal_final_compartment_summary.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 14. Overall Symbiodiniaceae sequence summaries

# Comparable IDs
meta$sample_uid <- as.character(meta$sample_uid)
rownames(postMED_counts) <- as.character(rownames(postMED_counts))

# Match metadata to postMED table
common_samples <- intersect(rownames(postMED_counts), meta$sample_uid)

postMED_sub <- postMED_counts[common_samples, , drop = FALSE]
meta_sub <- meta[match(common_samples, meta$sample_uid), ]

# Total number of Symbiodiniaceae post-MED sequences across all compartments
total_symbiodiniaceae_sequences <- sum(rowSums(postMED_sub))

# Average per sample depth of post-MED sequences across all compartments
average_postMED_depth_per_sample <- mean(rowSums(postMED_sub))

# Summary
overall_sequence_summary <- data.frame(
  metric = c(
    "Total Symbiodiniaceae post-MED sequences",
    "Average post-MED sequencing depth per sample"
  ),
  value = c(
    total_symbiodiniaceae_sequences,
    average_postMED_depth_per_sample
  )
)

overall_sequence_summary

# 15. Percentage of Symbiodiniaceae sequences by genus/clade

sym_cols <- grep("Clade", colnames(meta_sub), value = TRUE)

sym_cols

meta_sub[, sym_cols] <- lapply(meta_sub[, sym_cols, drop = FALSE], as.numeric)

# Total sequences per clade/genus
genus_sequence_totals <- colSums(meta_sub[, sym_cols, drop = FALSE], na.rm = TRUE)

# Remove zero-count clades
genus_sequence_totals <- genus_sequence_totals[genus_sequence_totals > 0]

# Convert to percentages
genus_sequence_percent <- 100 * genus_sequence_totals / sum(genus_sequence_totals)

genus_percentage_summary <- data.frame(
  genus_or_clade = gsub("noName\\.Clade\\.", "Clade ", names(genus_sequence_percent)),
  total_sequences = as.numeric(genus_sequence_totals),
  percentage = round(as.numeric(genus_sequence_percent), 2),
  stringsAsFactors = FALSE
)

genus_percentage_summary

write.table(
  overall_sequence_summary,
  file = file.path(table_dir, "symportal_overall_sequence_summary.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  genus_percentage_summary,
  file = file.path(table_dir, "symportal_genus_percentage_summary.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 16. Overall sequence + clade/genus percentages by compartment

# Comparable IDs
meta$sample_uid <- as.character(meta$sample_uid)
rownames(postMED_counts) <- as.character(rownames(postMED_counts))

# Match metadata to postMED table
common_samples <- intersect(rownames(postMED_counts), meta$sample_uid)

postMED_sub <- postMED_counts[common_samples, , drop = FALSE]
meta_sub <- meta[match(common_samples, meta$sample_uid), ]

# Sample_type as compartment
meta_sub$compartment <- meta_sub$sample_type

sym_cols <- grep("Clade", colnames(meta_sub), value = TRUE)

meta_sub[, sym_cols] <- lapply(meta_sub[, sym_cols, drop = FALSE], as.numeric)

# 16A. Total sequences and average depth per compartment

postMED_depth <- rowSums(postMED_sub)

depth_df <- data.frame(
  sample_uid = common_samples,
  compartment = meta_sub$compartment,
  postMED_depth = postMED_depth,
  stringsAsFactors = FALSE
)

depth_by_compartment <- aggregate(
  postMED_depth ~ compartment,
  data = depth_df,
  FUN = function(x) c(
    total_sequences = sum(x, na.rm = TRUE),
    mean_depth = mean(x, na.rm = TRUE),
    n_samples = length(x)
  )
)

depth_by_compartment <- do.call(data.frame, depth_by_compartment)

depth_by_compartment

# 16B. Percentage of sequences from each clade/genus for each compartment

genus_percent_by_compartment <- do.call(
  rbind,
  lapply(unique(meta_sub$compartment), function(comp) {
    
    meta_comp <- meta_sub[meta_sub$compartment == comp, ]
    
    genus_totals <- colSums(meta_comp[, sym_cols, drop = FALSE], na.rm = TRUE)
    genus_totals <- genus_totals[genus_totals > 0]
    
    genus_percent <- 100 * genus_totals / sum(genus_totals)
    
    data.frame(
      compartment = comp,
      genus_or_clade = gsub("noName\\.Clade\\.", "Clade ", names(genus_percent)),
      total_sequences = as.numeric(genus_totals),
      percentage = round(as.numeric(genus_percent), 2),
      stringsAsFactors = FALSE
    )
  })
)

genus_percent_by_compartment

write.table(
  depth_by_compartment,
  file = file.path(table_dir, "symportal_depth_by_compartment.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  genus_percent_by_compartment,
  file = file.path(table_dir, "symportal_genus_percentage_by_compartment.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 17. SymPortal tables
# One table per compartment
# One table per collection_site x compartment

# 17A. Comparable IDs

meta$sample_uid <- as.character(meta$sample_uid)
rownames(postMED_counts) <- as.character(rownames(postMED_counts))

# 17B. Match datasets

common_samples <- intersect(rownames(postMED_counts), meta$sample_uid)

postMED_sub <- postMED_counts[common_samples, , drop = FALSE]
meta_sub <- meta[match(common_samples, meta$sample_uid), ]

# 17C. Dominant postMED sequence per sample

meta_sub$majority_sequence <- apply(postMED_sub, 1, function(x) {
  colnames(postMED_sub)[which.max(x)]
})

# 17D. timepoint: "2024_4" = year_month

tp <- strsplit(as.character(meta_sub$timepoint), "_")

meta_sub$timepoint_year <- as.numeric(sapply(tp, `[`, 1))
meta_sub$timepoint_month <- as.numeric(sapply(tp, `[`, 2))

# 17E. ColonyID order

fix_colony_id <- function(x) {
  x <- as.character(x)
  
  out <- x
  
  has_number <- grepl("_(\\d+)$", x)
  
  number_part <- sub(".*_(\\d+)$", "\\1", x[has_number])
  prefix_part <- sub("_(\\d+)$", "", x[has_number])
  
  out[has_number] <- paste0(
    prefix_part,
    "_",
    sprintf("%02d", as.numeric(number_part))
  )
  
  out
}

meta_sub$colonyID_fixed <- fix_colony_id(meta_sub$colonyID)

# 17F. Ordered column labels

meta_sub$col_label <- paste0(
  meta_sub$colonyID_fixed,
  "_",
  meta_sub$timepoint_year,
  "_",
  sprintf("%02d", meta_sub$timepoint_month)
)


safe_name <- function(x) {
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^A-Za-z0-9_\\-]", "", x)
  x
}

# 17G. Function for table from subset

build_table <- function(df) {
  
  df <- df[
    order(
      df$colonyID_fixed,
      df$timepoint_year,
      df$timepoint_month
    ),
  ]
  
  mat <- rbind(
    dominant_postMED_sequence = df$majority_sequence,
    ITS2_type_profile = df$ITS2_profiles
  )
  
  colnames(mat) <- df$col_label
  
  return(mat)
}

# 17H. One table per compartment

compartments <- sort(unique(meta_sub$sample_type))

tables_by_compartment <- list()

for (comp in compartments) {
  
  df_comp <- meta_sub[meta_sub$sample_type == comp, ]
  
  tables_by_compartment[[comp]] <- build_table(df_comp)
  
  write.table(
    tables_by_compartment[[comp]],
    file = file.path(table_dir, paste0(
      "symportal_colony_timepoint_table_",
      safe_name(comp),
      ".txt"
    )),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

# 17I. One table per collection_site x compartment

site_compartment_tables <- list()

sites <- sort(unique(meta_sub$collection_site))

for (site in sites) {
  
  for (comp in compartments) {
    
    df_site_comp <- meta_sub[
      meta_sub$collection_site == site &
        meta_sub$sample_type == comp,
    ]
    
    if (nrow(df_site_comp) == 0) {
      next
    }
    
    table_name <- paste(site, comp, sep = "__")
    
    site_compartment_tables[[table_name]] <- build_table(df_site_comp)
    
    write.table(
      site_compartment_tables[[table_name]],
      file = file.path(table_dir, paste0(
        "symportal_colony_timepoint_table_",
        safe_name(site),
        "_",
        safe_name(comp),
        ".txt"
      )),
      sep = "\t",
      quote = FALSE,
      col.names = NA
    )
  }
}

length(tables_by_compartment)

length(site_compartment_tables)

names(site_compartment_tables)

list.files(pattern = "symportal_colony_timepoint_table_.*\\.txt$")

# 18. Grid tables per collection_site x compartment

meta$sample_uid <- as.character(meta$sample_uid)
rownames(postMED_counts) <- as.character(rownames(postMED_counts))

# 2. Match metadata with postMED table

common_samples <- intersect(rownames(postMED_counts), meta$sample_uid)

postMED_sub <- postMED_counts[common_samples, , drop = FALSE]
meta_sub <- meta[match(common_samples, meta$sample_uid), ]

# 3. Dominant postMED sequence

meta_sub$majority_sequence <- apply(postMED_sub, 1, function(x) {
  colnames(postMED_sub)[which.max(x)]
})

# 4. timepoint: "2024_4"

tp <- strsplit(as.character(meta_sub$timepoint), "_")

meta_sub$timepoint_year <- as.numeric(sapply(tp, `[`, 1))
meta_sub$timepoint_month <- as.numeric(sapply(tp, `[`, 2))

meta_sub$timepoint_label <- paste0(
  meta_sub$timepoint_year,
  "_",
  sprintf("%02d", meta_sub$timepoint_month)
)

# 5. colonyID ordering

fix_colony_id <- function(x) {
  x <- as.character(x)
  out <- x
  
  has_number <- grepl("_(\\d+)$", x)
  
  number_part <- sub(".*_(\\d+)$", "\\1", x[has_number])
  prefix_part <- sub("_(\\d+)$", "", x[has_number])
  
  out[has_number] <- paste0(
    prefix_part,
    "_",
    sprintf("%02d", as.numeric(number_part))
  )
  
  out
}

meta_sub$colonyID_fixed <- fix_colony_id(meta_sub$colonyID)

safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^A-Za-z0-9_\\-]", "", x)
  x
}

# 6. grid table function

build_grid_table <- function(df) {
  
  df <- df[
    order(
      df$timepoint_year,
      df$timepoint_month,
      df$colonyID_fixed
    ),
  ]
  
  df$cell_value <- paste0(
    df$majority_sequence,
    " | ",
    df$ITS2_profiles
  )
  
  timepoints <- unique(df$timepoint_label)
  colonies <- unique(df$colonyID_fixed)
  
  grid <- matrix(
    NA,
    nrow = length(timepoints),
    ncol = length(colonies),
    dimnames = list(timepoints, colonies)
  )
  
  for (i in seq_len(nrow(df))) {
    grid[df$timepoint_label[i], df$colonyID_fixed[i]] <- df$cell_value[i]
  }
  
  grid <- as.data.frame(grid, check.names = FALSE)
  
  grid <- data.frame(
    timepoint = rownames(grid),
    grid,
    check.names = FALSE
  )
  
  rownames(grid) <- NULL
  
  grid
}

# 7. export grids

grid_tables_site_compartment <- list()

sites <- sort(unique(meta_sub$collection_site))
compartments <- sort(unique(meta_sub$sample_type))

for (site in sites) {
  
  for (comp in compartments) {
    
    df_site_comp <- meta_sub[
      meta_sub$collection_site == site &
        meta_sub$sample_type == comp,
    ]
    
    if (nrow(df_site_comp) == 0) next
    
    table_name <- paste(site, comp, sep = "__")
    
    grid_tables_site_compartment[[table_name]] <- build_grid_table(df_site_comp)
    
    write.table(
      grid_tables_site_compartment[[table_name]],
      file = file.path(table_dir, paste0(
        "symportal_grid_",
        safe_name(site),
        "_",
        safe_name(comp),
        ".txt"
      )),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
}


length(grid_tables_site_compartment)
names(grid_tables_site_compartment)

list.files(table_dir, pattern = "symportal_grid_.*\\.txt$")

# Site x date table

coral_df <- meta_sub[meta_sub$sample_type == "coral_field", ]

coral_df$ITS2_profiles <- as.character(coral_df$ITS2_profiles)

coral_df <- coral_df[
  !is.na(coral_df$ITS2_profiles) &
    coral_df$ITS2_profiles != "",
]

tp <- strsplit(as.character(coral_df$timepoint), "_")

coral_df$timepoint_year <- as.numeric(sapply(tp, `[`, 1))
coral_df$timepoint_month <- as.numeric(sapply(tp, `[`, 2))

coral_df$timepoint_label <- paste0(
  coral_df$timepoint_year,
  "_",
  sprintf("%02d", coral_df$timepoint_month)
)

most_common_profile <- function(x) {
  
  # Remove NA/empty
  x <- x[!is.na(x) & x != ""]
  
  # Frequency table
  tab <- sort(table(x), decreasing = TRUE)
  
  # Remove "un_id"
  tab_no_unid <- tab[names(tab) != "un_id"]
  
  # If all are un_id
  if (length(tab_no_unid) == 0) {
    return("un_id")
  }
  
  # Return most common non-un_id profile
  names(tab_no_unid)[1]
}

# Ordered rows and columns
dates <- unique(
  coral_df$timepoint_label[
    order(coral_df$timepoint_year, coral_df$timepoint_month)
  ]
)

sites <- sort(unique(coral_df$collection_site))

# Empty grid
site_date_ITS2_profile <- matrix(
  NA,
  nrow = length(dates),
  ncol = length(sites),
  dimnames = list(dates, sites)
)

# Fill grid
for (d in dates) {
  for (s in sites) {
    
    x <- coral_df$ITS2_profiles[
      coral_df$timepoint_label == d &
        coral_df$collection_site == s
    ]
    
    x <- x[!is.na(x) & x != ""]
    
    if (length(x) > 0) {
      site_date_ITS2_profile[d, s] <- most_common_profile(x)
    }
  }
}

# Convert to data frame
site_date_ITS2_profile <- as.data.frame(
  site_date_ITS2_profile,
  check.names = FALSE
)

site_date_ITS2_profile <- data.frame(
  date = rownames(site_date_ITS2_profile),
  site_date_ITS2_profile,
  check.names = FALSE
)

rownames(site_date_ITS2_profile) <- NULL

site_date_ITS2_profile

# Export
write.table(
  site_date_ITS2_profile,
  file = file.path(table_dir, "symportal_coral_most_common_ITS2_profile_by_site_date.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# Current manuscript supplementary-table exports
write.table(final_summary, file = file.path(manuscript_table_dir, "Table_S1_library_statistics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(site_date_ITS2_profile, file = file.path(manuscript_table_dir, "Table_S2_dominant_ITS2_profiles.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nGeneral sequence statistics complete.\n")
