# psbA / ITS2 tanglegram comparison
#   1. combined_tanglegrams_vertical
#   2. entanglement_summary
#   3. mantel_summary
#   4. dominant_type_profiles_by_sample table

input_dir <- file.path("data", "marker_concordance")
figure_dir <- file.path("outputs", "figures")
table_dir <- file.path("outputs", "marker_concordance", "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# 1. LOAD LIBRARIES

library(Biostrings)
library(ape)
library(kmer)
library(GUniFrac)
library(vegan)
library(dendextend)
library(magick)

# 2. OUTPUT FOLDERS

temp_panel_dir <- file.path(tempdir(), "tanglegram_panels")
dir.create(temp_panel_dir, showWarnings = FALSE)

# 3. READ INPUT FILES

psba_F_seqs <- readDNAStringSet(file.path(input_dir, "psbAncr_forward_all.fasta"))
psba_R_seqs <- readDNAStringSet(file.path(input_dir, "psbAncr_reverse_all.fasta"))
its2_seqs   <- readDNAStringSet(file.path(input_dir, "ITS2_postMED_seqs.fasta"))

its2_tab <- read.table(
  file.path(input_dir, "ITS2_postMED_relative_abundances.txt"),
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

type_tab <- read.table(
  file.path(input_dir, "ITS2_type_profiles_relative_abundances.txt"),
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

# 4. CLEAN ITS2 TABLES


its2_tab <- its2_tab[, colSums(its2_tab) > 0, drop = FALSE]
its2_seqs <- its2_seqs[names(its2_seqs) %in% colnames(its2_tab)]

# Identify dominant profiles directly from their relative abundances.

# 5. FUNCTIONS

get_kmer_dist <- function(dna_set, k = 7) {
  
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
  
  attr(dist_mat, "Labels") <- names(dna_char)
  
  dist_mat
}

build_kmer_tree <- function(dna_set, k = 7) {
  
  dist_mat <- get_kmer_dist(dna_set, k = k)
  
  as.phylo(
    hclust(as.dist(dist_mat), method = "average")
  )
}

prepare_tanglegram_pair <- function(dend1, dend2) {
  
  common_labels <- intersect(labels(dend1), labels(dend2))
  
  dend1 <- dendextend::prune(
    dend1,
    setdiff(labels(dend1), common_labels)
  )
  
  dend2 <- dendextend::prune(
    dend2,
    setdiff(labels(dend2), common_labels)
  )
  
  dend_pair <- dendlist(dend1, dend2)
  dend_pair <- untangle(dend_pair, method = "step2side")
  
  list(
    dend1 = dend_pair[[1]],
    dend2 = dend_pair[[2]],
    entanglement = entanglement(dend_pair[[1]], dend_pair[[2]])
  )
}

export_single_tanglegram_panel <- function(tg,
                                           panel_letter,
                                           ent_value,
                                           file_path,
                                           width = 14,
                                           height = 6,
                                           dpi = 300,
                                           lab_cex = 1.7) {
  
  png(
    filename = file_path,
    width = width,
    height = height,
    units = "in",
    res = dpi
  )
  
  par(mar = c(4, 2, 2, 2))
  
  tanglegram(
    tg$dend1,
    tg$dend2,
    common_subtrees_color_lines = TRUE,
    highlight_distinct_edges = FALSE,
    margin_inner = 5,
    lwd = 0.8,
    lab.cex = lab_cex,
    columns_width = c(5, 3, 5),
    main = ""
  )
  
  mtext(
    panel_letter,
    side = 3,
    adj = 0,
    line = 0.2,
    font = 2,
    cex = 1.7
  )
  
  mtext(
    paste0("Entanglement = ", round(ent_value, 3)),
    side = 1,
    line = 2.2,
    cex = 1.7
  )
  
  dev.off()
}

# 6. SEQUENCE TREES

psba_F_tree <- build_kmer_tree(psba_F_seqs, k = 7)
psba_R_tree <- build_kmer_tree(psba_R_seqs, k = 7)
its2_tree   <- build_kmer_tree(its2_seqs, k = 7)

# 7. SAMPLE-LEVEL psbA TREES

psba_F_dist_matrix <- get_kmer_dist(psba_F_seqs, k = 7)
psba_R_dist_matrix <- get_kmer_dist(psba_R_seqs, k = 7)

psba_F_sample_tree <- as.phylo(
  hclust(as.dist(psba_F_dist_matrix), method = "average")
)

psba_R_sample_tree <- as.phylo(
  hclust(as.dist(psba_R_dist_matrix), method = "average")
)

# 8. SAMPLE-LEVEL ITS2 TREE USING GUniFrac

common_otus <- intersect(colnames(its2_tab), its2_tree$tip.label)

its2_tab_guni <- its2_tab[, common_otus, drop = FALSE]

its2_tree_guni <- drop.tip(
  its2_tree,
  setdiff(its2_tree$tip.label, common_otus)
)

its2_tab_guni <- its2_tab_guni[, its2_tree_guni$tip.label, drop = FALSE]

cat("ITS2 OTU match check:",
    all(colnames(its2_tab_guni) == its2_tree_guni$tip.label), "\n")

cat("Number of ITS2 OTUs retained for GUniFrac:",
    length(common_otus), "\n\n")

its2_guni <- GUniFrac(
  its2_tab_guni,
  its2_tree_guni,
  alpha = 0.5
)

its2_dist <- its2_guni$unifracs[, , "d_0.5"]

its2_sample_tree <- as.phylo(
  hclust(as.dist(its2_dist), method = "average")
)

# 9. TANGLEGRAMS

tg_a <- prepare_tanglegram_pair(
  as.dendrogram(psba_F_tree),
  as.dendrogram(psba_R_tree)
)

tg_b <- prepare_tanglegram_pair(
  as.dendrogram(its2_sample_tree),
  as.dendrogram(psba_F_sample_tree)
)

tg_c <- prepare_tanglegram_pair(
  as.dendrogram(its2_sample_tree),
  as.dendrogram(psba_R_sample_tree)
)

# 10. ENTANGLEMENT SUMMARY TABLE

entanglement_summary <- data.frame(
  Comparison = c(
    "psbA Forward vs Reverse",
    "ITS2 vs psbA Forward",
    "ITS2 vs psbA Reverse"
  ),
  Entanglement = c(
    tg_a$entanglement,
    tg_b$entanglement,
    tg_c$entanglement
  )
)

print(entanglement_summary)

write.table(
  entanglement_summary,
  file = file.path(table_dir, "entanglement_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# 11. MANTEL TESTS

its2_sample_dist <- as.dist(its2_dist)
psba_sample_dist_R <- as.dist(psba_R_dist_matrix)
psba_sample_dist_F <- as.dist(psba_F_dist_matrix)

common_samples <- Reduce(
  intersect,
  list(
    attr(its2_sample_dist, "Labels"),
    attr(psba_sample_dist_R, "Labels"),
    attr(psba_sample_dist_F, "Labels")
  )
)

its2_sample_dist_common <- as.dist(
  as.matrix(its2_sample_dist)[common_samples, common_samples]
)

psba_sample_dist_R_common <- as.dist(
  as.matrix(psba_sample_dist_R)[common_samples, common_samples]
)

psba_sample_dist_F_common <- as.dist(
  as.matrix(psba_sample_dist_F)[common_samples, common_samples]
)

set.seed(123)

mantel_R <- mantel(
  its2_sample_dist_common,
  psba_sample_dist_R_common,
  permutations = 9999,
  method = "pearson"
)

mantel_F <- mantel(
  its2_sample_dist_common,
  psba_sample_dist_F_common,
  permutations = 9999,
  method = "pearson"
)

mantel_summary <- data.frame(
  Comparison = c(
    "ITS2 vs psbA Reverse",
    "ITS2 vs psbA Forward"
  ),
  Mantel_r = c(
    mantel_R$statistic,
    mantel_F$statistic
  ),
  p_value = c(
    mantel_R$signif,
    mantel_F$signif
  ),
  Permutations = 9999
)

print(mantel_summary)

write.table(
  mantel_summary,
  file = file.path(table_dir, "mantel_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# 12. DOMINANT ITS2 TYPE PROFILE TABLE

common_samples_type <- intersect(
  rownames(type_tab),
  its2_sample_tree$tip.label
)

type_tab <- type_tab[common_samples_type, , drop = FALSE]

dominant_type <- colnames(type_tab)[
  max.col(type_tab, ties.method = "first")
]

names(dominant_type) <- rownames(type_tab)

type_labels <- dominant_type[
  match(its2_sample_tree$tip.label, names(dominant_type))
]

type_summary <- data.frame(
  Sample = its2_sample_tree$tip.label,
  Dominant_Type = type_labels
)

print(type_summary)

write.table(
  type_summary,
  file = file.path(table_dir, "dominant_type_profiles_by_sample.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# 13. EXPORT SINGLE TANGLEGRAM PANELS

panel_a_png <- file.path(temp_panel_dir, "panel_a_tanglegram.png")
panel_b_png <- file.path(temp_panel_dir, "panel_b_tanglegram.png")
panel_c_png <- file.path(temp_panel_dir, "panel_c_tanglegram.png")

export_single_tanglegram_panel(
  tg = tg_a,
  panel_letter = "(a)",
  ent_value = tg_a$entanglement,
  file_path = panel_a_png,
  lab_cex = 1.7
)

export_single_tanglegram_panel(
  tg = tg_b,
  panel_letter = "(b)",
  ent_value = tg_b$entanglement,
  file_path = panel_b_png,
  lab_cex = 1.7
)

export_single_tanglegram_panel(
  tg = tg_c,
  panel_letter = "(c)",
  ent_value = tg_c$entanglement,
  file_path = panel_c_png,
  lab_cex = 1.7
)

# 14. MERGED PANELS VERTICALLY

img_a <- image_read(panel_a_png)
img_b <- image_read(panel_b_png)
img_c <- image_read(panel_c_png)

combined_img <- image_append(
  c(img_a, img_b, img_c),
  stack = TRUE
)

image_write(
  combined_img,
  path = file.path(figure_dir, "Figure_S1_marker_tanglegrams.tiff"),
  format = "tiff",
  density = "600x600"
)

cat("\nFinal TIFF saved to:\n")
cat(file.path(figure_dir, "Figure_S1_marker_tanglegrams.tiff"), "\n\n")

cat("Tables saved to:\n")
cat(table_dir, "\n")
