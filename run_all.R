# Run all R analyses and manuscript figures from the repository root.

source(file.path("analysis", "general_statistics", "general_sequence_statistics.R"))
source(file.path("analysis", "thermal_density", "thermal_density_analysis.R"))
source(file.path("figures", "postMED_composition", "plot_postMED_composition.R"))
source(file.path("analysis", "community_multivariate", "community_multivariate_analysis.R"))

source(file.path("figures", "marker_concordance", "plot_marker_tanglegrams.R"))
source(file.path("figures", "gunifrac_tree", "plot_GUniFrac_type_profile_tree.R"))
source(file.path("figures", "phylogenetic_composition", "plot_phylogenetic_composition.R"))

cat("\nAll R workflows completed.\n")
