calc_ultrametric_stats <- function(taxon_name, family_name, focal_tree) {
  cat(taxon_name, family_name)
  treestats::gamma_statistic(focal_tree)
  cat(" gamma success ")
  treestats::nLTT_base(focal_tree)
  cat("nLTT success ")
  treestats::laplacian_spectrum(focal_tree)
  cat("Laplace success\n")
}


#setwd("/Users/thijsjanzen/Documents/GitHub/treestats-scripts/Figure_3/")

begin <- "/Users/thijsjanzen/Documents/GitHub/treestats-scripts/"

# tree collection files
tree_collection_files <- c("datasets/phylogenies/fracced/amphibia_fracced.rds",
                           "datasets/phylogenies/fracced/birds_fracced.rds",
                           "datasets/phylogenies/fracced/ferns_fracced.rds",
                           "datasets/phylogenies/fracced/mammals_fracced.rds",
                           "datasets/phylogenies/fracced/ray_finned_fish_fracced.rds",
                           "datasets/phylogenies/fracced/sharks_fracced.rds",
                           "datasets/phylogenies/fracced/vascular_plants_fracced.rds")

taxa_names <- c("Amphibians", "Birds", "Ferns", "Mammals", "Ray finned Fish", "Cartaliginous Fish", "Vascular Plants")



for (i in 1:length(tree_collection_files)) {

  file_name <- paste0(begin, tree_collection_files[i])

  tree_collection <- readRDS(file_name)

  families <- names(tree_collection)
  cat(taxa_names[i], "\n")
  for (j in 1:length(tree_collection)) {

    focal_tree <- tree_collection[[j]]

    # if (length(focal_tree$tip.label) >= 10) {
      calc_ultrametric_stats(taxa_names[i], families[j], focal_tree)
   #  }
  }
}

cat("done!\n")