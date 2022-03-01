#' function to simulate a birth-death tree, using a variety of packages
#' @param birth birth rate
#' @param death death rate
#' @param num_trees number of trees to simulate
#' @param max_t total time
#' @param max_lin total number of extant lineages at present
#' @param complete if true, returns a tree including extinct lineages
#' @param method chosen package used to simulate the tree, available are:
#' 'TreeSim', 'ape', 'TESS', 'DDD', 'geiger' and 'phytools'.
#' @return list of trees
#' @description when max_t is set, the tree is simulated conditional on
#' the maximum time. Similarly, when max_lin is set, the tree is simulated
#' conditional on the number of extant lineages. Although some packages (e.g.
#' TreeSim and TESS) provide functions to simulate conditional on both time
#' and number of lineages, this functionality is not universally available and
#' not supported here.
#' @export
sim_bd_tree <- function(birth, death, num_trees,
                        max_t = NULL, max_lin = NULL, complete = TRUE,
                        method = "TreeSim") {

  output <- list()

  if (method == "TreeSim") {
    if (!is.null(max_t)) {
      output <- TreeSim::sim.bd.age(age = max_t,
                                    numbsim = num_trees,
                                    lambda = birth,
                                    mu = death,
                                    mrca = TRUE,
                                    complete = complete)
    }
    if (!is.null(max_lin)) {
      output <- TreeSim::sim.bd.taxa(n = max_lin,
                                     numbsim = num_trees,
                                     lambda = birth,
                                     mu = death,
                                     complete = complete)
    }
  }

  if (method == "TESS") {
    if (!is.null(max_t)) {
      output <- TESS::tess.sim.age(n = num_trees,
                                   age = max_t,
                                   lambda = birth,
                                   mu = death)
    }
    if (!is.null(max_lin)) {
      output <- TESS::tess.sim.taxa(n = num_trees,
                                    nTaxa = max_lin,
                                    max = 100,
                                    lambda = birth,
                                    mu = death)
    }
  }

  if (method == "ape") {
     for (r in 1:num_trees) {
       output[[r]] <- ape::rphylo(n = max_lin,
                                  birth = birth,
                                  death = death,
                                  fossils = complete)
     }
  }

  if (method == "DDD") {
    for (r in 1:num_trees) {
      local_tree <- DDD::dd_sim(pars = c(birth, death, Inf),
                          age = max_t,
                          ddmodel = 1)
      if (complete) {
        output[[r]] <- local_tree$tas
      } else {
        output[[r]] <- local_tree$tes
      }
    }
  }

  if (method == "geiger") {
    for (r in 1:num_trees) {

      if (!is.null(max_t)) {
        local_tree <- geiger::sim.bdtree(b = birth,
                                       d = death,
                                       stop = "time",
                                       t = max_t)

      }
      if (!is.null(max_lin)) {
        local_tree <- geiger::sim.bdtree(b = birth,
                                         d = death,
                                         stop = "taxa",
                                         n = max_lin)
      }
      if (!complete) {
        local_tree <- geiger::drop.extinct(local_tree)
      }
      output[[r]] <- local_tree
    }
  }

  if (method == "phytools") {
    output <- phytools::pbtree(b = birth,
                               d = death,
                               n = max_lin,
                               t = max_t,
                               nsim = num_trees,
                               extant.only = complete)
  }


  return(output)
}
