\name{treestats-package}
\alias{treestats-package}
\alias{treestats}
\docType{package}
\title{
  Collection of phylogenetic tree statistics
}
\description{
  The 'treestats' package contains a collection of phylogenetic tree statistics,
  implemented in C++ to ensure high speed.
}
\details{
  Given a phylogenetic tree as a phylo object, the 'treestats' package provides
  a wide range of individual functions returning the relevant statistic. In
  addition, there are three functions available that calculate a collection
  of statistics at once: calc_all_statistics (which calculates all currently
  implemented statistics of treestats), calc_balance_stats, which calculates
  all (im)balance related statistics and calc_brts_stats, which calculates all
  branching times and branch length related statistics.
  Furthermore, there are a number of additional tools available that allow for
  phylogenetic tree manipulation: make_unbalanced_tree, which creates an
  imbalanced tree in a stepwise fashion. Then there are two functions related to
  conversion from and to an ltable, an alternative notation method used in some
  simulations. These are l_to_phylo which is a C++ based
  version of DDD::L2phylo, which converts an ltable to a phylo object, and
  phylo_to_l, which is a C+ based version of DDD::phylo2L, which converts a
  phylo object to an ltable. Lastly, the treestats package also includes a
  faster, C++ based, implementation of ape::branching.times (the
  function branching_times), which yields the same sequence of branching times,
  but omits the branching names in favour of speed.
}
\author{
Maintainer: Thijs Janzen <thijsjanzen@gmail.com>
}
\references{
  TBA
}

