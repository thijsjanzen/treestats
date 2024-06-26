% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beta.R
\name{beta_statistic}
\alias{beta_statistic}
\title{Aldous' beta statistic.}
\usage{
beta_statistic(
  phy,
  upper_lim = 10,
  algorithm = "COBYLA",
  abs_tol = 1e-04,
  rel_tol = 1e-06
)
}
\arguments{
\item{phy}{phylogeny or ltable}

\item{upper_lim}{Upper limit for beta parameter, default = 10.}

\item{algorithm}{optimization algorithm used, default is "COBYLA"
(Constrained Optimization BY Linear Approximations), also available are
"subplex" and "simplex". Subplex and Simplex seem to have difficulties with
unbalanced trees, e.g. if beta < 0.}

\item{abs_tol}{absolute stopping criterion of optimization. Default is 1e-4.}

\item{rel_tol}{relative stopping criterion of optimization. Default is 1e-6.}
}
\value{
Beta value
}
\description{
The Beta statistic fits a beta splitting model to each node,
assuming that the number of extant descendents of each daughter branch is
split following a beta distribution, such that the number of extant
descendentants x and y at a node follows \eqn{q(x, y) = s_n(beta)^-1
\frac{(gamma(x + 1 + beta)gamma(y + 1 + beta))}{gamma(x+1)gamma(y+1)}}, where
\eqn{s_n(beta)^-1} is a normalizing constant. When this model is fit to a
tree, different values of beta correspond to the expectation following from
different diversification models, such that a beta of 0 corresponds to a
Yule tree, a beta of -3/2 to a tree following from a PDA model. In general,
negative beta values correspond to trees more unbalanced than Yule trees, and
beta values larger than zero indicate trees more balanced than Yule trees.
The lower bound of the beta splitting parameter is -2.
}
\examples{
simulated_tree <- ape::rphylo(n = 100, birth = 1, death = 0)
balanced_tree <- treestats::create_fully_balanced_tree(simulated_tree)
unbalanced_tree <- treestats::create_fully_unbalanced_tree(simulated_tree)
beta_statistic(balanced_tree) # should be approximately 10
beta_statistic(simulated_tree) # should be near 0
beta_statistic(unbalanced_tree) # should be approximately -2
}
\references{
Aldous, David. "Probability distributions on cladograms." Random
discrete structures. Springer, New York, NY, 1996. 1-18.
Jones, Graham R. "Tree models for macroevolution and phylogenetic analysis."
Systematic biology 60.6 (2011): 735-746.
}
