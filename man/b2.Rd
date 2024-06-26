% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/b2.R
\name{b2}
\alias{b2}
\title{B2 metric}
\usage{
b2(phy, normalization = "none")
}
\arguments{
\item{phy}{phylo object or ltable}

\item{normalization}{"none" or "yule", when "yule" is chosen, the statistic
is divided by the Yule expectation, following from theorem 3.7 in Bienvenu
2020.}
}
\value{
Maximum depth (in number of edges)
}
\description{
Balance metric that uses the Shannon-Wiener statistic of
information content. The b2 measure is given by the sum over the depths of
all tips, divided by 2^depth: sum Ni / 2^Ni.
Although also defined on non-binary trees, the treestats package only
provides code for binary trees.
}
\references{
K.-T. Shao and R. R. Sokal. Tree Balance.
Systematic Zoology, 39(3):266, 1990. doi: 10.2307/2992186.

Bienvenu, François, Gabriel Cardona, and Celine Scornavacca.
"Revisiting Shao and Sokal’s $$ B_2 $$ B 2 index of phylogenetic balance."
Journal of Mathematical Biology 83.5 (2021): 1-43.
}
