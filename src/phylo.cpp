// Copyright 2022 - 2024 Thijs Janzen
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
//
//


#include <vector>
#include <array>
#include <Rcpp.h>

#include "util.h"        // NOLINT [build/include_subdir]
#include "phylo2L.h"     // NOLINT [build/include_subdir]
#include "L2newick.h"    // NOLINT [build/include_subdir]

#include "imbalance_steps.h"  // NOLINT [build/include_subdir]

//' Function to generate an ltable from a phy object.
//' @description This function is a C++
//' implementation of the function DDD::phylo2L. An L table summarises a
//' phylogeny in a table with four columns, being: 1) time at which a species
//' is born, 2) label of the parent of the species, where positive and negative
//' numbers indicate whether the species belongs to the left or right crown
//' lineage, 3) label of the daughter species itself (again positive or
//' negative depending on left or right crown lineage), and the last column 4)
//' indicates the time of extinction of a species, or -1 if the species is
//' extant.
//' @param phy phylo object
//' @return ltable (see description)
//' @export
//' @examples
//' simulated_tree <- ape::rphylo(n = 4, birth = 1, death = 0)
//' ltable <- phylo_to_l(simulated_tree)
//' reconstructed_tree <- DDD::L2phylo(ltable)
//' old_par <- par()
//' par(mfrow = c(1, 2))
//' # trees should be more or less similar, although labels may not match, and
//' # rotations might cause (initial) visual mismatches
//' plot(simulated_tree)
//' plot(reconstructed_tree)
//' par(old_par)
// [[Rcpp::export]]
Rcpp::NumericMatrix phylo_to_l(const Rcpp::List& phy) {
    const size_t ncol = 4;
    std::vector< std::array< double, ncol> > ltab = phylo_to_l_cpp(phy);
    size_t nrow = ltab.size();
    Rcpp::NumericMatrix out(nrow, ncol);
    for (size_t i = 0; i < ltab.size(); ++i) {
        for (size_t j = 0; j < ncol; ++j) {
            out(i, j) = ltab[i][j];
        }
    }
    return out;
}

// [[Rcpp::export]]
std::string l_to_newick(const Rcpp::NumericMatrix& ltable_R,
                        bool drop_extinct) {
    auto ltable_cpp = convert_to_ltable(ltable_R);
    auto newick_string = ltable_to_newick(ltable_cpp, drop_extinct);
    return newick_string;
}

// [[Rcpp::export]]
double imbalance_steps_cpp(const Rcpp::NumericMatrix& ltable_R,
                        bool normalization) {
    try {
       auto ltable_cpp = convert_to_ltable(ltable_R);
       return imbal_steps::number_of_steps(ltable_cpp, normalization);
    } catch(std::exception &ex) {
       forward_exception_to_r(ex);
    } catch (const char* msg) {
       Rcpp::Rcout << msg << std::endl;
    } catch(...) {
       ::Rf_error("c++ exception (unknown reason)");
    }
    return NA_REAL;
}
