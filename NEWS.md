# Version 1.70.7
- Fixed numerical issue in calculating mean pair distance for very large trees
- Now use armadillo to calculate Eigen values
- added the function 'sshape', which was previously called the blum statistic (
the blum function remains available as well). Thanks Sophie Kersting for 
pointing this out!
- polished the manual

# Version 1.70.6
- Added support for unrooted trees
- Added check if tree is unrooted, and aborting if statistic does not support
this
- Updated reference table in the README file, showing which statistics are
available for unrooted trees, and which statistics (appear) to be affected by
root position
- using treebalance functions when treestats functions do not suffice (e.g.
for rquartet for polytomous trees)
- updated Wiener test, to include distance matrix based calculation

# Version 1.70.5

Updated references.

# Version 1.70.4

Removed C++20 dependency

# Version 1.70.3

Version numbering has improved to include the number of statistics available.

# Version 1.1.3

-   Reduced dependencies, Matrix and RSpectra are no longer required (but availability will improve speed!)
-   calc_all_stats, calc_topology_stats and calc_brts_stats no longer return a named list, but return a named vector, for ease of rapid calculation across many trees

# Version 1.1.2

-   Removed dependencies to nodeSub by integrating functions to generate fully balanced and unbalanced into treestats.
-   Renamed 'calc_balance_stats' to 'calc_topology_stats', and included only statistics that take the topology (without the branch lengths) into account.

# Version 1.1.1

Added the following statistics:

-   minimum(\>0) eigenvalue of the Adjacency matrix

-   maximum eigenvalue of the Adjacency matrix

-   minimum(\>0) eigenvalue of the Laplacian matrix

-   maximum eigenvalue of the Laplacian matrix

-   double cherries

-   root imbalance

-   four prong

Renamed eigenvector to eigen_centrality (in line with the analogous function in igraph)

# Version 1.1.0

Add several new statistics:

-   Colless-corrected

-   Colless-Quadratic

-   Total path index

-   Total internal path index

-   Average vertex depth

-   Max Width over Max Depth

# Version 1.0.7

-   Added checks to each statistic verifying the statistic requires an ultrametric or binary tree - this should avoid some rare instances where memory access violations would pop up when providing a non-binary tree to a statistic assuming a binary input phylogeny. Many thanks to Fien Strijthaegen for pointing this out.

# Version 1.0.6

-   Added treeness statistic

# Version 1.0.5

-   Squashed a bug in imbalance_steps that would incorrectly calculate for smaller trees
-   Added bioRxiv link to DESCRIPTION
-   Added CITATION file

# Version 1.0.4

-   Corrected wording of variation in branch length statistics to correctly reflect variance, instead of variation

# Version 1.0.3

First release to CRAN
