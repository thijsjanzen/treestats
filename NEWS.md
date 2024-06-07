# Version 1.1.0

Add several new statistics:
  - Colless-corrected
  - Colless-Quadratic
  - Total path index
  - Total internal path index
  - Average vertex depth
  - Max Width over Max Depth

# Version 1.0.7
- Added checks to each statistic verifying the statistic requires an ultrametric
or binary tree - this should avoid some rare instances where memory access 
violations would pop up when providing a non-binary tree to a statistic 
assuming a binary input phylogeny. Many thanks to Fien Strijthaegen for 
pointing this out.

# Version 1.0.6
- Added treeness statistic

# Version 1.0.5
- Squashed a bug in imbalance_steps that would incorrectly calculate for smaller
trees
- Added bioRxiv link to DESCRIPTION
- Added CITATION file

# Version 1.0.4
- Corrected wording of variation in branch length statistics to correctly
reflect variance, instead of variation

# Version 1.0.3
First release to CRAN
