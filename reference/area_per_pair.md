# Area per pair index

The area per pair index calculates the sum of the number of edges on the
path between all two leaves. Instead, the area per pair index (APP) can
also be derived from the Sackin (S) and total cophenetic index (TC): \\
APP = \frac{2}{n}\cdot S - \frac{4}{n(n-1)}\cdot TC\\ \\APP = 2/n \* S -
4/(n(n-1)) \* TC\\

## Usage

``` r
area_per_pair(phy, normalization = "none")
```

## Arguments

- phy:

  phylo object or ltable

- normalization:

  "none" or "yule", in which case the acquired result is divided by the
  expectation for the Yule model.

## Value

Area per pair index

## References

T. Araújo Lima, F. M. D. Marquitti, and M. A. M. de Aguiar. Measuring
Tree Balance with Normalized Tree Area. arXiv e-prints, art.
arXiv:2008.12867, 2020.
