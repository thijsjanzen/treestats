# Root imbalance

Measures the distribution of tips over the two crown lineages, e.g. \\n1
/ (n1 + n2)\\, where n1 is the number of tips connected to crown lineage
1 and n2 is the number of tips connected to crown lineage 2. We always
take n1 \> n2, thus root imbalance is always in \[0.5, 1\].

## Usage

``` r
root_imbalance(phy)
```

## Arguments

- phy:

  phylo object or ltable

## Value

Root imbalance

## References

Guyer, Craig, and Joseph B. Slowinski. "Adaptive radiation and the
topology of large phylogenies." Evolution 47.1 (1993): 253-263.
