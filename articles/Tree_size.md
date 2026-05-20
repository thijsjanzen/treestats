# tree_size

## Tree Size

Most summary statistics for phylogenetic trees are sensitive to tree
size. Fortunately, some statistics provde corrections for tree size,
based on the expectation for the statistic under the Yule model. Let’s
explore to what extent these corrections are effective, and to what
extent the non-corrected statistics are influenced by tree size.

### Simulating trees

To simulate trees, we simulate trees under the Yule model, where we vary
the size of the tree, and we also vary the crown age, in order to sample
trees from a similar distribution. We keep the speciation rate constant.
Alternatively, we could also keep crown age constant and vary the
speciation rate.

``` r

found <- c()
for (tree_size in ceiling(10^seq(1, 3, length.out = 10))) {
  cat(tree_size, "\n")
  stats <- c()

  # more replicates for smaller trees:
  num_repl <- min(100, ceiling(1 + 1000 / tree_size))
  speciation <- 0.5

  # calculate the expected crown age for the relevant tree size:
  max_t <- (1 / speciation) * log(tree_size / 2)

  for (r in 1:num_repl) {
    tree <- TreeSim::sim.bd.taxa.age(n = tree_size,
                                     numbsim = 1,
                                     lambda = speciation,
                                     mu = 0,
                                     age = max_t,
                                     mrca = TRUE)[[1]]
    stats <- rbind(stats, treestats::calc_all_stats(tree, normalize = TRUE))
  }
  stats2 <- apply(stats, 2, median)
  found <- rbind(found, stats2)
}
```

    ## 10

    ## Loading required namespace: RSpectra

    ## 17 
    ## 28 
    ## 47 
    ## 78 
    ## 130 
    ## 216 
    ## 360 
    ## 600 
    ## 1000

``` r

found <- as_tibble(found)
```

### Yule trees

Before plotting relationships, we collect information on the type of
normalization applied to each statistic.

``` r

norm_info <- read.table("https://raw.githubusercontent.com/thijsjanzen/treestats-scripts/main/datasets/normalize.txt", header = TRUE) # nolint

found2 <- found |>
  tidyr::gather(key = "statistic", value = "value", -number_of_lineages)

found3 <- left_join(found2, norm_info)
```

    ## Joining with `by = join_by(statistic)`

Now, we can plot the relationship between tree size and the statistic:

``` r

found3 |>
  filter(normalize == TRUE & type == "Yule") |>
  ggplot(aes(x = number_of_lineages, y = value)) +
  geom_point(col = "blue", alpha = 0.3) +
  geom_line(col = "gold") +
  scale_x_log10() +
  facet_wrap(~statistic, scales = "free_y", ncol = 3) +
  theme_minimal() +
  xlab("Extant tree size") +
  ggtitle("Statistics normalized by the Yule expectation")
```

![](Tree_size_files/figure-html/yule_plots-1.png)

We see that variation in statistic values is reduced, but not zero.

### Normalization by tree size

For some statistics, no expectations based on the Yule model are
available, and instead, we can normalize by the number of tips. We then
find:

``` r

found3 |>
  filter(normalize == TRUE & type == "Tips") |>
  ggplot(aes(x = number_of_lineages, y = value)) +
  geom_point(col = "blue", alpha = 0.3) +
  geom_line(col = "gold") +
  scale_x_log10() +
  facet_wrap(~statistic, scales = "free_y", ncol = 4) +
  theme_minimal() +
  xlab("Extant tree size") +
  ggtitle("Statistics normalized by tree size")
```

![](Tree_size_files/figure-html/size_plots-1.png)

Clearly, this is not perfect either.

### No normalization

Lastly, we have those statistics for which no normalization is available

``` r

found3 |>
  filter(normalize == FALSE) |>
  ggplot(aes(x = number_of_lineages, y = value)) +
  geom_point(col = "blue", alpha = 0.3) +
  geom_line(col = "gold") +
  scale_x_log10() +
  facet_wrap(~statistic, scales = "free_y", ncol = 5) +
  theme_minimal() +
  xlab("Extant tree size") +
  ggtitle("Statistics not normalized")
```

![](Tree_size_files/figure-html/no_norm-1.png)

Here we now see clear patterns with tree size.
