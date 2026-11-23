# tree_size

### Tree Size

Most summary statistics for phylogenetic trees are sensitive to tree
size. Fortunately, some statistics provde corrections for tree size,
based on the expectation for the statistic under the Yule model. Let’s
explore to what extent these corrections are effective, and to what
extent the non-corrected statistics are influenced by tree size.

#### Simulating trees

To simulate trees, we simulate trees under the Yule model, where we vary
the size of the tree, and we also vary the crown age, in order to sample
trees from a similar distribution. We keep the speciation rate constant.
Alternatively, we could also keep crown age constant and vary the
speciation rate.

``` r

found <- c()
found_min <- c()
found_max <- c()
for (tree_size in ceiling(10^seq(1.75, 3, length.out = 10))) {
  cat(tree_size, "\n")
  stats <- c()

  # more replicates for smaller trees:
  num_repl <- min(100, ceiling(5 + 1000 / tree_size))
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
  stats2 <- apply(stats, 2, mean)

  min_stats <- apply(stats, 2,
                     function(x) return(quantile(x, na.rm = TRUE, 0.025)))
  max_stats <- apply(stats, 2,
                     function(x) return(quantile(x, na.rm = TRUE, 0.975)))

  found <- rbind(found, stats2)
  found_min <- rbind(found_min, min_stats)
  found_max <- rbind(found_max, max_stats)
}
```

    ## 57 
    ## 78 
    ## 107 
    ## 147 
    ## 203 
    ## 279 
    ## 384 
    ## 528 
    ## 727 
    ## 1000

``` r

found <- as_tibble(found)
found_min <- as_tibble(found_min)
found_max <- as_tibble(found_max)
```

#### Yule trees

Before plotting relationships, we collect information on the type of
normalization applied to each statistic.

``` r

norm_info <- read.table("https://raw.githubusercontent.com/thijsjanzen/treestats-scripts/main/datasets/families_2026.txt", header = TRUE) # nolint
colnames(norm_info)[colnames(norm_info) == "Statistic"] <- "statistic"

found2 <- found |>
  tidyr::gather(key = "statistic", value = "mean", -number_of_lineages)

found_min2 <- found_min |>
  tidyr::gather(key = "statistic", value = "min", -number_of_lineages)
found_max2 <- found_max |>
  tidyr::gather(key = "statistic", value = "max", -number_of_lineages)

found3 <- left_join(found2, norm_info,
                    by = join_by(statistic))
found3 <- left_join(found3, found_min2,
                    by = join_by(statistic, number_of_lineages))
found3 <- left_join(found3, found_max2,
                    by = join_by(statistic, number_of_lineages))
```

Now, we can plot the relationship between tree size and the statistic:

``` r

found3 |>
  filter(Normalization == "Yule") |>
  ggplot(aes(x = number_of_lineages, y = mean)) +
  geom_ribbon(aes(ymin = min, ymax = max), fill = "purple", alpha = 0.2) +
  geom_point(col = "blue") +
  geom_line(col = "gold") +
  scale_x_log10() +
  facet_wrap(~statistic, scales = "free_y", ncol = 5) +
  theme_minimal() +
  xlab("Extant tree size") +
  ggtitle("Statistics normalized by the Yule expectation")
```

![](Tree_size_files/figure-html/yule_plots-1.png)

We see that variation in statistic values is reduced, but not zero.

#### Normalization by tree size

For some statistics, no expectations based on the Yule model are
available, and instead, we can normalize by the number of tips. We then
find:

``` r

found3 |>
  filter(Normalization == "Tips") |>
  ggplot(aes(x = number_of_lineages, y = mean)) +
  geom_ribbon(aes(ymin = min, ymax = max), fill = "purple", alpha = 0.2) +
  geom_point(col = "blue") +
  geom_line(col = "gold") +
  scale_x_log10() +
  facet_wrap(~statistic, scales = "free_y", ncol = 5) +
  theme_minimal() +
  xlab("Extant tree size") +
  ggtitle("Statistics normalized by tree size")
```

![](Tree_size_files/figure-html/size_plots-1.png)

Clearly, this is not perfect either.

#### No normalization

Lastly, we have those statistics for which no normalization is available

``` r

found3 |>
  filter(Normalization == "None") |>
  ggplot(aes(x = number_of_lineages, y = mean)) +
  geom_ribbon(aes(ymin = min, ymax = max), fill = "purple", alpha = 0.2) +
  geom_point(col = "blue") +
  geom_line(col = "gold") +
  scale_x_log10() +
  facet_wrap(~statistic, scales = "free_y", ncol = 7) +
  theme_minimal() +
  xlab("Extant tree size") +
  ggtitle("Statistics not normalized")
```

![](Tree_size_files/figure-html/no_norm-1.png)

Here we now see clear patterns with tree size. This is worrying, as for
some statistics, patterns might arise that are simply a side-effect of
creating larger trees, rather than directly influencly the statistic of
choice. Overall, it appears there are few statistics insensitive to tree
size. When the Yule model correction is used, this tends to work when we
are looking at trees of varying size that generated under the Yule
model, but when the generating model is different, the correction no
longer works. See the original treestats paper (Janzen & Etienne 2024)
for a more thorough study on this.

## Normalization on non-Yule trees

When trees are not simulated using Yule, but using another model,
e.g. DDD, we see that results no longer hold up.

``` r

found <- c()
found_min <- c()
found_max <- c()
for (tree_size in ceiling(10^seq(1.75, 3, length.out = 10))) {
  cat(tree_size, "\n")

  # more replicates for smaller trees:
  num_repl <- min(100, ceiling(5 + 1000 / tree_size))
  speciation <- 0.5

  # calculate the expected crown age for the relevant tree size:
  max_t <- (1 / speciation) * log(tree_size / 2)

  K <- tree_size * 1.2

  stats <- c()
  for (r in 1:num_repl) {
    tree <- DDD::dd_sim(pars = c(speciation, 0, tree_size),
                        age = max_t,
                        ddmodel = 1)$tes

    stats <- rbind(stats, treestats::calc_all_stats(tree, normalize = TRUE))
  }
  stats2 <- apply(stats, 2, mean)

  min_stats <- apply(stats, 2,
                     function(x) return(quantile(x, na.rm = TRUE, 0.025)))
  max_stats <- apply(stats, 2,
                     function(x) return(quantile(x, na.rm = TRUE, 0.975)))

  stats2[names(stats2) == "number_of_lineages"]       <- tree_size
  min_stats[names(min_stats) == "number_of_lineages"] <- tree_size
  max_stats[names(max_stats) == "number_of_lineages"] <- tree_size

  found <- rbind(found, stats2)
  found_min <- rbind(found_min, min_stats)
  found_max <- rbind(found_max, max_stats)
}
```

    ## 57 
    ## 78

    ## Warning in calc_beta_cpp(phy, upper_lim, algorithm, abs_tol, rel_tol): Trees
    ## with only two tips have undefined beta

    ## Balance of a tree with two tips is undefined

    ## Warning in func(phylo): I statistic is only available for trees with at least 4
    ## tips.

    ## Warning in func(phylo, normalization = ifelse(normalize, norm_res[1],
    ## norm_res[2])): normalization not valid for trees of size 2

    ## 107 
    ## 147 
    ## 203 
    ## 279 
    ## 384 
    ## 528 
    ## 727 
    ## 1000

``` r

found <- as_tibble(found)
found_min <- as_tibble(found_min)
found_max <- as_tibble(found_max)

found2 <- found |>
  tidyr::gather(key = "statistic", value = "mean", -number_of_lineages)

found_min2 <- found_min |>
  tidyr::gather(key = "statistic", value = "min", -number_of_lineages)
found_max2 <- found_max |>
  tidyr::gather(key = "statistic", value = "max", -number_of_lineages)


found3 <- left_join(found2, found_min2,
                    by = join_by(statistic, number_of_lineages))
found3 <- left_join(found3, found_max2,
                    by = join_by(statistic, number_of_lineages))

found3 <- left_join(found3, norm_info,
                    by = join_by(statistic))
```

And now we can plot things.

``` r

found3 |>
  filter(Normalization == "Yule") |>
  ggplot(aes(x = number_of_lineages, y = mean)) +
  geom_ribbon(aes(ymin = min, ymax = max), fill = "purple", alpha = 0.2) +
  geom_point(col = "blue") +
  geom_line(col = "gold") +
  scale_x_log10() +
  facet_wrap(~statistic, scales = "free_y", ncol = 3) +
  theme_minimal() +
  xlab("Extant tree size") +
  ggtitle("Statistics normalized by the Yule expectation")
```

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](Tree_size_files/figure-html/yule_plots_ddd-1.png)

``` r

found3 |>
  filter(Normalization == "Tips") |>
  ggplot(aes(x = number_of_lineages, y = mean)) +
  geom_ribbon(aes(ymin = min, ymax = max), fill = "purple", alpha = 0.2) +
  geom_point(col = "blue") +
  geom_line(col = "gold") +
  scale_x_log10() +
  facet_wrap(~statistic, scales = "free_y", ncol = 5) +
  theme_minimal() +
  xlab("Extant tree size") +
  ggtitle("Statistics normalized by tree size")
```

    ## Warning: Removed 3 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](Tree_size_files/figure-html/yule_plots_ddd-2.png)

``` r

found3 |>
  filter(Normalization == "None") |>
  ggplot(aes(x = number_of_lineages, y = mean)) +
  geom_ribbon(aes(ymin = min, ymax = max), fill = "purple", alpha = 0.2) +
  geom_point(col = "blue") +
  geom_line(col = "gold") +
  scale_x_log10() +
  facet_wrap(~statistic, scales = "free_y", ncol = 7) +
  theme_minimal() +
  xlab("Extant tree size") +
  ggtitle("Statistics not normalized")
```

    ## Warning: Removed 8 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](Tree_size_files/figure-html/yule_plots_ddd-3.png)
