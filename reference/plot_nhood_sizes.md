# Plot the distribution of neighbourhood sizes

Visualises the distribution of cell counts per neighbourhood as a
histogram, following neighbourhood assignment with `create_communities`.

## Usage

``` r
plot_nhood_sizes(communities_summary)
```

## Arguments

- communities_summary:

  A data frame summarising neighbourhoods, as returned by
  `create_communities`, containing at minimum a numeric column `n_cells`
  giving the number of cells assigned to each neighbourhood.

## Value

A `ggplot` object showing a histogram of neighbourhood sizes (x-axis:
neighbourhood size; y-axis: count).

## See also

[`create_communities`](https://grace-hattersley.github.io/TEVRA/reference/create_communities.md)
to create communities and obtain communities summary
