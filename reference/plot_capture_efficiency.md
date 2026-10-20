# Plot VDJ contig capture efficiency across neighbourhoods

Generates two diagnostic plots summarising how efficiently VDJ contigs
were captured across transcriptional neighbourhoods: a jittered boxplot
of per-neighbourhood contig capture rate, and a scatterplot of contig
count against neighbourhood size (with a linear fit) to check whether
capture scales as expected with neighbourhood size.

## Usage

``` r
plot_capture_efficiency(nhood_diversity_df, group_by = NULL)
```

## Arguments

- nhood_diversity_df:

  The neighbourhood-level metadata data frame generated as the output of
  nhood_diversity, containing at minimum the columns
  `contig_capture_rate`, `n_contigs` and `nhood_size`.

- group_by:

  Character scalar giving the name of a column in `nhood_diversity_df`
  to colour/facet points by (e.g. cell type, timepoint, or tissue). If
  `NULL` (the default) or not found in `colnames(nhood_diversity_df)`,
  all neighbourhoods are treated as a single "All" group.

## Value

A named list of two `ggplot` nhood_diversity_dfs:

- capture_rate:

  Jittered boxplot of `contig_capture_rate` by `group_by`, with point
  size scaled to `nhood_size`.

- count_vs_size:

  Scatterplot of `n_contigs` against `nhood_size`, coloured by
  `group_by`, with a linear trend line.

## Details

If `group_by` is supplied but is not already a factor, it is coerced to
one internally; the input `nhood_diversity_df` itself is not modified.
This function is purely a visualisation wrapper and does not perform any
capture-rate calculation itself.

## See also

[`nhood_diversity`](https://grace-hattersley.github.io/TEVRA/reference/nhood_diversity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
plots <- plot_capture_efficiency(nhood_meta, group_by = "cell_type")
plots$capture_rate
plots$count_vs_size
} # }
```
