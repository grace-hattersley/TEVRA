# Plot mean community size across tested resolutions

Visualises the output of
[`test_resolution`](https://grace-hattersley.github.io/TEVRA/reference/test_resolution.md)
as a line plot of mean cells per community against resolution, with
optional error bars. Use this plot to identify a resolution that yields
communities of a suitable and stable size before running
[`create_communities`](https://grace-hattersley.github.io/TEVRA/reference/create_communities.md).

## Usage

``` r
plot_resolution_results(resolution_results, show_error_bars = TRUE)
```

## Arguments

- resolution_results:

  A `data.frame` as returned by
  [`test_resolution`](https://grace-hattersley.github.io/TEVRA/reference/test_resolution.md),
  containing at minimum the columns `resolution` and
  `mean_cells_per_community`. If `show_error_bars = TRUE`, the column
  `sd_cells_per_community` must also be present.

- show_error_bars:

  A logical value. If `TRUE` (the default), error bars showing ±1 SD
  (`sd_cells_per_community`) are added around each point.

## Value

A `ggplot` object showing mean cells per community against resolution.

## See also

[`test_resolution`](https://grace-hattersley.github.io/TEVRA/reference/test_resolution.md)
to generate the input `data.frame`;
[`create_communities`](https://grace-hattersley.github.io/TEVRA/reference/create_communities.md)
to assign cells to communities once a resolution has been chosen.

## Examples

``` r
if (FALSE) { # \dontrun{
# Test a range of resolutions
res_results <- test_resolution(sce, resolutions = c(0.5, 1, 2.5, 5))

# Plot the results with error bars
plot_resolution_results(res_results)

# Plot without error bars
plot_resolution_results(res_results, show_error_bars = FALSE)
} # }
```
