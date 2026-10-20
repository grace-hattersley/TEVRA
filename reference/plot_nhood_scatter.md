# Scatter plot of neighbourhood-level metrics

Plots two neighbourhood-level metrics against each other, with point
size optionally scaled by a neighbourhood-level column and colour mapped
to a grouping variable.

## Usage

``` r
plot_nhood_scatter(
  nhood_data,
  x_col = NULL,
  y_col = NULL,
  x_reliability_col = NULL,
  y_reliability_col = NULL,
  colour_by = NULL,
  size_by = NULL,
  size = NULL,
  alpha = 0.8,
  facet_by = NULL,
  facet_ncol = NULL,
  facet_scales = "free",
  add_density = FALSE,
  x_lab = NULL,
  y_lab = NULL
)
```

## Arguments

- nhood_data:

  Data frame containing all columns to plot.

- x_col:

  Column name (string) for the x-axis metric (default
  "shannon_diversity").

- y_col:

  Column name (string) for the y-axis metric (default "mean_shm").

- x_reliability_col:

  Column name of a logical reliability flag for the x-axis metric.
  Neighbourhoods where this column is FALSE or NA are dropped before
  plotting, regardless of their y-axis reliability. NULL skips this
  filter.

- y_reliability_col:

  Column name of a logical reliability flag for the y-axis metric.
  Neighbourhoods where this column is FALSE or NA are dropped before
  plotting, regardless of their x-axis reliability. NULL skips this
  filter.

- colour_by:

  Column name (string) for point colour (default "anat").

- size_by:

  Column name (string) to map point size to, or NULL (default) for a
  fixed point size.

- size:

  Controls point size and depends on `size_by`:

  - When `size_by` is `NULL` — a single numeric giving the fixed point
    size (default `2` when `size` is left `NULL`).

  - When `size_by` is a column name — a numeric vector of length 2
    giving the minimum and maximum point size of the size scale (default
    `c(1, 4)` when `size` is left `NULL`).

- alpha:

  Numeric between 0-1. Point transparency (default 0.8).

- facet_by:

  Column name to facet by, or NULL for no faceting.

- facet_ncol:

  Number of columns for facet_wrap (NULL = auto).

- facet_scales:

  Axis scale argument passed to facet_wrap: "free" (default), "fixed",
  "free_x", or "free_y".

- add_density:

  Logical. If TRUE, adds 2D density contours (default FALSE).

- x_lab:

  Optional custom x-axis label. NULL auto-derives from x_col.

- y_lab:

  Optional custom y-axis label. NULL auto-derives from y_col.

## Value

A ggplot object.

## See also

[`nhood_diversity`](https://grace-hattersley.github.io/TEVRA/reference/nhood_diversity.md),
[`nhood_shm`](https://grace-hattersley.github.io/TEVRA/reference/nhood_shm.md),
[`nhood_cdr3_convergence`](https://grace-hattersley.github.io/TEVRA/reference/nhood_cdr3_convergence.md),
[`plot_nhood_beeswarm`](https://grace-hattersley.github.io/TEVRA/reference/plot_nhood_beeswarm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fixed-size points
plot_nhood_scatter(
  nhood_data = nhood_data,
  x_col      = "shannon_diversity",
  y_col      = "mean_shm",
  colour_by  = "anat"
)

# Point size mapped to neighbourhood size, with a custom size range
plot_nhood_scatter(
  nhood_data = nhood_data,
  x_col      = "shannon_diversity",
  y_col      = "mean_shm",
  colour_by  = "anat",
  size_by    = "nhood_size",
  size       = c(2, 8)
)
} # }
```
