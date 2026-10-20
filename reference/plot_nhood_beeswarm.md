# Beeswarm plot for neighbourhood-level metrics

Accepts a data frame with neighbourhoods as rows.

## Usage

``` r
plot_nhood_beeswarm(
  nhood_data,
  plot_cols,
  reliability_col = NULL,
  group_by,
  summary_metric = c("global_mean", "none", "global_median", "group_mean",
    "group_median"),
  colour_by = NULL,
  p_val_col = "p_val",
  LFC_col = "LFC",
  sig_threshold = 0.1,
  size_by = NULL,
  size = NULL,
  alpha = 0.8,
  facet_by = NULL,
  facet_ncol = NULL,
  facet_scales = c("free", "fixed", "free_x", "free_y"),
  x_limits = NULL,
  x_lab = NULL,
  y_lab = NULL,
  title = "",
  subtitle = NULL
)
```

## Arguments

- nhood_data:

  Data frame containing all columns to plot.

- plot_cols:

  Character vector of column name(s) to use as the x-axis. Supply one
  name for a single-metric plot, or multiple names to get one panel row
  per metric (multi-column mode).

- reliability_col:

  Column name of a logical reliability flag. Rows where this column is
  FALSE or NA are dropped before plotting. NULL skips reliability
  filtering entirely.

- group_by:

  Column name (string) whose levels form the y-axis categories.

- summary_metric:

  Controls whether a summary line is drawn, which statistic is used, and
  whether it is computed globally or per group. One of:

  - `"none"` – no summary line is drawn.

  - `"global_mean"` – (default) one dashed vertical line per panel at
    the grand mean across all `group_by` levels.

  - `"global_median"` – as above but using the median.

  - `"group_mean"` – one short vertical segment per `group_by` level at
    the per-group mean. Segments span ±0.4 y-units so they sit within
    each beeswarm row.

  - `"group_median"` – as above but using the median.

- colour_by:

  Column name (string) for point colour. NULL defaults to group_by. Set
  to the value of LFC_col to activate p-value-filtered LFC colouring.

- p_val_col:

  Column name of p-values (used in LFC mode).

- LFC_col:

  Column name of log-fold-changes (used in LFC mode).

- sig_threshold:

  Significance threshold; points with p \< sig_threshold are coloured by
  LFC in LFC mode.

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

  Column name to facet by, or NULL. In multi-column mode this produces a
  facet_grid (metric rows x facet_by cols); in single-column mode it
  produces a facet_wrap.

- facet_ncol:

  Number of columns for facet_wrap (NULL = auto). Ignored when
  facet_grid is used.

- facet_scales:

  Axis scale argument passed to facet_wrap / facet_grid: "free"
  (default), "fixed", "free_x", or "free_y".

- x_limits:

  Numeric vector of length 2 for x-axis limits, or NULL for automatic
  limits. Example: c(0, 1) for proportions.

- x_lab:

  X-axis label string. NULL auto-derives from plot_cols: a single column
  name is title-cased; multiple columns produce "Value".

- y_lab:

  Y-axis label string. NULL auto-derives from group_by.

- title:

  Plot title string.

- subtitle:

  Plot subtitle string, or NULL for none.

## Value

A ggplot object.

## Details

Multi-column mode (plot_cols): when more than one column is supplied,
the data are pivoted to long format and each column becomes a row of
panels.

## See also

[`nhood_diversity`](https://grace-hattersley.github.io/TEVRA/reference/nhood_diversity.md),
[`nhood_shm`](https://grace-hattersley.github.io/TEVRA/reference/nhood_shm.md),
[`nhood_cdr3_convergence`](https://grace-hattersley.github.io/TEVRA/reference/nhood_cdr3_convergence.md),
[`plot_nhood_scatter`](https://grace-hattersley.github.io/TEVRA/reference/plot_nhood_scatter.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Single-metric plot
plot_nhood_beeswarm(
  nhood_data      = nhood_data,
  plot_cols       = "shannon_diversity",
  group_by        = "anat",
  reliability_col = "diversity_reliable",
  subtitle        = paste0("Reliable estimates (>=", attr(nhood_data, "min_bcr_tcr"), " BCR/TCR)")
)

# Multi-column plot
plot_nhood_beeswarm(
  nhood_data = ccall_df,
  plot_cols  = c("IgA", "IgG", "IgM"),
  group_by   = "anat",
  x_limits   = c(0, 1),
  x_lab      = "Proportion of Cells"
)

# Point size mapped to a column, with a custom size range
plot_nhood_beeswarm(
  nhood_data = nhood_data,
  plot_cols  = "shannon_diversity",
  group_by   = "anat",
  size_by    = "nhood_size",
  size       = c(1, 6)
)
} # }
```
