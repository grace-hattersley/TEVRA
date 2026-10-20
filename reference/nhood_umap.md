# Compute Neighbourhood UMAP Centroids

Aggregates per-cell UMAP coordinates into per-neighbourhood centroids

## Usage

``` r
nhood_umap(
  object,
  nhood_col = "community_id",
  aggregate_fun = "mean",
  result_name = "nhood_umap"
)
```

## Arguments

- object:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object with a UMAP stored in `reducedDims`. Seurat objects must be
  converted via
  [`as.SingleCellExperiment`](https://satijalab.org/seurat/reference/as.SingleCellExperiment.html)
  before use.

- nhood_col:

  Character string. Name of the `colData` column containing
  neighbourhood assignments. Default: `"community_id"`.

- aggregate_fun:

  Character string. Aggregation function used to compute neighbourhood
  centroids. Must be one of `"mean"` (arithmetic mean; default) or
  `"median"`. The choice affects both the computed coordinates and the
  names of the output columns (see Value).

- result_name:

  Character string. Key under which the results `data.frame` is stored
  in `metadata(object)`. If an entry with this name already exists a
  warning is emitted and the entry is overwritten. Default:
  `"nhood_umap"`.

## Value

The input `object`, with a `data.frame` of neighbourhood UMAP centroids
stored in `metadata(object)[[result_name]]`. The `data.frame` contains
the following columns:

- `nhood`:

  Neighbourhood identifier.

- `nhood_size`:

  Total number of cells in the neighbourhood.

- `mean_UMAP_1` / `median_UMAP_1`:

  UMAP dimension 1 centroid, named according to `aggregate_fun`.

- `mean_UMAP_2` / `median_UMAP_2`:

  UMAP dimension 2 centroid, named according to `aggregate_fun`.

## Details

The UMAP is accessed by name from `reducedDims(object)$UMAP` and is
assumed to have exactly two dimensions, which are renamed `UMAP_1` and
`UMAP_2` internally. If your UMAP is stored under a different name in
`reducedDims`, rename it before calling this function. Centroids are
computed using the function specified by `aggregate_fun` applied to all
cell coordinates within each neighbourhood. `"mean"` produces the
arithmetic centroid; `"median"` produces the geometric median and is
more robust to outlier cells at neighbourhood boundaries. If
`result_name` already exists in `metadata(object)` a warning is emitted
and the entry is overwritten.

## See also

[`plot_nhood_umap`](https://grace-hattersley.github.io/TEVRA/reference/plot_nhood_umap.md)
to visualise the neighbourhood umap

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Mean centroids (default) ---------------------------------------------
sce <- nhood_umap(
  object    = sce_object,
  nhood_col = "community_id"
)
umap_df <- S4Vectors::metadata(sce)[["nhood_umap"]]

# --- Median centroids -----------------------------------------------------
sce <- nhood_umap(
  object       = sce_object,
  nhood_col    = "community_id",
  aggregate_fun = "median"
)
umap_df <- S4Vectors::metadata(sce)[["nhood_umap"]]
} # }
```
