# Create cell communities using Leiden community detection

Assigns cells in a `SingleCellExperiment` object to communities using
the Leiden algorithm at a single chosen resolution. Community
assignments are added to `colData(sce)` and a per-community summary
table is returned alongside the updated object. Use
[`test_resolution`](https://grace-hattersley.github.io/TEVRA/reference/test_resolution.md)
first to identify a suitable resolution.

## Usage

``` r
create_communities(
  sce,
  graph.name = "RNA_snn",
  resolution = 0.5,
  min_cells = 5,
  feature.col = "nFeature_RNA",
  count.col = "nCount_RNA"
)
```

## Arguments

- sce:

  A `SingleCellExperiment` object. The object must contain a
  pre-computed nearest-neighbour graph stored in `metadata(sce)$graphs`
  under the name given by `graph.name`. See the **Graph preparation**
  section below for how to build and store the graph

- graph.name:

  A character string giving the name of the graph to retrieve from
  `metadata(sce)$graphs`. Defaults to `"RNA_snn"`.

- resolution:

  A single numeric value controlling the granularity of Leiden
  clustering. Higher values produce more, smaller communities. Use
  [`test_resolution`](https://grace-hattersley.github.io/TEVRA/reference/test_resolution.md)
  to explore candidate values before committing to one. Defaults to
  `0.5`.

- min_cells:

  A positive integer. Communities with fewer than this many cells are
  considered too small; their cells are reassigned to the nearest valid
  community. If all communities fall below this threshold the function
  stops with an informative error. Defaults to `5`.

- feature.col:

  A character string giving the name of the `colData(sce)` column
  holding the number of detected features (genes) per cell. Defaults to
  `"nFeature_RNA"`, the Seurat convention. Set this if your object uses
  a different column name (e.g. `"detected"`, as produced by
  [`scuttle::perCellQCMetrics`](https://rdrr.io/pkg/scuttle/man/perCellQCMetrics.html)).
  If the column is not found, `mean_detected` is returned as `NA` with a
  warning.

- count.col:

  A character string giving the name of the `colData(sce)` column
  holding total counts per cell. Defaults to `"nCount_RNA"`, the Seurat
  convention. Set this if your object uses a different column name (e.g.
  `"sum"`, as produced by
  [`scuttle::perCellQCMetrics`](https://rdrr.io/pkg/scuttle/man/perCellQCMetrics.html)).
  If the column is not found, `mean_total_counts` is returned as `NA`
  with a warning.

## Value

A named list with two elements:

- `sce`:

  The input `SingleCellExperiment` object with a new `community_id`
  column added to `colData`.

- `communities_summary`:

  A `data.frame` with one row per community containing:

  `community_id`

  :   Community identifier.

  `n_cells`

  :   Number of cells in the community.

  `mean_detected`

  :   Mean number of detected features per cell, taken from the
      `feature.col` column of `colData`.

  `mean_total_counts`

  :   Mean total counts per cell, taken from the `count.col` column of
      `colData`.

## Details

### Graph preparation:

This function does not build the nearest-neighbour graph itself. The
graph must be computed and stored in `metadata(sce)$graphs` before
calling this function. Two workflows are supported:

**1. From scratch with scran:**


    library(scran)

    # Build a shared nearest-neighbour graph on the PCA embedding
    g <- scran::buildSNNGraph(sce, use.dimred = "PCA", k = 20)

    # Store it under the name this function will look for
    metadata(sce)$graphs[["RNA_snn"]] <- g

**2. From an existing Seurat object:**


    library(Seurat)

    # Convert to SCE first
    sce <- as.SingleCellExperiment(seurat_obj)

    # Extract the adjacency matrix from the Seurat object, convert to igraph,
    # and store in the SCE metadata
    g <- igraph::graph_from_adjacency_matrix(
      seurat_obj@graphs[["RNA_snn"]],
      mode     = "max",
      weighted = TRUE
    )
    metadata(sce)$graphs[["RNA_snn"]] <- g

The value passed to `graph.name` must match the key used when storing
the graph in `metadata(sce)$graphs`.

## See also

[`test_resolution`](https://grace-hattersley.github.io/TEVRA/reference/test_resolution.md)
to evaluate resolution choices before running this function

## Examples

``` r
if (FALSE) { # \dontrun{
library(scran)

# Build and store the graph
g <- scran::buildSNNGraph(sce, use.dimred = "PCA", k = 20)
metadata(sce)$graphs[["RNA_snn"]] <- g

# Optionally explore resolutions first
test_resolution(sce, resolutions = c(0.5, 1, 2.5, 5))

# Create communities at the chosen resolution
out <- create_communities(sce, resolution = 1)
sce <- out$sce
summary_df <- out$communities_summary

# If sce was built with scuttle::perCellQCMetrics() rather than Seurat,
# specify the matching column names
out <- create_communities(sce, resolution = 1,
                           feature.col = "detected",
                           count.col = "sum")
} # }
```
