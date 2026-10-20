# Test Leiden clustering resolutions on a nearest-neighbour graph

Runs Leiden community detection across a range of resolution values and
returns a summary table of clustering metrics for each resolution. Use
this function to select an appropriate resolution before running
[`create_communities`](https://grace-hattersley.github.io/TEVRA/reference/create_communities.md).

## Usage

``` r
test_resolution(
  sce,
  graph.name = "RNA_snn",
  resolutions = c(0.5, 1, 2.5, 5, 10, 20),
  min_cells = 5
)
```

## Arguments

- sce:

  A `SingleCellExperiment` object. The object must contain a
  pre-computed nearest-neighbour graph stored in `metadata(sce)$graphs`
  under the name given by `graph.name`. See the **Graph preparation**
  section below for how to build and store the graph.

- graph.name:

  A character string giving the name of the graph to retrieve from
  `metadata(sce)$graphs`. Defaults to `"RNA_snn"`.

- resolutions:

  A numeric vector of resolution values to test. Higher values produce
  more, smaller communities. Defaults to `c(0.5, 1, 2.5, 5, 10, 20)`.

- min_cells:

  A positive integer that defines the minimum community size.
  Communities with fewer than this many cells are considered too small
  and each of the consistuent cells are reassigned to the valid
  community occupied by the majority of the cell's nearest neighbours/
  Defaults to `5`.

## Value

A `data.frame` with one row per tested resolution and the following
columns:

- `resolution`:

  The resolution value tested.

- `n_communities`:

  Number of communities after reassignment of small clusters.

- `mean_cells_per_community`:

  Mean community size.

- `sd_cells_per_community`:

  Standard deviation of community sizes.

- `median_cells_per_community`:

  Median community size.

- `min_cells_per_community`:

  Smallest community size.

- `max_cells_per_community`:

  Largest community size.

- `modularity`:

  Graph modularity of the partition.

Resolutions at which all communities fall below `min_cells` are skipped
with a warning and will not appear in the output.

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

[`create_communities`](https://grace-hattersley.github.io/TEVRA/reference/create_communities.md)
to assign cells to communities once a resolution has been chosen

## Examples

``` r
if (FALSE) { # \dontrun{
library(scran)

# Build and store the graph
g <- scran::buildSNNGraph(sce, use.dimred = "PCA", k = 20)
metadata(sce)$graphs[["RNA_snn"]] <- g

# Test a range of resolutions
res_results <- test_resolution(sce, resolutions = c(0.5, 1, 2.5, 5))
print(res_results)
} # }
```
