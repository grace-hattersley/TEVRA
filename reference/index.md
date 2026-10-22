# Package index

## Neighbourhood construction

Define transcriptional neighbourhoods and choose a clustering
resolution.

- [`test_resolution()`](https://grace-hattersley.github.io/TEVRA/reference/test_resolution.md)
  : Test Leiden clustering resolutions on a nearest-neighbour graph
- [`plot_resolution_results()`](https://grace-hattersley.github.io/TEVRA/reference/plot_resolution_results.md)
  : Plot mean community size across tested resolutions
- [`create_communities()`](https://grace-hattersley.github.io/TEVRA/reference/create_communities.md)
  : Create cell communities using Leiden community detection
- [`plot_nhood_sizes()`](https://grace-hattersley.github.io/TEVRA/reference/plot_nhood_sizes.md)
  : Plot the distribution of neighbourhood sizes
- [`nhood_umap()`](https://grace-hattersley.github.io/TEVRA/reference/nhood_umap.md)
  : Compute Neighbourhood UMAP Centroids

## Neighbourhood annotation

Assign cell type labels and composition to each neighbourhood.

- [`nhood_annotation()`](https://grace-hattersley.github.io/TEVRA/reference/nhood_annotation.md)
  : Assign Dominant Annotation to Neighbourhoods
- [`nhood_annotation_props()`](https://grace-hattersley.github.io/TEVRA/reference/nhood_annotation_props.md)
  : Calculate Cell Annotation Proportions per Neighbourhood

## Repertoire metrics

Compute receptor repertoire features for each neighbourhood.

- [`nhood_diversity()`](https://grace-hattersley.github.io/TEVRA/reference/nhood_diversity.md)
  : Calculate Neighbourhood Clonal Diversity Indices
- [`nhood_shm()`](https://grace-hattersley.github.io/TEVRA/reference/nhood_shm.md)
  : Calculate Neighbourhood Somatic Hypermutation Statistics
- [`nhood_cdr3_convergence()`](https://grace-hattersley.github.io/TEVRA/reference/nhood_cdr3_convergence.md)
  : Compute inter-clonal CDR3 amino acid convergence scores per
  neighbourhood
- [`nhood_cdr3_convergence_matrices()`](https://grace-hattersley.github.io/TEVRA/reference/nhood_cdr3_convergence_matrices.md)
  : Extract pairwise CDR3 convergence matrices for selected
  neighbourhoods
- [`merge_nhood_results()`](https://grace-hattersley.github.io/TEVRA/reference/merge_nhood_results.md)
  : Merge Neighbourhood Analysis Results

## Visualisation

Plot neighbourhood-level results.

- [`plot_nhood_umap()`](https://grace-hattersley.github.io/TEVRA/reference/plot_nhood_umap.md)
  : Plot neighbourhood UMAP with optional cell-level and edge overlays
- [`plot_nhood_beeswarm()`](https://grace-hattersley.github.io/TEVRA/reference/plot_nhood_beeswarm.md)
  : Beeswarm plot for neighbourhood-level metrics
- [`plot_nhood_scatter()`](https://grace-hattersley.github.io/TEVRA/reference/plot_nhood_scatter.md)
  : Scatter plot of neighbourhood-level metrics
- [`plot_capture_efficiency()`](https://grace-hattersley.github.io/TEVRA/reference/plot_capture_efficiency.md)
  : Plot VDJ contig capture efficiency across neighbourhoods

## Example data

- [`tevra_example`](https://grace-hattersley.github.io/TEVRA/reference/tevra_example.md)
  : Example dataset: murine splenic B cells
