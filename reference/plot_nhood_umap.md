# Plot neighbourhood UMAP with optional cell-level and edge overlays

Supports three composited layers:

1.  Cell-level UMAP scatter (optional)

2.  Neighbourhood-to-neighbourhood edges (optional)

3.  Neighbourhood centroids (always shown)

## Usage

``` r
plot_nhood_umap(
  nhood_umap,
  nhood_col = "community_id",
  x_col = "mean_UMAP_1",
  y_col = "mean_UMAP_2",
  fill_by = NULL,
  colour = NULL,
  size_by = NULL,
  size = NULL,
  show_edges = FALSE,
  nhood_connections = NULL,
  edge_col_1 = "mc_1",
  edge_col_2 = "mc_2",
  edge_colour = "grey70",
  edge_linewidth = 0.3,
  edge_alpha = 0.6,
  show_cells = FALSE,
  cell_umap_coords = NULL,
  cell_x_col = "UMAP_1",
  cell_y_col = "UMAP_2",
  cell_colour_by = NULL,
  cell_nhood_col = "nhood_id",
  cell_colour = NULL,
  cell_size = 0.3,
  cell_alpha = 0.5
)
```

## Arguments

- nhood_umap:

  A dataframe of neighbourhood centroids. Must contain neighbourhood
  IDs, UMAP coordinates, and any column used for fill_by or size_by.
  Typically SCE@metadata\$coords merged with
  SCE@metadata\$merged_nhood_analysis.

- nhood_col:

  Name of the neighbourhood ID column in nhood_umap (default
  "community_id").

- x_col:

  Name of the UMAP x-coordinate column in nhood_umap (default
  "mean_UMAP_1").

- y_col:

  Name of the UMAP y-coordinate column in nhood_umap (default
  "mean_UMAP_2").

- fill_by:

  Name of a column in nhood_umap to colour neighbourhood centroids by
  (e.g. "CD8" for assigned cell type). If NULL, uses colour as a single
  fixed colour.

- colour:

  Controls neighbourhood centroid colour and depends on `fill_by`:

  - When `fill_by` is `NULL` — a single hex colour string (default
    `"#3f4788"` when `colour` is left `NULL`).

  - When `fill_by` is a column name — a named character vector mapping
    values in fill_by to hex colour strings. Names must cover all levels
    present in `fill_by`. If left `NULL`, a default categorical palette
    is generated automatically.

- size_by:

  Name of a numeric column in nhood_umap to map centroid point size to,
  or NULL (default) for a fixed point size.

- size:

  Controls centroid point size and depends on `size_by`:

  - When `size_by` is `NULL` — a single numeric giving the fixed point
    size (default `2` when `size` is left `NULL`).

  - When `size_by` is a column name — a numeric vector of length 2
    giving the minimum and maximum point size of the size scale (default
    `c(1, 5)` when `size` is left `NULL`).

- show_edges:

  Logical. Whether to draw edges between connected neighbourhoods
  (default FALSE).

- nhood_connections:

  A dataframe with columns listing connected neighbourhood
  pairs.Self-loops are silently removed. Required when show_edges =
  TRUE.

- edge_col_1:

  Name of the first neighbourhood ID column in nhood_connections
  (default "mc_1").

- edge_col_2:

  Name of the second neighbourhood ID column in nhood_connections
  (default "mc_2").

- edge_colour:

  Colour for edges (default "grey70").

- edge_linewidth:

  Linewidth for edges (default 0.3).

- edge_alpha:

  Alpha for edges (default 0.6).

- show_cells:

  Logical. Whether to draw cell-level scatter as the bottom layer
  (default FALSE).

- cell_umap_coords:

  A dataframe with one row per cell containing UMAP coordinates and any
  cell-level metadata columns. Required when show_cells = TRUE.

- cell_x_col:

  Name of the UMAP x-coordinate column in cell_umap_coords (default
  "UMAP_1").

- cell_y_col:

  Name of the UMAP y-coordinate column in cell_umap_coords (default
  "UMAP_2").

- cell_colour_by:

  Either the name of a column in cell_umap_coords to colour cells by, or
  the string "nhood" to inherit the neighbourhood colour via colour and
  cell_nhood_col. Set NULL (default) for a fixed grey cell colour.

- cell_nhood_col:

  Name of the column in cell_umap_coords that holds each cell's
  neighbourhood assignment. Required when cell_colour_by = "nhood"
  (default "nhood_id").

- cell_colour:

  Controls cell point colour and depends on `cell_colour_by`:

  - When `cell_colour_by` is `NULL` — a single hex colour string
    (default `"grey80"` when `cell_colour` is left `NULL`).

  - When `cell_colour_by` is `"nhood"` — ignored; colour is inherited
    from `colour` via each cell's neighbourhood assignment.

  - When `cell_colour_by` is any other column name — a named character
    vector mapping values of cell_colour_by to hex colour strings. If
    left `NULL`, a default categorical palette is generated
    automatically.

- cell_size:

  Fixed point size for cells (default 0.3).

- cell_alpha:

  Alpha for cell points (default 0.5).

## Value

A `ggplot` object. Layers are composited in the order cells (optional)
-\> edges (optional) -\> neighbourhood centroids, so the result can be
extended with further ggplot2 layers, scales or themes.

## Details

Layer order (bottom to top): cells → edges → neighbourhood centroids.

## See also

[`create_communities`](https://grace-hattersley.github.io/TEVRA/reference/create_communities.md)
and
[`test_resolution`](https://grace-hattersley.github.io/TEVRA/reference/test_resolution.md)
for generating the neighbourhood assignments plotted here;
[`nhood_umap`](https://grace-hattersley.github.io/TEVRA/reference/nhood_umap.md)
and
[`test_resolution`](https://grace-hattersley.github.io/TEVRA/reference/test_resolution.md)
for generating the neighbourhood umap coordinates;
[`nhood_annotation`](https://grace-hattersley.github.io/TEVRA/reference/nhood_annotation.md)
for annotating neighbourhoods

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Neighbourhood centroids coloured by assigned cell type ---------------
plot_nhood_umap(
  nhood_umap = nhood_df,
  fill_by    = "CD8",
  colour     = celltype_colours
)

# --- Neighbourhood centroids with default palette -------------------------
plot_nhood_umap(
  nhood_umap = nhood_df,
  fill_by    = "CD8"
)

# --- Neighbourhood centroids scaled by neighbourhood size -----------------
plot_nhood_umap(
  nhood_umap = nhood_df,
  size_by    = "n_cells",
  size       = c(1, 6)
)

# --- Neighbourhood centroids, edges and cells -----------------------------
plot_nhood_umap(
  nhood_umap        = nhood_df,
  fill_by           = "CD8",
  colour            = celltype_colours,
  show_edges        = TRUE,
  nhood_connections = connections_df,
  show_cells        = TRUE,
  cell_umap_coords  = cell_df,
  cell_colour_by    = "nhood",
  cell_nhood_col    = "nhood_id"
)

# --- Neighbourhood centroids, edges and cells coloured independently ------
plot_nhood_umap(
  nhood_umap       = nhood_df,
  fill_by          = "CD8",
  colour           = celltype_colours,
  show_cells       = TRUE,
  cell_umap_coords = cell_df,
  cell_colour_by   = "timepoint",
  cell_colour      = c(Pre = "#8da0cb", On = "#fc8d62")
)
} # }
```
