# Calculate Cell Annotation Proportions per Neighbourhood

Computes the proportion of cells assigned to each neighbourhood that
carry each category of one or more annotation columns.

## Usage

``` r
nhood_annotation_props(
  object,
  nhood_col = NULL,
  annotation_cols = "celltype",
  result_name = "nhood_annotation_props"
)
```

## Arguments

- object:

  A [`Milo`](https://marionilab.github.io/miloR/reference/Milo.html) or
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object. Seurat objects must first be converted via
  [`as.SingleCellExperiment`](https://satijalab.org/seurat/reference/as.SingleCellExperiment.html).

- nhood_col:

  Character string. For plain `SingleCellExperiment` objects (i.e. not
  `Milo`), the name of the `colData` column containing neighbourhood
  assignments. Not required for `Milo` objects.

- annotation_cols:

  Character vector of `colData` column names containing per-cell
  annotation categories to summarise per neighbourhood (e.g. cell type,
  lineage, cluster label). All specified columns must be present in
  `colData`; missing columns raise an error. Default: `"celltype"`.

- result_name:

  Character string. Key under which the results `data.frame` is stored
  in `metadata(object)`. If an entry with this name already exists a
  warning is emitted and the entry is overwritten. Default:
  `"nhood_annotation_props"`.

## Value

The input `object`, with a `data.frame` stored in
`metadata(object)[[result_name]]`. The `data.frame` has one row per
neighbourhood and the following columns:

- `nhood`:

  Neighbourhood identifier.

- `nhood_size`:

  Total number of cells in the neighbourhood.

- `<annotation_col>_<category>_prop`:

  Proportion of cells carrying each annotation category. One column per
  category per entry in `annotation_cols` — e.g. if
  `annotation_cols = "celltype"` and the column contains `"CD4"` and
  `"CD8"`, the output columns will be `celltype_CD4_prop` and
  `celltype_CD8_prop`. Neighbourhoods where all cells have `NA` for a
  given annotation column have `NA` in all corresponding proportion
  columns.

## Details

### Object type handling

- **Milo objects**: Neighbourhood membership is read directly from the
  sparse `nhoods` matrix via
  [`nhoods`](https://marionilab.github.io/miloR/reference/methods.html).
  Row names of this matrix must match `colnames(object)`; a mismatch
  raises an error to prevent silent misassignment.

- **Plain SingleCellExperiment objects**: A binary neighbourhood
  membership matrix is constructed internally from the column specified
  by `nhood_col`. Cells with `NA` assignments are excluded with a
  warning.

### Proportion calculation

For each neighbourhood and each annotation column, per-cell labels are
tabulated across all valid (non-`NA`) cells. Proportions are computed
relative to the number of valid cells in the neighbourhood for that
annotation column (i.e. each row sums to 1 over the proportion columns
for a given `annotation_col`, excluding any `NA` cells). Cells with `NA`
in a given annotation column are excluded from that column's proportions
and their count is reported via a warning.

## See also

[`nhood_annotation`](https://grace-hattersley.github.io/TEVRA/reference/nhood_annotation.md),
[`nhood_diversity`](https://grace-hattersley.github.io/TEVRA/reference/nhood_diversity.md),
[`nhood_shm`](https://grace-hattersley.github.io/TEVRA/reference/nhood_shm.md),
[`nhood_cdr3_convergence`](https://grace-hattersley.github.io/TEVRA/reference/nhood_cdr3_convergence.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Milo object, single annotation column --------------------------------
milo <- nhood_annotation_props(
  object          = milo_object,
  annotation_cols = "celltype"
)
results <- S4Vectors::metadata(milo)[["nhood_annotation_props"]]
# Output columns: nhood, nhood_size, celltype_CD4_prop, celltype_CD8_prop, ...

# --- Multiple annotation columns ------------------------------------------
milo <- nhood_annotation_props(
  object          = milo_object,
  annotation_cols = c("celltype", "lineage"),
  result_name     = "nhood_annotation_props"
)

# --- Plain SingleCellExperiment object ------------------------------------
sce <- nhood_annotation_props(
  object          = sce_object,
  nhood_col       = "community_id",
  annotation_cols = "celltype"
)
} # }
```
