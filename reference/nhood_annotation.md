# Assign Dominant Annotation to Neighbourhoods

Computes per-neighbourhood the proportion of cells belonging to each
category of one or more annotation columns and assigns a dominant
identity based on a user-defined threshold. Neighbourhoods where no
single category exceeds the threshold are labelled `"mixed"`.

## Usage

``` r
nhood_annotation(
  object,
  nhood_col = NULL,
  annotation_cols = "celltype",
  threshold = 0.7,
  result_name = "nhood_annotation"
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

- threshold:

  Numeric value between 0-1. The minimum proportion a single annotation
  category must reach for a neighbourhood to be assigned that label.
  Neighbourhoods where the dominant category falls below this threshold
  are assigned `"mixed"`. Default: `0.7`.

- result_name:

  Character string. Key under which the results `data.frame` is stored
  in `metadata(object)`. If an entry with this name already exists a
  warning is emitted and the entry is overwritten. Default:
  `"nhood_annotation"`.

## Value

The input `object`, with a `data.frame` stored in
`metadata(object)[[result_name]]`. The `data.frame` has one row per
neighbourhood and the following columns:

- `nhood`:

  Neighbourhood identifier.

- `nhood_size`:

  Total number of cells in the neighbourhood.

- `<annotation_col>`:

  Dominant annotation label for this neighbourhood, or `"mixed"` if no
  single category reached `threshold`. One column per entry in
  `annotation_cols`.

- `<annotation_col>_dominant_prop`:

  Proportion of cells carrying the dominant annotation label. One column
  per entry in `annotation_cols`.

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

### Dominant label assignment

Each neighbourhood represents a group of transcriptionally similar
cells. For each neighbourhood and each annotation column, per-cell
labels are tabulated and converted to proportions. The category with the
highest proportion is the candidate dominant label. If that proportion
is \>= `threshold` the neighbourhood is assigned that label; otherwise
it is labelled `"mixed"`, indicating that no single annotation category
predominates. When multiple annotation columns are supplied the
procedure is applied independently for each.

## See also

[`nhood_annotation_props`](https://grace-hattersley.github.io/TEVRA/reference/nhood_annotation_props.md),
[`nhood_diversity`](https://grace-hattersley.github.io/TEVRA/reference/nhood_diversity.md),
[`nhood_shm`](https://grace-hattersley.github.io/TEVRA/reference/nhood_shm.md),
[`nhood_cdr3_convergence`](https://grace-hattersley.github.io/TEVRA/reference/nhood_cdr3_convergence.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Milo object, single annotation column --------------------------------
milo <- nhood_annotation(
  object          = milo_object,
  annotation_cols = "celltype",
  threshold       = 0.7
)
results <- S4Vectors::metadata(milo)[["nhood_annotation"]]

# --- Multiple annotation columns, relaxed threshold ----------------------
milo <- nhood_annotation(
  object          = milo_object,
  annotation_cols = c("celltype", "lineage"),
  threshold       = 0.6,
  result_name     = "nhood_annotation"
)

# --- Plain SingleCellExperiment object ------------------------------------
sce <- nhood_annotation(
  object          = sce_object,
  nhood_col       = "community_id",
  annotation_cols = "celltype"
)
} # }
```
