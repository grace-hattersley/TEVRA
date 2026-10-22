# Calculate Neighbourhood Somatic Hypermutation Statistics

Computes per-neighbourhood somatic hypermutation (SHM) summary
statistics from B cell receptor (BCR) data stored in `colData`. For each
neighbourhood, the function aggregates SHM metrics across cells that
carry BCR information, and flags neighbourhoods with insufficient BCR
coverage as unreliable.

## Usage

``` r
nhood_shm(
  object,
  nhood_col = NULL,
  shm_cols = c("mu_count_VDJ", "mu_count_VJ", "mu_freq_VDJ", "mu_freq_VJ", "mu_count",
    "mu_freq"),
  min_n_shm = 20,
  summary_metric = "mean",
  result_name = "nhood_shm"
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

- shm_cols:

  Character vector of column names in `colData(object)` containing
  per-cell SHM metrics to aggregate. Defaults to
  `c("mu_count_VDJ", "mu_count_VJ", "mu_freq_VDJ", "mu_freq_VJ", "mu_count", "mu_freq")`.
  All specified columns must be present in `colData`; missing columns
  will raise an error.

- min_n_shm:

  Integer. Minimum number of cells with SHM annotation required for a
  neighbourhood's SHM statistics to be considered reliable.
  Neighbourhoods with fewer than `min_n_shm` contigs are retained in the
  output but have `shm_reliable = FALSE` and `NA` values for all SHM
  columns. Default is `10`.

- summary_metric:

  Character string specifying the central tendency statistic to compute
  for each SHM metric within each neighbourhood. One of `"mean"`
  (default) or `"median"`. The corresponding spread metric is computed
  automatically: standard deviation (`sd`) for `"mean"`, interquartile
  range (`IQR`) for `"median"`.

- result_name:

  Character string. Name of the metadata slot in which to store the
  resulting SHM data frame, i.e. the results are accessible via
  `metadata(object)[[result_name]]`. If the slot already exists it will
  be overwritten with a warning. Default is `"nhood_shm"`.

## Value

The input `object` with an additional (or updated) entry in
`metadata(object)[[result_name]]`. The stored value is a `data.frame`
with one row per neighbourhood and the following columns:

- `nhood`:

  Neighbourhood identifier.

- `nhood_size`:

  Total number of cells in the neighbourhood.

- `n_shm`:

  Number of cells with at least one non-`NA` SHM value.

- `shm_capture_rate`:

  Proportion of neighbourhood cells with SHM data
  (`n_shm / nhood_size`).

- `shm_reliable`:

  Logical; `TRUE` if `n_shm >= min_n_shm`.

- `<shm_col>_<summary_metric>`:

  One column per entry in `shm_cols` containing the aggregated statistic
  (mean or median).

- `<shm_col>_<spread_metric>`:

  One column per entry in `shm_cols` containing the spread statistic
  (`sd` or `IQR`).

The data frame also carries the following attributes: `min_n_shm`,
`shm_cols`, `summary_metric`, and `spread_metric`.

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

### SHM reliability flag

The `shm_reliable` flag is intended to help downstream filtering.
Neighbourhoods below `min_n_shm` will have `NA` for all SHM statistics
and should generally be excluded from differential SHM analyses.

## See also

[`nhood_diversity`](https://grace-hattersley.github.io/TEVRA/reference/nhood_diversity.md),
[`nhood_annotation`](https://grace-hattersley.github.io/TEVRA/reference/nhood_annotation.md),
[`nhood_annotation_props`](https://grace-hattersley.github.io/TEVRA/reference/nhood_annotation_props.md),
[`nhood_cdr3_convergence`](https://grace-hattersley.github.io/TEVRA/reference/nhood_cdr3_convergence.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Milo object ---
milo_obj <- nhood_shm(
  object = milo_obj,
  shm_cols = c("mu_freq_VDJ", "mu_freq_VJ", "mu_freq"),
  min_n_shm = 15,
  summary_metric = "median",
  result_name = "nhood_shm"
)

# Retrieve results
shm_df <- metadata(milo_obj)[["nhood_shm"]]
head(shm_df)

# --- Plain SingleCellExperiment object ---
sce_obj <- nhood_shm(
  object = sce_obj,
  nhood_col = "neighbourhood_id",
  min_n_shm = 10,
  summary_metric = "mean"
)
} # }
```
