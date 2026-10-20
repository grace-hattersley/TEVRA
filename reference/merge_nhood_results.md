# Merge Neighbourhood Analysis Results

Merges any number of neighbourhood analysis dataframes stored in the
metadata slot of a neighbourhood object. Automatically detects which
dataframes exist, validates neighbourhood consistency across all
sources, and handles redundant columns.

## Usage

``` r
merge_nhood_results(
  object,
  dataframe_names = c("nhood_diversity", "nhood_shm"),
  merge_key = "nhood",
  size_col = "nhood_size",
  verbose = TRUE,
  result_name = "merged_nhood_analysis"
)
```

## Arguments

- object:

  A `SingleCellExperiment` or `Milo` object containing neighbourhood
  analysis results in its metadata slot.

- dataframe_names:

  Character vector of metadata slot names to merge. Defaults to
  `c("nhood_diversity", "nhood_shm")`. Must contain at least two names.

- merge_key:

  Column to merge on (default: `"nhood"`).

- size_col:

  Column name containing neighbourhood size information, used to
  validate that all dataframes describe the same neighbourhoods. Can be
  a single string (same column name expected in all dataframes) or a
  named character vector mapping each dataframe name to its size column,
  e.g. `c(nhood_data = "nhood_size", nhood_shm = "n_cells")`. Dataframes
  not present in a named vector fall back to `"nhood_size"`. Set to
  `NULL` to skip the size consistency check entirely.

- verbose:

  Print merge details (default: `TRUE`).

- result_name:

  Name for the merged output in metadata (default:
  `"merged_nhood_analysis"`).

## Value

Updated `object` with merged dataframe stored in
`metadata(object)[[result_name]]`.

## Details

### Dataframe availability:

1.  **Multiple dataframes found** — validates shared neighbourhoods then
    merges all sequentially on `merge_key` using inner joins, removing
    redundant columns at each step.

2.  **Only one dataframe found** — halts with an error (a merge requires
    at least two inputs)

3.  **No dataframes found** — halts with an error

4.  **Result already exists** — overwrites with a warning.

### Neighbourhood validation (two-stage):

1.  *merge_key overlap* — checks that the `merge_key` entries in every
    dataframe overlap with those in the reference (first) dataframe. If
    no overlap exists for any pair the function halts. If only partial
    overlap exists a warning is issued and only shared neighbourhoods
    proceed.

2.  *size consistency* — for neighbourhoods shared across dataframes,
    confirms that the size column contains identical values. Any
    mismatch halts the function as it indicates analyses were run on
    different data.

### Redundant column handling:

At each pairwise merge step, columns shared between the two dataframes
are compared value-by-value. Columns with identical values are dropped
from the right-hand dataframe. Columns with different values are
retained with `_1` / `_2` suffixes.

## See also

[`nhood_diversity`](https://grace-hattersley.github.io/TEVRA/reference/nhood_diversity.md),
[`nhood_cdr3_convergence`](https://grace-hattersley.github.io/TEVRA/reference/nhood_cdr3_convergence.md),
[`nhood_shm`](https://grace-hattersley.github.io/TEVRA/reference/nhood_shm.md),
[`nhood_annotation`](https://grace-hattersley.github.io/TEVRA/reference/nhood_annotation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Default: merge "nhood_diversity" and "nhood_shm"
object <- merge_nhood_results(object)

# Merge three custom dataframes
object <- merge_nhood_results(
  object,
  dataframe_names = c("diversity_results", "nhood_shm", "spatial_features")
)

# Dataframes use different column names for neighbourhood size
object <- merge_nhood_results(
  object,
  dataframe_names = c("diversity_results", "nhood_shm", "spatial_features"),
  size_col = c(diversity_results = "nhood_size",
               nhood_shm         = "n_cells",
               spatial_features  = "size")
)

# Custom merge key and output name
object <- merge_nhood_results(
  object,
  dataframe_names = c("results_a", "results_b"),
  result_name     = "combined_results",
  merge_key       = "NhoodID")
} # }
```
