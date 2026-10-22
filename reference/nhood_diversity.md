# Calculate Neighbourhood Clonal Diversity Indices

Computes diversity metrics for each neighbourhood in a `Milo` or
`SingleCellExperiment` object, using clonal identifiers (clone IDs or
CDR3 sequences) extracted from `colData`.

## Usage

``` r
nhood_diversity(
  object,
  nhood_col = NULL,
  clone_id_col = NULL,
  cdr3_col = "cdr3_aa_VDJ",
  min_n_contigs = 20,
  metrics = c("shannon", "inv.simpson", "gini.simpson", "chao1"),
  downsample_method = "rarefaction",
  n_boots = 100,
  feature_scale = FALSE,
  result_name = "nhood_diversity"
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

- clone_id_col:

  Character string. Name of the `colData` column containing clone
  identifiers. If `NULL` or not found, `cdr3_col` is used as a fallback.

- cdr3_col:

  Character string. Name of the `colData` column containing CDR3
  amino-acid sequences. Used when `clone_id_col` is unavailable.
  Default: `"cdr3_aa_VDJ"`.

- min_n_contigs:

  Integer. Minimum number of cells with contig availability required for
  a neighbourhood to be considered reliable. Neighbourhoods below this
  threshold are flagged with `diversity_reliable = FALSE` and their
  metric values are set to `NA`. Default: `20`.

- metrics:

  Character vector of diversity metrics to compute. Any combination of
  `"shannon"`, `"norm.entropy"`, `"inv.simpson"`, `"gini.simpson"`,
  `"chao1"`, `"ace"`, `"gini"` and `"standardised_n_clones"`.
  `"standardised_n_clones"` computes coverage-standardised observed
  clone richness via
  [`iNEXT::estimateD`](https://rdrr.io/pkg/iNEXT/man/estimateD.html)
  with Hill number order \\q = 0\\ and can only be requested if
  `downsample_method = "coverage"`. Default:
  `c("shannon", "inv.simpson", "gini.simpson", "chao1")`.

- downsample_method:

  Character string specifying how to standardise sequencing depth across
  neighbourhoods before computing metrics. One of:

  `"rarefaction"`

  :   Random subsampling to the median contig count across valid
      neighbourhoods (base R implementation). The rarefaction depth is
      reported at runtime. Neighbourhoods whose contig count is at or
      below the rarefaction depth are not downsampled; this can occur
      frequently when neighbourhood sizes vary widely. Rarefaction is
      preferred when communities are highly diverse or singleton-
      dominated, as it standardises sample size rather than coverage.
      `"chao1"` and `"ace"` are not rarefied; they are computed on
      observed counts and `<metric>_diversity_sd` is `NA` for these
      metrics.

  `"vegan_rarefaction"`

  :   As above, using
      [`vegan::rrarefy`](https://vegandevs.github.io/vegan/reference/rarefy.html).
      Requires the vegan package. The same guidance on rarefaction depth
      interpretation and `"chao1"`/ `"ace"` compatibility applies.

  `"coverage"`

  :   Coverage-based standardisation via
      [`iNEXT::estimateD`](https://rdrr.io/pkg/iNEXT/man/estimateD.html).
      Requires the iNEXT package. The target coverage is set to the
      median Good's coverage (\\1 - singletons / n\\) across valid
      neighbourhoods and is reported at runtime. Coverage-based
      standardisation is preferred when neighbourhood sizes are large
      and relatively uniform. If the reported target coverage is very
      low (near 0), rarefaction may be more appropriate. `"chao1"`,
      `"ace"`, `"norm.entropy"`, and `"gini"` are not coverage
      standardised; they are computed on observed counts and
      `<metric>_diversity_se` is `NA` for these metrics.

  `"none"`

  :   No downsampling; metrics are computed on the observed counts
      directly. No spread columns are added to the output.

  Default: `"rarefaction"`.

- n_boots:

  Integer. Number of bootstrap replicates used when `downsample_method`
  is `"rarefaction"` or `"vegan_rarefaction"`. Default: `100`. For
  reproducible results, call
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before invoking
  this function.

- feature_scale:

  Logical. If `TRUE`, a min-max scaled version of each diversity metric
  (range 0-1) is appended to the output as `<metric>_diversity_scaled`.
  Default: `FALSE`.

- result_name:

  Character string. Key under which the results `data.frame` is stored
  in `metadata(object)`. Default: `"nhood_diversity"`.

## Value

The input `object`, with a `data.frame` of diversity results stored in
`metadata(object)[[result_name]]`. The `data.frame` contains the
following columns:

- `nhood`:

  Neighbourhood identifier.

- `nhood_size`:

  Total number of cells in the neighbourhood.

- `n_contigs`:

  Number of cells with a contig.

- `contig_capture_rate`:

  Fraction of cells with a contig (`n_contigs / nhood_size`).

- `n_clones`:

  Observed number of unique clones.

- `clone_per_contig`:

  Ratio of unique clones to cells with a valid contig
  (`n_clones / n_contigs`). Values near 1 indicate high clonal
  diversity; values near 0 indicate strong clonal expansion.

- `coverage`:

  Good's coverage estimate (\\1 - singletons / n\\).

- `diversity_reliable`:

  Logical; `TRUE` if `n_contigs >= min_n_contigs`.

- `<metric>_diversity`:

  Estimated diversity for each requested metric.

- `<metric>_diversity_sd`:

  Standard deviation of the metric across bootstrap iterations. Present
  only when `downsample_method` is `"rarefaction"` or
  `"vegan_rarefaction"`. `NA` for `"chao1"` and `"ace"`, which are not
  bootstrapped.

- `<metric>_diversity_se`:

  Analytical standard error derived from `iNEXT` confidence intervals.
  Present only when `downsample_method = "coverage"`. `NA` for
  `"chao1"`, `"ace"`, `"norm.entropy"`, and `"gini"`, which are not
  coverage-standardised.

- `<metric>_diversity_scaled`:

  Min-max scaled diversity (only present when `feature_scale = TRUE`).

The following attributes are also set on the returned `data.frame`:
`min_n_contigs` and `method` are always present. `rarefy_depth` is set
only when `downsample_method` is `"rarefaction"` or
`"vegan_rarefaction"`. `target_coverage` is set only when
`downsample_method` is `"coverage"`.

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

### Clone identity resolution

`clone_id_col` is tried first; if it is `NULL` or absent from `colData`,
the function falls back to `cdr3_col` with a warning. Cells whose
clone/CDR3 value is `NA`, `""`, `"No_contig"`, or `"None"` are excluded
from diversity calculations but are still counted in `nhood_size`.

### Metric selection

The bias-corrected Chao1 formula (\\S\_{obs} + F_1(F_1-1) /
(2(F_2+1))\\) is used rather than the classicformula (\\S\_{obs} + F_1^2
/ (2 F_2)\\).

### Downsampling and metric compatibility

Rarefaction (`"rarefaction"` and `"vegan_rarefaction"`) is not applied
to `"chao1"` or `"ace"`, as subsampling removes the singleton and
doubleton structure these estimators rely on. Both are always computed
on observed counts under rarefaction; `_diversity_sd` is `NA` for these
metrics.

Under `"coverage"`, only `"shannon"`, `"inv.simpson"`, `"gini.simpson"`,
and `"standardised_n_clones"` are coverage-standardised via `iNEXT`.
`"chao1"`, `"ace"`, `"norm.entropy"`, and `"gini"` fall back to
observed-count calculation with a runtime warning; `_diversity_se` is
`NA` for these metrics.

`"standardised_n_clones"` is only valid with
`downsample_method = "coverage"`; requesting it under any other method
raises an error.

### Reproducibility

Bootstrap-based methods (`"rarefaction"`, `"vegan_rarefaction"`) use
[`base::sample`](https://rdrr.io/r/base/sample.html) internally. For
reproducible results, call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before invoking this
function.

## Note

The `vegan` package is required when
`downsample_method = "vegan_rarefaction"` and `iNEXT` is required when
`downsample_method = "coverage"`.

## See also

[`nhood_shm`](https://grace-hattersley.github.io/TEVRA/reference/nhood_shm.md),
[`nhood_cdr3_convergence`](https://grace-hattersley.github.io/TEVRA/reference/nhood_cdr3_convergence.md),
[`nhood_annotation`](https://grace-hattersley.github.io/TEVRA/reference/nhood_annotation.md),
[`nhood_annotation_props`](https://grace-hattersley.github.io/TEVRA/reference/nhood_annotation_props.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Milo object, default metrics and rarefaction -------------------------
milo <- nhood_diversity(
  object       = milo_object,
  clone_id_col = "clone_id",
  cdr3_col     = "cdr3_aa_VDJ"
)
results <- S4Vectors::metadata(milo)[["nhood_diversity"]]

# --- Coverage standardisation with standardised_n_clones ----------------
milo <- nhood_diversity(
  object            = milo_object,
  clone_id_col      = "clone_id",
  metrics           = c("shannon", "inv.simpson", "gini.simpson",
                        "standardised_n_clones"),
  downsample_method = "coverage",
  feature_scale     = TRUE,
  result_name       = "nhood_diversity_coverage"
)

# --- Plain SingleCellExperiment object ------------------------------------
sce <- nhood_diversity(
  object       = sce_object,
  clone_id_col = "clone_id",
  nhood_col    = "neighbourhood_id"
)
} # }
```
