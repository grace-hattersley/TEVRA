# Compute inter-clonal CDR3 amino acid convergence scores per neighbourhood

Computes pairwise complementary determining region 3 (CDR3) amino acid
sequence similarities between different B cell receptor (BCR) or T cell
receptor (TCR) clones within a neighbourhood, and summarises this score
per neighbourhood.

## Usage

``` r
nhood_cdr3_convergence(
  object,
  nhood_col = NULL,
  receptor_type = c("BCR", "TCR"),
  chain = c("VDJ", "VJ", "both"),
  aa_col,
  clone_id_col = NULL,
  min_n_contigs = 20L,
  min_clones = NULL,
  substitution_weighting = c("levenshtein", "blosum62", "custom"),
  substitution_matrix = NULL,
  summary_metric = "converge_score",
  knn_k = NULL,
  downsample = TRUE,
  n_boots = 100,
  seed = NULL,
  BPPARAM = BiocParallel::SerialParam(),
  result_name = "nhood_cdr3_convergence"
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

- receptor_type:

  Character string; one of `"BCR"` (default) or `"TCR"`. Controls
  clone-level representative selection: for BCRs, CDR3 sequences
  belonging to the same clone may have different sequences due to
  somatic hypermutation (SHM), hence a centroid sequence is chosen per
  clone (the sequence with the lowest mean pairwise distance to all
  other members of the same clone). For TCRs, unique sequences can be
  used directly as representatives. Ignored if `clone_id_col` is `NULL`.

- chain:

  Character string; one of `"VDJ"` (default), `"VJ"`, or `"both"`.
  Determines which receptor chain(s) are analysed. When
  `chain = "both"`, `aa_col` must be a named character vector with
  elements `"VDJ"` and `"VJ"` (see `aa_col`).

- aa_col:

  Character string or named character vector specifying the `colData`
  column(s) containing CDR3 amino acid sequences. When `chain` is
  `"VDJ"` or `"VJ"`, supply a single column name. When `chain = "both"`,
  supply a named vector with elements `"VDJ"` and `"VJ"`, e.g.:
  `aa_col = c(VDJ = "cdr3_aa_heavy", VJ = "cdr3_aa_light")`. Cells with
  `NA`, empty string, `"No_contig"`, or `"None"` in the relevant column
  are excluded from analysis.

- clone_id_col:

  Character string; name of the `colData` column containing clone
  identifiers. When supplied and `receptor_type = "BCR"`, one centroid
  sequence is selected per clone per neighbourhood. When `NULL`
  (default), all unique CDR3 sequences are treated as individual
  representatives regardless of `receptor_type`. If the column is not
  found in `colData`, a warning is issued and the function falls back to
  unique-sequence mode.

- min_n_contigs:

  Integer; minimum number of cells with valid contig annotation required
  for a neighbourhood to be included in the analysis. Neighbourhoods
  below this threshold are assigned `NA` for all distance metrics and
  `converge_score_reliable = FALSE`. Default is `20`. When `min_clones`
  is also supplied, setting `min_n_contigs` below that value gives it no
  independent gating effect; a warning is issued in this case.

- min_clones:

  Integer or `NULL` (default). When supplied, neighbourhoods with fewer
  representative sequences than this threshold are excluded: all metrics
  are set to `NA` and `converge_score_reliable = FALSE`. When `NULL`,
  all neighbourhoods that pass `min_n_contigs` are processed regardless
  of representative count. Must be at least equal to `min_n_contigs` to
  have any independent gating effect; a warning is issued otherwise.
  `converge_score` may be unstable when `n_clones` is small relative to
  \\k\\: when the adaptive \\k\\ approaches `n_clones`, each sequence's
  nearest neighbours span nearly the entire neighbourhood and the kNN
  score converges toward the global mean. If this is a concern, set
  `min_clones` to a value comfortably above the adaptive \\k\\ floor
  (e.g. `20`) to exclude such neighbourhoods, or inspect
  `<chain>_n_clones` in the output to identify small neighbourhoods
  post-hoc.

- substitution_weighting:

  Character string; one of `"levenshtein"` (default), `"blosum62"`, or
  `"custom"`.

  `"levenshtein"`

  :   Flat substitution cost: 0 for identical amino acids, 1 for any
      mismatch.

  `"blosum62"`

  :   Substitution cost derived from BLOSUM62 via \\c(i,j) = (s(i,i) +
      s(j,j)) / 2 - s(i,j)\\, normalised to \\\[0, 1\]\\ by the maximum
      off-diagonal cost. Physicochemically similar amino acids incur
      lower cost than dissimilar ones.

  `"custom"`

  :   User-supplied substitution weightings passed as a matrix to
      `substitution_matrix`.

- substitution_matrix:

  A square numeric matrix with single-letter amino acid codes as named
  rows and columns. The matrix must be bounded -1-0, where higher values
  indicate more similar amino acids. Used when
  `substitution_weighting = "custom"`.

- summary_metric:

  Character vector; one or more of `"converge_score"` (default),
  `"converge_score_global"`, and `"converge_score_global_median"`.
  Controls which summary convergence metrics are calculated:

  `"converge_score"`

  :   Local convergence score computed via a kNN-based approach. `NA`
      when `n_clones <= resolved_k`, or when `min_clones` is supplied
      and not met.

  `"converge_score_global"`

  :   Mean pairwise CDR3 convergence score across the upper triangle of
      the representative distance matrix.

  `"converge_score_global_median"`

  :   Median pairwise CDR3 convergence score across the upper triangle
      of the representative distance matrix.

  Any combination of the three values may be supplied. The default is
  `"converge_score"`.

- knn_k:

  Integer or `NULL` (default). Number of nearest neighbours to use when
  computing `converge_score`. When `NULL`, \\k\\ is set adaptively to
  \\\lfloor 0.25 \times n\_{\text{reps}} \rfloor\\ per neighbourhood. 3
  is the minimum permitted adaptive \\k\\ to ensure the metric captures
  more than just the nearest pairwise comparisons (not applied to a
  fixed \\k\\). Within each bootstrap subsample \\k\\ is additionally
  capped at \\n - 1\\ to guard against subsamples smaller than \\k\\.

- downsample:

  Logical; whether to downsample large neighbourhoods to a common
  representative depth before computing distances. The downsampling
  depth is set to the median number of unique representatives across all
  passing neighbourhoods (excluding singletons), computed separately per
  chain. Downsampling is performed via bootstrapping (`n_boots`
  iterations). Default is `TRUE`. When `FALSE`, the full representative
  matrix is used directly and `_sd` columns are not reported.

- n_boots:

  Integer; number of bootstrap iterations used during downsampling. Each
  iteration subsamples `downsample_depth` representatives without
  replacement and computes distance summaries; reported values are the
  mean across iterations and the corresponding `_sd` column reports the
  standard deviation across iterations for each requested metric.
  Default is `100`. Ignored when `downsample = FALSE` or when a
  neighbourhood has fewer representatives than the downsampling depth.
  When `downsample = TRUE` and `"converge_score"` is in
  `summary_metric`, the effective sample size for kNN computation within
  each bootstrap iteration is the downsampling depth rather than the
  full representative count.

- seed:

  Integer or `NULL` (default). Random seed for reproducibility of the
  bootstrap downsampling procedure. When supplied, a
  neighbourhood-specific seed is derived as `seed + i` (where `i` is the
  neighbourhood index), ensuring that each neighbourhood's bootstrap is
  independently reproducible while the full run remains deterministic.
  Set this parameter whenever results will be reported or compared
  across runs. Has no effect when `downsample = FALSE`.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  instance controlling parallelisation. The default,
  [`SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html),
  runs sequentially with per-neighbourhood progress messages. Parallel
  backends (e.g.
  [`MulticoreParam`](https://rdrr.io/pkg/BiocParallel/man/MulticoreParam-class.html),
  [`SnowParam`](https://rdrr.io/pkg/BiocParallel/man/SnowParam-class.html))
  suppress per-neighbourhood messages but respect all other arguments.
  Note that when using process-based parallel backends, all data objects
  are serialised to each worker, which may increase peak memory usage
  substantially for large datasets.

- result_name:

  Character string; the key under which results are stored in
  `metadata(object)`. If an entry with this name already exists a
  warning is emitted and the entry is overwritten. Default is
  `"nhood_cdr3_convergence"`.

## Value

The input `object` with a `data.frame` added to
`metadata(object)[[result_name]]`. The data.frame has one row per
neighbourhood and the following columns:

- `nhood`:

  Neighbourhood identifier.

- `nhood_size`:

  Total number of cells in the neighbourhood, including cells without a
  contig. Used as the denominator for `contig_capture_rate`.

- `<chain>_n_contigs`:

  Number of cells with a contig for this chain.

- `<chain>_contig_capture_rate`:

  Proportion of neighbourhood cells with a contig
  (`n_contigs / nhood_size`). `NA` if `nhood_size` is zero.

- `<chain>_n_clones`:

  Number of representative sequences used for distance computation
  (unique sequences for contigs without clone labels; one centroid per
  clone for BCR with clone labels).

- `<chain>_converge_score_reliable`:

  Logical; `FALSE` when `n_contigs < min_n_contigs` (all metrics `NA`),
  or when `min_clones` is supplied and `n_clones` does not meet it (all
  metrics `NA`). `TRUE` for all other neighbourhoods, including those
  where `converge_score` is `NA` due to insufficient sequences for kNN.

- `CDR3_converge_score_<chain>`:

  kNN convergence score per neighbourhood, averaged across all
  representatives. Captures local sequence density. When downsampling,
  this is the mean across bootstrap iterations. `NA` when
  `n_clones <= resolved_k` (insufficient sequences for kNN), or when
  `min_clones` is supplied and `n_clones` does not meet it. Only present
  when `"converge_score"` is in `summary_metric`. Example:
  `CDR3_converge_score_VDJ`.

- `CDR3_converge_score_sd_<chain>`:

  Standard deviation of `CDR3_converge_score_<chain>` across bootstrap
  iterations. Only present when `downsample = TRUE` and
  `"converge_score"` is in `summary_metric`. `NA` when `converge_score`
  is `NA` (kNN could not be computed). Example:
  `CDR3_converge_score_sd_VDJ`.

- `CDR3_converge_score_global_<chain>`:

  Mean pairwise CDR3 convergence score across the upper triangle of the
  representative distance matrix. When downsampling, this is the mean
  across bootstrap iterations. Only present when
  `"converge_score_global"` is in `summary_metric`. Example:
  `CDR3_converge_score_global_VDJ`.

- `CDR3_converge_score_global_sd_<chain>`:

  Standard deviation of `CDR3_converge_score_global_<chain>` across
  bootstrap iterations. Only present when `downsample = TRUE` and
  `"converge_score_global"` is in `summary_metric`. Example:
  `CDR3_converge_score_global_sd_VDJ`.

- `CDR3_converge_score_global_median_<chain>`:

  Median pairwise CDR3 convergence score across the upper triangle of
  the representative distance matrix. When downsampling, this is the
  mean of bootstrap-iteration medians. Only present when
  `"converge_score_global_median"` is in `summary_metric`. Example:
  `CDR3_converge_score_global_median_VJ`.

- `CDR3_converge_score_global_median_sd_<chain>`:

  Standard deviation of `CDR3_converge_score_global_median_<chain>`
  across bootstrap iterations. Only present when `downsample = TRUE` and
  `"converge_score_global_median"` is in `summary_metric`. Example:
  `CDR3_converge_score_global_median_sd_VDJ`.

Where `chain = "both"`, columns are present for both `VDJ` and `VJ`
chains.

The returned `data.frame` also carries the following attributes
recording the analysis parameters: `receptor_type`,
`substitution_weighting`, `chains`, `downsample_depths`,
`min_n_contigs`, `min_clones`, `downsample`, `n_boots`,
`summary_metric`.

## Details

### Summary of steps

For each neighbourhood the following steps are performed:

1.  A single pairwise distance matrix is computed over all unique CDR3
    sequences per neighbourhood.

2.  For BCR data with clone labels, one centroid sequence is selected
    per clone as the sequence with the lowest mean pairwise distance to
    all other unique sequences within that clone in the neighbourhood.
    The raw pairwise distance matrix is filtered to these centroids.

3.  If `downsample = TRUE` and the number of representatives exceeds the
    target depth, bootstrapped subsampling is performed over the
    pre-computed distance matrix (no re-alignment), and summary
    statistics are averaged across iterations.

4.  The distance matrix is summarised per neighbourhood using the
    requested metric

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

### Distance calculation methods

Pairwise distances are computed using
[`pwalign::stringDist()`](https://bioconductor.org/packages/release/bioc/manuals/pwalign/man/pwalign.pdf)
across all unique CDR3 amino acid sequences per neighbourhood requested
exceeding min_n_contigs and optionally min_clones. CDR3 sequences can
differ via insertions, deletions or substitutions. In all cases,
insertions and deletions are assigned a cost of 1. When
`substitution_weighting = "levenshtein"`, all substitutions also incur a
cost of 1. Substitutions can also be weighted according to their
physicochemical significance when `substitution_weighting = "blosum62"`.
The raw BLOSUM substitution matrix is converted to a cost matrix via
\\c(i,j) = (s(i,i) + s(j,j))/2 - s(i,j)\\. This cost matrix is then
normalised to 0-1 by dividing each cost by the maximum off-diagonal
cost. The normalised cost matrix is then negated such that the most
costly amino acid substitutions score -1 and the least significant
substitutions score 0. If the user wishes to use an alternative method
of substitution weighting via `substitution_weighting = "custom"`, they
must provide their own substitution matrix which matches the same logic
(bounded -1-0, more similar amino acids have a higher score). Overall,
independent of the substitution scoring used, the maximum distance that
can be achieved is the length of the longest sequence (all amino acids
substituted incurrring the maximum possible cost if blosum-weighted or
custom substitution scores have been used, insertion of the extra amino
acids necessary to account for sequence length variation). All distances
are subsequently normalised to 0-1 by dividing by the length of the
longer sequence.

### Selecting a representative sequence for BCRs

Pairwise distances are computed between all unique CDR3 amino acid
sequences per aggregate using the above method. However, for BCRs, some
of these unique sequences belong to the same clone due to the effects of
somatic hypermutation. In this scenario, the CDR3 sequence with the
lowest average pairwise distance to all other sequences in the clone is
selected. All other sequences are discarded from the per aggregate
distance matrix, such that it contains only one sequence per clone.

### Downsampling depth

Downsampling is recommended where neighbourhood sizes differ
substantially to ensure distances are comparable across
neighbourhoods.The downsampling depth is the median number of unique
representatives across all neighbourhoods that pass `min_n_contigs`,
excluding singleton- representative neighbourhoods. This is computed
separately per chain when `chain = "both"`, as VDJ and VJ clone counts
can differ. The depth is reported in the progress messages emitted to
the console.

### Summarising the distance matrix and obtaining a convergence score

The distance matrix is then summarised to obtain a singular convergence
score representing the similarity of CDR3 sequences of different clones
per aggregate. This can be a global convergence score, defined as the
mean or median pairwise distance across all representative CDR3
sequences, or a local score. Local convergence scores are computed as
the mean distance of each CDR3 sequence to its k nearest neighbours,
averaged across all representatives in the aggregate. This metric is
more sensitive to repertoire structure than a global average. For
example, a neighbourhood containing clones that respond to two distinct
epitopes would contain two clusters of highly similar CDR3s that are
dissimilar from one another. Here, global convergence would be high due
to the difference between these two antigens, whilst local convergence
would reflect sequence similarity within these two groups. By default, k
is set adaptively at n/4,where n is the number of representatives, with
a minimum k of 3; a fixed k can also be supplied. Local convergence is
not calculated when the number of representatives does not exceed k. All
distance-based metrics are converted to convergence scores via 1 -
distance.

### Parallelisation

Each neighbourhood is processed independently and dispatched as a
separate task via
[`bplapply`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html). With
[`SerialParam`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
(default), progress is reported at every 10\\ progress is suppressed; a
completion summary is printed after collation.

### Neighbourhood size and capture rate

`nhood_size` and `contig_capture_rate` are based on all cells assigned
to the neighbourhood in the binary nhoods matrix, regardless of contig
annotation status. This reflects the overall cellular composition of the
neighbourhood rather than immune cell abundance alone, and is the
intended denominator for capture rate.

## See also

[`pwalign::stringDist()`](https://bioconductor.org/packages/release/bioc/manuals/pwalign/man/pwalign.pdf)
for the underlying alignment engine.
[`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
for parallelisation options.
[`nhood_diversity`](https://grace-hattersley.github.io/TEVRA/reference/nhood_diversity.md),
[`nhood_shm`](https://grace-hattersley.github.io/TEVRA/reference/nhood_shm.md),
[`nhood_annotation`](https://grace-hattersley.github.io/TEVRA/reference/nhood_annotation.md),
[`nhood_annotation_props`](https://grace-hattersley.github.io/TEVRA/reference/nhood_annotation_props.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Milo object, BCR, VDJ chain only, Levenshtein distance -------------
milo <- nhood_cdr3_convergence(
  object        = milo_object,
  receptor_type = "BCR",
  chain         = "VDJ",
  aa_col        = "cdr3_aa_heavy",
  clone_id_col  = "clone_id",
  seed          = 42
)
results <- S4Vectors::metadata(milo)$nhood_cdr3_convergence

# --- Milo object, TCR, both chains, BLOSUM62, all metrics ----------------
milo <- nhood_cdr3_convergence(
  object          = milo_object,
  receptor_type   = "TCR",
  chain           = "both",
  aa_col          = c(VDJ = "cdr3_aa_TRB", VJ = "cdr3_aa_TRA"),
  substitution_weighting = "blosum62",
  summary_metric  = c("converge_score", "converge_score_global", "converge_score_global_median"),
  min_n_contigs   = 5,
  seed            = 42,
  BPPARAM         = BiocParallel::MulticoreParam(workers = 4)
)

# --- SingleCellExperiment object, BCR, custom substitution matrix --------
custom_mat <- matrix(...)  # 20x20 numeric, named by single-letter AA code
sce <- nhood_cdr3_convergence(
  object              = sce_object,
  receptor_type       = "BCR",
  chain               = "VDJ",
  aa_col              = "cdr3_aa_heavy",
  clone_id_col        = "clone_id",
  nhood_col           = "nhood_id",
  substitution_weighting  = "custom",
  substitution_matrix = custom_mat,
  seed                = 42
)
} # }

```
