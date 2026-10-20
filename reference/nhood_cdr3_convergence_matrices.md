# Extract pairwise CDR3 convergence matrices for selected neighbourhoods

Computes the full pairwise normalised CDR3 amino acid convergence matrix
per requested neighbourhood over representative sequences (clone
centroids for BCR, unique sequences for TCR). Convergence scores are
defined as \\1 - \text{normalised distance}\\, so higher values indicate
greater sequence similarity.

## Usage

``` r
nhood_cdr3_convergence_matrices(
  object,
  nhood_col = NULL,
  nhood_ids,
  receptor_type = c("BCR", "TCR"),
  chain = c("VDJ", "VJ", "both"),
  aa_col,
  clone_id_col = NULL,
  min_n_contigs = 20L,
  min_clones = NULL,
  distance_method = c("levenshtein", "blosum62", "custom"),
  substitution_matrix = NULL
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

- nhood_ids:

  Character or integer vector of neighbourhood identifiers to process.
  Values are matched against `colnames(nhoods(object))` for Milo
  objects, or against the values in `nhood_col` for SingleCellExperiment
  objects. An error is raised for any identifier not found.

- receptor_type:

  Character string; one of `"BCR"` (default) or `"TCR"`. Controls
  clone-level representative selection: for BCR, a centroid sequence is
  chosen per clone (the sequence with the lowest mean pairwise distance
  to all other members of the same clone). For TCR, unique sequences are
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
  identifiers (e.g. output from a clonotyping tool). When supplied and
  `receptor_type = "BCR"`, one centroid sequence is selected per clone
  per neighbourhood. When `NULL` (default), all unique CDR3 sequences
  are treated as individual representatives regardless of
  `receptor_type`. If the column is not found in `colData`, a warning is
  issued and the function falls back to unique-sequence mode.

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
  post-hoc. `converge_score` is reported as `NA` only when
  `n_clones <= resolved_k`, since kNN cannot be computed when there are
  insufficient sequences to form a neighbourhood.

- distance_method:

  Character string; one of `"levenshtein"` (default), `"blosum62"`, or
  `"custom"`.

  `"levenshtein"`

  :   Flat substitution cost: 0 for identical amino acids, 1 for any
      mismatch, normalised by the length of the longer sequence.
      Equivalent to normalised edit distance.

  `"blosum62"`

  :   Substitution cost derived from BLOSUM62 via \\c(i,j) = (s(i,i) +
      s(j,j)) / 2 - s(i,j)\\, normalised to \\\[0, 1\]\\ by the maximum
      off-diagonal cost. Biochemically similar amino acids incur lower
      cost than dissimilar ones.

  `"custom"`

  :   User-supplied substitution matrix passed via
      `substitution_matrix`. Must follow the same sign convention as
      [`pwalign::stringDist`](https://rdrr.io/pkg/pwalign/man/stringDist.html)
      (higher score = more similar). See `substitution_matrix`.

- substitution_matrix:

  A square numeric matrix with named rows and columns (single-letter
  amino acid codes) to use when `distance_method = "custom"`. The matrix
  must follow the
  [`pwalign::stringDist`](https://rdrr.io/pkg/pwalign/man/stringDist.html)
  sign convention: higher values indicate more similar amino acids (i.e.
  a score matrix, not a cost matrix). Ignored unless
  `distance_method = "custom"`.

## Value

A named list with one element per requested neighbourhood ID. Each
element is itself a named list with one entry per chain (`"VDJ"`,
`"VJ"`, or both). Each chain entry is one of:

- A symmetric numeric matrix:

  Rows and columns are labelled by representative amino acid sequences.
  Values are normalised pairwise CDR3 convergence scores in \\\[0,
  1\]\\, where 1 indicates identical sequences and 0 indicates maximum
  dissimilarity; the diagonal is 0.

- `NULL`:

  The neighbourhood did not have enough cells or representatives to
  compute convergence scores (a warning is emitted).

## Details

This function is a targeted companion to
[`nhood_cdr3_convergence`](https://grace-hattersley.github.io/TEVRA/reference/nhood_cdr3_convergence.md),
intended for downstream visualisation (e.g. heatmaps) of selected
neighbourhoods. It reuses all the same internal alignment helpers and
applies the same filtering, centroid-selection, and normalisation logic,
but returns the convergence matrix per neighbourhood instead of
summarising the convergence scores.

## See also

[`nhood_cdr3_convergence`](https://grace-hattersley.github.io/TEVRA/reference/nhood_cdr3_convergence.md)
for the full neighbourhood-wide summary analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract convergence matrices for three neighbourhoods
conv_mats <- nhood_cdr3_convergence_matrices(
  object        = milo_object,
  nhood_ids     = c("5", "12", "27"),
  receptor_type = "BCR",
  chain         = "VDJ",
  aa_col        = "cdr3_aa_heavy",
  clone_id_col  = "clone_id"
)

# Plot a heatmap for neighbourhood 5, VDJ chain
mat <- conv_mats[["5"]][["VDJ"]]
pheatmap::pheatmap(mat, main = "Nhood 5 CDR3 convergence")
} # }
```
