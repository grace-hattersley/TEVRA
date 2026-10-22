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

- substitution_matrix:

  A square numeric matrix with single-letter amino acid codes as named
  rows and columns. The matrix must be bounded -1-0, where higher values
  indicate more similar amino acids. Used when
  `substitution_weighting = "custom"`.

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

### Summary of steps

For each neighbourhood the following steps are performed:

1.  A single pairwise distance matrix is computed over all unique CDR3
    sequences per neighbourhood.

2.  For BCR data with clone labels, one centroid sequence is selected
    per clone as the sequence with the lowest mean pairwise distance to
    all other unique sequences within that clone in the neighbourhood.
    The raw pairwise distance matrix is filtered to these centroids.

3.  The distance matrix is converted to a convergence matrix and
    returned without summarisation

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
longer sequence and converted to convergence scores via 1 - distance.

### Selecting a representative sequence for BCRs

Pairwise distances are computed between all unique CDR3 amino acid
sequences per aggregate using the above method. However, for BCRs, some
of these unique sequences belong to the same clone due to the effects of
somatic hypermutation. In this scenario, the CDR3 sequence with the
lowest average pairwise distance to all other sequences in the clone is
selected. All other sequences are discarded from the per aggregate
distance matrix, such that it contains only one sequence per clone.

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
