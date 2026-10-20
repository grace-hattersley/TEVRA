#' Extract pairwise CDR3 convergence matrices for selected neighbourhoods
#'
#' Computes the full pairwise normalised CDR3 amino acid convergence matrix per
#' requested neighbourhood over representative sequences (clone centroids for
#' BCR, unique sequences for TCR). Convergence scores are defined as
#' \eqn{1 - \text{normalised distance}}, so higher values indicate greater
#' sequence similarity.
#'
#' This function is a targeted companion to
#' \code{\link{nhood_cdr3_convergence}}, intended for downstream
#' visualisation (e.g. heatmaps) of selected neighbourhoods. It reuses all
#' the same internal alignment helpers and applies the same filtering,
#' centroid-selection, and normalisation logic, but returns the convergence matrix
#' per neighbourhood instead of summarising the convergence scores.
#'
#' @param object A \code{\link[miloR]{Milo}} or
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}} object. Seurat
#'   objects must first be converted via
#'   \code{\link[Seurat]{as.SingleCellExperiment}}.
#' @param nhood_col Character string. For plain \code{SingleCellExperiment}
#'   objects (i.e. not \code{Milo}), the name of the \code{colData} column
#'   containing neighbourhood assignments. Not required for \code{Milo}
#'   objects.
#' @param nhood_ids Character or integer vector of neighbourhood identifiers
#'   to process. Values are matched against \code{colnames(nhoods(object))}
#'   for Milo objects, or against the values in \code{nhood_col} for
#'   SingleCellExperiment objects. An error is raised for any identifier not
#'   found.
#' @param receptor_type Character string; one of \code{"BCR"} (default) or
#'   \code{"TCR"}. Controls clone-level representative selection: for BCR,
#'   a centroid sequence is chosen per clone (the sequence with the lowest mean
#'   pairwise distance to all other members of the same clone). For TCR,
#'   unique sequences are used directly as representatives. Ignored if
#'   \code{clone_id_col} is \code{NULL}.
#' @param chain Character string; one of \code{"VDJ"} (default), \code{"VJ"},
#'   or \code{"both"}. Determines which receptor chain(s) are analysed. When
#'   \code{chain = "both"}, \code{aa_col} must be a named character vector with
#'   elements \code{"VDJ"} and \code{"VJ"} (see \code{aa_col}).
#' @param aa_col Character string or named character vector specifying the
#'   \code{colData} column(s) containing CDR3 amino acid sequences. When
#'   \code{chain} is \code{"VDJ"} or \code{"VJ"}, supply a single column name.
#'   When \code{chain = "both"}, supply a named vector with elements \code{"VDJ"}
#'   and \code{"VJ"}, e.g.:
#'   \code{aa_col = c(VDJ = "cdr3_aa_heavy", VJ = "cdr3_aa_light")}.
#'   Cells with \code{NA}, empty string, \code{"No_contig"}, or \code{"None"}
#'   in the relevant column are excluded from analysis.
#' @param clone_id_col Character string; name of the \code{colData} column
#'   containing clone identifiers (e.g. output from a clonotyping tool).
#'   When supplied and \code{receptor_type = "BCR"}, one centroid sequence is
#'   selected per clone per neighbourhood. When \code{NULL} (default), all
#'   unique CDR3 sequences are treated as individual representatives regardless
#'   of \code{receptor_type}. If the column is not found in \code{colData}, a
#'   warning is issued and the function falls back to unique-sequence mode.
#' @param min_n_contigs Integer; minimum number of cells with valid contig
#'   annotation required for a neighbourhood to be included in the analysis.
#'   Neighbourhoods below this threshold are assigned \code{NA} for all
#'   distance metrics and \code{converge_score_reliable = FALSE}. Default is
#'   \code{20}. When \code{min_clones} is also supplied, setting
#'   \code{min_n_contigs} below that value gives it no independent gating effect;
#'   a warning is issued in this case.
#' @param min_clones Integer or \code{NULL} (default). When supplied,
#'   neighbourhoods with fewer representative sequences than this threshold are
#'   excluded: all metrics are set to \code{NA} and \code{converge_score_reliable =
#'   FALSE}. When \code{NULL}, all neighbourhoods that pass \code{min_n_contigs}
#'   are processed regardless of representative count. Must be at least equal
#'   to \code{min_n_contigs} to have any independent gating effect; a warning is
#'   issued otherwise.
#'   \code{converge_score} may be unstable when \code{n_clones}
#'   is small relative to \eqn{k}: when the adaptive \eqn{k} approaches \code{n_clones}, each sequence's
#'   nearest neighbours span nearly the entire neighbourhood and the kNN score
#'   converges toward the global mean. If this is a concern, set
#'   \code{min_clones} to a value comfortably above the adaptive
#'   \eqn{k} floor (e.g. \code{20}) to exclude such neighbourhoods, or inspect
#'   \code{<chain>_n_clones} in the output to identify small
#'   neighbourhoods post-hoc. \code{converge_score} is reported as \code{NA} only
#'   when \code{n_clones <= resolved_k}, since kNN cannot be computed
#'   when there are insufficient sequences to form a neighbourhood.
#' @param distance_method Character string; one of \code{"levenshtein"}
#'   (default), \code{"blosum62"}, or \code{"custom"}.
#'   \describe{
#'     \item{\code{"levenshtein"}}{Flat substitution cost: 0 for identical
#'       amino acids, 1 for any mismatch, normalised by the length of the
#'       longer sequence. Equivalent to normalised edit distance.}
#'     \item{\code{"blosum62"}}{Substitution cost derived from BLOSUM62 via
#'       \eqn{c(i,j) = (s(i,i) + s(j,j)) / 2 - s(i,j)}, normalised to
#'       \eqn{[0, 1]} by the maximum off-diagonal cost. Biochemically similar
#'       amino acids incur lower cost than dissimilar ones.}
#'     \item{\code{"custom"}}{User-supplied substitution matrix passed via
#'       \code{substitution_matrix}. Must follow the same sign convention as
#'       \code{pwalign::stringDist} (higher score = more similar). See
#'       \code{substitution_matrix}.}
#'   }
#' @param substitution_matrix A square numeric matrix with named rows and
#'   columns (single-letter amino acid codes) to use when
#'   \code{distance_method = "custom"}. The matrix must follow the
#'   \code{pwalign::stringDist} sign convention: higher values indicate more
#'   similar amino acids (i.e. a score matrix, not a cost matrix). Ignored
#'   unless \code{distance_method = "custom"}.
#'
#' @return A named list with one element per requested neighbourhood ID. Each
#'   element is itself a named list with one entry per chain
#'   (\code{"VDJ"}, \code{"VJ"}, or both). Each chain entry is one of:
#'   \describe{
#'     \item{A symmetric numeric matrix}{Rows and columns are labelled by
#'       representative amino acid sequences. Values are normalised pairwise
#'       CDR3 convergence scores in \eqn{[0, 1]}, where 1 indicates identical
#'       sequences and 0 indicates maximum dissimilarity; the diagonal is 0.}
#'     \item{\code{NULL}}{The neighbourhood did not have enough cells or
#'       representatives to compute convergence scores (a warning is emitted).}
#'   }
#'
#' @seealso \code{\link{nhood_cdr3_convergence}} for the full
#'   neighbourhood-wide summary analysis.
#'
#' @examples
#' \dontrun{
#' # Extract convergence matrices for three neighbourhoods
#' conv_mats <- nhood_cdr3_convergence_matrices(
#'   object        = milo_object,
#'   nhood_ids     = c("5", "12", "27"),
#'   receptor_type = "BCR",
#'   chain         = "VDJ",
#'   aa_col        = "cdr3_aa_heavy",
#'   clone_id_col  = "clone_id"
#' )
#'
#' # Plot a heatmap for neighbourhood 5, VDJ chain
#' mat <- conv_mats[["5"]][["VDJ"]]
#' pheatmap::pheatmap(mat, main = "Nhood 5 CDR3 convergence")
#' }
#'
#' @export
#'
nhood_cdr3_convergence_matrices <- function(object,
                                            nhood_col           = NULL,
                                            nhood_ids,
                                            receptor_type       = c("BCR", "TCR"),
                                            chain               = c("VDJ", "VJ", "both"),
                                            aa_col,
                                            clone_id_col        = NULL,
                                            min_n_contigs       = 20L,
                                            min_clones          = NULL,
                                            distance_method     = c("levenshtein", "blosum62", "custom"),
                                            substitution_matrix = NULL) {

  # ── 0. Dependencies ────────────────────────────────────────────────────────
  for (pkg in c("pwalign", "Biostrings")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf(
        "Package '%s' is required. Install with: BiocManager::install('%s')", pkg, pkg
      ))
    }
  }

  # ── 1. Input validation ────────────────────────────────────────────────────
  receptor_type   <- match.arg(receptor_type)
  chain           <- match.arg(chain)
  distance_method <- match.arg(distance_method)

  if (chain == "both") {
    if (!is.character(aa_col) || is.null(names(aa_col)) ||
        !all(c("VDJ", "VJ") %in% names(aa_col))) {
      stop("When chain = 'both', aa_col must be a named character vector with elements 'VDJ' and 'VJ'.\n",
           "  Example: aa_col = c(VDJ = 'cdr3_aa_VDJ', VJ = 'cdr3_aa_VJ')")
    }
    chains_to_run <- c("VDJ", "VJ")
    aa_col_map    <- aa_col
  } else {
    if (!is.character(aa_col) || length(aa_col) != 1L) {
      stop("When chain = '", chain, "', aa_col must be a single column name string.")
    }
    chains_to_run <- chain
    aa_col_map    <- stats::setNames(aa_col, chain)
  }

  if (distance_method == "custom") {
    if (is.null(substitution_matrix)) {
      stop("'substitution_matrix' must be provided when distance_method = 'custom'.")
    }
    if (!is.matrix(substitution_matrix) || !is.numeric(substitution_matrix)) {
      stop("'substitution_matrix' must be a numeric matrix.")
    }
    if (is.null(rownames(substitution_matrix)) || is.null(colnames(substitution_matrix))) {
      stop("'substitution_matrix' must have named rows and columns (amino acid single-letter codes).")
    }
  }

  if (!is.null(min_clones)) {
    if (!is.numeric(min_clones) || length(min_clones) != 1L ||
        min_clones < 1L) {
      stop("'min_clones' must be a single positive integer or NULL.")
    }
    min_clones <- as.integer(min_clones)
    if (min_clones < min_n_contigs) {
      warning(sprintf(
        paste0(
          "'min_clones' (%d) is less than 'min_n_contigs' (%d). ",
          "This gives 'min_clones' no independent gating effect since ",
          "neighbourhoods failing 'min_n_contigs' are already excluded."
        ),
        min_clones, min_n_contigs
      ))
    }
  }

  score_mat <- .build_substitution_matrix(distance_method, substitution_matrix)

  # ── 2. Object type determination ───────────────────────────────────────────
  is_milo     <- inherits(object, "Milo") && inherits(object, "SingleCellExperiment")
  is_sce_only <- !is_milo && inherits(object, "SingleCellExperiment")

  if (!is_milo && !is_sce_only) {
    stop("'object' must be a Milo or SingleCellExperiment object.")
  }

  # ── 3. Extract colData and neighbourhood matrix ────────────────────────────
  col_data <- SummarizedExperiment::colData(object)

  missing_cols <- setdiff(unname(aa_col_map), colnames(col_data))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Column(s) not found in colData: %s",
                 paste(missing_cols, collapse = ", ")))
  }
  if (!is.null(clone_id_col) && !clone_id_col %in% colnames(col_data)) {
    warning(sprintf("clone_id_col '%s' not found in colData. Falling back to unique CDR3 sequences.",
                    clone_id_col))
    clone_id_col <- NULL
  }

  if (is_milo) {
    if (!requireNamespace("miloR", quietly = TRUE)) {
      stop("Package 'miloR' is required for Milo objects. Install with: BiocManager::install('miloR')")
    }
    nhoods        <- miloR::nhoods(object)
    all_nhood_ids <- colnames(nhoods)
  } else {
    if (is.null(nhood_col)) {
      stop("For SingleCellExperiment objects, 'nhood_col' must be specified.")
    }
    if (!nhood_col %in% colnames(col_data)) {
      stop(sprintf("Column '%s' not found in colData.", nhood_col))
    }
    nhood_assignments <- col_data[[nhood_col]]
    valid_cells       <- !is.na(nhood_assignments)
    if (any(!valid_cells)) {
      warning(sprintf("%d cells have NA neighbourhood assignments and will be excluded.",
                      sum(!valid_cells)))
    }
    unique_nhoods_raw <- unique(nhood_assignments[valid_cells])
    numeric_ids       <- suppressWarnings(as.numeric(unique_nhoods_raw))
    all_nhood_ids     <- if (!anyNA(numeric_ids)) {
      as.character(sort(numeric_ids))
    } else {
      as.character(sort(unique_nhoods_raw))
    }
    n_nhoods_total <- length(all_nhood_ids)
    nhoods <- matrix(0L,
                     nrow = ncol(object),
                     ncol = n_nhoods_total,
                     dimnames = list(colnames(object), all_nhood_ids))
    for (idx in seq_along(all_nhood_ids)) {
      cells_in <- as.character(nhood_assignments) == all_nhood_ids[idx] & valid_cells
      nhoods[cells_in, idx] <- 1L
    }
    if (requireNamespace("Matrix", quietly = TRUE)) {
      nhoods <- Matrix::Matrix(nhoods, sparse = TRUE)
    }
  }

  # ── 4. Validate requested nhood_ids ───────────────────────────────────────
  nhood_ids <- as.character(nhood_ids)
  missing_nhoods <- setdiff(nhood_ids, all_nhood_ids)
  if (length(missing_nhoods) > 0L) {
    stop(sprintf(
      "The following nhood_ids were not found: %s\nAvailable IDs: %s",
      paste(missing_nhoods, collapse = ", "),
      paste(utils::head(all_nhood_ids, 20L), collapse = ", ")
    ))
  }

  # ── 5. Pre-filter colData to cells with valid CDR3 annotation ─────────────
  col_data_filtered <- lapply(chains_to_run, function(ch) {
    aa_col_name <- aa_col_map[[ch]]
    valid_aa    <- !is.na(col_data[[aa_col_name]]) &
      !col_data[[aa_col_name]] %in% c("No_contig", "None", "")
    if (!is.null(clone_id_col)) {
      valid_aa <- valid_aa & !is.na(col_data[[clone_id_col]])
    }
    col_data[valid_aa, , drop = FALSE]
  })
  names(col_data_filtered) <- chains_to_run

  cell_row_index <- lapply(chains_to_run, function(ch) {
    cd <- col_data_filtered[[ch]]
    stats::setNames(seq_len(nrow(cd)), rownames(cd))
  })
  names(cell_row_index) <- chains_to_run

  nhood_barcodes <- rownames(nhoods)

  # ── 6. Process each requested neighbourhood ────────────────────────────────
  out <- lapply(nhood_ids, function(nid) {

    col_idx        <- match(nid, colnames(nhoods))
    cell_indices   <- which(nhoods[, col_idx] == 1L)
    nhood_barcodes_i <- nhood_barcodes[cell_indices]

    chain_out <- lapply(chains_to_run, function(ch) {

      aa_col_name <- aa_col_map[[ch]]
      idx_lookup  <- cell_row_index[[ch]]
      cd          <- col_data_filtered[[ch]]

      rows       <- idx_lookup[nhood_barcodes_i]
      rows       <- rows[!is.na(rows)]
      nhood_data <- cd[rows, , drop = FALSE]
      n_contigs  <- nrow(nhood_data)

      if (n_contigs < min_n_contigs) {
        warning(sprintf(
          "Neighbourhood '%s', chain '%s': only %d cells with valid contig annotation (min_n_contigs = %d). Returning NULL.",
          nid, ch, n_contigs, min_n_contigs
        ))
        return(NULL)
      }

      aa_seqs     <- nhood_data[[aa_col_name]]
      unique_seqs <- unique(aa_seqs)

      # Compute full pairwise distance matrix over all unique sequences,
      # then convert to convergence scores via 1 - normalised distance.
      full_dist_mat <- .run_pairwise_alignment(unique_seqs, score_mat)

      # Select representatives: centroid per clone for BCR, unique seqs for TCR
      if (receptor_type == "BCR" &&
          !is.null(clone_id_col) &&
          clone_id_col %in% colnames(nhood_data)) {

        clone_ids     <- nhood_data[[clone_id_col]]
        unique_clones <- unique(clone_ids)

        centroid_indices <- vapply(unique_clones, function(cid) {
          member_seqs <- unique(aa_seqs[clone_ids == cid])
          if (length(member_seqs) == 1L) {
            return(match(member_seqs, unique_seqs))
          }
          clone_idx  <- match(member_seqs, unique_seqs)
          clone_sub  <- full_dist_mat[clone_idx, clone_idx, drop = FALSE]
          mean_dists <- rowMeans(clone_sub, na.rm = TRUE)
          clone_idx[which.min(mean_dists)]
        }, FUN.VALUE = integer(1))

        rep_dist_mat <- full_dist_mat[centroid_indices, centroid_indices, drop = FALSE]

      } else {
        rep_dist_mat <- full_dist_mat
      }

      if (nrow(rep_dist_mat) < 2L) {
        warning(sprintf(
          "Neighbourhood '%s', chain '%s': fewer than 2 representatives after centroid selection. Returning NULL.",
          nid, ch
        ))
        return(NULL)
      }

      # Optional minimum representative gate, applied on top of the
      # unconditional 2-representative floor above. Mirrors min_clones
      # behaviour in nhood_cdr3_convergence.
      n_reps <- nrow(rep_dist_mat)
      if (!is.null(min_clones) && n_reps < min_clones) {
        warning(sprintf(
          "Neighbourhood '%s', chain '%s': only %d representative(s) after centroid selection (min_clones = %d). Returning NULL.",
          nid, ch, n_reps, min_clones
        ))
        return(NULL)
      }

      # Convert distance scores to convergence scores: 1 - distance.
      # Higher values indicate greater sequence similarity.
      # Reset diagonal to 0 (self-comparison) after inversion.
      diag(rep_dist_mat) <- 0
      conv_mat <- 1 - rep_dist_mat
      conv_mat
    })
    names(chain_out) <- chains_to_run
    chain_out
  })
  names(out) <- nhood_ids

  out
}
