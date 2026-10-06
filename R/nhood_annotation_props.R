#' Calculate Cell Annotation Proportions per Neighbourhood
#'
#' Computes the proportion of cells assigned to each neighbourhood that carry
#' each category of one or more annotation columns.
#'
#' @param object A \code{\link[miloR]{Milo}} or
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}} object. Seurat
#'   objects must first be converted via
#'   \code{\link[Seurat]{as.SingleCellExperiment}}.
#' @param annotation_cols Character vector of \code{colData} column names containing
#'   per-cell annotation categories to summarise per neighbourhood (e.g.
#'   cell type, lineage, cluster label). All specified columns must be present
#'   in \code{colData}; missing columns raise an error. Default:
#'   \code{"celltype"}.
#' @param nhood_col Character string. For plain \code{SingleCellExperiment}
#'   objects (i.e. not \code{Milo}), the name of the \code{colData} column
#'   containing neighbourhood assignments. Not required for \code{Milo}
#'   objects. Default: \code{"nhood"}.
#' @param result_name Character string. Key under which the results
#'   \code{data.frame} is stored in \code{metadata(object)}. If an entry with
#'   this name already exists a warning is emitted and the entry is
#'   overwritten. Default: \code{"nhood_annotation_props"}.
#'
#' @return The input \code{object}, with a \code{data.frame} stored in
#'   \code{metadata(object)[[result_name]]}. The \code{data.frame} has one row
#'   per neighbourhood and the following columns:
#'   \describe{
#'     \item{\code{nhood}}{Neighbourhood identifier.}
#'     \item{\code{nhood_size}}{Total number of cells in the neighbourhood.}
#'     \item{\code{<annotation_col>_<category>_prop}}{Proportion of cells carrying
#'       each annotation category. One column per category per entry in
#'       \code{annotation_cols} — e.g. if \code{annotation_cols = "celltype"} and the
#'       column contains \code{"CD4"} and \code{"CD8"}, the output columns will
#'       be \code{celltype_CD4_prop} and \code{celltype_CD8_prop}.
#'       Neighbourhoods where all cells have \code{NA} for a given annotation
#'       column have \code{NA} in all corresponding proportion columns.}
#'   }
#'
#' @details
#' ## Object type handling
#' \itemize{
#'   \item \strong{Milo objects}: Neighbourhood membership is read directly
#'     from the sparse \code{nhoods} matrix via \code{\link[miloR]{nhoods}}.
#'     Row names of this matrix must match \code{colnames(object)}; a mismatch
#'     raises an error to prevent silent misassignment.
#'   \item \strong{Plain SingleCellExperiment objects}: A binary neighbourhood
#'     membership matrix is constructed internally from the column specified by
#'     \code{nhood_col}. Cells with \code{NA} assignments are excluded with a
#'     warning.
#' }
#'
#' ## Proportion calculation
#' For each neighbourhood and each annotation column, per-cell labels are
#' tabulated across all valid (non-\code{NA}) cells. Proportions are computed
#' relative to the number of valid cells in the neighbourhood for that
#' annotation column (i.e. each row sums to 1 over the proportion columns for
#' a given \code{annotation_col}, excluding any \code{NA} cells). Cells with
#' \code{NA} in a given annotation column are excluded from that column's
#' proportions and their count is reported via a warning.
#'
#' @seealso
#' \code{\link{nhood_annotation}}, \code{\link{nhood_diversity}},
#' \code{\link{nhood_shm}}, \code{\link{nhood_cdr3_convergence}}
#'
#' @examples
#' \dontrun{
#' # --- Milo object, single annotation column --------------------------------
#' milo <- nhood_annotation_props(
#'   object          = milo_object,
#'   annotation_cols = "celltype"
#' )
#' results <- S4Vectors::metadata(milo)[["nhood_annotation_props"]]
#' # Output columns: nhood, nhood_size, celltype_CD4_prop, celltype_CD8_prop, ...
#'
#' # --- Multiple annotation columns ------------------------------------------
#' milo <- nhood_annotation_props(
#'   object          = milo_object,
#'   annotation_cols = c("celltype", "lineage"),
#'   result_name     = "nhood_annotation_props"
#' )
#'
#' # --- Plain SingleCellExperiment object ------------------------------------
#' sce <- nhood_annotation_props(
#'   object          = sce_object,
#'   annotation_cols = "celltype",
#'   nhood_col       = "neighbourhood_id"
#' )
#' }
#'
#' @export
#'
nhood_annotation_props <- function(object,
                                   annotation_cols = "celltype",
                                   nhood_col       = "nhood",
                                   result_name     = "nhood_annotation_props") {

  # ── 1. Determine object type ────────────────────────────────────────────────
  is_milo     <- inherits(object, "Milo") && inherits(object, "SingleCellExperiment")
  is_sce_only <- !is_milo && inherits(object, "SingleCellExperiment")

  if (!is_milo && !is_sce_only) {
    stop("'object' must be a Milo or SingleCellExperiment object. If using Seurat, convert via as.SingleCellExperiment")
  }

  # ── 2. Extract colData ───────────────────────────────────────────────────────
  cd <- as.data.frame(SummarizedExperiment::colData(object))

  # ── 3. Validate annotation columns ──────────────────────────────────────────
  missing_cols <- annotation_cols[!annotation_cols %in% colnames(cd)]
  if (length(missing_cols) > 0) {
    stop("Annotation column(s) not found in colData: ", paste(missing_cols, collapse = ", "))
  }

  # ── 4. Build neighbourhood membership structure ─────────────────────────────
  if (is_milo) {
    if (!requireNamespace("miloR", quietly = TRUE)) {
      stop("Package 'miloR' is required for Milo objects. Please install it with: BiocManager::install('miloR')")
    }
    if (is.null(miloR::nhoods(object))) {
      stop("nhoods slot not found in Milo object")
    }

    nhood_mat <- miloR::nhoods(object)

    if (nrow(nhood_mat) != ncol(object)) {
      stop(sprintf(
        "nhoods matrix has %d rows but the object has %d cells. These must match.",
        nrow(nhood_mat), ncol(object)
      ))
    }

    if (!is.null(rownames(nhood_mat)) && !is.null(colnames(object))) {
      if (!identical(rownames(nhood_mat), colnames(object))) {
        stop(
          "rownames(nhoods(object)) do not match colnames(object). ",
          "Cell ordering is inconsistent; cannot safely assign annotation values."
        )
      }
    } else {
      message(
        "Note: rownames(nhoods(object)) or colnames(object) are NULL. ",
        "Assuming cell order in nhoods matrix matches colData row order."
      )
    }

    nhood_ids <- colnames(nhood_mat)
    n_nhoods  <- ncol(nhood_mat)

    # Pre-compute cell indices per neighbourhood (reused across all annotation_cols)
    nhood_cell_idx <- lapply(seq_len(n_nhoods), function(j) which(nhood_mat[, j] != 0))

  } else {

    if (!nhood_col %in% colnames(cd)) {
      stop(sprintf("Column '%s' not found in colData. Check 'nhood_col'.", nhood_col))
    }

    nhood_assignments <- cd[[nhood_col]]

    if (any(is.na(nhood_assignments))) {
      warning(sprintf(
        "%d cells have NA neighbourhood assignments and will be excluded.",
        sum(is.na(nhood_assignments))
      ))
    }

    nhood_ids <- sort(unique(nhood_assignments[!is.na(nhood_assignments)]))

    if (length(nhood_ids) == 0) {
      stop(sprintf("No valid neighbourhood assignments found in colData column '%s'.", nhood_col))
    }

    n_nhoods <- length(nhood_ids)

    # Pre-compute cell indices per neighbourhood
    nhood_cell_idx <- lapply(nhood_ids, function(nh) which(nhood_assignments == nh))
  }

  # ── 5. Initialise results data frame ────────────────────────────────────────
  out_df <- data.frame(
    nhood      = nhood_ids,
    nhood_size = vapply(nhood_cell_idx, length, integer(1)),
    stringsAsFactors = FALSE
  )

  # ── 6. Calculate proportions for each annotation column ─────────────────────
  for (annotation_col in annotation_cols) {

    raw_labels <- cd[[annotation_col]]

    n_na <- sum(is.na(raw_labels))
    if (n_na > 0) {
      warning(sprintf(
        "%d cell(s) have NA in annotation column '%s' and will be excluded from proportions.",
        n_na, annotation_col
      ))
    }

    valid_mask      <- !is.na(raw_labels)
    all_categories  <- sort(unique(raw_labels[valid_mask]))

    if (length(all_categories) == 0) {
      warning(sprintf(
        "No valid (non-NA) values found in annotation column '%s'. Skipping.",
        annotation_col
      ))
      next
    }

    prop_col_names <- paste0(annotation_col, "_", all_categories, "_prop")

    # Initialise proportion matrix for this annotation column
    prop_mat <- matrix(
      NA_real_,
      nrow     = n_nhoods,
      ncol     = length(all_categories),
      dimnames = list(NULL, prop_col_names)
    )

    for (i in seq_len(n_nhoods)) {
      cell_idx       <- nhood_cell_idx[[i]]
      valid_in_nhood <- cell_idx[valid_mask[cell_idx]]

      if (length(valid_in_nhood) == 0) {
        next  # leave row as NA
      }

      counts          <- table(factor(raw_labels[valid_in_nhood], levels = all_categories))
      prop_mat[i, ]   <- as.numeric(counts) / sum(counts)
    }

    out_df <- cbind(out_df, as.data.frame(prop_mat, stringsAsFactors = FALSE))
  }

  # ── 7. Store results and return ─────────────────────────────────────────────
  if (result_name %in% names(S4Vectors::metadata(object))) {
    warning(sprintf(
      "metadata(object)[['%s']] already exists and will be overwritten.",
      result_name
    ))
  }

  message(sprintf(
    "Stored annotation proportions (%d neighbourhood(s) x %d annotation column(s)) in metadata(object)[['%s']].",
    n_nhoods, length(annotation_cols), result_name
  ))

  S4Vectors::metadata(object)[[result_name]] <- out_df

  return(object)
}
