#' Assign Dominant Annotation to Neighbourhoods
#'
#' Computes per-neighbourhood the proportion of cells belonging to each
#' category of one or more annotation columns and assigns a dominant identity
#' based on a user-defined threshold. Neighbourhoods where no single category
#' exceeds the threshold are labelled \code{"mixed"}.
#'
#' @param object A \code{\link[miloR]{Milo}} or
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}} object. Seurat
#'   objects must first be converted via
#'   \code{\link[Seurat]{as.SingleCellExperiment}}.
#' @param nhood_col Character string. For plain \code{SingleCellExperiment}
#'   objects (i.e. not \code{Milo}), the name of the \code{colData} column
#'   containing neighbourhood assignments. Not required for \code{Milo}
#'   objects.
#' @param annotation_cols Character vector of \code{colData} column names containing
#'   per-cell annotation categories to summarise per neighbourhood (e.g.
#'   cell type, lineage, cluster label). All specified columns must be present
#'   in \code{colData}; missing columns raise an error. Default:
#'   \code{"celltype"}.
#' @param threshold Numeric value between 0-1. The minimum proportion a single
#'   annotation category must reach for a neighbourhood to be assigned that
#'   label. Neighbourhoods where the dominant category falls below this
#'   threshold are assigned \code{"mixed"}. Default: \code{0.7}.
#' @param result_name Character string. Key under which the results
#'   \code{data.frame} is stored in \code{metadata(object)}. If an entry with
#'   this name already exists a warning is emitted and the entry is
#'   overwritten. Default: \code{"nhood_annotation"}.
#'
#' @return The input \code{object}, with a \code{data.frame} stored in
#'   \code{metadata(object)[[result_name]]}. The \code{data.frame} has one
#'   row per neighbourhood and the following columns:
#'   \describe{
#'     \item{\code{nhood}}{Neighbourhood identifier.}
#'     \item{\code{nhood_size}}{Total number of cells in the neighbourhood.}
#'     \item{\code{<annotation_col>}}{Dominant annotation
#'       label for this neighbourhood, or \code{"mixed"} if no single category
#'       reached \code{threshold}. One column per entry in \code{annotation_cols}.}
#'     \item{\code{<annotation_col>_dominant_prop}}{Proportion of cells
#'       carrying the dominant annotation label. One column per entry in
#'       \code{annotation_cols}.}
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
#' ## Dominant label assignment
#' Each neighbourhood represents a group of transcriptionally similar cells.
#' For each neighbourhood and each annotation column, per-cell labels are
#' tabulated and converted to proportions. The category with the highest
#' proportion is the candidate dominant label. If that proportion is >=
#' \code{threshold} the neighbourhood is assigned that label; otherwise it is
#' labelled \code{"mixed"}, indicating that no single annotation category
#' predominates. When multiple annotation columns are supplied the procedure
#' is applied independently for each.
#'
#' @seealso
#' \code{\link{nhood_annotation_props}}, \code{\link{nhood_diversity}},
#' \code{\link{nhood_shm}}, \code{\link{nhood_cdr3_convergence}}
#'
#' @examples
#' \dontrun{
#' # --- Milo object, single annotation column --------------------------------
#' milo <- nhood_annotation(
#'   object          = milo_object,
#'   annotation_cols = "celltype",
#'   threshold       = 0.7
#' )
#' results <- S4Vectors::metadata(milo)[["nhood_annotation"]]
#'
#' # --- Multiple annotation columns, relaxed threshold ----------------------
#' milo <- nhood_annotation(
#'   object          = milo_object,
#'   annotation_cols = c("celltype", "lineage"),
#'   threshold       = 0.6,
#'   result_name     = "nhood_annotation"
#' )
#'
#' # --- Plain SingleCellExperiment object ------------------------------------
#' sce <- nhood_annotation(
#'   object          = sce_object,
#'   nhood_col       = "community_id",
#'   annotation_cols = "celltype"
#' )
#' }
#'
#' @export
#'
nhood_annotation <- function(object,
                             nhood_col = NULL,
                             annotation_cols = "celltype",
                             threshold = 0.7,
                             result_name = "nhood_annotation") {

  # ── 1. Determine object type ────────────────────────────────────────────────
  is_milo     <- inherits(object, "Milo") && inherits(object, "SingleCellExperiment")
  is_sce_only <- !is_milo && inherits(object, "SingleCellExperiment")

  if (!is_milo && !is_sce_only) {
    stop("'object' must be a Milo or SingleCellExperiment object. If using Seurat, convert via as.SingleCellExperiment")
  }

  # ── 2. Extract colData ───────────────────────────────────────────────────────
  col_data <- SummarizedExperiment::colData(object)

  # ── 3. Validate annotation columns ──────────────────────────────────────────
  missing_cols <- annotation_cols[!annotation_cols %in% colnames(col_data)]
  if (length(missing_cols) > 0) {
    stop("Annotation column(s) not found in colData: ", paste(missing_cols, collapse = ", "))
  }

  # ── 4. Build neighbourhood matrix ───────────────────────────────────────────
  if (is_milo) {
    if (!requireNamespace("miloR", quietly = TRUE)) {
      stop("Package 'miloR' is required for Milo objects. Please install it with: BiocManager::install('miloR')")
    }
    if (is.null(miloR::nhoods(object))) {
      stop("nhoods slot not found in Milo object")
    }
    nhoods <- miloR::nhoods(object)
    unique_nhoods <- colnames(nhoods)
    n_nhoods <- ncol(nhoods)

    # Confirm row count of nhood matrix matches number of cells
    if (nrow(nhoods) != ncol(object)) {
      stop(sprintf(
        "nhoods matrix has %d rows but the object has %d cells. These must match.",
        nrow(nhoods), ncol(object)
      ))
    }

    # Confirm row order: rownames of nhood matrix should match colnames of object
    if (!is.null(rownames(nhoods)) && !is.null(colnames(object))) {
      if (!identical(rownames(nhoods), colnames(object))) {
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

  } else if (is_sce_only) {
    if (is.null(nhood_col)) {
      stop("For SingleCellExperiment objects, 'nhood_col' must be specified")
    }
    if (!nhood_col %in% colnames(col_data)) {
      stop("Neighbourhood column '", nhood_col, "' not found in colData")
    }

    nhood_assignments <- col_data[[nhood_col]]

    if (any(is.na(nhood_assignments))) {
      warning(sprintf("%d cells have NA neighbourhood assignments and will be excluded",
                      sum(is.na(nhood_assignments))))
      valid_cells <- !is.na(nhood_assignments)
    } else {
      valid_cells <- rep(TRUE, length(nhood_assignments))
    }

    unique_nhoods <- unique(nhood_assignments[valid_cells])
    unique_nhoods <- as.character(sort(as.numeric(unique_nhoods)))
    n_nhoods <- length(unique_nhoods)

    nhoods <- matrix(0,
                     nrow = ncol(object),
                     ncol = n_nhoods,
                     dimnames = list(colnames(object),
                                     paste0("nhood_", unique_nhoods)))
    for (i in seq_along(unique_nhoods)) {
      cells_in_nhood <- nhood_assignments == unique_nhoods[i] & valid_cells
      nhoods[cells_in_nhood, i] <- 1
    }
    if (requireNamespace("Matrix", quietly = TRUE)) {
      nhoods <- Matrix::Matrix(nhoods, sparse = TRUE)
    }
  }

  # ── 5. Initialise results dataframe ────────────────────────────────────────
  annotation_results <- data.frame(
    nhood      = unique_nhoods,
    nhood_size = integer(n_nhoods),
    stringsAsFactors = FALSE
  )

  # Pre-compute cell IDs per neighbourhood (used across all annotation_cols)
  nhood_cell_ids <- lapply(seq_len(n_nhoods), function(i) {
    rownames(nhoods)[which(nhoods[, i] > 0)]
  })

  # Fill nhood_size once
  annotation_results[["nhood_size"]] <- vapply(nhood_cell_ids, length, integer(1))

  # ── 6. Calculate annotation proportions for each neighbourhood ──────────────
  for (annotation_col in annotation_cols) {
    dom_prop_col <- paste0(annotation_col, "_dominant_prop")

    annotation_results[[annotation_col]] <- character(n_nhoods)
    annotation_results[[dom_prop_col]]   <- numeric(n_nhoods)

    for (i in seq_len(n_nhoods)) {
      nhood_annotation_types <- col_data[[annotation_col]][rownames(col_data) %in% nhood_cell_ids[[i]]]

      annotation_type_props  <- table(nhood_annotation_types) / length(nhood_annotation_types)
      dominant_annotation    <- names(which.max(annotation_type_props))
      dominant_prop          <- max(annotation_type_props)

      annotation_results[i, dom_prop_col]      <- dominant_prop
      annotation_results[i, annotation_col]    <- if (dominant_prop >= threshold) {
        dominant_annotation
      } else {
        "mixed"
      }
    }
  }

  # ── 7. Summary ──────────────────────────────────────────────────────────────
  for (annotation_col in annotation_cols) {
    dom_prop_col <- paste0(annotation_col, "_dominant_prop")
    assignments  <- annotation_results[[annotation_col]]
    n_mixed      <- sum(assignments == "mixed")
    n_assigned   <- n_nhoods - n_mixed

    message(sprintf("Annotation column '%s':", annotation_col))
    message(sprintf(
      "  %d / %d neighbourhoods assigned a dominant identity (threshold = %.2f).",
      n_assigned, n_nhoods, threshold
    ))
    message(sprintf(
      "  %d / %d neighbourhoods labelled 'mixed' (no identity exceeded threshold = %.2f).",
      n_mixed, n_nhoods, threshold
    ))

    if (n_assigned > 0) {
      identity_counts <- sort(table(assignments[assignments != "mixed"]), decreasing = TRUE)
      message("  Dominant identity breakdown:")
      for (nm in names(identity_counts)) {
        message(sprintf("    %s: %d neighbourhood(s)", nm, identity_counts[[nm]]))
      }
    }
  }

  # ── 8. Store results in metadata and return ─────────────────────────────────
  if (result_name %in% names(S4Vectors::metadata(object))) {
    warning(sprintf(
      "metadata(object)[['%s']] already exists and will be overwritten.",
      result_name
    ))
  }

  S4Vectors::metadata(object)[[result_name]] <- annotation_results

  message(sprintf(
    "Stored annotation results (%d neighbourhoods x %d annotation column(s)) in metadata(object)[['%s']].",
    n_nhoods, length(annotation_cols), result_name
  ))

  return(object)
}
