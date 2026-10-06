#' Calculate Neighbourhood Somatic Hypermutation Statistics
#'
#' Computes per-neighbourhood somatic hypermutation (SHM) summary statistics
#' from BCR (B cell receptor) data stored in \code{colData}. For each
#' neighbourhood, the function aggregates SHM metrics across cells that carry
#' BCR information, and flags neighbourhoods with insufficient BCR coverage as
#' unreliable.
#'
#' @param object A \code{\link[miloR]{Milo}} or
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}} object. Seurat
#'   objects must first be converted via
#'   \code{\link[Seurat]{as.SingleCellExperiment}}.
#' @param nhood_col Character string naming a \code{colData} column that holds
#'   pre-computed neighbourhood assignments (as a numeric or character vector).
#'   Required when \code{object} is a plain \code{SingleCellExperiment};
#'   ignored for \code{Milo} objects, which carry neighbourhood information in
#'   the \code{nhoods} slot.
#' @param shm_cols Character vector of column names in \code{colData(object)}
#'   containing per-cell SHM metrics to aggregate. Defaults to
#'   \code{c("mu_count_VDJ", "mu_count_VJ", "mu_freq_VDJ", "mu_freq_VJ",
#'   "mu_count", "mu_freq")}. All specified columns must be present in
#'   \code{colData}; missing columns will raise an error.
#' @param min_n_shm Integer. Minimum number of cells with SHM annotation required for a
#'   neighbourhood's SHM statistics to be considered reliable. Neighbourhoods
#'   with fewer than \code{min_n_shm} contigs are retained in the
#'   output but have \code{shm_reliable = FALSE} and \code{NA} values for all
#'   SHM columns. Default is \code{10}.
#' @param summary_metric Character string specifying the central tendency
#'   statistic to compute for each SHM metric within each neighbourhood. One
#'   of \code{"mean"} (default) or \code{"median"}. The corresponding spread
#'   metric is computed automatically: standard deviation (\code{sd}) for
#'   \code{"mean"}, interquartile range (\code{IQR}) for \code{"median"}.
#' @param result_name Character string. Name of the metadata slot in which to
#'   store the resulting SHM data frame, i.e. the results are accessible via
#'   \code{metadata(object)[[result_name]]}. If the slot already exists it will
#'   be overwritten with a warning. Default is \code{"nhood_shm"}.
#'
#' @return The input \code{object} with an additional (or updated) entry in
#'   \code{metadata(object)[[result_name]]}. The stored value is a
#'   \code{data.frame} with one row per neighbourhood and the following columns:
#'   \describe{
#'     \item{\code{nhood}}{Neighbourhood identifier.}
#'     \item{\code{nhood_size}}{Total number of cells in the neighbourhood.}
#'     \item{\code{n_shm}}{Number of cells with at least one non-\code{NA} SHM
#'       value.}
#'     \item{\code{shm_capture_rate}}{Proportion of neighbourhood cells with
#'       SHM data (\code{n_shm / nhood_size}).}
#'     \item{\code{shm_reliable}}{Logical; \code{TRUE} if
#'       \code{n_shm >= min_n_shm}.}
#'     \item{\code{<shm_col>_<summary_metric>}}{One column per entry in
#'       \code{shm_cols} containing the aggregated statistic (mean or median).}
#'     \item{\code{<shm_col>_<spread_metric>}}{One column per entry in
#'       \code{shm_cols} containing the spread statistic (\code{sd} or
#'       \code{IQR}).}
#'   }
#'   The data frame also carries the following attributes:
#'   \code{min_n_shm}, \code{shm_cols}, \code{summary_metric}, and
#'   \code{spread_metric}.
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
#' ## SHM reliability flag
#' The \code{shm_reliable} flag is intended to help downstream filtering.
#' Neighbourhoods below \code{min_n_shm} will have \code{NA} for all SHM
#' statistics and should generally be excluded from differential SHM analyses.
#'
#' @seealso
#' \code{\link{nhood_diversity}}, \code{\link{nhood_annotation}},
#' \code{\link{nhood_annotation_props}}, \code{\link{nhood_cdr3_convergence}}
#'
#' @examples
#' \dontrun{
#' # --- Milo object ---
#' milo_obj <- nhood_shm(
#'   object = milo_obj,
#'   shm_cols = c("mu_freq_VDJ", "mu_freq_VJ", "mu_freq"),
#'   min_n_shm = 15,
#'   summary_metric = "median",
#'   result_name = "nhood_shm"
#' )
#'
#' # Retrieve results
#' shm_df <- metadata(milo_obj)[["nhood_shm"]]
#' head(shm_df)
#'
#' # --- Plain SingleCellExperiment object ---
#' sce_obj <- nhood_shm(
#'   object = sce_obj,
#'   nhood_col = "neighbourhood_id",
#'   min_n_shm = 10,
#'   summary_metric = "mean"
#' )
#' }
#'
#' @export
#'
nhood_shm <- function(object,
                      nhood_col = NULL,
                      shm_cols = c('mu_count_VDJ', 'mu_count_VJ', 'mu_freq_VDJ', 'mu_freq_VJ', 'mu_count', 'mu_freq'),
                      min_n_shm = 20,
                      summary_metric = "mean",
                      result_name = 'nhood_shm') {

  # ── 0. Validate summary_metric ────────────────────────────────────────
  if (!is.character(summary_metric) || length(summary_metric) != 1L ||
      is.na(summary_metric) || !summary_metric %in% c("mean", "median")) {
    stop(sprintf(
      "'summary_metric' must be either 'mean' or 'median'. Received: %s",
      paste(format(summary_metric), collapse = ", ")
    ), call. = FALSE)
  }

  # ── 1. Determine object type ────────────────────────────────────────────────
  is_milo    <- inherits(object, "Milo") && inherits(object, "SingleCellExperiment")
  is_sce_only <- !is_milo && inherits(object, "SingleCellExperiment")

  if (!is_milo && !is_sce_only) {
    stop("'object' must be a Milo or SingleCellExperiment object. If using Seurat, convert via as.SingleCellExperiment")
  }

  # ── 2. Extract colData and neighbourhood matrix ─────────────────────────────
  col_data <- SummarizedExperiment::colData(object)

  if (is_milo) {
    if (!requireNamespace("miloR", quietly = TRUE)) {
      stop("Package 'miloR' is required for Milo objects. Please install it with: BiocManager::install('miloR')")
    }
    if (is.null(miloR::nhoods(object))) {
      stop("nhoods slot not found in Milo object")
    }
    nhoods <- miloR::nhoods(object)
    unique_nhoods <- colnames(nhoods)

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
          "Cell ordering is inconsistent; cannot safely assign SHM values."
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
      stop(sprintf("Column '%s' not found in colData. Available columns: %s",
                   nhood_col,
                   paste(colnames(col_data), collapse = ", ")))
    }

    nhood_assignments <- col_data[[nhood_col]]

    if (any(is.na(nhood_assignments))) {
      n_missing <- sum(is.na(nhood_assignments))
      warning(sprintf("%d cells have NA neighbourhood assignments and will be excluded from analysis",
                      n_missing))
      valid_cells <- !is.na(nhood_assignments)
    } else {
      valid_cells <- rep(TRUE, length(nhood_assignments))
    }

    unique_nhoods <- unique(nhood_assignments[valid_cells])
    unique_nhoods <- as.character(sort(as.numeric(unique_nhoods)))
    n_nhoods_total <- length(unique_nhoods)

    nhoods <- matrix(0,
                     nrow = ncol(object),
                     ncol = n_nhoods_total,
                     dimnames = list(colnames(object),
                                     paste0("nhood_", unique_nhoods)))

    for (i in seq_along(unique_nhoods)) {
      nhood_id <- unique_nhoods[i]
      cells_in_nhood <- nhood_assignments == nhood_id & valid_cells
      nhoods[cells_in_nhood, i] <- 1
    }

    if (requireNamespace("Matrix", quietly = TRUE)) {
      nhoods <- Matrix::Matrix(nhoods, sparse = TRUE)
    }
  }

  n_nhoods <- ncol(nhoods)

  # ── 3. Validate SHM columns ─────────────────────────────────────────────────
  missing_cols <- setdiff(shm_cols, colnames(col_data))
  if (length(missing_cols) > 0) {
    stop("Missing SHM columns in colData: ", paste(missing_cols, collapse = ", "))
  }

  # Determine spread metric based on aggregation method
  spread_metric <- if (summary_metric == "mean") "sd" else "IQR"

  # ── 4. Initialise results dataframe ────────────────────────────────────────
  shm_results <- data.frame(
    nhood = unique_nhoods,
    nhood_size = integer(n_nhoods),
    n_shm = integer(n_nhoods),
    shm_capture_rate = numeric(n_nhoods),
    shm_reliable = logical(n_nhoods),
    stringsAsFactors = FALSE
  )

  # Pre-allocate SHM statistic columns (both aggregation and spread)
  for (shm_col in shm_cols) {
    # Aggregation column (mean or median)
    col_name_agg <- paste0(shm_col, "_", summary_metric)
    shm_results[[col_name_agg]] <- numeric(n_nhoods)

    # Spread column (sd or IQR)
    col_name_spread <- paste0(shm_col, "_", spread_metric)
    shm_results[[col_name_spread]] <- numeric(n_nhoods)
  }

  # ── 5. Calculate SHM statistics for each neighbourhood ──────────────────────
  for (i in seq_len(n_nhoods)) {
    # Get cells in neighbourhood
    cells_in_nhood <- which(nhoods[, i] == 1)
    cell_ids <- rownames(nhoods)[cells_in_nhood]

    # Filter for cells with BCR data
    nhood_col_data <- col_data[rownames(col_data) %in% cell_ids, ]

    # Filter out rows where all SHM columns are NA
    valid_cells <- rowSums(!is.na(nhood_col_data[, shm_cols, drop = FALSE])) > 0
    valid_nhood_col_data <- nhood_col_data[valid_cells, ]

    n_shm <- nrow(valid_nhood_col_data)

    # Store basic metrics
    shm_results[i, "nhood_size"] <- length(cells_in_nhood)
    shm_results[i, "n_shm"] <- n_shm
    shm_results[i, "shm_capture_rate"] <- n_shm / length(cells_in_nhood)

    if (n_shm < min_n_shm) {
      # Insufficient data for reliable SHM calculation
      shm_results[i, "shm_reliable"] <- FALSE

      for (shm_col in shm_cols) {
        col_name_agg <- paste0(shm_col, "_", summary_metric)
        col_name_spread <- paste0(shm_col, "_", spread_metric)

        shm_results[i, col_name_agg] <- NA
        shm_results[i, col_name_spread] <- NA
      }

    } else {
      # Sufficient data - calculate SHM statistics
      shm_results[i, "shm_reliable"] <- TRUE

      # Calculate aggregated SHM statistics and spread metrics
      for (shm_col in shm_cols) {
        col_name_agg <- paste0(shm_col, "_", summary_metric)
        col_name_spread <- paste0(shm_col, "_", spread_metric)
        shm_data <- valid_nhood_col_data[[shm_col]]

        if (length(shm_data) > 0) {
          # Calculate aggregation statistic
          if (summary_metric == "mean") {
            shm_results[i, col_name_agg] <- mean(shm_data, na.rm = TRUE)
            # Calculate standard deviation
            shm_results[i, col_name_spread] <- stats::sd(shm_data, na.rm = TRUE)
          } else if (summary_metric == "median") {
            shm_results[i, col_name_agg] <- stats::median(shm_data, na.rm = TRUE)
            # Calculate IQR
            shm_results[i, col_name_spread] <- stats::IQR(shm_data, na.rm = TRUE)
          }
        } else {
          shm_results[i, col_name_agg] <- NA
          shm_results[i, col_name_spread] <- NA
        }
      }
    }
  }

  # ── 6. Summary ──────────────────────────────────────────────────────────────
  total_cells <- ncol(object)
  total_shm   <- sum(rowSums(!is.na(as.data.frame(col_data)[, shm_cols, drop = FALSE])) > 0)
  message("SHM analysis summary:")
  message("Cells with SHM data: ", total_shm, " / ", total_cells)
  message("Neighbourhoods with >= ", min_n_shm, " cells with SHM: ",
          sum(shm_results$n_shm >= min_n_shm), " / ", n_nhoods)
  message("Spread metric used: ", spread_metric)

  # ── 7. Add attributes ───────────────────────────────────────────────────────
  attr(shm_results, "min_n_shm") <- min_n_shm
  attr(shm_results, "shm_cols") <- shm_cols
  attr(shm_results, "summary_metric") <- summary_metric
  attr(shm_results, "spread_metric") <- spread_metric

  # ── 8. Store results in metadata and return ─────────────────────────────────
  if (result_name %in% names(S4Vectors::metadata(object))) {
    warning(sprintf(
      "metadata(object)[['%s']] already exists and will be overwritten.",
      result_name
    ))
  }

  S4Vectors::metadata(object)[[result_name]] <- shm_results

  message(sprintf(
    "Stored SHM results (%d neighbourhoods x %d SHM columns) in metadata(object)[['%s']].",
    n_nhoods, length(shm_cols), result_name
  ))

  return(object)
}
