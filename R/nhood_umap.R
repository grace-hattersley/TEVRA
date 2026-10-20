#' Compute Neighbourhood UMAP Centroids
#'
#' Aggregates per-cell UMAP coordinates into per-neighbourhood centroids
#'
#' @param object A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   object with a UMAP stored in \code{reducedDims}. Seurat objects must be
#'   converted via \code{\link[Seurat]{as.SingleCellExperiment}}
#'   before use.
#' @param nhood_col Character string. Name of the \code{colData} column
#'   containing neighbourhood assignments. Default: \code{"community_id"}.
#' @param aggregate_fun Character string. Aggregation function used to compute
#'   neighbourhood centroids. Must be one of \code{"mean"} (arithmetic mean;
#'   default) or \code{"median"}. The choice affects both the computed
#'   coordinates and the names of the output columns (see Value).
#' @param result_name Character string. Key under which the results
#'   \code{data.frame} is stored in \code{metadata(object)}. If an entry with
#'   this name already exists a warning is emitted and the entry is
#'   overwritten. Default: \code{"nhood_umap"}.
#'
#' @return The input \code{object}, with a \code{data.frame} of neighbourhood
#'   UMAP centroids stored in \code{metadata(object)[[result_name]]}. The
#'   \code{data.frame} contains the following columns:
#'   \describe{
#'     \item{\code{nhood}}{Neighbourhood identifier.}
#'     \item{\code{nhood_size}}{Total number of cells in the neighbourhood.}
#'     \item{\code{mean_UMAP_1} / \code{median_UMAP_1}}{UMAP dimension 1
#'       centroid, named according to \code{aggregate_fun}.}
#'     \item{\code{mean_UMAP_2} / \code{median_UMAP_2}}{UMAP dimension 2
#'       centroid, named according to \code{aggregate_fun}.}
#'   }
#'
#' @details
#' The UMAP is accessed by name from \code{reducedDims(object)$UMAP} and is
#' assumed to have exactly two dimensions, which are renamed \code{UMAP_1}
#' and \code{UMAP_2} internally. If your UMAP is stored under a different
#' name in \code{reducedDims}, rename it before calling this function.
#' Centroids are computed using the function specified by \code{aggregate_fun}
#' applied to all cell coordinates within each neighbourhood.
#' \code{"mean"} produces the arithmetic centroid; \code{"median"} produces
#' the geometric median and is more robust to outlier cells at neighbourhood
#' boundaries.
#' If \code{result_name} already exists in \code{metadata(object)} a warning
#' is emitted and the entry is overwritten.
#'
#' @seealso
#' \code{\link{plot_nhood_umap}} to visualise the neighbourhood umap
#'
#' @examples
#' \dontrun{
#' # --- Mean centroids (default) ---------------------------------------------
#' sce <- nhood_umap(
#'   object    = sce_object,
#'   nhood_col = "community_id"
#' )
#' umap_df <- S4Vectors::metadata(sce)[["nhood_umap"]]
#'
#' # --- Median centroids -----------------------------------------------------
#' sce <- nhood_umap(
#'   object       = sce_object,
#'   nhood_col    = "community_id",
#'   aggregate_fun = "median"
#' )
#' umap_df <- S4Vectors::metadata(sce)[["nhood_umap"]]
#' }
#'
#' @export
#'
nhood_umap <- function(object,
                       nhood_col     = "community_id",
                       aggregate_fun = "mean",
                       result_name   = "nhood_umap") {

  # ── 1. Input validation ────────────────────────────────────────────────────
  if (!inherits(object, "SingleCellExperiment")) {
    stop("'object' must be a SingleCellExperiment object. If using Seurat, convert via as.SingleCellExperiment")
  }

  if (is.null(SingleCellExperiment::reducedDims(object)$UMAP)) {
    stop("No UMAP found in reducedDims(object). Run UMAP reduction first.")
  }

  if (!nhood_col %in% colnames(SummarizedExperiment::colData(object))) {
    stop(sprintf("Column '%s' not found in colData. Available columns: %s",
                 nhood_col,
                 paste(colnames(SummarizedExperiment::colData(object)), collapse = ", ")))
  }

  aggregate_fun <- match.arg(aggregate_fun, choices = c("mean", "median"))

  # ── 2. Extract UMAP coordinates and cell metadata ──────────────────────────
  umap_coords <- as.data.frame(SingleCellExperiment::reducedDims(object)$UMAP)
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  umap_coords_cells <- merge(
    umap_coords,
    as.data.frame(SummarizedExperiment::colData(object)),
    by = "row.names"
  )

  # ── 3. Aggregate to neighbourhood centroids ────────────────────────────────
  agg_fn <- switch(aggregate_fun, mean = mean, median = stats::median)

  neighbourhood_centroids <- stats::aggregate(
    cbind(UMAP_1, UMAP_2) ~ .,
    data = umap_coords_cells[, c("UMAP_1", "UMAP_2", nhood_col)],
    FUN  = agg_fn
  )
  colnames(neighbourhood_centroids)[2:3] <- paste0(aggregate_fun, c("_UMAP_1", "_UMAP_2"))

  nhood_size <- stats::aggregate(
    UMAP_1 ~ .,
    data = umap_coords_cells[, c("UMAP_1", nhood_col)],
    FUN  = length
  )
  colnames(nhood_size)[2] <- "nhood_size"

  neighbourhood_centroids <- merge(neighbourhood_centroids, nhood_size, by = nhood_col)
  colnames(neighbourhood_centroids)[colnames(neighbourhood_centroids) == nhood_col] <- "nhood"
  neighbourhood_centroids <- neighbourhood_centroids[, c("nhood", "nhood_size",
                                                         paste0(aggregate_fun, "_UMAP_1"),
                                                         paste0(aggregate_fun, "_UMAP_2"))]

  # ── 4. Store result in metadata and return object ──────────────────────────
  if (result_name %in% names(S4Vectors::metadata(object))) {
    warning(sprintf(
      "metadata(object)[['%s']] already exists and will be overwritten.",
      result_name
    ))
  }

  S4Vectors::metadata(object)[[result_name]] <- neighbourhood_centroids

  message(sprintf(
    "Neighbourhood UMAP centroids (%s) for %d neighbourhoods stored in metadata(object)[['%s']].",
    aggregate_fun,
    nrow(neighbourhood_centroids),
    result_name
  ))

  return(object)
}
