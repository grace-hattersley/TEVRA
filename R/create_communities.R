#' Create cell communities using Leiden community detection
#'
#' Assigns cells in a \code{SingleCellExperiment} object to communities
#' using the Leiden algorithm at a single chosen resolution. Community
#' assignments are added to \code{colData(sce)} and a per-community
#' summary table is returned alongside the updated object. Use
#' \code{\link{test_resolution}} first to identify a suitable resolution.
#'
#' @param sce A \code{SingleCellExperiment} object. The object must contain a
#'   pre-computed nearest-neighbour graph stored in
#'   \code{metadata(sce)$graphs} under the name given by \code{graph.name}.
#'   See the \strong{Graph preparation} section below for how to build and
#'   store the graph
#' @param graph.name A character string giving the name of the graph to
#'   retrieve from \code{metadata(sce)$graphs}. Defaults to \code{"RNA_snn"}.
#' @param resolution A single numeric value controlling the granularity of
#'   Leiden clustering. Higher values produce more, smaller communities.
#'   Use \code{\link{test_resolution}} to explore candidate values before
#'   committing to one. Defaults to \code{0.5}.
#' @param min_cells A positive integer. Communities with fewer than this
#'   many cells are considered too small; their cells are reassigned to the
#'   nearest valid community. If all communities fall below this
#'   threshold the function stops with an informative error. Defaults to
#'   \code{5}.
#' @param feature.col A character string giving the name of the
#'   \code{colData(sce)} column holding the number of detected features
#'   (genes) per cell. Defaults to \code{"nFeature_RNA"}, the Seurat
#'   convention. Set this if your object uses a different column name (e.g.
#'   \code{"detected"}, as produced by \code{scuttle::perCellQCMetrics}).
#'   If the column is not found, \code{mean_detected} is returned as
#'   \code{NA} with a warning.
#' @param count.col A character string giving the name of the
#'   \code{colData(sce)} column holding total counts per cell. Defaults to
#'   \code{"nCount_RNA"}, the Seurat convention. Set this if your object
#'   uses a different column name (e.g. \code{"sum"}, as produced by
#'   \code{scuttle::perCellQCMetrics}). If the column is not found,
#'   \code{mean_total_counts} is returned as \code{NA} with a warning.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{sce}}{The input \code{SingleCellExperiment} object with a
#'       new \code{community_id} column added to \code{colData}.}
#'     \item{\code{communities_summary}}{A \code{data.frame} with one row
#'       per community containing:
#'       \describe{
#'         \item{\code{community_id}}{Community identifier.}
#'         \item{\code{n_cells}}{Number of cells in the community.}
#'         \item{\code{mean_detected}}{Mean number of detected features per
#'           cell, taken from the \code{feature.col} column of
#'           \code{colData}.}
#'         \item{\code{mean_total_counts}}{Mean total counts per cell, taken
#'           from the \code{count.col} column of \code{colData}.}
#'       }
#'     }
#'   }
#'
#' @details
#' ## Graph preparation:
#' This function does not build the nearest-neighbour graph itself.
#' The graph must be computed and stored in \code{metadata(sce)$graphs} before
#' calling this function. Two workflows are supported:
#'
#' \strong{1. From scratch with scran:}
#' \preformatted{
#' library(scran)
#'
#' # Build a shared nearest-neighbour graph on the PCA embedding
#' g <- scran::buildSNNGraph(sce, use.dimred = "PCA", k = 20)
#'
#' # Store it under the name this function will look for
#' metadata(sce)$graphs[["RNA_snn"]] <- g
#' }
#'
#' \strong{2. From an existing Seurat object:}
#' \preformatted{
#' library(Seurat)
#'
#' # Convert to SCE first
#' sce <- as.SingleCellExperiment(seurat_obj)
#'
#' # Extract the adjacency matrix from the Seurat object, convert to igraph,
#' # and store in the SCE metadata
#' g <- igraph::graph_from_adjacency_matrix(
#'   seurat_obj@graphs[["RNA_snn"]],
#'   mode     = "max",
#'   weighted = TRUE
#' )
#' metadata(sce)$graphs[["RNA_snn"]] <- g
#' }
#'
#' The value passed to \code{graph.name} must match the key used when storing
#' the graph in \code{metadata(sce)$graphs}.
#'
#' @seealso \code{\link{test_resolution}} to evaluate resolution choices before
#' running this function
#'
#' @examples
#' \dontrun{
#' library(scran)
#'
#' # Build and store the graph
#' g <- scran::buildSNNGraph(sce, use.dimred = "PCA", k = 20)
#' metadata(sce)$graphs[["RNA_snn"]] <- g
#'
#' # Optionally explore resolutions first
#' test_resolution(sce, resolutions = c(0.5, 1, 2.5, 5))
#'
#' # Create communities at the chosen resolution
#' out <- create_communities(sce, resolution = 1)
#' sce <- out$sce
#' summary_df <- out$communities_summary
#'
#' # If sce was built with scuttle::perCellQCMetrics() rather than Seurat,
#' # specify the matching column names
#' out <- create_communities(sce, resolution = 1,
#'                            feature.col = "detected",
#'                            count.col = "sum")
#' }
#'
#' @importFrom rlang .data
#'
#' @export
#'
create_communities <- function(sce,
                               graph.name = "RNA_snn",
                               resolution = 0.5,
                               min_cells = 5,
                               feature.col = "nFeature_RNA",
                               count.col = "nCount_RNA") {

  # ── 1. Determine object type ───────────────────────────────────────────────
  if (inherits(sce, "Seurat")) {
    stop(
      "Input is a Seurat object. This function expects a SingleCellExperiment object.\n",
      "Please convert first with: sce <- as.SingleCellExperiment(seurat_obj)"
    )
  }
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("Input must be a SingleCellExperiment object.")
  }

  # ── 2. Extract kNN graph ───────────────────────────────────────────────────
  # Check there is a kNN graph stored in the SCE metadata
  graphs <- S4Vectors::metadata(sce)[["graphs"]]
  if (is.null(graphs) || is.null(graphs[[graph.name]])) {
    available_graphs <- names(graphs)
    if (length(available_graphs) > 0) {
      stop(paste0(
        "Graph '", graph.name, "' not found in metadata(sce)$graphs. ",
        "Available graphs: ", paste(available_graphs, collapse = ", "),
        ". Run scran::buildSNNGraph() and store the result in metadata(sce)$graphs, ",
        "or specify the correct graph.name."
      ))
    } else {
      stop(
        "No graphs found in metadata(sce)$graphs. ",
        "Run scran::buildSNNGraph() and store the result in metadata(sce)$graphs first."
      )
    }
  }

  # Extract the kNN graph (expected to be a sparse matrix or igraph object)
  knn_graph <- graphs[[graph.name]]

  # Convert to igraph if not already
  if (inherits(knn_graph, "igraph")) {
    igraph_obj <- knn_graph
    knn_matrix <- igraph::as_adjacency_matrix(igraph_obj, attr = "weight", sparse = TRUE)
  } else {
    knn_matrix <- knn_graph
    igraph_obj <- igraph::graph_from_adjacency_matrix(knn_matrix,
                                                      mode = "max",
                                                      weighted = TRUE)
  }

  # ── 4. Compute communities ─────────────────────────────────────────────────
  # Run leiden
  communities <- igraph::cluster_leiden(igraph_obj,
                                        resolution = resolution,
                                        objective_function = "modularity")

  community_id <- communities$membership
  names(community_id) <- rownames(knn_matrix)

  # Filter out small communities
  community_counts <- table(community_id)
  valid_communities <- names(community_counts)[community_counts >= min_cells]

  # Check if any valid communities exist
  if (length(valid_communities) == 0) {
    stop(paste0(
      "All communities have fewer than ", min_cells,
      " cells. Try lowering min_cells or run test_resolution() to select an appropriate resolution."
    ))
  }

  # Reassign cells in small communities
  cells_to_reassign <- names(community_id)[!community_id %in% valid_communities]
  for (cell in cells_to_reassign) {
    # Find nearest community of this cell
    nearest <- which(knn_matrix[cell, ] > 0)
    nearest_communities <- community_id[names(nearest)]
    # Assign to most common valid nearest community
    valid_nearest_communities <- nearest_communities[nearest_communities %in% valid_communities]
    if (length(valid_nearest_communities) > 0) {
      community_id[cell] <- names(sort(table(valid_nearest_communities), decreasing = TRUE))[1]
    }
  }

  # ── 5. Add community assignments to colData ────────────────────────────────
  sce$community_id <- community_id[colnames(sce)]

  # ── 6. Summary ─────────────────────────────────────────────────────────────
  col_data <- as.data.frame(SummarizedExperiment::colData(sce))

  # Check that the user-specified (or default) QC columns are present
  available_cols <- names(col_data)

  if (!feature.col %in% available_cols) {
    warning(paste0(
      "Column '", feature.col, "' (feature.col) not found in colData(sce). ",
      "mean_detected will be returned as NA. "
    ))
  }
  if (!count.col %in% available_cols) {
    warning(paste0(
      "Column '", count.col, "' (count.col) not found in colData(sce). ",
      "mean_total_counts will be returned as NA. "
    ))
  }

  communities_summary <- col_data |>
    dplyr::group_by(community_id) |>
    dplyr::summarise(
      n_cells = dplyr::n(),
      mean_detected  = if (feature.col %in% available_cols) mean(.data[[feature.col]], na.rm = TRUE) else NA_real_,
      mean_total_counts = if (count.col %in% available_cols)  mean(.data[[count.col]],   na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )

  print(paste("Created", length(unique(community_id)), "communities"))
  print(paste("Mean cells per community:", round(mean(communities_summary$n_cells), 2)))

  return(list(sce = sce,
              communities_summary = communities_summary))
}
