#' Test Leiden clustering resolutions on a nearest-neighbour graph
#'
#' Runs Leiden community detection across a range of resolution values and
#' returns a summary table of clustering metrics for each resolution. Use this
#' function to select an appropriate resolution before running
#' \code{\link{create_communities}}.
#'
#' @param sce A \code{SingleCellExperiment} object. The object must contain a
#'   pre-computed nearest-neighbour graph stored in
#'   \code{metadata(sce)$graphs} under the name given by \code{graph.name}.
#'   See the \strong{Graph preparation} section below for how to build and
#'   store the graph. Supplying a Seurat object will raise an error.
#' @param graph.name A character string giving the name of the graph to
#'   retrieve from \code{metadata(sce)$graphs}. Defaults to \code{"RNA_snn"}.
#' @param resolutions A numeric vector of resolution values to test. Higher
#'   values produce more, smaller communities. Defaults to
#'   \code{c(0.5, 1, 2.5, 5, 10, 20)}.
#' @param min_cells A positive integer. Communities with fewer than this
#'   many cells are considered too small; their cells are reassigned to the
#'   nearest valid community. Defaults to \code{5}.
#'
#' @return A \code{data.frame} with one row per tested resolution and the
#'   following columns:
#'   \describe{
#'     \item{\code{resolution}}{The resolution value tested.}
#'     \item{\code{n_communities}}{Number of communities after
#'       reassignment of small clusters.}
#'     \item{\code{mean_cells_per_community}}{Mean community size.}
#'     \item{\code{sd_cells_per_community}}{Standard deviation of
#'       community sizes.}
#'     \item{\code{median_cells_per_community}}{Median community size.}
#'     \item{\code{min_cells_per_community}}{Smallest community size.}
#'     \item{\code{max_cells_per_community}}{Largest community size.}
#'     \item{\code{modularity}}{Graph modularity of the partition.}
#'   }
#'   Resolutions at which all communities fall below \code{min_cells} are
#'   skipped with a warning and will not appear in the output.
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
#' @seealso \code{\link{create_communities}} to assign cells to communities
#' once a resolution has been chosen
#'
#' @examples
#' \dontrun{
#' library(scran)
#'
#' # Build and store the graph
#' g <- scran::buildSNNGraph(sce, use.dimred = "PCA", k = 20)
#' metadata(sce)$graphs[["RNA_snn"]] <- g
#'
#' # Test a range of resolutions
#' res_results <- test_resolution(sce, resolutions = c(0.5, 1, 2.5, 5))
#' print(res_results)
#' }
#'
#' @export
#'
test_resolution <- function(sce,
                            graph.name = "RNA_snn",
                            resolutions = c(0.5,1,2.5,5,10,20),
                            min_cells = 5) {

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
        ". Run scran::buildSNNGraph() and store the result in metadata(sce)$graphs first, ",
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
    # Reconstruct the adjacency matrix for cell-lookup operations below
    knn_matrix <- igraph::as_adjacency_matrix(igraph_obj, attr = "weight", sparse = TRUE)
  } else {
    knn_matrix <- knn_graph
    igraph_obj <- igraph::graph_from_adjacency_matrix(knn_matrix,
                                                      mode = "max",
                                                      weighted = TRUE)
  }

  # ── 4. Compute communities at resolutions for testing ──────────────────────
  results <- data.frame()

  cat("Testing resolutions...\n")

  for (res in resolutions) {
    cat(paste0("Resolution: ", res, "\n"))

    # Run leiden
    communities <- igraph::cluster_leiden(igraph_obj,
                                          resolution = res,
                                          objective_function = "modularity")
    community_id <- communities$membership
    names(community_id) <- rownames(knn_matrix)

    # Filter out small communities
    community_counts <- table(community_id)
    valid_communities <- names(community_counts)[community_counts >= min_cells]

    # Check if any valid communities exist
    if (length(valid_communities) == 0) {
      warning(paste0("Resolution ", res, ": All communities below min_cells threshold. Skipping."))
      next
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

    # Calculate metrics from single run
    community_sizes <- table(community_id)

    results <- rbind(results, data.frame(
      resolution = res,
      n_communities = length(unique(community_id)),
      mean_cells_per_community = mean(community_sizes),
      sd_cells_per_community = stats::sd(community_sizes),
      median_cells_per_community = stats::median(community_sizes),
      min_cells_per_community = min(community_sizes),
      max_cells_per_community = max(community_sizes),
      modularity = igraph::modularity(igraph_obj, as.numeric(factor(community_id)))
    ))
  }

  cat("\nResolution testing complete!\n")

  return(results)
}
