#' Calculate Neighbourhood Clonal Diversity Indices
#'
#' Computes diversity metrics for each neighbourhood in a \code{Milo} or
#' \code{SingleCellExperiment} object, using clonal identifiers (clone IDs or
#' CDR3 sequences) extracted from \code{colData}.
#'
#' @param object A \code{Milo} or \code{SingleCellExperiment} object.
#'   Seurat objects must be converted via \code{as.SingleCellExperiment()}
#'   before use.
#' @param nhood_col Character string. For plain \code{SingleCellExperiment}
#'   objects (i.e. not \code{Milo}), the name of the \code{colData} column
#'   that assigns each cell to a neighbourhood.  Not required for
#'   \code{Milo} objects.
#' @param clone_id_col Character string. Name of the \code{colData} column
#'   containing clone identifiers.  If \code{NULL} or not found,
#'   \code{cdr3_col} is used as a fallback.
#' @param cdr3_col Character string. Name of the \code{colData} column
#'   containing CDR3 amino-acid sequences.  Used when \code{clone_id_col}
#'   is unavailable.  Default: \code{"cdr3_aa_VDJ"}.
#' @param min_n_contigs Integer. Minimum number of cells with contig
#'   availability required for a neighbourhood to be considered reliable.
#'   Neighbourhoods below this threshold are flagged with
#'   \code{diversity_reliable = FALSE} and their metric values are set to
#'   \code{NA}.  Default: \code{20}.
#' @param metrics Character vector of diversity metrics to compute.  Any
#'   combination of \code{"shannon"}, \code{"norm.entropy"},
#'   \code{"inv.simpson"}, \code{"gini.simpson"}, \code{"chao1"}, \code{"ace"},
#'   \code{"gini"} and \code{"standardised_n_clones"}.
#'   Default: \code{c("shannon", "inv.simpson", "gini.simpson", "chao1")}.
#'   \code{"standardised_n_clones"} can only be requested if
#'   \code{downsample_method = "coverage"}. This computes coverage-standardised
#'   observed clone richness via \code{iNEXT::estimateD} with Hill number order
#'   \eqn{q = 0}, and is not available under any other downsampling method.
#'   Note that \code{"chao1"}, \code{"ace"}, \code{"norm.entropy"}, and
#'   \code{"gini"} are not coverage-standardised when
#'   \code{downsample_method = "coverage"}; they are computed on observed
#'   counts directly.
#' @param downsample_method Character string specifying how to standardise
#'   sequencing depth across neighbourhoods before computing metrics.
#'   One of:
#'   \describe{
#'     \item{\code{"rarefaction"}}{Random subsampling to the median contig
#'       count across valid neighbourhoods (base R implementation). The
#'       rarefaction depth is reported at runtime. Neighbourhoods whose contig
#'       count is at or below the rarefaction depth are not downsampled and are
#'       used as-is; this can occur frequently when neighbourhood sizes vary
#'       widely.  Rarefaction is preferred when communities are highly diverse
#'       or singleton-dominated, as it standardises sample size rather than
#'       coverage. \code{"chao1"} and \code{"ace"} are not rarefied; they are
#'       computed on observed counts and \code{<metric>_diversity_sd} is
#'       \code{NA} for these metrics.}
#'     \item{\code{"vegan_rarefaction"}}{As above, using
#'       \code{vegan::rrarefy}.  Requires the \pkg{vegan} package.  The same
#'       guidance on rarefaction depth interpretation and \code{"chao1"}/
#'       \code{"ace"} compatibility applies.}
#'     \item{\code{"coverage"}}{Coverage-based standardisation via
#'       \code{iNEXT::estimateD}. Requires the \pkg{iNEXT} package. The target
#'       coverage is set to the median Good's coverage (1 - singletons / n)
#'       across valid neighbourhoods and is reported at runtime. Coverage-based
#'       standardisation is preferred when neighbourhood sizes are large and
#'       relatively uniform.  If the reported target coverage is very low (near
#'       0), rarefaction may be more appropriate. \code{"chao1"}, \code{"ace"},
#'       \code{"norm.entropy"}, and \code{"gini"} are not coverage standardised;
#'       they are computed on observed counts and \code{<metric>_diversity_se}
#'       is \code{NA} for these metrics.}
#'     \item{\code{"none"}}{No downsampling; metrics are computed on the
#'       observed counts directly. No spread columns are added to the output.}
#'   }
#'   Default: \code{"rarefaction"}.
#' @param n_boots Integer. Number of bootstrap replicates used when
#'   \code{downsample_method} is \code{"rarefaction"} or
#'   \code{"vegan_rarefaction"}.  Default: \code{100}.
#'   For reproducible results, call \code{set.seed()} before invoking this
#'   function.
#' @param feature_scale Logical.  If \code{TRUE}, a min-max scaled version
#'   of each diversity metric (range 0-1) is appended to the output as
#'   \code{<metric>_diversity_scaled}.  Default: \code{FALSE}.
#' @param result_name Character string.  Key under which the results
#'   \code{data.frame} is stored in \code{metadata(object)}.
#'   Default: \code{"nhood_diversity"}.
#'
#' @return The input \code{object}, with a \code{data.frame} of diversity
#'   results stored in \code{metadata(object)[[result_name]]}.  The
#'   \code{data.frame} contains the following columns:
#'   \describe{
#'     \item{\code{nhood}}{Neighbourhood identifier.}
#'     \item{\code{nhood_size}}{Total number of cells in the neighbourhood.}
#'     \item{\code{n_contigs}}{Number of cells with a contig.}
#'     \item{\code{contig_capture_rate}}{Fraction of cells with a contig
#'       (\code{n_contigs / nhood_size}).}
#'     \item{\code{n_clones}}{Observed number of unique clones.}
#'     \item{\code{clone_per_contig}}{Ratio of unique clones to cells with a
#'       valid contig (\code{n_clones / n_contigs}). Values near 1 indicate
#'       high clonal diversity; values near 0 indicate strong clonal
#'       expansion.}
#'     \item{\code{coverage}}{Good's coverage estimate
#'       (1 - singletons / n_contigs).}
#'     \item{\code{diversity_reliable}}{Logical; \code{TRUE} if
#'       \code{n_contigs >= min_n_contigs}.}
#'     \item{\code{<metric>_diversity}}{Estimated diversity for each
#'       requested metric.}
#'     \item{\code{<metric>_diversity_sd}}{Standard deviation of the metric
#'       across bootstrap iterations.  Present only when
#'       \code{downsample_method} is \code{"rarefaction"} or
#'       \code{"vegan_rarefaction"}.  \code{NA} for \code{"chao1"} and
#'       \code{"ace"}, which are not bootstrapped.}
#'     \item{\code{<metric>_diversity_se}}{Analytical standard error derived
#'       from \code{iNEXT} confidence intervals.  Present only when
#'       \code{downsample_method = "coverage"}.  \code{NA} for
#'       \code{"chao1"}, \code{"ace"}, \code{"norm.entropy"}, and
#'       \code{"gini"}, which are not coverage-standardised.}
#'     \item{\code{<metric>_diversity_scaled}}{Min-max scaled diversity
#'       (only present when \code{feature_scale = TRUE}).}
#'   }
#'   The following attributes are also set on the returned \code{data.frame}:
#'   \code{min_n_contigs} and \code{method} are always present.
#'   \code{rarefy_depth} is set only when \code{downsample_method} is
#'   \code{"rarefaction"} or \code{"vegan_rarefaction"}.
#'   \code{target_coverage} is set only when \code{downsample_method} is
#'   \code{"coverage"}.
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
#' ## Clone identity resolution
#' \code{clone_id_col} is tried first; if it is \code{NULL} or absent from
#' \code{colData}, the function falls back to \code{cdr3_col} with a warning.
#' Cells whose clone/CDR3 value is \code{NA}, \code{""}, \code{"No_contig"},
#' or \code{"None"} are excluded from diversity calculations but are still
#' counted in \code{nhood_size}.
#'
#' ## Metric selection
#' The \code{.chao1_richness} helper uses the bias-corrected Chao1 formula
#' (\eqn{S_{obs} + F_1(F_1-1) / (2(F_2+1))}) rather than the classic formula
#' (\eqn{S_{obs} + F_1^2 / (2 F_2)}).
#'
#' ## Downsampling and metric compatibility
#' Rarefaction (\code{"rarefaction"} and \code{"vegan_rarefaction"}) is not
#' applied to \code{"chao1"} or \code{"ace"}, as subsampling removes the
#' singleton and doubleton structure these estimators rely on.  Both are
#' always computed on observed counts under rarefaction; \code{_diversity_sd}
#' is \code{NA} for these metrics.
#'
#' Under \code{"coverage"}, only \code{"shannon"}, \code{"inv.simpson"},
#' \code{"gini.simpson"}, and \code{"standardised_n_clones"} are
#' coverage-standardised via \code{iNEXT}.  \code{"chao1"}, \code{"ace"},
#' \code{"norm.entropy"}, and \code{"gini"} fall back to observed-count
#' calculation with a runtime warning; \code{_diversity_se} is \code{NA}
#' for these metrics.
#'
#' \code{"standardised_n_clones"} is only valid with
#' \code{downsample_method = "coverage"}; requesting it under any other
#' method raises an error.
#'
#' ## Reproducibility
#' Bootstrap-based methods (\code{"rarefaction"}, \code{"vegan_rarefaction"})
#' use \code{base::sample} internally. For reproducible results, call
#' \code{set.seed()} before invoking this function.
#'
#' @note
#' The \code{vegan} package is required when \code{downsample_method =
#' "vegan_rarefaction"} and \code{iNEXT} is required when
#' \code{downsample_method = "coverage"}.
#'
#' @seealso
#' \code{\link{nhood_shm}}, \code{\link{nhood_cdr3_convergence}},
#' \code{\link{nhood_annotation}}, \code{\link{nhood_annotation_props}}
#'
#' @examples
#' \dontrun{
#' # --- Milo object, default metrics and rarefaction -------------------------
#' milo <- nhood_diversity(
#'   object       = milo_object,
#'   clone_id_col = "clone_id",
#'   cdr3_col     = "cdr3_aa_VDJ"
#' )
#' results <- S4Vectors::metadata(milo)[["nhood_diversity"]]
#'
#' # --- Coverage standardisation with standardised_n_clones ----------------
#' milo <- nhood_diversity(
#'   object            = milo_object,
#'   clone_id_col      = "clone_id",
#'   metrics           = c("shannon", "inv.simpson", "gini.simpson",
#'                         "standardised_n_clones"),
#'   downsample_method = "coverage",
#'   feature_scale     = TRUE,
#'   result_name       = "nhood_diversity_coverage"
#' )
#'
#' # --- Plain SingleCellExperiment object ------------------------------------
#' sce <- nhood_diversity(
#'   object       = sce_object,
#'   clone_id_col = "clone_id",
#'   nhood_col    = "neighbourhood_id"
#' )
#' }
#'
#' @export
#'
nhood_diversity <- function(object,
                            nhood_col = NULL,
                            clone_id_col = NULL,
                            cdr3_col = "cdr3_aa_VDJ",
                            min_n_contigs = 20,
                            metrics = c("shannon", "inv.simpson", "gini.simpson", "chao1"),
                            downsample_method = "rarefaction",
                            n_boots = 100,
                            feature_scale = FALSE,
                            result_name = "nhood_diversity") {

  # ── 0. Validate downsample_method and check dependencies ───────────────────
  valid_downsample_methods <- c("rarefaction", "vegan_rarefaction", "coverage", "none")
  if (!is.character(downsample_method) || length(downsample_method) != 1L ||
      is.na(downsample_method) || !downsample_method %in% valid_downsample_methods) {
    stop(sprintf(
      "'downsample_method' must be one of: %s. Received: %s",
      paste(sprintf("'%s'", valid_downsample_methods), collapse = ", "),
      paste(format(downsample_method), collapse = ", ")
    ), call. = FALSE)
  }

  if (downsample_method == "vegan_rarefaction" && !requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required for rarefaction. Please install it with: install.packages('vegan')")
  }

  if (downsample_method == "coverage" && !requireNamespace("iNEXT", quietly = TRUE)) {
    stop("Package 'iNEXT' is required for iNEXT method. Please install it with: install.packages('iNEXT')")
  }

  # ── 1. Validate metrics against downsample_method ──────────────────────────
  valid_metrics <- c("shannon", "inv.simpson", "gini.simpson", "chao1",
                     "norm.entropy", "gini", "ace")
  coverage_only_metrics <- "standardised_n_clones"

  if ("standardised_n_clones" %in% metrics && downsample_method != "coverage") {
    stop("'standardised_n_clones' is only valid when downsample_method = 'coverage'.")
  }

  invalid_metrics <- setdiff(metrics, c(valid_metrics, coverage_only_metrics))
  if (length(invalid_metrics) > 0) {
    stop(sprintf("Unknown metric(s): %s. Valid metrics are: %s",
                 paste(invalid_metrics, collapse = ", "),
                 paste(c(valid_metrics, coverage_only_metrics), collapse = ", ")))
  }

  if (any(metrics %in% c("chao1", "ace")) &&
      downsample_method %in% c("rarefaction", "vegan_rarefaction")) {
    warning(paste0("Rarefaction is not compatible with chao1 or ace estimators. ",
                   "These metrics will be computed on observed counts without rarefaction. ",
                   "_diversity_sd will be NA for these metrics."))
  }

  if (any(metrics %in% c("chao1", "ace", "norm.entropy", "gini")) &&
      downsample_method == "coverage") {
    warning(paste0("Coverage is not compatible with chao1 or ace estimators, normalised entropy and gini coefficient. ",
                   "These metrics will be computed on observed counts. ",
                   "_diversity_se will be NA for these metrics."))
  }

  # ── 2. Determine object type ───────────────────────────────────────────────
  is_milo     <- inherits(object, "Milo") && inherits(object, "SingleCellExperiment")
  is_sce_only <- !is_milo && inherits(object, "SingleCellExperiment")

  if (!is_milo && !is_sce_only) {
    stop("'object' must be a Milo or SingleCellExperiment object. If using Seurat, convert via as.SingleCellExperiment")
  }

  # ── 3. Extract colData and neighbourhood matrix ────────────────────────────
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
          "Cell ordering is inconsistent; cannot safely assign clonal identity values."
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

    unique_nhoods <- sort(as.character(unique(nhood_assignments[valid_cells])))
    n_nhoods_total <- length(unique_nhoods)

    nhoods <- matrix(0,
                     nrow = ncol(object),
                     ncol = n_nhoods_total,
                     dimnames = list(colnames(object),
                                     paste0("Nhood_", unique_nhoods)))

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

  # ── 4. Initialise results dataframe and pre-allocate metric vectors ────────
  diversity_results <- data.frame(
    nhood                = unique_nhoods,
    nhood_size           = integer(n_nhoods),
    n_contigs            = integer(n_nhoods),
    contig_capture_rate  = numeric(n_nhoods),
    n_clones             = integer(n_nhoods),
    clone_per_contig     = numeric(n_nhoods),
    coverage             = numeric(n_nhoods),
    diversity_reliable   = logical(n_nhoods),
    stringsAsFactors     = FALSE
  )

  # Pre-allocate per-metric result vectors.
  use_sd <- downsample_method %in% c("rarefaction", "vegan_rarefaction")
  use_se <- downsample_method == "coverage"

  spread_suffix <- if (use_sd) "_diversity_sd" else if (use_se) "_diversity_se" else NULL

  metric_col_names <- paste0(metrics, "_diversity")
  if (!is.null(spread_suffix)) {
    metric_col_names <- as.vector(rbind(
      paste0(metrics, "_diversity"),
      paste0(metrics, spread_suffix)
    ))
  }

  metric_cols <- vector("list", length(metric_col_names))
  names(metric_cols) <- metric_col_names
  for (nm in names(metric_cols)) {
    metric_cols[[nm]] <- rep(NA_real_, n_nhoods)
  }

  # ── 5. Single pass: cache clone counts, contig counts, and coverages ────────
  if (is.null(clone_id_col) || !(clone_id_col %in% colnames(col_data))) {
    if (!is.null(cdr3_col) && cdr3_col %in% colnames(col_data)) {
      warning(sprintf("'clone_id_col' (%s) is not found in colData. Falling back to CDR3 column '%s'.",
                      if (is.null(clone_id_col)) "NULL" else clone_id_col,
                      cdr3_col))
    } else {
      stop(sprintf(
        "Neither 'clone_id_col' (%s) nor 'cdr3_col' (%s) were found in colData. Specify a valid clone_id_col or cdr3_col.",
        if (is.null(clone_id_col)) "NULL" else clone_id_col,
        if (is.null(cdr3_col)) "NULL" else cdr3_col
      ))
    }
  }

  cached_clone_counts <- vector("list", n_nhoods)
  cached_clone_ids <- vector("list", n_nhoods)
  all_n_contigs <- integer(n_nhoods)
  all_coverages <- numeric(n_nhoods)

  for (i in seq_len(n_nhoods)) {
    cells_in_nhood <- which(nhoods[, i] == 1)
    cell_ids       <- rownames(nhoods)[cells_in_nhood]
    nhood_col_data <- col_data[cell_ids, , drop = FALSE]

    if (!is.null(clone_id_col) && clone_id_col %in% colnames(nhood_col_data)) {
      valid_cells <- !is.na(nhood_col_data[[clone_id_col]]) &
        nhood_col_data[[clone_id_col]] != "No_contig" &
        nhood_col_data[[clone_id_col]] != "None" &
        nhood_col_data[[clone_id_col]] != ""
      clone_identifiers <- nhood_col_data[[clone_id_col]][valid_cells]
    } else {
      valid_cells <- !is.na(nhood_col_data[[cdr3_col]]) &
        nhood_col_data[[cdr3_col]] != "No_contig" &
        nhood_col_data[[cdr3_col]] != "None" &
        nhood_col_data[[cdr3_col]] != ""
      clone_identifiers <- nhood_col_data[[cdr3_col]][valid_cells]
    }

    n_contigs <- length(clone_identifiers)
    all_n_contigs[i] <- n_contigs
    cached_clone_ids[[i]] <- clone_identifiers

    if (n_contigs >= min_n_contigs) {
      clone_counts          <- table(clone_identifiers)
      cached_clone_counts[[i]] <- clone_counts
      f1                    <- sum(clone_counts == 1)
      all_coverages[i]      <- 1 - (f1 / n_contigs)
    } else {
      cached_clone_counts[[i]] <- NULL
      all_coverages[i]         <- NA_real_
    }
  }

  valid_contigs <- all_n_contigs[all_n_contigs >= min_n_contigs]

  rarefy_depth   <- NULL
  target_coverage <- NULL

  if (downsample_method %in% c("rarefaction", "vegan_rarefaction")) {
    rarefy_depth <- if (length(valid_contigs) > 0) {
      round(stats::quantile(valid_contigs, 0.5, na.rm = TRUE))
    } else {
      min_n_contigs
    }
  } else if (downsample_method == "coverage") {
    target_coverage <- stats::quantile(all_coverages, 0.5, na.rm = TRUE)
  }

  # ── 6. Calculate diversity for each neighbourhood ──────────────────────────
  for (i in seq_len(n_nhoods)) {
    cells_in_nhood  <- which(nhoods[, i] == 1)
    n_contigs       <- all_n_contigs[i]
    clone_identifiers <- cached_clone_ids[[i]]

    diversity_results[i, "nhood_size"]          <- length(cells_in_nhood)
    diversity_results[i, "n_contigs"]           <- n_contigs
    diversity_results[i, "contig_capture_rate"] <- n_contigs / length(cells_in_nhood)

    if (n_contigs < min_n_contigs) {
      diversity_results[i, "n_clones"]         <- length(unique(clone_identifiers))
      diversity_results[i, "clone_per_contig"] <- if (n_contigs > 0) length(unique(clone_identifiers)) / n_contigs else NA_real_
      diversity_results[i, "diversity_reliable"] <- FALSE
      diversity_results[i, "coverage"]           <- NA

    } else {
      clone_counts <- cached_clone_counts[[i]]
      diversity_results[i, "n_clones"]         <- length(clone_counts)
      diversity_results[i, "clone_per_contig"] <- length(clone_counts) / n_contigs
      diversity_results[i, "diversity_reliable"] <- TRUE

      f1 <- sum(clone_counts == 1)
      diversity_results[i, "coverage"] <- 1 - (f1 / n_contigs)

      metric_result <- .calculate_diversity_metrics(
        clone_counts      = clone_counts,
        clone_identifiers = clone_identifiers,
        metrics           = metrics,
        downsample_method = downsample_method,
        n_boots           = n_boots,
        rarefy_depth      = rarefy_depth,
        target_coverage   = target_coverage
      )

      for (nm in names(metric_cols)) {
        metric_cols[[nm]][i] <- metric_result[[nm]]
      }
    }
  }

  # ── 7. Collate metric vectors into the results frame ───────────────────────
  for (nm in names(metric_cols)) {
    diversity_results[[nm]] <- metric_cols[[nm]]
  }

  # ── 8. Summary ─────────────────────────────────────────────────────────────
  total_cells    <- ncol(object)
  clone_col_used <- if (!is.null(clone_id_col) && clone_id_col %in% colnames(col_data)) clone_id_col else cdr3_col
  total_n_contigs  <- sum(!is.na(col_data[[clone_col_used]]) &
                            col_data[[clone_col_used]] != "No_contig" &
                            col_data[[clone_col_used]] != "None" &
                            col_data[[clone_col_used]] != "")
  message("Contig summary:")
  message("Cells with contig: ", total_n_contigs, " / ", total_cells)
  message("Neighbourhoods with >= ", min_n_contigs, " cells with contig: ",
          sum(all_n_contigs >= min_n_contigs), " / ", n_nhoods)
  if (downsample_method %in% c("rarefaction", "vegan_rarefaction")) {
    message("Rarefaction depth: ", rarefy_depth)

    rarefiable_metrics <- setdiff(metrics, c("chao1", "ace"))
    if (length(rarefiable_metrics) > 0) {
      sd_col <- paste0(rarefiable_metrics[1], "_diversity_sd")
      reliable_rows <- diversity_results[["diversity_reliable"]]
      n_rarefied <- sum(reliable_rows & !is.na(diversity_results[[sd_col]]))
      n_reliable <- sum(reliable_rows)
      message("Neighbourhoods successfully rarefied (n_contigs > rarefaction depth): ",
              n_rarefied, " / ", n_reliable)
    }

  } else if (downsample_method == "coverage") {
    message("Target coverage: ", round(target_coverage, 3))

    coverage_unsupported <- c("chao1", "ace", "norm.entropy", "gini")
    coverage_eligible_metrics <- setdiff(metrics, coverage_unsupported)
    reliable_rows <- diversity_results[["diversity_reliable"]]
    n_reliable <- sum(reliable_rows)

    for (metric in coverage_eligible_metrics) {
      se_col <- paste0(metric, "_diversity_se")
      n_achieved <- sum(reliable_rows & !is.na(diversity_results[[se_col]]))
      message("Neighbourhoods where target coverage was achieved for '", metric, "': ",
              n_achieved, " / ", n_reliable)
    }
  }

  # ── 9. Add attributes ──────────────────────────────────────────────────────
  if (downsample_method %in% c("rarefaction", "vegan_rarefaction")) {
    attr(diversity_results, "rarefy_depth") <- rarefy_depth
  } else if (downsample_method == "coverage") {
    attr(diversity_results, "target_coverage") <- target_coverage
  }
  attr(diversity_results, "min_n_contigs") <- min_n_contigs
  attr(diversity_results, "method") <- downsample_method

  # ── 10. Scale output ───────────────────────────────────────────────────────
  if (feature_scale) {
    for (metric in metrics) {
      col <- paste0(metric, "_diversity")
      scaled_col <- paste0(metric, "_diversity_scaled")
      vals <- diversity_results[[col]]
      min_val <- min(vals, na.rm = TRUE)
      max_val <- max(vals, na.rm = TRUE)
      if (!is.finite(min_val) || !is.finite(max_val) || min_val == max_val) {
        diversity_results[[scaled_col]] <- NA_real_
      } else {
        diversity_results[[scaled_col]] <- (vals - min_val) / (max_val - min_val)
      }
      col_order <- colnames(diversity_results)
      col_order <- col_order[col_order != scaled_col]
      spread_col <- if (use_sd) paste0(metric, "_diversity_sd") else
        if (use_se) paste0(metric, "_diversity_se") else NULL
      anchor_col <- if (!is.null(spread_col) && spread_col %in% col_order) {
        spread_col
      } else {
        paste0(metric, "_diversity")
      }
      anchor_pos <- which(col_order == anchor_col)
      diversity_results <- diversity_results[, append(col_order, scaled_col, after = anchor_pos)]
    }
  }

  # ── 11. Store results in metadata and return ───────────────────────────────
  if (result_name %in% names(S4Vectors::metadata(object))) {
    warning(sprintf(
      "metadata(object)[['%s']] already exists and will be overwritten.",
      result_name
    ))
  }

  S4Vectors::metadata(object)[[result_name]] <- diversity_results

  message(sprintf(
    "Stored diversity results (%d neighbourhoods x %d metrics) in metadata(object)[['%s']].",
    n_nhoods, length(metrics), result_name
  ))

  return(object)
}

#' Compute Diversity Metrics for a Single Neighbourhood
#'
#' Internal helper called once per neighbourhood inside
#' \code{nhood_diversity}.  Applies the requested downsampling strategy and
#' computes each metric via \code{.calculate_single_metric}, returning point
#' estimates and spread values where applicable.
#'
#' @param clone_counts A named integer \code{table} (or equivalent) of
#'   clone counts for the neighbourhood, with zeros already removed.
#' @param clone_identifiers Character vector of per-cell clone labels for
#'   the neighbourhood (length equal to \code{n_contigs} for that
#'   neighbourhood).
#' @param metrics Character vector of metric names to compute.  See
#'   \code{nhood_diversity} for valid values.
#' @param downsample_method Character string; the downsampling strategy.
#'   See \code{nhood_diversity} for details.
#' @param n_boots Integer; number of bootstrap replicates.
#' @param rarefy_depth Integer; target sample size for rarefaction,
#'   pre-computed as the median contig count across valid neighbourhoods.
#' @param target_coverage Numeric between 0-1; target Good's coverage for
#'   the \code{"coverage"} method, pre-computed as the median coverage
#'   across valid neighbourhoods.
#'
#' @return A named list.  Each metric contributes one or two elements:
#'   \code{<metric>_diversity} (point estimate, always present) and either
#'   \code{<metric>_diversity_sd} (SD across bootstrap iterations, only
#'   under rarefaction methods) or \code{<metric>_diversity_se} (analytical
#'   SE from iNEXT, only under coverage).  No spread element is returned
#'   under \code{"none"}.
#'
#' @noRd
.calculate_diversity_metrics <- function(clone_counts, clone_identifiers,
                                         metrics, downsample_method, n_boots,
                                         rarefy_depth, target_coverage) {

  n_contigs  <- length(clone_identifiers)
  use_sd     <- downsample_method %in% c("rarefaction", "vegan_rarefaction")
  use_se     <- downsample_method == "coverage"
  out        <- list()

  for (metric in metrics) {

    div_val    <- NA_real_
    div_spread <- NA_real_

    # No downsampling: compute on observed counts; no spread column produced.
    if (downsample_method == "none") {
      div_val <- .calculate_single_metric(metric, clone_counts)

      # Rarefaction
    } else if (downsample_method %in% c("rarefaction", "vegan_rarefaction")) {

      # chao1 and ace bypass bootstrapping
      if (metric %in% c("chao1", "ace")) {
        div_val <- .calculate_single_metric(metric, clone_counts)
      }

      # Neighbourhood at or below rarefaction depth use observed counts
      else if (n_contigs <= rarefy_depth) {
        div_val <- .calculate_single_metric(metric, clone_counts)
      }

      # Base R rarefaction
      else if (downsample_method == "rarefaction") {
        expanded_clones <- rep(names(clone_counts), clone_counts)
        boot_values <- replicate(n_boots, {
          rarefied_sample <- sample(expanded_clones, size = rarefy_depth, replace = FALSE)
          rarefied_counts <- table(rarefied_sample)
          .calculate_single_metric(metric, rarefied_counts)
        })
        div_val    <- mean(boot_values, na.rm = TRUE)
        div_spread <- stats::sd(boot_values,   na.rm = TRUE)
      }

      # Vegan rarefaction
      else {
        if (!requireNamespace("vegan", quietly = TRUE)) {
          stop("vegan package required for vegan_rarefaction method")
        }
        comm_matrix <- t(as.matrix(clone_counts))
        boot_values <- replicate(n_boots, {
          rarefied_comm <- vegan::rrarefy(comm_matrix, sample = rarefy_depth)
          .calculate_single_metric(metric, as.integer(rarefied_comm[rarefied_comm > 0]))
        })
        div_val    <- mean(boot_values, na.rm = TRUE)
        div_spread <- stats::sd(boot_values,   na.rm = TRUE)
      }

      # Coverage
    } else if (downsample_method == "coverage") {
      if (!requireNamespace("iNEXT", quietly = TRUE)) {
        stop("iNEXT package required for coverage downsampling method. Install with: install.packages('iNEXT')")
      }

      # Metrics not supported by iNEXT: fall back to observed counts with warning.
      if (metric %in% c("chao1", "ace", "norm.entropy", "gini")) {
        div_val <- .calculate_single_metric(metric, clone_counts)
      }

      # iNEXT-supported metrics
      else {
        abundance_vector <- as.vector(clone_counts)
        coverage_result <- tryCatch({

          # Shannon
          if (metric == "shannon") {
            inext_result <- iNEXT::estimateD(abundance_vector, q = 1,
                                             datatype = "abundance",
                                             base = "coverage",
                                             level = target_coverage,
                                             conf = 0.95)
            qD    <- inext_result$qD[1]
            qD_se <- (inext_result$qD.UCL[1] - inext_result$qD.LCL[1]) / (2 * 1.96)
            # iNEXT returns the Hill number (exp(H')); convert to Shannon scale
            # and propagate SE via the delta method: SE(log(x)) = SE(x) / x
            list(val = log(qD), se = qD_se / qD)

            # Inverse simpson
          } else if (metric == "inv.simpson") {
            inext_result <- iNEXT::estimateD(abundance_vector, q = 2,
                                             datatype = "abundance",
                                             base = "coverage",
                                             level = target_coverage,
                                             conf = 0.95)
            list(val = inext_result$qD[1],
                 se  = (inext_result$qD.UCL[1] - inext_result$qD.LCL[1]) / (2 * 1.96))

            # Gini simpson
          } else if (metric == "gini.simpson") {
            inext_result <- iNEXT::estimateD(abundance_vector, q = 2,
                                             datatype = "abundance",
                                             base = "coverage",
                                             level = target_coverage,
                                             conf = 0.95)
            inv_simp    <- inext_result$qD[1]
            inv_simp_se <- (inext_result$qD.UCL[1] - inext_result$qD.LCL[1]) / (2 * 1.96)
            list(val = 1 - (1 / inv_simp),
                 se  = inv_simp_se / (inv_simp^2))

            # Coverage-standardised observed richness via Hill q=0
          } else if (metric == "standardised_n_clones") {
            inext_result <- iNEXT::estimateD(abundance_vector, q = 0,
                                             datatype = "abundance",
                                             base = "coverage",
                                             level = target_coverage,
                                             conf = 0.95)
            list(val = inext_result$qD[1],
                 se  = (inext_result$qD.UCL[1] - inext_result$qD.LCL[1]) / (2 * 1.96))
          }
        }, error = function(e) {
          list(val = tryCatch(.calculate_single_metric(metric, clone_counts),
                              error = function(e2) NA_real_),
               se  = NA_real_)
        })
        div_val    <- coverage_result$val
        div_spread <- coverage_result$se
      }
    }

    # Store point estimate; add spread element only when the method requires it.
    out[[paste0(metric, "_diversity")]] <- div_val
    if (use_sd) out[[paste0(metric, "_diversity_sd")]] <- div_spread
    if (use_se) out[[paste0(metric, "_diversity_se")]] <- div_spread
  }

  out
}

#' Dispatch a Single Diversity Metric Calculation
#'
#' Thin \code{switch}-based dispatcher that routes a named metric to the
#' appropriate diversity function.  Used by
#' \code{.calculate_diversity_metrics} to keep per-metric dispatch logic
#' in one place.
#'
#' @param metric Character string; one of \code{"shannon"},
#'   \code{"inv.simpson"}, \code{"gini.simpson"}, \code{"chao1"},
#'   \code{"norm.entropy"}, \code{"gini"}, or \code{"ace"}.
#' @param counts Numeric vector of non-negative counts (zeros are handled
#'   by the underlying metric functions).
#'
#' @return A single numeric value, or \code{NA} if \code{metric} is not
#'   recognised.
#'
#' @noRd
.calculate_single_metric <- function(metric, counts) {
  switch(metric,
         "shannon"             = .shannon_entropy(counts),
         "inv.simpson"         = .inv_simpson(counts),
         "gini.simpson"        = .gini_simpson(counts),
         "chao1"               = .chao1_richness(counts),
         "norm.entropy"        = .norm_entropy(counts),
         "gini"                = .gini_coef(counts),
         "ace"                 = .ace_richness(counts),
         "standardised_n_clones" = NA_real_,
         NA
  )
}

##### Diversity calculation functions #####
## Adapted from immApex (BorchLab), R/diversity.R, v1.7.2 (devel):
## https://github.com/BorchLab/immApex/blob/devel/R/diversity.R
## Copyright (c) Nick Borcherding; immApex is MIT-licensed. Function bodies are
## unchanged apart from being renamed with a leading '.'. Only the functions
## used by nhood_diversity() are retained, and roxygen documentation has been
## edited. See LICENSE.note / inst/COPYRIGHTS.

#' Shannon Diversity Index (Entropy)
#'
#' Calculates Shannon’s information entropy (often denoted *H*) for a set
#' of clone or sequence counts.
#'
#' \deqn{H' \;=\; -\sum_{i = 1}^{S} p_i \ln p_i}
#'
#' where *p*<sub>*i*</sub> = *n*<sub>*i*</sub> / *N* are the relative
#' frequencies (proportions) of each of the *S* distinct categories.
#'
#' @param cnt Numeric vector of non-negative counts (one entry per clone/
#'   residue/OTU).  Zero counts are ignored.
#'
#' @return A single numeric value (>= 0).  When `cnt` contains exactly one
#'   positive entry the function returns `0`.
#'
#' @references
#' Shannon, C. E. (1948). *A mathematical theory of communication*.
#' Bell System Technical Journal, 27(3), 379-423.
#'
#' @seealso \code{\link{.norm_entropy}}, \code{\link{.inv_simpson}}
#' @examples
#' counts <- c(A = 12, B = 4, C = 4)
#' .shannon_entropy(counts)
#' @noRd
.shannon_entropy <- function(cnt) {
  cnt <- cnt[cnt > 0]
  p   <- cnt / sum(cnt)
  -sum(p * log(p))
}

#' Inverse Simpson Diversity
#'
#' Computes the inverse of Simpson’s concentration index, sometimes written
#' as *1/D*.  This metric emphasizes dominant categories.
#'
#' \deqn{1/D \;=\; \frac{1}{\sum_{i} p_i^{\,2}}}
#'
#' @param cnt Numeric vector of non-negative counts. Zero counts are ignored.
#' @return Numeric value >= 1.  Equals 1 when all observations belong to a
#'   single category.
#'
#' @references
#' Simpson, E. H. (1949). *Measurement of diversity*. Nature, 163, 688.
#'
#' @seealso \code{\link{.gini_simpson}}, \code{\link{.shannon_entropy}}
#' @examples
#' .inv_simpson(c(10, 5, 1))
#' @noRd
.inv_simpson <- function(cnt) {
  cnt <- cnt[cnt > 0]
  p   <- cnt / sum(cnt)
  1 / sum(p * p)
}

#' Gini–Simpson Diversity
#'
#' Computes the complement of Simpson’s index (also called the
#' Gini–Simpson index or probability of interspecific encounter):
#'
#' \deqn{1 - \lambda = 1 - \sum_{i} p_i^{\,2}}
#'
#' @param cnt Numeric vector of non-negative counts. Zero counts are ignored.
#' @return Value between 0-1.  Higher numbers indicate greater heterogeneity.
#'
#' @references
#' Simpson, E. H. (1949). *Measurement of diversity*. Nature, 163, 688.
#'
#' @seealso \code{\link{.inv_simpson}}, \code{\link{.gini_coef}}
#' @examples
#' .gini_simpson(c(10, 5, 5))
#' @noRd
.gini_simpson <- function(cnt) 1 - 1 / .inv_simpson(cnt)

#' Normalised Shannon Entropy
#'
#' Shannon entropy scaled to 0-1 by its maximum possible value given *S*
#' observed categories:
#'
#' \deqn{H^* = \frac{H'}{\ln S}}
#'
#' (also known as “Shannon evenness”).
#'
#' @param cnt Numeric vector of non-negative counts. Zero counts are ignored.
#' @return Numeric value between 0-1; `0` when all observations are in a
#'   single category.
#'
#' @references
#' Pielou, E. C. (1966). *The measurement of diversity in different types of
#' biological collections*. Journal of Theoretical Biology, 13, 131-144.
#'
#' @seealso \code{\link{.shannon_entropy}}
#' @examples
#' .norm_entropy(c(40, 10, 10, 10))
#' @noRd
.norm_entropy <- function(cnt) {
  cnt <- cnt[cnt > 0]
  S   <- length(cnt)
  if (S == 1L) 0 else .shannon_entropy(cnt) / log(S)
}

#' ACE Richness Estimator
#'
#' Calculates the Abundance-based Coverage Estimator (ACE) of species richness.
#' This metric is particularly useful for datasets with a large number of rare
#' species.
#'
#' \deqn{S_{ace} = S_{abund} + \frac{S_{rare}}{C_{ace}} + \frac{F_1}{C_{ace}} \gamma^2_{ace}}
#'
#' where the classification of rare and abundant species is based on a
#' threshold of 10 individuals, *F*1 is the count of singletons,
#' *S*rare is the number of rare species, and *C*ace is the
#' sample coverage for rare species.
#'
#' @param cnt Numeric vector of non-negative counts. Zero counts are ignored.
#' @return A single numeric value representing the estimated total number of
#'   species. The estimate is constrained to be at least the number of
#'   observed species.
#'
#' @references
#' Chao, A., & Lee, S.-M. (1992). *Estimating the number of classes via sample coverage*.
#' Journal of the American Statistical Association, 87(417), 210-217.
#'
#' @seealso \code{\link{.chao1_richness}}, \code{\link{.shannon_entropy}}
#' @examples
#' counts <- rpois(50, lambda=1.5)
#' .ace_richness(counts)
#' @noRd
.ace_richness <- function(cnt) {
  cnt <- cnt[cnt > 0]
  S_obs <- length(cnt)
  if (S_obs == 0) return(0)

  rare_thresh <- 10
  is_rare <- cnt <= rare_thresh
  S_rare <- sum(is_rare)
  S_abund <- sum(!is_rare)

  if (S_rare == 0) return(S_obs)

  cnt_rare <- cnt[is_rare]
  N_rare <- sum(cnt_rare)
  F1 <- sum(cnt_rare == 1)

  if (N_rare == 0) return(S_abund)

  C_ace <- 1 - F1 / N_rare
  if (C_ace == 0) return(S_obs) # Avoid division by zero

  sum_gamma <- sum(cnt_rare * (cnt_rare - 1))

  gamma_sq_num <- S_rare * sum_gamma
  gamma_sq_den <- C_ace * N_rare * (N_rare - 1)

  gamma_sq <- if (gamma_sq_den == 0) {
    0
  } else {
    max(0, (gamma_sq_num / gamma_sq_den) - 1)
  }

  est <- S_abund + (S_rare / C_ace) + (F1 / C_ace) * gamma_sq

  max(S_obs, est)
}

#' Gini Coefficient of Abundance Inequality
#'
#' Calculates the Gini coefficient, a measure of inequality, for a vector
#' of clone/sequence counts. It ranges from 0 (perfect equality) to nearly 1
#' (maximal inequality).
#'
#' \deqn{G = \frac{\sum_{i=1}^{S} (2i - S - 1) n_i}{S \sum_{i=1}^{S} n_i}}
#'
#' where *n*i are the counts of each of the *S* categories,
#' sorted in non-decreasing order.
#'
#' @param cnt Numeric vector of non-negative counts. Zero counts are ignored.
#' @return A numeric value between 0-1. Returns `0` if there is only one
#'   category.
#'
#' @seealso \code{\link{.gini_simpson}}
#' @examples
#' # High inequality
#' .gini_coef(c(100, 1, 1, 1))
#' # Perfect equality
#' .gini_coef(c(10, 10, 10, 10))
#' @noRd
.gini_coef <- function(cnt) {
  cnt <- cnt[cnt > 0]
  S <- length(cnt)
  if (S <= 1) return(0)

  cnt_sorted <- sort(cnt)
  i <- 1:S

  num <- sum((2 * i - S - 1) * cnt_sorted)
  den <- S * sum(cnt_sorted)

  num / den
}

#' Chao1 Richness Estimator
#'
#' Calculates the Chao1 non-parametric estimator of species richness.
#'
#' The bias-corrected formula is used:
#' \deqn{S_{chao1} = S_{obs} + \frac{F_1(F_1 - 1)}{2(F_2 + 1)}}
#' where *S*obs is the number of observed species, *F*1
#' is the count of singletons, and *F*2 is the count of doubletons.
#'
#' If the conditions for the formula are not met (*F*1 <= 1 or
#' *F*2 = 0), the function returns the observed richness (*S*obs).
#'
#' @param cnt Numeric vector of non-negative counts. Zero counts are ignored.
#' @return A single numeric value representing the estimated total number of
#'   species.
#'
#' @references
#' Chao, A. (1984). *Nonparametric estimation of the number of classes in a
#' population*. Scandinavian Journal of Statistics, 11(4), 265-270.
#'
#' @seealso \code{\link{.ace_richness}}, \code{\link{.shannon_entropy}}
#' @examples
#' # Sample with singletons and doubletons
#' counts <- c(rep(1, 10), rep(2, 5), 5, 8, 12)
#' .chao1_richness(counts)
#'
#' # Sample without doubletons returns observed richness
#' .chao1_richness(c(rep(1, 5), 3, 4, 5))
#' @noRd
.chao1_richness <- function(cnt) {
  cnt <- cnt[cnt > 0]
  S_obs <- length(cnt)
  if (S_obs == 0) return(0)

  F1 <- sum(cnt == 1)
  F2 <- sum(cnt == 2)

  # Use bias-corrected Chao1 formula
  if (F1 > 1 && F2 > 0) {
    S_obs + (F1 * (F1 - 1)) / (2 * (F2 + 1))
  } else {
    # Fallback to observed richness if conditions not met
    S_obs
  }
}
