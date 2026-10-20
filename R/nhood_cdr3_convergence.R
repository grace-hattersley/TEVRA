#' Compute CDR3 amino acid convergence scores across neighbourhoods
#'
#' Computes a per neighbourhood summary convergence statistic of pairwise CDR3
#' amino acid sequence distances between B cell receptor (BCR) or T cell
#' receptor (TCR) clones.
#'
#' @param object A \code{\link[miloR]{Milo}} or
#'   \code{\link[SingleCellExperiment]{SingleCellExperiment}} object. Seurat
#'   objects must first be converted via
#'   \code{\link[Seurat]{as.SingleCellExperiment}}.
#' @param nhood_col Character string. For plain \code{SingleCellExperiment}
#'   objects (i.e. not \code{Milo}), the name of the \code{colData} column
#'   containing neighbourhood assignments. Not required for \code{Milo}
#'   objects.
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
#' @param summary_metric Character vector; one or more of \code{"converge_score"}
#'   (default), \code{"converge_score_global"}, and \code{"converge_score_global_median"}.
#'   Controls which summary convergence metrics are calculated and returned:
#'   \describe{
#'     \item{\code{"converge_score"}}{kNN convergence score per neighbourhood.
#'       Captures local sequence density. \code{NA} when
#'       \code{n_clones <= resolved_k} (kNN cannot be computed), or
#'       when \code{min_clones} is supplied and not met.}
#'     \item{\code{"converge_score_global"}}{Mean pairwise CDR3 convergence score
#'       across the upper triangle of the representative distance matrix.}
#'     \item{\code{"converge_score_global_median"}}{Median pairwise CDR3 convergence
#'       score across the upper triangle of the representative distance matrix.}
#'   }
#'   Any combination of the three values may be supplied. The default is
#'   \code{"converge_score"}.
#' @param knn_k Integer or \code{NULL} (default). Number of nearest neighbours
#'   to use when computing \code{converge_score}. When \code{NULL}, \eqn{k} is
#'   set adaptively to \eqn{\lfloor 0.25 \times n_{\text{reps}} \rfloor} per
#'   neighbourhood. A floor of 3 ensures the metric
#'   captures more than just the nearest pairwise comparisons. When a fixed integer
#'   is supplied it is used as-is with no clamping applied. Within each bootstrap subsample \eqn{k} is additionally
#'   capped at \eqn{n - 1} to guard against subsamples smaller than \eqn{k}.
#' @param downsample Logical; whether to downsample large neighbourhoods to a
#'   common representative depth before computing distances. The downsampling
#'   depth is set to the median number of unique representatives across all
#'   passing neighbourhoods (excluding singletons), computed separately per
#'   chain. Downsampling is performed via bootstrapping (\code{n_boots}
#'   iterations). Default is \code{TRUE}. When \code{FALSE}, the full
#'   representative matrix is used directly and \code{_sd} columns are not
#'   reported.
#' @param n_boots Integer; number of bootstrap iterations used during
#'   downsampling. Each iteration subsamples \code{downsample_depth}
#'   representatives without replacement and computes distance summaries;
#'   reported values are the mean across iterations and the corresponding
#'   \code{_sd} column reports the standard deviation across iterations for
#'   each requested metric. Default is \code{100}. Ignored when
#'   \code{downsample = FALSE} or when a neighbourhood has fewer
#'   representatives than the downsampling depth.
#'   When \code{downsample = TRUE} and \code{"converge_score"} is in
#'   \code{summary_metric}, the effective sample size for kNN computation
#'   within each bootstrap iteration is the downsampling depth rather than the
#'   full representative count. If this depth is \eqn{\leq} \eqn{k},
#'   \code{converge_score} will be \code{NA} for all downsampled neighbourhoods
#'   and a warning is emitted. In this scenario, set \code{downsample = FALSE}
#'   or use \code{summary_metric = "converge_score_global"} to obtain values for
#'   all neighbourhoods.
#' @param seed Integer or \code{NULL} (default). Random seed for
#'   reproducibility of the bootstrap downsampling procedure. When supplied, a
#'   neighbourhood-specific seed is derived as \code{seed + i} (where \code{i}
#'   is the neighbourhood index), ensuring that each neighbourhood's bootstrap
#'   is independently reproducible while the full run remains deterministic.
#'   Set this parameter whenever results will be reported or compared across
#'   runs. Has no effect when \code{downsample = FALSE}.
#' @param result_name Character string; the key under which results are stored
#'   in \code{metadata(object)}. If an entry with this name already exists a
#'   warning is emitted and the entry is overwritten. Default is
#'   \code{"nhood_cdr3_convergence"}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} instance
#'   controlling parallelisation. The default,
#'   \code{\link[BiocParallel]{SerialParam}()}, runs sequentially with
#'   per-neighbourhood progress messages. Parallel backends (e.g.
#'   \code{\link[BiocParallel]{MulticoreParam}},
#'   \code{\link[BiocParallel]{SnowParam}}) suppress per-neighbourhood messages
#'   but respect all other arguments. Note that when using process-based
#'   parallel backends, all data objects are serialised to each worker, which
#'   may increase peak memory usage substantially for large datasets.
#'
#' @return The input \code{object} with a \code{data.frame} added to
#'   \code{metadata(object)[[result_name]]}. The data.frame has one row per
#'   neighbourhood and the following columns:
#'   \describe{
#'     \item{\code{nhood}}{Neighbourhood identifier.}
#'     \item{\code{nhood_size}}{Total number of cells in the neighbourhood,
#'       including cells without a contig. Used as the denominator
#'       for \code{contig_capture_rate}.}
#'     \item{\code{<chain>_n_contigs}}{Number of cells with a contig for
#'       this chain.}
#'     \item{\code{<chain>_contig_capture_rate}}{Proportion of neighbourhood
#'       cells with a contig (\code{n_contigs / nhood_size}).
#'       \code{NA} if \code{nhood_size} is zero.}
#'     \item{\code{<chain>_n_clones}}{Number of representative sequences used
#'       for distance computation (unique sequences for contigs without
#'       clone labels; one centroid per clone for BCR with clone labels).}
#'     \item{\code{<chain>_converge_score_reliable}}{Logical; \code{FALSE} when
#'       \code{n_contigs < min_n_contigs} (all metrics \code{NA}), or when
#'       \code{min_clones} is supplied and \code{n_clones} does not meet it
#'       (all metrics \code{NA}). \code{TRUE} for all other neighbourhoods,
#'       including those where \code{converge_score} is \code{NA} due to
#'       insufficient sequences for kNN.}
#'     \item{\code{CDR3_converge_score_<chain>}}{kNN convergence score per
#'       neighbourhood, averaged across all representatives. Captures local
#'       sequence density. When downsampling, this is the mean across bootstrap
#'       iterations. \code{NA} when \code{n_clones <= resolved_k}
#'       (insufficient sequences for kNN), or when \code{min_clones}
#'       is supplied and \code{n_clones} does not meet it.
#'       Only present when \code{"converge_score"} is in \code{summary_metric}.
#'       Example: \code{CDR3_converge_score_VDJ}.}
#'     \item{\code{CDR3_converge_score_sd_<chain>}}{Standard deviation of
#'       \code{CDR3_converge_score_<chain>} across bootstrap iterations.
#'       Only present when \code{downsample = TRUE} and \code{"converge_score"}
#'       is in \code{summary_metric}. \code{NA} when \code{converge_score} is
#'       \code{NA} (kNN could not be computed).
#'       Example: \code{CDR3_converge_score_sd_VDJ}.}
#'     \item{\code{CDR3_converge_score_global_<chain>}}{Mean pairwise CDR3
#'       convergence score across the upper triangle of the representative
#'       distance matrix. When downsampling, this is the mean across bootstrap
#'       iterations. Only present when \code{"converge_score_global"} is in
#'       \code{summary_metric}. Example: \code{CDR3_converge_score_global_VDJ}.}
#'     \item{\code{CDR3_converge_score_global_sd_<chain>}}{Standard deviation of
#'       \code{CDR3_converge_score_global_<chain>} across bootstrap iterations.
#'       Only present when \code{downsample = TRUE} and
#'       \code{"converge_score_global"} is in \code{summary_metric}.
#'       Example: \code{CDR3_converge_score_global_sd_VDJ}.}
#'     \item{\code{CDR3_converge_score_global_median_<chain>}}{Median pairwise
#'       CDR3 convergence score across the upper triangle of the representative
#'       distance matrix. When downsampling, this is the mean of
#'       bootstrap-iteration medians. Only present when
#'       \code{"converge_score_global_median"} is in \code{summary_metric}.
#'       Example: \code{CDR3_converge_score_global_median_VJ}.}
#'     \item{\code{CDR3_converge_score_global_median_sd_<chain>}}{Standard
#'       deviation of \code{CDR3_converge_score_global_median_<chain>} across
#'       bootstrap iterations. Only present when \code{downsample = TRUE} and
#'       \code{"converge_score_global_median"} is in \code{summary_metric}.
#'       Example: \code{CDR3_converge_score_global_median_sd_VDJ}.}
#'   }
#'   Where \code{chain = "both"}, columns are present for both \code{VDJ} and
#'   \code{VJ} chains.
#'
#'   The returned \code{data.frame} also carries the following attributes
#'   recording the analysis parameters: \code{receptor_type},
#'   \code{distance_method}, \code{chains}, \code{downsample_depths},
#'   \code{min_n_contigs}, \code{min_clones}, \code{downsample},
#'   \code{n_boots}, \code{summary_metric}.
#'
#' @details
#' ## Distance computation pipeline
#' For each neighbourhood the following steps are performed:
#' \enumerate{
#'   \item Cells are filtered to those with valid CDR3 annotation.
#'   \item A single pairwise distance matrix is computed over all unique CDR3
#'     sequences. Distances are normalised to \eqn{[0, 1]} by dividing each raw
#'     pairwise distance by the length of the longer sequence in that pair.
#'   \item For BCR data with clone labels, one centroid sequence is selected
#'     per clone as the sequence with the lowest mean pairwise distance to all
#'     other unique sequences within that clone in the neighbourhood. The raw
#'     pairwise distance matrix is filtered to these centroids.
#'   \item If \code{downsample = TRUE} and the number of representatives
#'     exceeds the target depth, bootstrapped subsampling is performed over the
#'     pre-computed distance matrix (no re-alignment), and summary statistics
#'     are averaged across iterations.
#' }
#'
#' ## Supported distance methods
#' The function supports three distance methods: Levenshtein edit distance,
#' BLOSUM62-derived substitution cost, and a user-supplied custom matrix, all
#' normalised to 0-1 by the length of the longer sequence in each
#' pair.
#'
#' ## Downsampling depth
#' Downsampling is recommended where neighbourhood sizes differ substantially
#' to ensure distances are comparable across neighbourhoods.
#' The downsampling depth is the median number of unique representatives across
#' all neighbourhoods that pass \code{min_n_contigs}, excluding singleton-
#' representative neighbourhoods. This is computed separately per chain when
#' \code{chain = "both"}, as VDJ and VJ clone counts can differ. The depth is
#' reported in the progress messages emitted to the console.
#'
#' ## Neighbourhood size and capture rate
#' \code{nhood_size} and \code{contig_capture_rate} are based on all cells
#' assigned to the neighbourhood in the binary nhoods matrix, regardless of
#' contig annotation status. This reflects the overall cellular composition of
#' the neighbourhood rather than immune cell abundance alone, and is the
#' intended denominator for capture rate.
#'
#' ## Parallelisation
#' Each neighbourhood is processed independently and dispatched as a separate
#' task via \code{\link[BiocParallel]{bplapply}}. With
#' \code{\link[BiocParallel]{SerialParam}} (default), progress is reported at
#' every 10\% of neighbourhoods. With parallel backends, per-neighbourhood
#' progress is suppressed; a completion summary is printed after collation.
#'
#' @seealso
#' \href{https://bioconductor.org/packages/release/bioc/manuals/pwalign/man/pwalign.pdf}{\code{pwalign::stringDist()}}
#' for the underlying alignment engine.
#' \code{\link[BiocParallel]{BiocParallelParam}} for parallelisation options.
#' \code{\link{nhood_diversity}}, \code{\link{nhood_shm}},
#' \code{\link{nhood_annotation}}, \code{\link{nhood_annotation_props}}
#'
#' @examples
#' \dontrun{
#' # --- Milo object, BCR, VDJ chain only, Levenshtein distance -------------
#' milo <- nhood_cdr3_convergence(
#'   object        = milo_object,
#'   receptor_type = "BCR",
#'   chain         = "VDJ",
#'   aa_col        = "cdr3_aa_heavy",
#'   clone_id_col  = "clone_id",
#'   seed          = 42
#' )
#' results <- S4Vectors::metadata(milo)$nhood_cdr3_convergence
#'
#' # --- Milo object, TCR, both chains, BLOSUM62, all metrics ----------------
#' milo <- nhood_cdr3_convergence(
#'   object          = milo_object,
#'   receptor_type   = "TCR",
#'   chain           = "both",
#'   aa_col          = c(VDJ = "cdr3_aa_TRB", VJ = "cdr3_aa_TRA"),
#'   distance_method = "blosum62",
#'   summary_metric  = c("converge_score", "converge_score_global", "converge_score_global_median"),
#'   min_n_contigs     = 5,
#'   seed            = 42,
#'   BPPARAM         = BiocParallel::MulticoreParam(workers = 4)
#' )
#'
#' # --- SingleCellExperiment object, BCR, custom substitution matrix --------
#' custom_mat <- matrix(...)  # 20x20 numeric, named by single-letter AA code
#' sce <- nhood_cdr3_convergence(
#'   object              = sce_object,
#'   receptor_type       = "BCR",
#'   chain               = "VDJ",
#'   aa_col              = "cdr3_aa_heavy",
#'   clone_id_col        = "clone_id",
#'   nhood_col           = "nhood_id",
#'   distance_method     = "custom",
#'   substitution_matrix = custom_mat,
#'   seed                = 42
#' )
#' }
#'
#'
#' @export
#'
nhood_cdr3_convergence <- function(object,
                                   nhood_col = NULL,
                                   receptor_type = c("BCR", "TCR"),
                                   chain = c("VDJ", "VJ", "both"),
                                   aa_col,
                                   clone_id_col = NULL,
                                   min_n_contigs = 20L,
                                   min_clones = NULL,
                                   distance_method = c("levenshtein", "blosum62", "custom"),
                                   substitution_matrix = NULL,
                                   summary_metric = "converge_score",
                                   knn_k = NULL,
                                   downsample = TRUE,
                                   n_boots = 100,
                                   seed = NULL,
                                   BPPARAM = BiocParallel::SerialParam(),
                                   result_name = "nhood_cdr3_convergence") {

  # ── 0. Dependencies ──────────────────────────────────────────────────────────
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("Package 'BiocParallel' is required. Install with: BiocManager::install('BiocParallel')")
  }
  if (!requireNamespace("pwalign", quietly = TRUE)) {
    stop("Package 'pwalign' is required. Install with: BiocManager::install('pwalign')")
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required. Install with: BiocManager::install('Biostrings')")
  }

  # ── 1. Input validation ───────────────────────────────────────────────────
  receptor_type   <- match.arg(receptor_type)
  chain           <- match.arg(chain)
  distance_method <- match.arg(distance_method)

  valid_metrics <- c("converge_score", "converge_score_global", "converge_score_global_median")
  if (!is.character(summary_metric) || length(summary_metric) == 0L ||
      !all(summary_metric %in% valid_metrics)) {
    stop(sprintf(
      "'summary_metric' must be a character vector containing one or more of: %s.\n  Received: %s",
      paste(valid_metrics, collapse = ", "),
      paste(summary_metric, collapse = ", ")
    ))
  }
  summary_metric <- unique(summary_metric)

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed))) {
    stop("'seed' must be a single finite numeric value or NULL.")
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

  # Validate aa_col against chain selection
  if (chain == "both") {
    if (!is.character(aa_col) || is.null(names(aa_col)) ||
        !all(c("VDJ", "VJ") %in% names(aa_col))) {
      stop("When chain = 'both', aa_col must be a named character vector with elements 'VDJ' and 'VJ'.\n",
           "  Example: aa_col = c(VDJ = 'cdr3_aa_VDJ', VJ = 'cdr3_aa_VJ')")
    }
    chains_to_run <- c("VDJ", "VJ")
    aa_col_map    <- aa_col                    # named vector, e.g. c(VDJ="col1", VJ="col2")
  } else {
    if (!is.character(aa_col) || length(aa_col) != 1) {
      stop("When chain = '", chain, "', aa_col must be a single column name string.")
    }
    chains_to_run <- chain
    aa_col_map    <- stats::setNames(aa_col, chain)
  }

  # ── 2. Object type determination ──────────────────────────────────────────
  is_milo     <- inherits(object, "Milo") && inherits(object, "SingleCellExperiment")
  is_sce_only <- !is_milo && inherits(object, "SingleCellExperiment")

  if (!is_milo && !is_sce_only) {
    stop("'object' must be a Milo or SingleCellExperiment object. If using Seurat, convert via as.SingleCellExperiment")
  }

  # ── 3. Build substitution cost matrix ────────────────────────────────────
  # If distance_method = custom, the caller-supplied matrix is used directly.
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

  score_mat <- .build_substitution_matrix(distance_method, substitution_matrix)

  # ── 4. Extract colData and neighbourhood matrix ───────────────────────────
  message("Extracting neighbourhood membership...")
  col_data <- SummarizedExperiment::colData(object)

  # Validate aa columns exist in colData
  missing_cols <- setdiff(aa_col_map, colnames(col_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Column(s) not found in colData: %s\nAvailable columns: %s",
                 paste(missing_cols, collapse = ", "),
                 paste(colnames(col_data), collapse = ", ")))
  }

  if (!is.null(clone_id_col) && !clone_id_col %in% colnames(col_data)) {
    warning(sprintf("clone_id_col '%s' not found in colData. Falling back to unique CDR3 sequences as clones.",
                    clone_id_col))
    clone_id_col <- NULL
  }

  # Extract neighbourhood matrix
  if (is_milo) {
    if (!requireNamespace("miloR", quietly = TRUE)) {
      stop("Package 'miloR' is required for Milo objects. Install with: BiocManager::install('miloR')")
    }
    if (is.null(miloR::nhoods(object))) {
      stop("nhoods slot not found in Milo object")
    }
    nhoods        <- miloR::nhoods(object)
    unique_nhoods <- colnames(nhoods)

    # Warn when the user has also supplied nhood_col
    if (!is.null(nhood_col) && nrow(nhoods) > 0L) {
      message(
        "'nhood_col' is ignored for Milo objects whose nhoods matrix is ",
        "already populated. Neighbourhood membership is read directly from ",
        "nhoods(object). To use 'nhood_col', supply a plain SingleCellExperiment."
      )
    }

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
        "rownames(nhoods(object)) or colnames(object) are NULL. ",
        "Assuming cell order in nhoods matrix matches colData row order."
      )
    }

  } else {
    if (is.null(nhood_col)) {
      stop("For SingleCellExperiment objects, 'nhood_col' must be specified")
    }
    if (!nhood_col %in% colnames(col_data)) {
      stop(sprintf("Column '%s' not found in colData. Available columns: %s",
                   nhood_col,
                   paste(colnames(col_data), collapse = ", ")))
    }

    nhood_assignments <- col_data[[nhood_col]]
    valid_cells       <- !is.na(nhood_assignments)

    if (any(!valid_cells)) {
      warning(sprintf("%d cells have NA neighbourhood assignments and will be excluded",
                      sum(!valid_cells)))
    }

    unique_nhoods_raw <- unique(nhood_assignments[valid_cells])
    # Sort numerically if all labels are numeric, otherwise lexicographically.
    numeric_ids <- suppressWarnings(as.numeric(unique_nhoods_raw))
    if (!anyNA(numeric_ids)) {
      unique_nhoods <- as.character(sort(numeric_ids))
    } else {
      unique_nhoods <- as.character(sort(unique_nhoods_raw))
    }
    n_nhoods_total <- length(unique_nhoods)

    nhoods <- matrix(0,
                     nrow = ncol(object),
                     ncol = n_nhoods_total,
                     dimnames = list(colnames(object), unique_nhoods))

    for (idx in seq_along(unique_nhoods)) {
      # Compare against the original nhood_assignments values
      cells_in <- as.character(nhood_assignments) == unique_nhoods[idx] & valid_cells
      nhoods[cells_in, idx] <- 1
    }

    if (requireNamespace("Matrix", quietly = TRUE)) {
      nhoods <- Matrix::Matrix(nhoods, sparse = TRUE)
    }
  }

  n_nhoods      <- ncol(nhoods)
  unique_nhoods <- colnames(nhoods)

  # ── 5. Pre-filter colData to valid cells only, per chain ──────────────────
  # Retain only cells that have both a valid aa sequence and a non-NA clone ID
  # (when clone_id_col is provided) for each chain.
  message("Filtering to cells with complete contig annotation...")

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

  # Named integer lookup: cell barcode -> row index in col_data_filtered[[ch]].
  cell_row_index <- lapply(chains_to_run, function(ch) {
    cd <- col_data_filtered[[ch]]
    stats::setNames(seq_len(nrow(cd)), rownames(cd))
  })
  names(cell_row_index) <- chains_to_run

  # Precompute neighbourhood membership once as a list of integer index vectors.
  nhood_cell_indices <- lapply(seq_len(n_nhoods), function(i) which(nhoods[, i] == 1))

  # nhood_size is a property of the full neighbourhood, computed from the
  # unfiltered nhoods matrix
  nhood_sizes <- lengths(nhood_cell_indices)

  # ── 6. Pre-pass: count valid cells and unique clones per neighbourhood ─────
  # Both counts are computed in a single pass per chain using the pre-filtered
  # colData and precomputed membership indices
  message("Pre-pass: counting clones per neighbourhood...")

  pre_pass <- lapply(chains_to_run, function(ch) {
    aa_col_name    <- aa_col_map[[ch]]
    cd             <- col_data_filtered[[ch]]
    idx_lookup     <- cell_row_index[[ch]]
    nhood_barcodes <- rownames(nhoods)   # all cell barcodes in nhoods row order

    .count_nhood_stats(
      n_nhoods           = n_nhoods,
      nhood_cell_indices = nhood_cell_indices,
      nhood_barcodes     = nhood_barcodes,
      cd                 = cd,
      idx_lookup         = idx_lookup,
      aa_col_name        = aa_col_name,
      clone_id_col       = clone_id_col
    )
  })
  names(pre_pass) <- chains_to_run

  bcr_tcr_counts       <- lapply(pre_pass, `[[`, "n_valid_cells")
  rep_counts_per_chain <- lapply(pre_pass, `[[`, "n_unique_clones")

  # Downsampling depth: median unique clone count across neighbourhoods that
  # pass the quality gate, excluding singleton-clone neighbourhoods.
  # Computed per chain as clone counts can differ between VDJ and VJ.
  downsample_depths <- lapply(chains_to_run, function(ch) {
    valid_rep_counts <- rep_counts_per_chain[[ch]][bcr_tcr_counts[[ch]] >= min_n_contigs]
    if (length(valid_rep_counts) > 0 && any(valid_rep_counts > 1)) {
      as.integer(round(stats::median(valid_rep_counts[valid_rep_counts > 1], na.rm = TRUE)))
    } else {
      2L
    }
  })
  names(downsample_depths) <- chains_to_run

  for (ch in chains_to_run) {
    n_valid_nhoods <- sum(bcr_tcr_counts[[ch]] >= min_n_contigs)
    n_skip_nhoods  <- n_nhoods - n_valid_nhoods
    message(sprintf(
      "  %s: %d / %d neighbourhoods will be processed (downsampling depth = %d representatives) | %d skipped (< %d cells with contig annotation)",
      ch, n_valid_nhoods, n_nhoods, downsample_depths[[ch]], n_skip_nhoods, min_n_contigs
    ))
    min_k_for_chain <- if (!is.null(knn_k)) knn_k else 3L
    if (downsample && "converge_score" %in% summary_metric &&
        downsample_depths[[ch]] <= min_k_for_chain) {
      warning(sprintf(paste0(
        "Chain %s: downsampling depth (%d representatives) is <= k (%d). ",
        "converge_score will be NA for all downsampled neighbourhoods since kNN ",
        "cannot be computed when n_clones <= k. ",
        "To obtain converge_score values, either set downsample = FALSE or use ",
        "summary_metric = 'converge_score_global'."
      ), ch, downsample_depths[[ch]], min_k_for_chain))
    }
  }

  # ── 7. Initialise results dataframe ──────────────────────────────────────
  # nhood_size counts all cells assigned to this neighbourhood (i.e. all cells
  # in the nhoods binary matrix with a 1 entry), regardless of contig
  # annotation status. contig_capture_rate is therefore expressed as a
  # fraction of total neighbourhood size, not of immune-annotated cells only.
  distance_results <- data.frame(
    nhood            = unique_nhoods,
    nhood_size       = nhood_sizes,
    stringsAsFactors = FALSE
  )

  # ── 8. Main loop over neighbourhoods (parallelised via BiocParallel) ──────
  # Progress messages are emitted per neighbourhood when using SerialParam.
  # With parallel backends, per-neighbourhood messages are suppressed since
  # workers run in separate processes; the summary printed after will confirm
  # completion.
  is_serial    <- inherits(BPPARAM, "SerialParam")
  n_to_process <- max(vapply(chains_to_run, function(ch) {
    sum(bcr_tcr_counts[[ch]] >= min_n_contigs)
  }, integer(1)))
  message(sprintf(
    "Computing CDR3 distances for %d valid neighbourhoods (of %d total)%s...",
    n_to_process, n_nhoods,
    if (is_serial) "" else sprintf(" using %d workers", BiocParallel::bpworkers(BPPARAM))
  ))

  nhood_results <- BiocParallel::bplapply(
    seq_len(n_nhoods),
    function(i) {
      if (is_serial && i %% max(1L, n_nhoods %/% 10L) == 0L) {
        message(sprintf("  Processing neighbourhood %d / %d (%.0f%%)...",
                        i, n_nhoods, 100 * i / n_nhoods))
      }
      .process_single_nhood(
        i                              = i,
        nhood_cell_indices             = nhood_cell_indices,
        nhood_barcodes                 = rownames(nhoods),
        nhood_sizes                    = nhood_sizes,
        col_data_filtered              = col_data_filtered,
        cell_row_index                 = cell_row_index,
        aa_col_map                     = aa_col_map,
        chains_to_run                  = chains_to_run,
        clone_id_col                   = clone_id_col,
        receptor_type                  = receptor_type,
        score_mat                      = score_mat,
        min_n_contigs                    = min_n_contigs,
        min_clones            = min_clones,
        downsample                     = downsample,
        n_boots                        = n_boots,
        downsample_depths              = downsample_depths,
        knn_k                          = knn_k,
        seed                           = seed,
        summary_metric                 = summary_metric
      )
    },
    BPPARAM = BPPARAM
  )

  # ── 9. Collate results into dataframe ─────────────────────────────────────
  # Pre-allocate one vector per output column and fill by neighbourhood index,
  # then assign each completed vector to the dataframe in a single operation.
  message("Collating results...")

  for (ch in chains_to_run) {
    col_n_contigs                 <- integer(n_nhoods)
    col_capture_rate              <- numeric(n_nhoods)
    col_n_clones                  <- integer(n_nhoods)
    col_converge_score_reliable      <- logical(n_nhoods)
    col_converge_score               <- if ("converge_score" %in% summary_metric) numeric(n_nhoods) else NULL
    col_converge_score_sd            <- if (downsample && "converge_score" %in% summary_metric) numeric(n_nhoods) else NULL
    col_converge_score_global        <- if ("converge_score_global" %in% summary_metric) numeric(n_nhoods) else NULL
    col_converge_score_global_sd     <- if (downsample && "converge_score_global" %in% summary_metric) numeric(n_nhoods) else NULL
    col_converge_score_global_median <- if ("converge_score_global_median" %in% summary_metric) numeric(n_nhoods) else NULL
    col_converge_score_global_median_sd <- if (downsample && "converge_score_global_median" %in% summary_metric) numeric(n_nhoods) else NULL

    for (i in seq_len(n_nhoods)) {
      res                  <- nhood_results[[i]][[ch]]
      col_n_contigs[i]     <- res$n_contigs
      col_capture_rate[i]  <- res$contig_capture_rate
      col_n_clones[i]      <- res$n_clones
      col_converge_score_reliable[i] <- res$converge_score_reliable
      if ("converge_score" %in% summary_metric) col_converge_score[i] <- res$converge_score
      if (downsample && "converge_score" %in% summary_metric) col_converge_score_sd[i]  <- res$converge_score_sd
      if ("converge_score_global" %in% summary_metric)  col_converge_score_global[i] <- res$converge_score_global
      if (downsample && "converge_score_global" %in% summary_metric) col_converge_score_global_sd[i] <- res$converge_score_global_sd
      if ("converge_score_global_median" %in% summary_metric) col_converge_score_global_median[i] <- res$converge_score_global_median
      if (downsample && "converge_score_global_median" %in% summary_metric) col_converge_score_global_median_sd[i] <- res$converge_score_global_median_sd
    }

    # Convert distance scores to similarity (convergence) scores via 1 - distance
    if ("converge_score" %in% summary_metric) col_converge_score  <- 1 - col_converge_score
    if ("converge_score_global" %in% summary_metric) col_converge_score_global <- 1 - col_converge_score_global
    if ("converge_score_global_median" %in% summary_metric) col_converge_score_global_median <- 1 - col_converge_score_global_median

    # Assign to results frame
    distance_results[[paste0(ch, "_n_contigs")]]           <- col_n_contigs
    distance_results[[paste0(ch, "_contig_capture_rate")]] <- col_capture_rate
    distance_results[[paste0(ch, "_n_clones")]]            <- col_n_clones
    distance_results[[paste0(ch, "_converge_score_reliable")]] <- col_converge_score_reliable
    if ("converge_score" %in% summary_metric) {
      distance_results[[paste0("CDR3_converge_score_", ch)]]    <- col_converge_score
      if (downsample) distance_results[[paste0("CDR3_converge_score_sd_", ch)]] <- col_converge_score_sd
    }
    if ("converge_score_global" %in% summary_metric) {
      distance_results[[paste0("CDR3_converge_score_global_", ch)]]    <- col_converge_score_global
      if (downsample) distance_results[[paste0("CDR3_converge_score_global_sd_", ch)]] <- col_converge_score_global_sd
    }
    if ("converge_score_global_median" %in% summary_metric) {
      distance_results[[paste0("CDR3_converge_score_global_median_", ch)]]    <- col_converge_score_global_median
      if (downsample) distance_results[[paste0("CDR3_converge_score_global_median_sd_", ch)]] <- col_converge_score_global_median_sd
    }
  }

  # ── 10. Summary ───────────────────────────────────────────────────────────
  message("\nCDR3 convergence score computation complete.")
  message(sprintf("  Receptor type  : %s", receptor_type))
  message(sprintf("  Distance method: %s", distance_method))
  message(sprintf("  Chain(s)       : %s", paste(chains_to_run, collapse = ", ")))
  message(sprintf("  Summary metric(s): %s", paste(summary_metric, collapse = ", ")))
  for (ch in chains_to_run) {
    n_reliable <- sum(distance_results[[paste0(ch, "_converge_score_reliable")]], na.rm = TRUE)
    if ("converge_score" %in% summary_metric) {
      n_knn_na <- sum(
        is.na(distance_results[[paste0("CDR3_converge_score_", ch)]]) &
          distance_results[[paste0(ch, "_converge_score_reliable")]],
        na.rm = TRUE
      )
      message(sprintf(
        "  %s: %d / %d neighbourhoods with reliable convergence estimates (%d with NA converge_score due to n_clones <= k)",
        ch, n_reliable, n_nhoods, n_knn_na
      ))
    } else {
      message(sprintf("  %s: %d / %d neighbourhoods with reliable convergence estimates",
                      ch, n_reliable, n_nhoods))
    }
  }

  # ── 11. Store attributes and write to metadata ────────────────────────────
  attr(distance_results, "receptor_type")       <- receptor_type
  attr(distance_results, "distance_method")      <- distance_method
  attr(distance_results, "chains")               <- chains_to_run
  attr(distance_results, "downsample_depths")    <- downsample_depths
  attr(distance_results, "min_n_contigs")          <- min_n_contigs
  attr(distance_results, "min_clones")  <- min_clones
  attr(distance_results, "downsample")           <- downsample
  attr(distance_results, "n_boots")              <- n_boots
  attr(distance_results, "summary_metric")       <- summary_metric

  if (!is.null(S4Vectors::metadata(object)[[result_name]])) {
    warning(sprintf("Overwriting existing metadata entry '%s'.", result_name))
  }

  S4Vectors::metadata(object)[[result_name]] <- distance_results

  message(sprintf(
    "Stored CDR3 convergence results (%d neighbourhoods x %d columns) in metadata(object)[['%s']].",
    nrow(distance_results), ncol(distance_results), result_name
  ))

  return(object)
}

#' Build a substitution score matrix for pairwise CDR3 alignment
#'
#' Constructs a 20x20 amino acid substitution score matrix suitable for use
#' with \code{pwalign::stringDist}. The matrix is determined by
#' \code{distance_method}: a flat Levenshtein-style cost matrix, a
#' BLOSUM62-derived cost matrix normalised to \eqn{[0, 1]}, or a
#' caller-supplied custom matrix.
#'
#' @param distance_method Character string; one of \code{"levenshtein"},
#'   \code{"blosum62"}, or \code{"custom"}.
#' @param substitution_matrix A numeric matrix; required only when
#'   \code{distance_method = "custom"}, otherwise ignored.
#'
#' @return A named numeric matrix of dimensions 20x20, with single-letter
#'   amino acid codes as row and column names. The \code{distance_method}
#'   attribute is set on the returned matrix for downstream dispatch. The
#'   sign convention follows \code{pwalign::stringDist}: higher values
#'   indicate more similar amino acids. Diagonal entries are 0 (exact match
#'   incurs no cost).
#'
#'   For \code{"levenshtein"}: off-diagonal entries are \code{-1}.
#'
#'   For \code{"blosum62"}: entries are the negated, normalised BLOSUM62 cost
#'   matrix, where cost is computed as
#'   \eqn{c(i,j) = (s(i,i) + s(j,j)) / 2 - s(i,j)} and then divided by the
#'   maximum off-diagonal cost to scale to \eqn{[0, 1]} before negation.
#'
#'   For \code{"custom"}: the caller-supplied matrix is returned unchanged
#'   except for the addition of the \code{distance_method} attribute.
#'
#' @noRd


.build_substitution_matrix <- function(distance_method, substitution_matrix = NULL) {

  aa_alphabet <- c("A","C","D","E","F","G","H","I","K","L",
                   "M","N","P","Q","R","S","T","V","W","Y")

  if (distance_method == "levenshtein") {
    # Flat score matrix: 0 for identical amino acids, -1 for any substitution.
    # Used with pwalign::stringDist (substitutionMatrix argument).
    sub_mat       <- matrix(-1, nrow = 20L, ncol = 20L,
                            dimnames = list(aa_alphabet, aa_alphabet))
    diag(sub_mat) <- 0L
    attr(sub_mat, "distance_method") <- "levenshtein"
    return(sub_mat)
  }

  if (distance_method == "blosum62") {
    # BLOSUM62 cost via formula: c(i,j) = s(i,i)/2 + s(j,j)/2 - s(i,j)
    # Normalised to [0, 1] by dividing by the maximum off-diagonal cost.
    # Used with pwalign::stringDist (substitutionMatrix argument).
    # Load BLOSUM62 into an isolated environment to avoid polluting the local
    # frame. The name is verified explicitly.
    e <- new.env(parent = emptyenv())
    utils::data("BLOSUM62", package = "pwalign", envir = e)
    if (is.null(e$BLOSUM62)) {
      stop("Could not load BLOSUM62 from pwalign - the dataset may have been renamed in this version of the package.")
    }
    blosum62_raw <- e$BLOSUM62
    blosum62     <- blosum62_raw[aa_alphabet, aa_alphabet]
    self_scores  <- diag(blosum62)

    cost_mat       <- outer(self_scores, self_scores, FUN = `+`) / 2 - blosum62
    d_max          <- max(cost_mat)
    cost_mat       <- cost_mat / d_max
    diag(cost_mat) <- 0   # guard against floating-point drift on diagonal

    # Negate to convert cost matrix to score matrix for pwalign convention
    score_mat <- -cost_mat
    attr(score_mat, "distance_method") <- "blosum62"
    return(score_mat)
  }

  if (distance_method == "custom") {
    # Caller is responsible for constructing a valid cost matrix.
    # Attach the method attribute for dispatch in .run_pairwise_alignment.
    attr(substitution_matrix, "distance_method") <- "custom"
    return(substitution_matrix)
  }
}

#' Count contig-annotated cells and unique clones per neighbourhood
#'
#' A single-pass pre-computation step that counts, for each neighbourhood,
#' the number of cells with valid contig annotation and the number of unique
#' clones (or unique sequences, when no clone column is supplied). Results are
#' used to determine the downsampling depth and to skip under-represented
#' neighbourhoods before the main alignment loop.
#'
#' Uses a pre-built barcode-to-row-index hash (\code{idx_lookup}) to subset
#' pre-filtered \code{colData} in O(neighbourhood_size) rather than
#' O(n_cells), avoiding repeated full-object subsetting.
#'
#' @param n_nhoods Integer; total number of neighbourhoods.
#' @param nhood_cell_indices List of integer vectors; element \code{i} gives
#'   the row indices (into the nhoods matrix) of cells belonging to
#'   neighbourhood \code{i}.
#' @param nhood_barcodes Character vector; cell barcodes in nhoods matrix row
#'   order (\code{rownames(nhoods)}).
#' @param cd A \code{DataFrame} or \code{data.frame}; pre-filtered
#'   \code{colData} containing only cells with valid CDR3 annotation for the
#'   chain being processed.
#' @param idx_lookup Named integer vector; maps cell barcodes to row indices
#'   in \code{cd}.
#' @param aa_col_name Character string; name of the CDR3 amino acid column in
#'   \code{cd}, used to count unique sequences when \code{clone_id_col} is
#'   \code{NULL}.
#' @param clone_id_col Character string or \code{NULL}; name of the clone
#'   identifier column in \code{cd}. When \code{NULL}, unique sequences are
#'   counted instead.
#'
#' @return A named list with two integer vectors, each of length
#'   \code{n_nhoods}:
#'   \describe{
#'     \item{\code{n_valid_cells}}{Number of cells with valid contig
#'       annotation per neighbourhood.}
#'     \item{\code{n_unique_clones}}{Number of unique clones (or unique CDR3
#'       sequences) per neighbourhood.}
#'   }
#'
#' @noRd
.count_nhood_stats <- function(n_nhoods, nhood_cell_indices, nhood_barcodes,
                               cd, idx_lookup, aa_col_name, clone_id_col) {

  n_valid_cells   <- integer(n_nhoods)
  n_unique_clones <- integer(n_nhoods)

  for (i in seq_len(n_nhoods)) {
    barcodes   <- nhood_barcodes[nhood_cell_indices[[i]]]
    rows       <- idx_lookup[barcodes]
    rows       <- rows[!is.na(rows)]
    nhood_data <- cd[rows, , drop = FALSE]

    n_valid_cells[i] <- nrow(nhood_data)

    n_unique_clones[i] <- if (
      !is.null(clone_id_col) && clone_id_col %in% colnames(nhood_data)
    ) {
      length(unique(nhood_data[[clone_id_col]]))
    } else {
      length(unique(nhood_data[[aa_col_name]]))
    }
  }

  list(n_valid_cells = n_valid_cells, n_unique_clones = n_unique_clones)
}

#' Process a single neighbourhood: alignment, centroid selection, and
#' distance summarisation
#'
#' Performs the full per-neighbourhood CDR3 distance analysis for all
#' requested chains. Called once per neighbourhood via
#' \code{\link[BiocParallel]{bplapply}}.
#'
#' The function:
#' \enumerate{
#'   \item Subsets pre-filtered \code{colData} to neighbourhood cells using a
#'     hash lookup.
#'   \item Computes one pairwise alignment over all unique CDR3 sequences in
#'     the neighbourhood.
#'   \item For BCR data with clone labels, selects a centroid sequence per
#'     clone (lowest mean distance to clone members) from the pre-computed
#'     matrix.
#'   \item Applies bootstrapped downsampling if required, operating on index
#'     subsets of the pre-computed matrix (no re-alignment).
#'   \item Returns per-chain summary statistics.
#' }
#'
#' @param i Integer; 1-based index of the neighbourhood being processed.
#' @param nhood_cell_indices List of integer vectors; precomputed cell
#'   membership indices per neighbourhood (output of \code{which(nhoods[,i] == 1)}).
#' @param nhood_barcodes Character vector; cell barcodes in nhoods matrix row
#'   order.
#' @param nhood_sizes Integer vector; total number of cells per neighbourhood,
#'   used as the denominator for \code{contig_capture_rate}.
#' @param col_data_filtered Named list of pre-filtered \code{colData}
#'   \code{DataFrame}s, one per chain.
#' @param cell_row_index Named list of barcode-to-row-index lookup vectors,
#'   one per chain.
#' @param aa_col_map Named character vector mapping chain name to CDR3 amino
#'   acid column name.
#' @param chains_to_run Character vector; chains to process (\code{"VDJ"},
#'   \code{"VJ"}, or both).
#' @param clone_id_col Character string or \code{NULL}; clone identifier
#'   column name.
#' @param receptor_type Character string; \code{"BCR"} or \code{"TCR"}.
#' @param score_mat Numeric matrix; substitution score matrix from
#'   \code{.build_substitution_matrix}.
#' @param min_n_contigs Integer; minimum annotated cells required to process a
#'   neighbourhood. When not met, all metrics are \code{NA} and
#'   \code{converge_score_reliable = FALSE}.
#' @param min_clones Integer or \code{NULL}; when supplied, minimum
#'   representative sequences required to process a neighbourhood. When not
#'   met, all metrics are \code{NA} and \code{converge_score_reliable = FALSE}.
#' @param downsample Logical; whether bootstrapped downsampling is applied.
#' @param n_boots Integer; number of bootstrap iterations.
#' @param downsample_depths Named integer list; target representative count per
#'   chain, precomputed from the global median.
#' @param knn_k Integer or \code{NULL}; user-supplied fixed \eqn{k}. When
#'   \code{NULL}, adaptive k is used with floor 3
#' @param seed Integer or \code{NULL}; user-supplied base seed.
#'   Neighbourhood-specific seed is \code{seed + i}.
#' @param summary_metric Character vector; which convergence metrics to compute.
#'   One or more of \code{"converge_score"}, \code{"converge_score_global"},
#'   \code{"converge_score_global_median"}.
#'
#' @return A named list with one element per chain in \code{chains_to_run}.
#'   Each element is a named list with the following fields:
#'   \describe{
#'     \item{\code{n_contigs}}{Integer; number of cells with a contig.}
#'     \item{\code{contig_capture_rate}}{Numeric; fraction of neighbourhood
#'       cells with a contig. \code{NA} if \code{nhood_size} is 0.}
#'     \item{\code{n_clones}}{Integer; number of representative sequences.}
#'     \item{\code{converge_score_reliable}}{Logical; \code{FALSE} when
#'       \code{n_contigs < min_n_contigs} or when \code{min_clones} is
#'       supplied and not met.}
#'     \item{\code{converge_score}}{Numeric; kNN convergence score. \code{NA}
#'       when \code{n_clones <= resolved_k} or when \code{min_clones} is
#'       not met.}
#'     \item{\code{converge_score_sd}}{Numeric; SD of \code{converge_score} across
#'       bootstrap iterations. \code{NA} if not downsampled or kNN could not
#'       be computed.}
#'     \item{\code{converge_score_global}}{Numeric; mean pairwise convergence
#'       score (\code{NA} if not requested).}
#'     \item{\code{converge_score_global_sd}}{Numeric; SD of
#'       \code{converge_score_global} across bootstrap iterations. \code{NA} if
#'       not downsampled.}
#'     \item{\code{converge_score_global_median}}{Numeric; median pairwise
#'       convergence score (\code{NA} if not requested).}
#'     \item{\code{converge_score_global_median_sd}}{Numeric; SD of
#'       \code{converge_score_global_median} across bootstrap iterations.
#'       \code{NA} if not downsampled.}
#'   }
#'
#' @noRd
.process_single_nhood <- function(i, nhood_cell_indices, nhood_barcodes,
                                  nhood_sizes, col_data_filtered, cell_row_index,
                                  aa_col_map, chains_to_run, clone_id_col,
                                  receptor_type, score_mat, min_n_contigs,
                                  min_clones,
                                  downsample, n_boots,
                                  downsample_depths, knn_k, seed,
                                  summary_metric) {

  nhood_size     <- nhood_sizes[[i]]
  # Barcodes of all cells in this neighbourhood
  nhood_barcodes_i <- nhood_barcodes[nhood_cell_indices[[i]]]

  result <- list()

  for (ch in chains_to_run) {
    aa_col_name      <- aa_col_map[[ch]]
    downsample_depth <- downsample_depths[[ch]]
    idx_lookup       <- cell_row_index[[ch]]
    cd               <- col_data_filtered[[ch]]

    # Subset pre-filtered colData to cells in this neighbourhood
    rows       <- idx_lookup[nhood_barcodes_i]
    rows       <- rows[!is.na(rows)]
    nhood_data <- cd[rows, , drop = FALSE]
    n_contigs  <- nrow(nhood_data)

    # Insufficient contigs for this chain - return NA
    if (n_contigs < min_n_contigs) {
      # Count unique clones directly from the filtered subset
      n_clones <- if (!is.null(clone_id_col) && clone_id_col %in% colnames(nhood_data)) {
        length(unique(nhood_data[[clone_id_col]]))
      } else {
        length(unique(nhood_data[[aa_col_name]]))
      }
      result[[ch]] <- list(
        n_contigs                    = n_contigs,
        contig_capture_rate          = if (nhood_size > 0L) n_contigs / nhood_size else NA_real_,
        n_clones                     = n_clones,
        converge_score_reliable         = FALSE,
        converge_score                  = NA_real_,
        converge_score_sd               = NA_real_,
        converge_score_global           = NA_real_,
        converge_score_global_sd        = NA_real_,
        converge_score_global_median    = NA_real_,
        converge_score_global_median_sd = NA_real_
      )
      next
    }

    # Extract aa sequences
    aa_seqs <- nhood_data[[aa_col_name]]

    # Build one distance matrix over all unique sequences in this neighbourhood
    unique_seqs   <- unique(aa_seqs)
    full_dist_mat <- .run_pairwise_alignment(unique_seqs, score_mat)

    # Identify centroid indices and derive representative submatrix
    if (receptor_type == "BCR" &&
        !is.null(clone_id_col) &&
        clone_id_col %in% colnames(nhood_data)) {

      # BCR: for each clone, subset the neighbourhood-wide matrix to clone
      # members and pick the centroid as the sequence with the lowest mean
      # distance to all other members
      clone_ids     <- nhood_data[[clone_id_col]]
      unique_clones <- unique(clone_ids)

      centroid_indices <- vapply(unique_clones, function(cid) {
        member_seqs <- unique(aa_seqs[clone_ids == cid])

        if (length(member_seqs) == 1L) {
          # Singleton or all-identical clone: centroid is the only sequence
          return(match(member_seqs, unique_seqs))
        }

        # Subset the pre-computed neighbourhood matrix to clone members only
        clone_idx  <- match(member_seqs, unique_seqs)
        clone_sub  <- full_dist_mat[clone_idx, clone_idx, drop = FALSE]
        mean_dists <- rowMeans(clone_sub, na.rm = TRUE)
        clone_idx[which.min(mean_dists)]

      }, FUN.VALUE = integer(1))

      # Subset the neighbourhood matrix to centroid rows/cols only.
      rep_dist_mat <- full_dist_mat[centroid_indices, centroid_indices,
                                    drop = FALSE]

    } else {
      # TCR (or BCR without clone_id_col): use unique sequences
      rep_dist_mat <- full_dist_mat
      centroid_indices <- seq_len(nrow(full_dist_mat))
    }

    n_reps <- nrow(rep_dist_mat)

    # Optional minimum representative gate. When supplied and not met, all
    # metrics are NA and converge_score_reliable = FALSE
    if (!is.null(min_clones) && n_reps < min_clones) {
      result[[ch]] <- list(
        n_contigs                    = n_contigs,
        contig_capture_rate          = if (nhood_size > 0L) n_contigs / nhood_size else NA_real_,
        n_clones                     = n_reps,
        converge_score_reliable         = FALSE,
        converge_score                  = NA_real_,
        converge_score_sd               = NA_real_,
        converge_score_global           = NA_real_,
        converge_score_global_sd        = NA_real_,
        converge_score_global_median    = NA_real_,
        converge_score_global_median_sd = NA_real_
      )
      next
    }

    # Single representative: distance is 0
    if (n_reps == 1L) {
      result[[ch]] <- list(
        n_contigs                    = n_contigs,
        contig_capture_rate          = if (nhood_size > 0L) n_contigs / nhood_size else NA_real_,
        n_clones                     = 1L,
        converge_score_reliable         = TRUE,
        converge_score                  = if ("converge_score" %in% summary_metric)               0 else NA_real_,
        converge_score_sd               = NA_real_,
        converge_score_global           = if ("converge_score_global" %in% summary_metric)        0 else NA_real_,
        converge_score_global_sd        = NA_real_,
        converge_score_global_median    = if ("converge_score_global_median" %in% summary_metric) 0 else NA_real_,
        converge_score_global_median_sd = NA_real_
      )
      next
    }

    # Resolve knn_k for this neighbourhood.
    effective_n <- if (downsample && n_reps > downsample_depth) downsample_depth else n_reps
    resolved_k <- if (!is.null(knn_k)) {
      knn_k
    } else {
      max(3L, floor(0.25 * effective_n))
    }

    # converge_score requires n_clones > k to form a valid kNN graph.
    # When this is not met, converge_score is NA.
    can_compute_knn <- effective_n > resolved_k

    if (!downsample || n_reps <= downsample_depth) {
      # No downsampling needed, summarise centroid matrix directly.
      dist_summary <- .summarise_dist_mat(rep_dist_mat, n_reps = n_reps, knn_k = resolved_k,
                                          summary_metric = summary_metric)

      result[[ch]] <- list(
        n_contigs                    = n_contigs,
        contig_capture_rate          = if (nhood_size > 0L) n_contigs / nhood_size else NA_real_,
        n_clones                     = n_reps,
        converge_score_reliable         = TRUE,
        converge_score                  = if ("converge_score" %in% summary_metric) {
          if (can_compute_knn) dist_summary$converge_score else NA_real_
        } else NA_real_,
        converge_score_sd               = NA_real_,
        converge_score_global           = if ("converge_score_global" %in% summary_metric)        dist_summary$mean   else NA_real_,
        converge_score_global_sd        = NA_real_,
        converge_score_global_median    = if ("converge_score_global_median" %in% summary_metric) dist_summary$median else NA_real_,
        converge_score_global_median_sd = NA_real_
      )
    } else {
      # Bootstrapped downsampling: subsample row/col indices into rep_dist_mat
      if (!is.null(seed)) set.seed(seed + i)
      boot_converge_score               <- if ("converge_score" %in% summary_metric && can_compute_knn) numeric(n_boots) else NULL
      boot_converge_score_global        <- if ("converge_score_global" %in% summary_metric)             numeric(n_boots) else NULL
      boot_converge_score_global_median <- if ("converge_score_global_median" %in% summary_metric)      numeric(n_boots) else NULL

      for (b in seq_len(n_boots)) {
        idx          <- sample(n_reps, size = downsample_depth, replace = FALSE)
        sub_mat      <- rep_dist_mat[idx, idx, drop = FALSE]
        dist_summary <- .summarise_dist_mat(sub_mat, n_reps = downsample_depth, knn_k = resolved_k,
                                            summary_metric = summary_metric)
        if ("converge_score" %in% summary_metric && can_compute_knn)    boot_converge_score[b]               <- dist_summary$converge_score
        if ("converge_score_global" %in% summary_metric)                boot_converge_score_global[b]        <- dist_summary$mean
        if ("converge_score_global_median" %in% summary_metric)         boot_converge_score_global_median[b] <- dist_summary$median
      }

      result[[ch]] <- list(
        n_contigs                    = n_contigs,
        contig_capture_rate          = if (nhood_size > 0L) n_contigs / nhood_size else NA_real_,
        n_clones                     = n_reps,
        converge_score_reliable         = TRUE,
        converge_score                  = if ("converge_score" %in% summary_metric) {
          if (can_compute_knn) mean(boot_converge_score, na.rm = TRUE) else NA_real_
        } else NA_real_,
        converge_score_sd               = if ("converge_score" %in% summary_metric && can_compute_knn) stats::sd(boot_converge_score,               na.rm = TRUE) else NA_real_,
        converge_score_global           = if ("converge_score_global" %in% summary_metric)        mean(boot_converge_score_global,        na.rm = TRUE) else NA_real_,
        converge_score_global_sd        = if ("converge_score_global" %in% summary_metric)        stats::sd(boot_converge_score_global,          na.rm = TRUE) else NA_real_,
        converge_score_global_median    = if ("converge_score_global_median" %in% summary_metric) mean(boot_converge_score_global_median, na.rm = TRUE) else NA_real_,
        converge_score_global_median_sd = if ("converge_score_global_median" %in% summary_metric) stats::sd(boot_converge_score_global_median,  na.rm = TRUE) else NA_real_
      )
    }
  }

  return(result)
}

#' Compute a normalised pairwise CDR3 distance matrix
#'
#' Computes all pairwise amino acid sequence distances for a set of CDR3
#' sequences using \code{\link[pwalign]{stringDist}}, then normalises each
#' raw distance by the length of the longer sequence in that pair, yielding
#' values in \eqn{[0, 1]}.
#'
#' All three distance methods (Levenshtein, BLOSUM62, custom) are routed
#' through this function using the substitution matrix produced by
#' \code{.build_substitution_matrix}.
#'
#' @param sequences Character vector of unique CDR3 amino acid sequences.
#' @param score_mat Numeric matrix; the substitution score matrix returned by
#'   \code{.build_substitution_matrix}.
#' @param gap_open Numeric; gap opening penalty passed to
#'   \code{pwalign::stringDist}. Default \code{0}.
#' @param gap_extend Numeric; gap extension penalty passed to
#'   \code{pwalign::stringDist}. Default \code{1}.
#'
#' @return A symmetric numeric matrix of dimensions
#'   \eqn{n \times n} (where \eqn{n = \code{length(sequences)}}), with
#'   sequences as both row and column names. Values are in \eqn{[0, 1]};
#'   the diagonal is 0.
#'
#' @noRd
.run_pairwise_alignment <- function(sequences, score_mat,
                                    gap_open = 0, gap_extend = 1) {
  n           <- length(sequences)
  seq_lengths <- nchar(sequences)

  seq_set <- Biostrings::AAStringSet(sequences)

  dist_obj <- pwalign::stringDist(
    seq_set,
    method             = "substitutionMatrix",
    substitutionMatrix = score_mat,
    gapOpening         = gap_open,
    gapExtension       = gap_extend,
    diag = TRUE,
    upper = TRUE
  )

  raw_mat            <- as.matrix(dist_obj)
  dimnames(raw_mat)  <- list(sequences, sequences)

  # Normalise by length of longer sequence in each pair
  len_mat        <- outer(seq_lengths, seq_lengths, FUN = pmax)
  norm_mat       <- raw_mat / len_mat
  diag(norm_mat) <- 0

  norm_mat
}

#' Summarise a pre-computed pairwise distance matrix
#'
#' Computes converge_score_global (mean), converge_score_global_median (median),
#' and/or converge_score (kNN convergence score) from the upper triangle of a
#' symmetric pairwise distance matrix. Called on the full representative matrix
#' when downsampling is not applied, and on each bootstrap subsample when it is.
#' The matrix is never recomputed here; only index-level operations are
#' performed. Only the metrics requested via \code{summary_metric} are computed.
#'
#' @param dist_mat Numeric matrix; symmetric pairwise distance matrix with
#'   values in \eqn{[0, 1]} and zeros on the diagonal.
#' @param n_reps Integer; number of rows/columns in \code{dist_mat} (i.e.
#'   the number of representatives).
#' @param knn_k Integer or \code{NULL}; number of nearest neighbours for
#'   converge_score computation.
#' @param summary_metric Character vector; which metrics to compute. One or
#'   more of \code{"converge_score"}, \code{"converge_score_global"},
#'   \code{"converge_score_global_median"}.
#'
#' @return A named list with up to three elements (those not requested are
#'   \code{NULL}):
#'   \describe{
#'     \item{\code{mean}}{Mean of upper-triangle distances (converge_score_global).}
#'     \item{\code{median}}{Median of upper-triangle distances (converge_score_global_median).}
#'     \item{\code{converge_score}}{kNN convergence score from
#'       \code{.compute_converge_score}.}
#'   }
#'
#' @noRd
.summarise_dist_mat <- function(dist_mat, n_reps, knn_k,
                                summary_metric = c("converge_score", "converge_score_global", "converge_score_global_median")) {

  need_upper <- any(c("converge_score_global", "converge_score_global_median") %in% summary_metric)
  upper_vals <- if (need_upper) {
    uv <- dist_mat[upper.tri(dist_mat)]
    uv[!is.na(uv)]
  } else {
    NULL
  }

  list(
    mean       = if ("converge_score_global" %in% summary_metric)        mean(upper_vals,                na.rm = TRUE) else NULL,
    median     = if ("converge_score_global_median" %in% summary_metric) stats::median(upper_vals,              na.rm = TRUE) else NULL,
    converge_score = if ("converge_score" %in% summary_metric)              .compute_converge_score(dist_mat, k = knn_k) else NULL
  )
}

#' Compute converge_score via k-nearest neighbours
#'
#' For each sequence (row), identifies its \eqn{k} nearest neighbours among
#' all other sequences by taking the \eqn{k} smallest distances in that row
#' (excluding the self-comparison on the diagonal). Computes the mean distance
#' to those neighbours per sequence, then averages across all sequences.
#'
#' This captures local sequence density: a neighbourhood of tightly clustered
#' sequences will have a low converge_score even if some outliers inflate the
#' global convergence score. Ties at position \eqn{k} are broken by index
#' order. \eqn{k} is capped at \eqn{n - 1} to guard against small bootstrap
#' subsamples.
#'
#' @param dist_mat Numeric matrix; symmetric pairwise distance matrix.
#' @param k Integer; number of nearest neighbours. Capped internally at
#'   \eqn{n - 1}.
#'
#' @return A single numeric value: the mean across all sequences of their mean
#'   distance to their \eqn{k} nearest neighbours, rounded to 4 decimal
#'   places.
#'
#' @noRd
.compute_converge_score <- function(dist_mat, k) {
  n <- nrow(dist_mat)

  # Cap k so it never exceeds the number of available neighbours
  k <- min(k, n - 1L)

  per_seq_means <- vapply(seq_len(n), function(i) {
    # All distances from sequence i, excluding self (diagonal)
    dists <- dist_mat[i, -i]
    # k nearest neighbours: partial sort, take first k (ties broken by index)
    knn_dists <- sort(dists, partial = k)[seq_len(k)]
    mean(knn_dists, na.rm = TRUE)
  }, FUN.VALUE = numeric(1))

  round(mean(per_seq_means, na.rm = TRUE), 4)
}
