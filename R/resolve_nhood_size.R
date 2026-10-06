#' Resolve consistent point-size aesthetics for plot_nhood_* functions
#'
#' Internal helper shared by \code{plot_nhood_beeswarm}, \code{plot_nhood_scatter}
#' and \code{plot_nhood_umap} to resolve the \code{size_by} / \code{size}
#' parameter pair into a single consistent internal representation.
#'
#' Contract:
#' \itemize{
#'   \item \code{size_by = NULL} (fixed-size mode) — \code{size} must be
#'     \code{NULL} (use \code{default_fixed_size}) or a single numeric value.
#'   \item \code{size_by = "col_name"} (mapped mode) — \code{size} must be
#'     \code{NULL} (use \code{default_size_range}) or a numeric vector of
#'     length 2 giving the min/max of the size scale.
#' }
#'
#' @param data      Data frame in which \code{size_by} (if given) is looked up.
#' @param size_by   Column name (string) to map point size to, or NULL for a
#'                  fixed size.
#' @param size      NULL, a single numeric (fixed size; only when
#'                  \code{size_by} is NULL), or a length-2 numeric (size
#'                  range; only when \code{size_by} is set).
#' @param default_fixed_size Default fixed size used when \code{size_by} is
#'                  NULL and \code{size} is NULL.
#' @param default_size_range Default size range used when \code{size_by} is
#'                  set and \code{size} is NULL.
#' @param arg_name  Name used for \code{size_by} in error messages
#'
#' @return A list with \code{size_col} (column name or NULL),
#'   \code{fixed_size} (numeric or NULL) and \code{size_range} (length-2
#'   numeric or NULL). Exactly one of \code{fixed_size} / \code{size_range}
#'   is non-NULL.
#'
#' @noRd
#'
.resolve_nhood_size <- function(data,
                                size_by,
                                size,
                                default_fixed_size = 2,
                                default_size_range = c(1, 4),
                                arg_name = "size_by") {

  if (!is.null(size_by)) {
    if (!is.character(size_by) || length(size_by) != 1)
      stop(sprintf("'%s' must be a single column name (string) or NULL.", arg_name))
    if (!size_by %in% colnames(data))
      stop(sprintf(
        "'%s' column '%s' not found in data. Available columns: %s",
        arg_name, size_by, paste(colnames(data), collapse = ", ")
      ))
  }

  if (is.null(size_by)) {
    # ── Fixed-size mode ──────────────────────────────────────────────────
    if (is.null(size)) {
      fixed_size <- default_fixed_size
    } else if (is.numeric(size) && length(size) == 1) {
      fixed_size <- size
    } else {
      stop(sprintf(
        "'size' must be a single numeric value (fixed point size) when '%s' is NULL.",
        arg_name
      ))
    }
    list(size_col = NULL, fixed_size = fixed_size, size_range = NULL)

  } else {
    # ── Mapped mode ──────────────────────────────────────────────────────
    if (is.null(size)) {
      size_range <- default_size_range
    } else if (is.numeric(size) && length(size) == 2) {
      size_range <- size
    } else {
      stop(sprintf(
        "'size' must be a numeric vector of length 2 (min, max) when '%s' is set.",
        arg_name
      ))
    }
    list(size_col = size_by, fixed_size = NULL, size_range = size_range)
  }
}
