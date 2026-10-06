#' Scatter plot of neighbourhood-level metrics
#'
#' Plots two neighbourhood-level metrics against each other, with point size
#' optionally scaled by a neighbourhood-level column and colour mapped to a
#' grouping variable.
#'
#' @param nhood_data      Data frame containing all columns to plot.
#' @param x_col           Column name (string) for the x-axis metric
#'                        (default "shannon_diversity").
#' @param y_col           Column name (string) for the y-axis metric
#'                        (default "mean_shm").
#' @param x_reliability_col Column name of a logical reliability flag for the
#'                        x-axis metric. Neighbourhoods where this column is
#'                        FALSE or NA are dropped before plotting, regardless
#'                        of their y-axis reliability. NULL skips this filter.
#' @param y_reliability_col Column name of a logical reliability flag for the
#'                        y-axis metric. Neighbourhoods where this column is
#'                        FALSE or NA are dropped before plotting, regardless
#'                        of their x-axis reliability. NULL skips this filter.
#' @param colour_by      Column name (string) for point colour
#'                        (default "anat").
#' @param size_by         Column name (string) to map point size to, or NULL
#'                        (default) for a fixed point size.
#' @param size            Controls point size and depends on \code{size_by}:
#'                        \itemize{
#'                          \item When \code{size_by} is \code{NULL} — a
#'                            single numeric giving the fixed point size
#'                            (default \code{2} when \code{size} is left
#'                            \code{NULL}).
#'                          \item When \code{size_by} is a column name — a
#'                            numeric vector of length 2 giving the minimum
#'                            and maximum point size of the size scale
#'                            (default \code{c(1, 4)} when \code{size} is
#'                            left \code{NULL}).
#'                        }
#' @param alpha           Numeric between 0-1. Point transparency (default 0.8).
#' @param facet_by        Column name to facet by, or NULL for no faceting.
#' @param facet_ncol      Number of columns for facet_wrap (NULL = auto).
#' @param facet_scales    Axis scale argument passed to facet_wrap:
#'                        "free" (default), "fixed", "free_x", or "free_y".
#' @param add_density     Logical. If TRUE, adds 2D density contours
#'                        (default FALSE).
#' @param x_lab         Optional custom x-axis label. NULL auto-derives from
#'                        x_col.
#' @param y_lab         Optional custom y-axis label. NULL auto-derives from
#'                        y_col.
#'
#' @return A ggplot object.
#'
#' @seealso \code{\link{nhood_diversity}}, \code{\link{nhood_shm}},
#' \code{\link{nhood_cdr3_convergence}}, \code{\link{plot_nhood_beeswarm}}
#'
#' @examples
#' \dontrun{
#' # Fixed-size points
#' plot_nhood_scatter(
#'   nhood_data = nhood_data,
#'   x_col      = "shannon_diversity",
#'   y_col      = "mean_shm",
#'   colour_by  = "anat"
#' )
#'
#' # Point size mapped to neighbourhood size, with a custom size range
#' plot_nhood_scatter(
#'   nhood_data = nhood_data,
#'   x_col      = "shannon_diversity",
#'   y_col      = "mean_shm",
#'   colour_by  = "anat",
#'   size_by    = "nhood_size",
#'   size       = c(2, 8)
#' )
#' }
#'
#' @importFrom rlang .data
#'
#' @export
#'
plot_nhood_scatter <- function(nhood_data,
                               x_col            = NULL,
                               y_col            = NULL,
                               x_reliability_col= NULL,
                               y_reliability_col= NULL,
                               colour_by        = NULL,
                               size_by          = NULL,
                               size             = NULL,
                               alpha            = 0.8,
                               facet_by         = NULL,
                               facet_ncol       = NULL,
                               facet_scales     = "free",
                               add_density      = FALSE,
                               x_lab            = NULL,
                               y_lab            = NULL) {

  # ── 1. Input validation ───────────────────────────────────────────────────────
  required_cols <- c(x_col, y_col, colour_by)
  missing_cols  <- setdiff(required_cols, colnames(nhood_data))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.null(facet_by) && !facet_by %in% colnames(nhood_data)) {
    stop("Facet column '", facet_by, "' not found in data")
  }

  if (!is.null(x_reliability_col) && !x_reliability_col %in% colnames(nhood_data))
    stop(sprintf("x_reliability_col '%s' not found in 'nhood_data'.", x_reliability_col))

  if (!is.null(y_reliability_col) && !y_reliability_col %in% colnames(nhood_data))
    stop(sprintf("y_reliability_col '%s' not found in 'nhood_data'.", y_reliability_col))

  # Resolve size_by/size into the two internal flags used throughout:
  #   size_col   : column name to map to size aes (NULL when fixed)
  #   fixed_size : numeric size passed to geom directly (NULL when col mode)
  #   size_range : length-2 numeric for scale_size_continuous (NULL when fixed)
  size_info  <- .resolve_nhood_size(nhood_data, size_by, size)
  size_col   <- size_info$size_col
  fixed_size <- size_info$fixed_size
  size_range <- size_info$size_range

  # ── 2. Reliability filtering ──────────────────────────────────────────────────
  # A neighbourhood is dropped if *either* reliability flag is FALSE or NA.
  # This is done on the full nhood_data before column selection so that
  # reliability columns do not need to be retained in the plot data.
  if (!is.null(x_reliability_col)) {
    nhood_data <- nhood_data[
      !is.na(nhood_data[[x_reliability_col]]) & nhood_data[[x_reliability_col]] == TRUE, ]
  }
  if (!is.null(y_reliability_col)) {
    nhood_data <- nhood_data[
      !is.na(nhood_data[[y_reliability_col]]) & nhood_data[[y_reliability_col]] == TRUE, ]
  }

  # ── 3. Prepare data ───────────────────────────────────────────────────────────
  filter_cols <- c(x_col, y_col, colour_by, size_col)
  if (!is.null(facet_by)) {
    filter_cols <- c(filter_cols, facet_by)
  }

  plot_data <- nhood_data |>
    dplyr::select(dplyr::all_of(filter_cols)) |>
    dplyr::filter(!is.na(.data[[x_col]]) & !is.na(.data[[y_col]]))

  if (!is.null(facet_by)) {
    plot_data <- dplyr::filter(plot_data, !is.na(.data[[facet_by]]))
  }

  # ── 4. Detect colour scale type ───────────────────────────────────────────────
  is_continuous <- is.numeric(plot_data[[colour_by]]) &&
    !is.factor(plot_data[[colour_by]])

  # ── 5. Axis labels ────────────────────────────────────────────────────────────
  fmt_label    <- function(col) gsub("_", " ", tools::toTitleCase(col))
  x_axis_label <- if (!is.null(x_lab)) x_lab else fmt_label(x_col)
  y_axis_label <- if (!is.null(y_lab)) y_lab else fmt_label(y_col)

  # ── 6. Base plot ──────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x     = .data[[x_col]],
                                    y     = .data[[y_col]],
                                    color = .data[[colour_by]]))

  if (!is.null(size_col)) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(size = .data[[size_col]]), alpha = alpha) +
      ggplot2::scale_size_continuous(range = size_range,
                                     name  = fmt_label(size_col))
  } else {
    p <- p + ggplot2::geom_point(size = fixed_size, alpha = alpha)
  }

  p <- p + ggplot2::labs(x = x_axis_label, y = y_axis_label)

  # ── 7. Colour scale ───────────────────────────────────────────────────────────
  if (is_continuous) {
    p <- p +
      ggplot2::scale_color_viridis_c(
        option   = "viridis",
        na.value = "grey80",
        name     = fmt_label(colour_by)
      )
  } else {
    p <- p +
      ggplot2::scale_color_discrete(name = fmt_label(colour_by))
  }

  # ── 8. Density contours ───────────────────────────────────────────────────────
  if (add_density) {
    p <- p +
      ggplot2::geom_density_2d(
        ggplot2::aes(x = .data[[x_col]],
                     y = .data[[y_col]]),
        color        = "black",
        alpha        = 0.3,
        linewidth    = 0.5,
        inherit.aes  = FALSE,
        data         = plot_data
      )
  }

  # ── 9. Faceting ───────────────────────────────────────────────────────────────
  if (!is.null(facet_by)) {
    p <- p + ggplot2::facet_wrap(
      ggplot2::vars(!!rlang::sym(facet_by)),
      ncol   = facet_ncol,
      scales = facet_scales
    )
  }

  # ── 10. Theme ─────────────────────────────────────────────────────────────────
  p <- p +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle    = ggplot2::element_text(size = 11, color = "gray30"),
      legend.position  = "right",
      legend.box       = "vertical",
      strip.background = ggplot2::element_rect(fill = "white"),
      strip.text       = ggplot2::element_text(face = "bold")
    )

  return(p)
}
