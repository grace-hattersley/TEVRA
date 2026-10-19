#' Plot mean community size across tested resolutions
#'
#' Visualises the output of \code{\link{test_resolution}} as a line plot of
#' mean cells per community against resolution, with optional error bars.
#' Use this plot to identify a resolution that yields communities of a
#' suitable and stable size before running \code{\link{create_communities}}.
#'
#' @param resolution_results A \code{data.frame} as returned by
#'   \code{\link{test_resolution}}, containing at minimum the columns
#'   \code{resolution} and \code{mean_cells_per_community}. If
#'   \code{show_error_bars = TRUE}, the column \code{sd_cells_per_community}
#'   must also be present.
#' @param show_error_bars A logical value. If \code{TRUE} (the default),
#'   error bars showing ±1 SD (\code{sd_cells_per_community}) are
#'   added around each point.
#'
#' @return A \code{ggplot} object showing mean cells per community against
#'   resolution.
#'
#' @seealso \code{\link{test_resolution}} to generate the input
#' \code{data.frame}; \code{\link{create_communities}} to assign cells to
#' communities once a resolution has been chosen.
#'
#' @examples
#' \dontrun{
#' # Test a range of resolutions
#' res_results <- test_resolution(sce, resolutions = c(0.5, 1, 2.5, 5))
#'
#' # Plot the results with error bars
#' plot_resolution_results(res_results)
#'
#' # Plot without error bars
#' plot_resolution_results(res_results, show_error_bars = FALSE)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_errorbar labs
#'                     theme_classic theme element_text
#' @importFrom rlang .data
#'
#' @export
#'
plot_resolution_results <- function(resolution_results, show_error_bars = TRUE) {

  # ── 1. Extract mean and SD columns ─────────────────────────────────────────
  metric <- "mean_cells_per_community"
  error_col <- "sd_cells_per_community"
  y_label <- "Mean Cells per Community"

  # Check that resolution_results contains the necessary columns
  required_cols <- if (show_error_bars) c(metric, error_col) else metric
  missing_cols <- setdiff(required_cols, colnames(resolution_results))
  if (length(missing_cols) > 0) {
    stop(paste0(
      "resolution_results is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      ". Run test_resolution() first."
    ))
  }

  # ── 2. Create base plot ────────────────────────────────────────────────────
  p <- ggplot(resolution_results, aes(x = .data$resolution, y = .data[[metric]])) +
    geom_line(color = "#2E86AB", linewidth = 1) +
    geom_point(color = "#2E86AB", size = 3) +
    labs(
      title = paste(y_label, "vs Resolution"),
      x = "Resolution Parameter",
      y = y_label
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )

  # ── 3. Add error bars (optional) ───────────────────────────────────────────
  if (show_error_bars) {
    p <- p + geom_errorbar(
      aes(
        ymin = .data[[metric]] - .data[[error_col]],
        ymax = .data[[metric]] + .data[[error_col]]
      ),
      width = 0.02,
      color = "#2E86AB",
      alpha = 0.7
    )
  }

  return(p)
}
