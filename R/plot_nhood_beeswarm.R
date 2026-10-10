#' Beeswarm plot for neighbourhood-level metrics
#'
#' Accepts a data frame with neighbourhoods as rows.
#'
#' Multi-column mode (plot_cols): when more than one column is supplied,
#' the data are pivoted to long format and each column becomes a row of panels.
#'
#' @param nhood_data    Data frame containing all columns to plot.
#' @param plot_cols     Character vector of column name(s) to use as the x-axis.
#'                      Supply one name for a single-metric plot, or multiple
#'                      names to get one panel row per metric (multi-column mode).
#' @param reliability_col  Column name of a logical reliability flag. Rows where
#'                      this column is FALSE or NA are dropped before plotting.
#'                      NULL skips reliability filtering entirely.
#' @param group_by      Column name (string) whose levels form the y-axis categories.
#' @param summary_metric Controls whether a summary line is drawn, which statistic
#'                      is used, and whether it is computed globally or per group.
#'                      One of:
#'                      \itemize{
#'                        \item \code{"none"}          – no summary line is drawn.
#'                        \item \code{"global_mean"}   – (default) one dashed vertical
#'                          line per panel at the grand mean across all
#'                          \code{group_by} levels.
#'                        \item \code{"global_median"} – as above but using the median.
#'                        \item \code{"group_mean"}    – one short vertical segment per
#'                          \code{group_by} level at the per-group mean. Segments span
#'                          ±0.4 y-units so they sit within each beeswarm row.
#'                        \item \code{"group_median"}  – as above but using the median.
#'                      }
#' @param colour_by     Column name (string) for point colour. NULL defaults to
#'                      group_by. Set to the value of LFC_col to activate
#'                      p-value-filtered LFC colouring.
#' @param p_val_col     Column name of p-values (used in LFC mode).
#' @param LFC_col       Column name of log-fold-changes (used in LFC mode).
#' @param sig_threshold Significance threshold; points with p < sig_threshold
#'                      are coloured by LFC in LFC mode.
#' @param size_by       Column name (string) to map point size to, or NULL
#'                      (default) for a fixed point size.
#' @param size          Controls point size and depends on \code{size_by}:
#'                      \itemize{
#'                        \item When \code{size_by} is \code{NULL} — a single
#'                          numeric giving the fixed point size (default
#'                          \code{2} when \code{size} is left \code{NULL}).
#'                        \item When \code{size_by} is a column name — a
#'                          numeric vector of length 2 giving the minimum and
#'                          maximum point size of the size scale (default
#'                          \code{c(1, 4)} when \code{size} is left \code{NULL}).
#'                      }
#' @param alpha         Numeric between 0-1. Point transparency (default 0.8).
#' @param facet_by      Column name to facet by, or NULL. In multi-column mode
#'                      this produces a facet_grid (metric rows x facet_by cols);
#'                      in single-column mode it produces a facet_wrap.
#' @param facet_ncol    Number of columns for facet_wrap (NULL = auto).
#'                      Ignored when facet_grid is used.
#' @param facet_scales  Axis scale argument passed to facet_wrap / facet_grid:
#'                      "free" (default), "fixed", "free_x", or "free_y".
#' @param x_limits      Numeric vector of length 2 for x-axis limits, or NULL
#'                      for automatic limits. Example: c(0, 1) for proportions.
#' @param x_lab         X-axis label string. NULL auto-derives from plot_cols:
#'                      a single column name is title-cased; multiple columns
#'                      produce "Value".
#' @param y_lab         Y-axis label string. NULL auto-derives from group_by.
#' @param title         Plot title string.
#' @param subtitle      Plot subtitle string, or NULL for none.

#'
#' @return A ggplot object.
#'
#' @seealso \code{\link{nhood_diversity}}, \code{\link{nhood_shm}},
#' \code{\link{nhood_cdr3_convergence}}, \code{\link{plot_nhood_scatter}}
#'
#' @examples
#' \dontrun{
#' # Single-metric plot
#' plot_nhood_beeswarm(
#'   nhood_data      = nhood_data,
#'   plot_cols       = "shannon_diversity",
#'   group_by        = "anat",
#'   reliability_col = "diversity_reliable",
#'   subtitle        = paste0("Reliable estimates (>=", attr(nhood_data, "min_bcr_tcr"), " BCR/TCR)")
#' )
#'
#' # Multi-column plot
#' plot_nhood_beeswarm(
#'   nhood_data = ccall_df,
#'   plot_cols  = c("IgA", "IgG", "IgM"),
#'   group_by   = "anat",
#'   x_limits   = c(0, 1),
#'   x_lab      = "Proportion of Cells"
#' )
#'
#' # Point size mapped to a column, with a custom size range
#' plot_nhood_beeswarm(
#'   nhood_data = nhood_data,
#'   plot_cols  = "shannon_diversity",
#'   group_by   = "anat",
#'   size_by    = "nhood_size",
#'   size       = c(1, 6)
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes scale_y_discrete scale_x_continuous
#'                     geom_vline geom_segment scale_colour_gradient2 scale_colour_discrete
#'                     scale_colour_viridis_c scale_size_continuous facet_grid facet_wrap vars labs
#'                     theme_classic theme element_text
#' @importFrom rlang sym .data
#'
#' @export
#'
plot_nhood_beeswarm <- function(nhood_data,
                                plot_cols,
                                reliability_col = NULL,
                                group_by,
                                summary_metric  = c("global_mean", "none", "global_median",
                                                    "group_mean", "group_median"),
                                colour_by       = NULL,
                                p_val_col       = "p_val",
                                LFC_col         = "LFC",
                                sig_threshold   = 0.1,
                                size_by         = NULL,
                                size            = NULL,
                                alpha           = 0.8,
                                facet_by        = NULL,
                                facet_ncol      = NULL,
                                facet_scales    = c("free", "fixed", "free_x", "free_y"),
                                x_limits        = NULL,
                                x_lab           = NULL,
                                y_lab           = NULL,
                                title           = "",
                                subtitle        = NULL) {

  # ── 0. Dependencies ───────────────────────────────────────────────────────────
  if (!requireNamespace("ggbeeswarm", quietly = TRUE))
    stop("Package 'ggbeeswarm' required. Install with: install.packages('ggbeeswarm')")

  # ── 1. Basic input validation ─────────────────────────────────────────────────
  summary_metric <- match.arg(summary_metric)
  facet_scales   <- match.arg(facet_scales)

  if (!is.data.frame(nhood_data))
    stop("'nhood_data' must be a data frame.")

  missing_cols <- setdiff(
    c(plot_cols, group_by, facet_by, colour_by, reliability_col),
    c(colnames(nhood_data), NULL)
  )
  if (length(missing_cols) > 0)
    stop("The following columns were not found in 'nhood_data': ",
         paste(missing_cols, collapse = ", "))

  # Resolve size_by/size into the two internal flags used throughout:
  #   size_col   : column name to map to size aes (NULL when fixed)
  #   fixed_size : numeric size passed to geom directly (NULL when col mode)
  #   size_range : length-2 numeric for scale_size_continuous (NULL when fixed)
  size_info  <- .resolve_nhood_size(nhood_data, size_by, size)
  size_col   <- size_info$size_col
  fixed_size <- size_info$fixed_size
  size_range <- size_info$size_range

  # ── 2. Reliability filtering ──────────────────────────────────────────────────
  if (!is.null(reliability_col)) {
    if (!reliability_col %in% colnames(nhood_data))
      stop(sprintf("reliability_col '%s' not found in 'nhood_data'.", reliability_col))
    # Capture factor levels before subsetting — subsetting drops unused levels
    # by default, which destroys any user-supplied ordering on group_by.
    orig_group_levels <- if (is.factor(nhood_data[[group_by]])) levels(nhood_data[[group_by]]) else NULL
    nhood_data <- nhood_data[!is.na(nhood_data[[reliability_col]]) & nhood_data[[reliability_col]] == TRUE, ]
    # Re-apply the original levels so the y-axis order is preserved.
    if (!is.null(orig_group_levels))
      nhood_data[[group_by]] <- factor(nhood_data[[group_by]], levels = orig_group_levels)
  }

  # ── 3. Coerce plot columns to numeric ─────────────────────────────────────────
  for (col in plot_cols) {
    nhood_data[[col]] <- as.numeric(nhood_data[[col]])
  }

  # ── 4. Determine mode and resolve colour ──────────────────────────────────────
  multi_col_mode <- length(plot_cols) > 1
  colour_col      <- if (is.null(colour_by)) group_by else colour_by
  colour_by_LFC   <- !is.null(colour_by) && colour_by == LFC_col

  if (colour_by_LFC) {
    if (is.null(p_val_col))
      stop("'p_val_col' must be specified when colour_by == LFC_col.")
    if (is.null(LFC_col))
      stop("'LFC_col' must be specified when colour_by == LFC_col.")
    if (!p_val_col %in% colnames(nhood_data))
      stop(sprintf("p-value column '%s' not found in 'nhood_data'. Check 'p_val_col'.",
                   p_val_col))
    if (!LFC_col %in% colnames(nhood_data))
      stop(sprintf("LFC column '%s' not found in 'nhood_data'. Check 'LFC_col'.",
                   LFC_col))
  } else if (!is.null(colour_by) && !colour_by %in% colnames(nhood_data)) {
    stop(sprintf("'colour_by' column '%s' not found in 'nhood_data'.", colour_by))
  }

  # ── 5. Factor conversions ─────────────────────────────────────────────────────
  if (!is.factor(nhood_data[[group_by]])) {
    nhood_data[[group_by]] <- factor(nhood_data[[group_by]],
                                     levels = sort(unique(as.character(nhood_data[[group_by]]))))
  }

  if (!is.null(facet_by) && !is.factor(nhood_data[[facet_by]]))
    nhood_data[[facet_by]] <- factor(nhood_data[[facet_by]],
                                     levels = sort(unique(as.character(nhood_data[[facet_by]]))))

  is_discrete_colour <- !colour_by_LFC &&
    (is.character(nhood_data[[colour_col]]) || is.factor(nhood_data[[colour_col]]))

  if (is_discrete_colour && !is.factor(nhood_data[[colour_col]]))
    nhood_data[[colour_col]] <- factor(nhood_data[[colour_col]],
                                       levels = sort(unique(as.character(nhood_data[[colour_col]]))))

  # ── 6. LFC mode: flag significance and derive colour variable ─────────────────
  if (colour_by_LFC) {
    nhood_data$is_significant <- !is.na(nhood_data[[p_val_col]]) &
      nhood_data[[p_val_col]] < sig_threshold
    nhood_data$colour_var <- ifelse(nhood_data$is_significant, nhood_data[[LFC_col]], NA_real_)
  }

  # ── 7. Multi-column mode: pivot to long format ────────────────────────────────
  if (multi_col_mode) {
    meta_cols <- unique(c(group_by, facet_by, colour_col,
                          if (colour_by_LFC) c(p_val_col, LFC_col,
                                               "is_significant", "colour_var"),
                          size_col))
    meta_cols <- intersect(meta_cols, colnames(nhood_data))

    plot_df <- do.call(rbind, lapply(plot_cols, function(col) {
      sub <- nhood_data[, unique(c(meta_cols, col)), drop = FALSE]
      sub$.metric       <- col
      sub$.metric_value <- sub[[col]]
      sub[[col]]        <- NULL
      sub
    }))

    # Preserve user-supplied column order as factor levels for panel ordering
    plot_df$.metric <- factor(plot_df$.metric, levels = plot_cols)
    x_col <- ".metric_value"
  } else {
    plot_df <- nhood_data
    x_col   <- plot_cols
  }

  # ── 8. Compute summary line data (controlled by summary_metric) ───────────────
  #
  #  Parsed into two orthogonal decisions:
  #    scope  : "global" → one marker per panel (ignores group_by)
  #             "group"  → one short segment per group_by row (×panel ×facet)
  #    stat_fn: mean or median
  #
  #  "none" → summary_layers stays an empty list; nothing is drawn.

  summary_layers <- list()

  if (summary_metric != "none") {

    # Derive scope and stat function from the combined option string
    use_group_scope <- grepl("^group_", summary_metric)
    stat_fn         <- if (grepl("_median$", summary_metric)) stats::median else mean

    # Panel-level grouping columns (shared by both scopes)
    panel_cols <- if (multi_col_mode) {
      if (is.null(facet_by)) ".metric" else c(".metric", facet_by)
    } else {
      if (is.null(facet_by)) character(0) else facet_by
    }

    if (!use_group_scope) {
      # ── Global scope: one dashed vline per panel ────────────────────────────
      if (length(panel_cols) == 0) {
        vline_data <- data.frame(
          xintercept = stat_fn(plot_df[[x_col]], na.rm = TRUE)
        )
      } else {
        vline_data <- stats::aggregate(
          plot_df[[x_col]],
          by  = lapply(panel_cols, function(col) plot_df[[col]]),
          FUN = function(x) stat_fn(x, na.rm = TRUE)
        )
        colnames(vline_data) <- c(panel_cols, "xintercept")
      }

      summary_layers <- list(
        geom_vline(
          data     = vline_data,
          mapping  = aes(xintercept = .data$xintercept),
          linetype = "dashed", colour = "grey50", alpha = 0.7
        )
      )

    } else {
      # ── Group scope: short segment per group_by level (×panel ×facet) ──────
      group_cols <- c(panel_cols, group_by)

      seg_data <- stats::aggregate(
        plot_df[[x_col]],
        by  = lapply(group_cols, function(col) plot_df[[col]]),
        FUN = function(x) stat_fn(x, na.rm = TRUE)
      )
      colnames(seg_data) <- c(group_cols, "xintercept")

      # Map group_by levels to numeric y positions (ggplot2 discrete scale:
      # levels numbered 1..n bottom to top) so segments sit within each row.
      grp_levels        <- levels(plot_df[[group_by]])
      grp_y_map         <- stats::setNames(seq_along(grp_levels), grp_levels)
      seg_data$y_centre <- grp_y_map[as.character(seg_data[[group_by]])]
      seg_data$ymin     <- seg_data$y_centre - 0.4
      seg_data$ymax     <- seg_data$y_centre + 0.4

      summary_layers <- list(
        geom_segment(
          data    = seg_data,
          mapping = aes(x    = .data$xintercept, xend = .data$xintercept,
                        y    = .data$ymin,       yend = .data$ymax),
          colour   = "grey30", linewidth = 0.8, linetype = "solid"
        )
      )
    }
  }

  # ── 9. Base plot ──────────────────────────────────────────────────────────────
  p <- ggplot()

  # Lock the y-axis order to the factor levels of group_by
  p <- p + scale_y_discrete(limits = levels(plot_df[[group_by]]))

  # ── 10. Beeswarm layers ───────────────────────────────────────────────────────
  if (colour_by_LFC) {
    nonsig <- plot_df[!plot_df$is_significant, ]
    sig    <- plot_df[plot_df$is_significant, ]

    # Preserve factor levels in subsets
    nonsig[[group_by]] <- factor(nonsig[[group_by]], levels = levels(plot_df[[group_by]]))
    sig[[group_by]]    <- factor(sig[[group_by]],    levels = levels(plot_df[[group_by]]))

    if (!is.null(size_col)) {
      p <- p +
        ggbeeswarm::geom_quasirandom(
          data    = nonsig,
          mapping = aes(x    = !!sym(x_col),
                        y    = !!sym(group_by),
                        size = !!sym(size_col)),
          alpha = alpha, method = "quasirandom", colour = "grey", stroke = 0) +
        ggbeeswarm::geom_quasirandom(
          data    = sig,
          mapping = aes(x     = !!sym(x_col),
                        y     = !!sym(group_by),
                        colour = .data$colour_var,
                        size  = !!sym(size_col)),
          alpha = alpha, method = "quasirandom", stroke = 0)
    } else {
      p <- p +
        ggbeeswarm::geom_quasirandom(
          data    = nonsig,
          mapping = aes(x = !!sym(x_col),
                        y = !!sym(group_by)),
          alpha = alpha, method = "quasirandom", size = fixed_size, colour = "grey", stroke = 0) +
        ggbeeswarm::geom_quasirandom(
          data    = sig,
          mapping = aes(x     = !!sym(x_col),
                        y     = !!sym(group_by),
                        colour = .data$colour_var),
          alpha = alpha, method = "quasirandom", size = fixed_size, stroke = 0)
    }

    p <- p + scale_colour_gradient2(
      low      = "blue", mid = "white", high = "red",
      midpoint = 0, name = "LFC", na.value = "grey"
    )

  } else {
    if (!is.null(size_col)) {
      p <- p + ggbeeswarm::geom_quasirandom(
        data    = plot_df,
        mapping = aes(x     = !!sym(x_col),
                      y     = !!sym(group_by),
                      colour = !!sym(colour_col),
                      size  = !!sym(size_col)),
        alpha = alpha, method = "quasirandom", stroke = 0)
    } else {
      p <- p + ggbeeswarm::geom_quasirandom(
        data    = plot_df,
        mapping = aes(x     = !!sym(x_col),
                      y     = !!sym(group_by),
                      colour = !!sym(colour_col)),
        alpha = alpha, method = "quasirandom", size = fixed_size, stroke = 0)
    }

    legend_name <- gsub("_", " ", tools::toTitleCase(colour_col))
    if (is_discrete_colour) {
      p <- p + scale_colour_discrete(
        name = if (colour_col == group_by) "Cell Type\nCluster" else legend_name
      )
    } else {
      p <- p + scale_colour_viridis_c(name = legend_name, option = "viridis")
    }
  }

  # Add summary line and label on top of all point layers
  p <- p + summary_layers

  # ── 11. Faceting ──────────────────────────────────────────────────────────────
  if (multi_col_mode) {
    if (!is.null(facet_by)) {
      # Two-way grid: metric rows x facet_by columns
      p <- p + facet_grid(
        rows   = vars(!!sym(".metric")),
        cols   = vars(!!sym(facet_by)),
        scales = facet_scales
      )
    } else if (length(plot_cols) > 1) {
      # Metric panels stacked as rows
      p <- p + facet_wrap(
        vars(!!sym(".metric")),
        ncol   = if (is.null(facet_ncol)) 1L else facet_ncol,
        scales = facet_scales
      )
    }
  } else if (!is.null(facet_by)) {
    p <- p + facet_wrap(
      vars(!!sym(facet_by)),
      ncol   = facet_ncol,
      scales = facet_scales
    )
  }

  # ── 12. Axis limits ───────────────────────────────────────────────────────────
  x_limit_val <- if (!is.null(x_limits)) {
    x_limits
  } else if (facet_scales == "fixed") {
    NULL   # let ggplot2 auto-range; pass NULL explicitly
  } else {
    NULL
  }
  p <- p + scale_x_continuous(limits = x_limit_val)

  # ── 13. Axis labels ───────────────────────────────────────────────────────────
  resolved_x_lab <- if (!is.null(x_lab)) {
    x_lab
  } else if (!multi_col_mode) {
    gsub("_", " ", tools::toTitleCase(plot_cols))
  } else {
    "Value"
  }

  resolved_y_lab <- if (!is.null(y_lab)) {
    y_lab
  } else {
    gsub("_", " ", tools::toTitleCase(group_by))
  }

  # ── 14. Labels and theme ──────────────────────────────────────────────────────
  p <- p +
    labs(
      x        = resolved_x_lab,
      y        = resolved_y_lab,
      title    = title,
      subtitle = subtitle
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title       = element_text(hjust = 0.5),
      plot.subtitle    = element_text(hjust = 0.5, size = 10),
      legend.position  = "right",
      axis.text.y      = element_text(size = 10)
    )

  if (!is.null(size_col))
    p <- p + scale_size_continuous(name  = gsub("_", " ", tools::toTitleCase(size_col)),
                                   range = size_range)

  return(p)
}
