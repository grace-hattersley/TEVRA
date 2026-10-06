#' Plot VDJ contig capture efficiency across neighbourhoods
#'
#' Generates two diagnostic plots summarising how efficiently VDJ contigs
#' were captured across transcriptional neighbourhoods: a jittered
#' boxplot of per-neighbourhood contig capture rate, and a scatterplot of
#' contig count against neighbourhood size (with a linear fit) to check
#' whether capture scales as expected with neighbourhood size.
#'
#' @param object The neighbourhood-level metadata data frame generated as
#'   the output of nhood_diversity, containing at minimum the columns
#'   \code{contig_capture_rate}, \code{n_contigs} and
#'   \code{nhood_size}.
#' @param group_by Character scalar giving the name of a column in
#'   \code{object} to colour/facet points by (e.g. cell type, timepoint, or
#'   tissue). If \code{NULL} (the default) or not found in
#'   \code{colnames(object)}, all neighbourhoods are treated as a single
#'   "All" group.
#'
#' @return A named list of two \code{ggplot} objects:
#'   \describe{
#'     \item{capture_rate}{Jittered boxplot of \code{contig_capture_rate} by
#'       \code{group_by}, with point size scaled to \code{nhood_size}.}
#'     \item{count_vs_size}{Scatterplot of \code{n_contigs} against
#'       \code{nhood_size}, coloured by \code{group_by}, with a linear
#'       trend line.}
#'   }
#'
#' @details
#' If \code{group_by} is supplied but is not already a factor, it is
#' coerced to one internally; the input \code{object} itself is not
#' modified. This function is purely a visualisation wrapper and does not
#' perform any capture-rate calculation itself.
#'
#' @seealso
#' \code{\link{nhood_diversity}}
#'
#' @examples
#' \dontrun{
#' plots <- plot_capture_efficiency(nhood_meta, group_by = "cell_type")
#' plots$capture_rate
#' plots$count_vs_size
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_jitter geom_boxplot geom_point
#'                     geom_smooth labs theme_classic theme element_text
#' @importFrom rlang sym .data
#'
#' @export
#'
plot_capture_efficiency <- function(object, group_by = NULL) {

  plot_data <- object

  if (!is.null(group_by) && group_by %in% colnames(plot_data)) {
    if (!is.factor(plot_data[[group_by]])) {
      plot_data[[group_by]] <- factor(plot_data[[group_by]])
    }
    group_var <- group_by
  } else {
    plot_data$group_by <- "All"
    group_var <- "group_by"
  }

  p1 <- ggplot(plot_data, aes(x = !!sym(group_var), y = .data$contig_capture_rate)) +
    geom_jitter(aes(size = .data$nhood_size), alpha = 0.3, width = 0.2) +
    geom_boxplot(outlier.alpha = 0, alpha = 0, color = "red") +
    labs(y = "Contig Capture Rate", title = "Contig Capture Efficiency by Neighbourhood") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p2 <- ggplot(plot_data, aes(x = .data$nhood_size, y = .data$n_contigs)) +
    geom_point(aes(color = !!sym(group_var)), alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
    labs(x = "Neighbourhood Size", y = "Number of Contigs",
         title = "Contig Count vs Neighbourhood Size") +
    theme_classic()

  return(list(capture_rate = p1, count_vs_size = p2))
}

