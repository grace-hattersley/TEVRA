#' Plot the distribution of neighbourhood sizes
#'
#' Visualises the distribution of cell counts per neighbourhood as a
#' histogram, following neighbourhood assignment with \code{create_communities}.
#'
#' @param communities_summary A data frame summarising neighbourhoods, as
#'   returned by \code{create_communities}, containing at minimum a numeric
#'   column \code{n_cells} giving the number of cells assigned to each
#'   neighbourhood.
#'
#' @return A \code{ggplot} object showing a histogram of neighbourhood sizes
#'   (x-axis: neighbourhood size; y-axis: count).
#'
#' @seealso \code{\link{create_communities}} to create communities and obtain
#' communities summary
#'
#' @importFrom ggplot2 ggplot aes geom_histogram labs theme_classic
#' @importFrom rlang .data
#'
#' @export
#'
plot_nhood_sizes <- function(communities_summary){
  ggplot(communities_summary, aes(x = .data$n_cells)) +
    geom_histogram(bins = 50) +
    labs(x = 'Neighbourhood size', y = 'count') +
    theme_classic()}
