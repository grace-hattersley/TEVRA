#' Plot neighbourhood UMAP with optional cell-level and edge overlays
#'
#' Supports three composited layers:
#'   1. Cell-level UMAP scatter (optional)
#'   2. Neighbourhood-to-neighbourhood edges (optional)
#'   3. Neighbourhood centroids (always shown)
#'
#' Layer order (bottom to top): cells → edges → neighbourhood centroids.
#'
#' @param nhood_umap          A dataframe of neighbourhood centroids. Must
#'                            contain neighbourhood IDs, UMAP coordinates, and
#'                            any column used for fill_by or size_by.
#'                            Typically SCE@metadata$coords merged with
#'                            SCE@metadata$merged_nhood_analysis.
#' @param nhood_col           Name of the neighbourhood ID column in nhood_umap
#'                            (default "community_id").
#' @param x_col               Name of the UMAP x-coordinate column in
#'                            nhood_umap (default "mean_UMAP_1").
#' @param y_col               Name of the UMAP y-coordinate column in
#'                            nhood_umap (default "mean_UMAP_2").
#' @param fill_by             Name of a column in nhood_umap to colour
#'                            neighbourhood centroids by (e.g. "CD8" for
#'                            assigned cell type). If NULL, uses colour as a
#'                            single fixed colour.
#' @param colour              Controls neighbourhood centroid colour and
#'                            depends on \code{fill_by}:
#'                            \itemize{
#'                              \item When \code{fill_by} is \code{NULL} — a
#'                                single hex colour string (default
#'                                \code{"#3f4788"} when \code{colour} is left
#'                                \code{NULL}).
#'                              \item When \code{fill_by} is a column name —
#'                                a named character vector mapping values in
#'                                fill_by to hex colour strings. Names must
#'                                cover all levels present in \code{fill_by}.
#'                                If left \code{NULL}, a default categorical
#'                                palette is generated automatically.
#'                            }
#' @param size_by             Name of a numeric column in nhood_umap to map
#'                            centroid point size to, or NULL (default) for a
#'                            fixed point size.
#' @param size                Controls centroid point size and depends on
#'                            \code{size_by}:
#'                            \itemize{
#'                              \item When \code{size_by} is \code{NULL} — a
#'                                single numeric giving the fixed point size
#'                                (default \code{2} when \code{size} is left
#'                                \code{NULL}).
#'                              \item When \code{size_by} is a column name —
#'                                a numeric vector of length 2 giving the
#'                                minimum and maximum point size of the size
#'                                scale (default \code{c(1, 5)} when
#'                                \code{size} is left \code{NULL}).
#'                            }
#' @param show_edges          Logical. Whether to draw edges between connected
#'                            neighbourhoods (default FALSE).
#' @param nhood_connections   A dataframe with columns listing connected
#'                            neighbourhood pairs.Self-loops are silently
#'                            removed. Required when show_edges = TRUE.
#' @param edge_col_1          Name of the first neighbourhood ID column in
#'                            nhood_connections (default "mc_1").
#' @param edge_col_2          Name of the second neighbourhood ID column in
#'                            nhood_connections (default "mc_2").
#' @param edge_colour         Colour for edges (default "grey70").
#' @param edge_linewidth      Linewidth for edges (default 0.3).
#' @param edge_alpha          Alpha for edges (default 0.6).
#' @param show_cells          Logical. Whether to draw cell-level scatter as
#'                            the bottom layer (default FALSE).
#' @param cell_umap_coords    A dataframe with one row per cell containing UMAP
#'                            coordinates and any cell-level metadata columns.
#'                            Required when show_cells = TRUE.
#' @param cell_x_col          Name of the UMAP x-coordinate column in
#'                            cell_umap_coords (default "UMAP_1").
#' @param cell_y_col          Name of the UMAP y-coordinate column in
#'                            cell_umap_coords (default "UMAP_2").
#' @param cell_colour_by      Either the name of a column in cell_umap_coords
#'                            to colour cells by, or the string "nhood" to
#'                            inherit the neighbourhood colour via colour
#'                            and cell_nhood_col. Set NULL (default) for a
#'                            fixed grey cell colour.
#' @param cell_nhood_col      Name of the column in cell_umap_coords that
#'                            holds each cell's neighbourhood assignment.
#'                            Required when cell_colour_by = "nhood"
#'                            (default "nhood_id").
#' @param cell_colour         Controls cell point colour and depends on
#'                            \code{cell_colour_by}:
#'                            \itemize{
#'                              \item When \code{cell_colour_by} is
#'                                \code{NULL} — a single hex colour string
#'                                (default \code{"grey80"} when
#'                                \code{cell_colour} is left \code{NULL}).
#'                              \item When \code{cell_colour_by} is
#'                                \code{"nhood"} — ignored; colour is
#'                                inherited from \code{colour} via each
#'                                cell's neighbourhood assignment.
#'                              \item When \code{cell_colour_by} is any other
#'                                column name — a named character vector
#'                                mapping values of cell_colour_by to hex
#'                                colour strings. If left \code{NULL}, a
#'                                default categorical palette is generated
#'                                automatically.
#'                            }
#' @param cell_size           Fixed point size for cells (default 0.3).
#' @param cell_alpha          Alpha for cell points (default 0.5).
#'
#' @return A \code{ggplot} object. Layers are composited in the order cells
#'   (optional) -> edges (optional) -> neighbourhood centroids, so the result
#'   can be extended with further ggplot2 layers, scales or themes.
#'
#' @seealso
#' \code{\link{create_communities}} and \code{\link{test_resolution}} for
#' generating the neighbourhood assignments plotted here;
#' \code{\link{nhood_umap}} and \code{\link{test_resolution}} for
#' generating the neighbourhood umap coordinates;
#' \code{\link{nhood_annotation}} for annotating neighbourhoods
#'
#' @examples
#' \dontrun{
#' # --- Neighbourhood centroids coloured by assigned cell type ---------------
#' plot_nhood_umap(
#'   nhood_umap = nhood_df,
#'   fill_by    = "CD8",
#'   colour     = celltype_colours
#' )
#'
#' # --- Neighbourhood centroids with default palette -------------------------
#' plot_nhood_umap(
#'   nhood_umap = nhood_df,
#'   fill_by    = "CD8"
#' )
#'
#' # --- Neighbourhood centroids scaled by neighbourhood size -----------------
#' plot_nhood_umap(
#'   nhood_umap = nhood_df,
#'   size_by    = "n_cells",
#'   size       = c(1, 6)
#' )
#'
#' # --- Neighbourhood centroids, edges and cells -----------------------------
#' plot_nhood_umap(
#'   nhood_umap        = nhood_df,
#'   fill_by           = "CD8",
#'   colour            = celltype_colours,
#'   show_edges        = TRUE,
#'   nhood_connections = connections_df,
#'   show_cells        = TRUE,
#'   cell_umap_coords  = cell_df,
#'   cell_colour_by    = "nhood",
#'   cell_nhood_col    = "nhood_id"
#' )
#'
#' # --- Neighbourhood centroids, edges and cells coloured independently ------
#' plot_nhood_umap(
#'   nhood_umap       = nhood_df,
#'   fill_by          = "CD8",
#'   colour           = celltype_colours,
#'   show_cells       = TRUE,
#'   cell_umap_coords = cell_df,
#'   cell_colour_by   = "timepoint",
#'   cell_colour      = c(Pre = "#8da0cb", On = "#fc8d62")
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_segment
#'                     scale_colour_manual scale_fill_manual scale_size_continuous
#'                     theme_classic theme element_blank element_text unit labs
#' @importFrom rlang .data sym
#' @importFrom scales hue_pal
#' @importFrom stats setNames
#'
#' @export
#'
plot_nhood_umap <- function(
    # ── Neighbourhood centroid arguments ─────────────────────────────────────
    nhood_umap,
    nhood_col           = "community_id",
    x_col               = "mean_UMAP_1",
    y_col               = "mean_UMAP_2",
    fill_by             = NULL,
    colour              = NULL,
    size_by             = NULL,
    size                = NULL,
    # ── Edge arguments ───────────────────────────────────────────────────────
    show_edges          = FALSE,
    nhood_connections   = NULL,
    edge_col_1          = "mc_1",
    edge_col_2          = "mc_2",
    edge_colour         = "grey70",
    edge_linewidth      = 0.3,
    edge_alpha          = 0.6,
    # ── Cell-level arguments ─────────────────────────────────────────────────
    show_cells          = FALSE,
    cell_umap_coords    = NULL,
    cell_x_col          = "UMAP_1",
    cell_y_col          = "UMAP_2",
    cell_colour_by      = NULL,
    cell_nhood_col      = "nhood_id",
    cell_colour         = NULL,
    cell_size           = 0.3,
    cell_alpha          = 0.5
) {

  # ── 1. Validate neighbourhood dataframe ────────────────────────────────────
  if (!is.data.frame(nhood_umap)) {
    stop(
      "'nhood_umap' must be a dataframe. ",
      "Pass SCE@metadata$coords merged with annotation columns."
    )
  }

  required_nhood_cols <- c(nhood_col, x_col, y_col)
  missing_nhood_cols  <- setdiff(required_nhood_cols, colnames(nhood_umap))
  if (length(missing_nhood_cols) > 0) {
    stop(sprintf(
      "nhood_umap is missing required column(s): %s",
      paste(missing_nhood_cols, collapse = ", ")
    ))
  }

  # Resolve size_by/size into the two internal flags used throughout:
  #   size_col   : column name to map to size aes (NULL when fixed)
  #   fixed_size : numeric size passed to geom directly (NULL when col mode)
  #   size_range : length-2 numeric for scale_size_continuous (NULL when fixed)
  size_info  <- .resolve_nhood_size(nhood_umap, size_by, size,
                                    default_size_range = c(1, 5))
  size_col   <- size_info$size_col
  fixed_size <- size_info$fixed_size
  size_range <- size_info$size_range

  if (!is.null(fill_by)) {
    if (!fill_by %in% colnames(nhood_umap)) {
      stop(sprintf(
        "'fill_by' (\"%s\") not found in nhood_umap. Available columns: %s",
        fill_by, paste(colnames(nhood_umap), collapse = ", ")
      ))
    }

    fill_levels <- unique(nhood_umap[[fill_by]])

    if (is.null(colour)) {
      # No mapping supplied — fall back to ggplot2's default discrete palette
      fill_levels_nona <- fill_levels[!is.na(fill_levels)]
      colour <- stats::setNames(
        scales::hue_pal()(length(fill_levels_nona)),
        as.character(fill_levels_nona)
      )
    } else {
      if (is.null(names(colour))) {
        stop(
          "'colour' must be a named character vector mapping each level ",
          "of 'fill_by' to a colour when 'fill_by' is not NULL (a single ",
          "unnamed colour was supplied instead)."
        )
      }
      missing_levels <- setdiff(fill_levels, names(colour))
      if (length(missing_levels) > 0) {
        stop(sprintf(
          "'colour' is missing entries for: %s",
          paste(missing_levels, collapse = ", ")
        ))
      }
    }
  } else {
    if (is.null(colour)) {
      colour <- "#3f4788"
    } else if (!is.null(names(colour))) {
      stop(
        "'colour' must be a single (unnamed) colour string when 'fill_by' ",
        "is NULL (a named vector was supplied instead, implying a mapping)."
      )
    }
  }

  # ── 2. Validate and prepare edges ──────────────────────────────────────────
  if (show_edges) {
    if (is.null(nhood_connections)) {
      stop("'nhood_connections' must be provided when show_edges = TRUE.")
    }
    if (!is.data.frame(nhood_connections)) {
      stop("'nhood_connections' must be a dataframe.")
    }
    missing_edge_cols <- setdiff(c(edge_col_1, edge_col_2), colnames(nhood_connections))
    if (length(missing_edge_cols) > 0) {
      stop(sprintf(
        "nhood_connections is missing column(s): %s",
        paste(missing_edge_cols, collapse = ", ")
      ))
    }

    # Remove self-loops and coerce IDs to character to match nhood_umap
    edges <- nhood_connections[
      as.character(nhood_connections[[edge_col_1]]) !=
        as.character(nhood_connections[[edge_col_2]]), ,
      drop = FALSE
    ]
    edges[[edge_col_1]] <- as.character(edges[[edge_col_1]])
    edges[[edge_col_2]] <- as.character(edges[[edge_col_2]])

    # Ensure nhood IDs are character for joining
    nhood_coords_for_edges <- nhood_umap
    nhood_coords_for_edges[[nhood_col]] <- as.character(
      nhood_coords_for_edges[[nhood_col]]
    )

    # Join coordinates for both endpoints
    edges <- merge(
      edges,
      nhood_coords_for_edges[, c(nhood_col, x_col, y_col)],
      by.x = edge_col_1,
      by.y = nhood_col,
      all.x = TRUE
    )
    colnames(edges)[colnames(edges) == x_col] <- "x_start"
    colnames(edges)[colnames(edges) == y_col] <- "y_start"

    edges <- merge(
      edges,
      nhood_coords_for_edges[, c(nhood_col, x_col, y_col)],
      by.x = edge_col_2,
      by.y = nhood_col,
      all.x = TRUE
    )
    colnames(edges)[colnames(edges) == x_col] <- "x_end"
    colnames(edges)[colnames(edges) == y_col] <- "y_end"

    # Drop edges where either endpoint has no coordinates (unmatched IDs)
    n_before <- nrow(edges)
    edges <- edges[
      is.finite(edges$x_start) & is.finite(edges$y_start) &
        is.finite(edges$x_end)   & is.finite(edges$y_end), ,
      drop = FALSE
    ]
    n_dropped <- n_before - nrow(edges)
    if (n_dropped > 0) {
      warning(sprintf(
        "%d edge(s) dropped: one or both endpoint IDs not found in nhood_umap.",
        n_dropped
      ))
    }
  }

  # ── 3. Validate and prepare cell data ──────────────────────────────────────
  if (show_cells) {
    if (is.null(cell_umap_coords)) {
      stop("'cell_umap_coords' must be provided when show_cells = TRUE.")
    }
    if (!is.data.frame(cell_umap_coords)) {
      stop("'cell_umap_coords' must be a dataframe.")
    }
    missing_cell_cols <- setdiff(c(cell_x_col, cell_y_col), colnames(cell_umap_coords))
    if (length(missing_cell_cols) > 0) {
      stop(sprintf(
        "cell_umap_coords is missing column(s): %s",
        paste(missing_cell_cols, collapse = ", ")
      ))
    }

    if (!is.null(cell_colour_by)) {
      if (cell_colour_by == "nhood") {
        # Cells inherit colour from their neighbourhood's fill_by value
        if (is.null(fill_by)) {
          stop(
            "cell_colour_by = \"nhood\" requires fill_by to be set, ",
            "so that neighbourhood cell types can be looked up."
          )
        }
        if (!cell_nhood_col %in% colnames(cell_umap_coords)) {
          stop(sprintf(
            "'cell_nhood_col' (\"%s\") not found in cell_umap_coords.",
            cell_nhood_col
          ))
        }

        # Map each cell's neighbourhood ID → fill_by value → colour
        nhood_lookup <- nhood_umap[, c(nhood_col, fill_by)]
        nhood_lookup[[nhood_col]] <- as.character(nhood_lookup[[nhood_col]])
        cell_umap_coords[[cell_nhood_col]] <- as.character(
          cell_umap_coords[[cell_nhood_col]]
        )

        cell_umap_coords <- merge(
          cell_umap_coords,
          nhood_lookup,
          by.x  = cell_nhood_col,
          by.y  = nhood_col,
          all.x = TRUE
        )
        # Cells with no neighbourhood match will have NA in fill_by column.
        # Colour itself is inherited from the already-resolved 'colour'
        # (validated/defaulted above against fill_by) — cell_colour is not
        # used in this mode.
        cell_colour_by_col <- fill_by

      } else {
        # cell_colour_by is a column name in cell_umap_coords
        if (!cell_colour_by %in% colnames(cell_umap_coords)) {
          stop(sprintf(
            "'cell_colour_by' (\"%s\") not found in cell_umap_coords. Available columns: %s",
            cell_colour_by,
            paste(colnames(cell_umap_coords), collapse = ", ")
          ))
        }

        cell_levels <- unique(cell_umap_coords[[cell_colour_by]])

        if (is.null(cell_colour)) {
          # No mapping supplied — fall back to ggplot2's default discrete palette
          cell_levels_nona <- cell_levels[!is.na(cell_levels)]
          cell_colour <- stats::setNames(
            scales::hue_pal()(length(cell_levels_nona)),
            as.character(cell_levels_nona)
          )
        } else {
          if (is.null(names(cell_colour))) {
            stop(
              "'cell_colour' must be a named character vector mapping each ",
              "level of 'cell_colour_by' to a colour when cell_colour_by is ",
              "a column name (a single unnamed colour was supplied instead)."
            )
          }
          missing_cell_levels <- setdiff(cell_levels, names(cell_colour))
          if (length(missing_cell_levels) > 0) {
            stop(sprintf(
              "'cell_colour' is missing entries for: %s",
              paste(missing_cell_levels, collapse = ", ")
            ))
          }
        }
        cell_colour_by_col <- cell_colour_by
      }
    } else {
      # cell_colour_by is NULL — fixed single colour for all cells
      if (is.null(cell_colour)) {
        cell_colour <- "grey80"
      } else if (!is.null(names(cell_colour))) {
        stop(
          "'cell_colour' must be a single (unnamed) colour string when ",
          "'cell_colour_by' is NULL (a named vector was supplied instead, ",
          "implying a mapping)."
        )
      }
    }
  }

  # ── 4. Build plot ──────────────────────────────────────────────────────────

  # Initialise on the neighbourhood dataframe
  p <- ggplot2::ggplot() +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      axis.line       = ggplot2::element_blank(),
      axis.ticks      = ggplot2::element_blank(),
      axis.text       = ggplot2::element_blank(),
      legend.position = "right",
      legend.key.size = ggplot2::unit(0.4, "cm"),
      legend.text     = ggplot2::element_text(size = 8),
      legend.title    = ggplot2::element_text(size = 9)
    ) +
    ggplot2::labs(x = "UMAP1", y = "UMAP2")

  # ── 4.1 Cell layer ─────────────────────────────────────────────────────────
  if (show_cells) {
    if (!is.null(cell_colour_by)) {
      # Resolve colour map: "nhood" mode inherits from colour; column mode uses cell_colour
      active_cell_colour_map <- if (cell_colour_by == "nhood") colour else cell_colour

      p <- p +
        ggplot2::geom_point(
          data    = cell_umap_coords,
          mapping = ggplot2::aes(
            x     = .data[[cell_x_col]],
            y     = .data[[cell_y_col]],
            colour = .data[[cell_colour_by_col]]
          ),
          size  = cell_size,
          alpha = cell_alpha,
          shape = 16
        ) +
        ggplot2::scale_colour_manual(
          values   = active_cell_colour_map,
          na.value = "grey85",
          guide    = "none"
        )
    } else {
      p <- p +
        ggplot2::geom_point(
          data    = cell_umap_coords,
          mapping = ggplot2::aes(
            x = .data[[cell_x_col]],
            y = .data[[cell_y_col]]
          ),
          colour = cell_colour,
          # cell_colour is already resolved to a single hex string above
          size  = cell_size,
          alpha = cell_alpha,
          shape = 16
        )
    }
  }

  # ── 4.2 Edges ──────────────────────────────────────────────────────────────
  if (show_edges && nrow(edges) > 0) {
    p <- p +
      ggplot2::geom_segment(
        data    = edges,
        mapping = ggplot2::aes(
          x    = .data$x_start,
          y    = .data$y_start,
          xend = .data$x_end,
          yend = .data$y_end
        ),
        colour    = edge_colour,
        linewidth = edge_linewidth,
        alpha     = edge_alpha
      )
  }

  # ── 4.3 Neighbourhood layer ────────────────────────────────────────────────
  nhood_aes_args <- list(
    x = rlang::sym(x_col),
    y = rlang::sym(y_col)
  )
  if (!is.null(size_col)) nhood_aes_args$size <- rlang::sym(size_col)
  if (!is.null(fill_by))  nhood_aes_args$fill <- rlang::sym(fill_by)
  nhood_aes <- do.call(ggplot2::aes, nhood_aes_args)

  # Fixed aesthetics — size omitted when mapped via aes()
  nhood_fixed <- list(shape = 21, colour = "black", stroke = 0.4)
  if (is.null(size_col)) nhood_fixed$size <- fixed_size

  if (!is.null(fill_by)) {
    p <- p +
      do.call(
        ggplot2::geom_point,
        c(list(data = nhood_umap, mapping = nhood_aes), nhood_fixed)
      ) +
      ggplot2::scale_fill_manual(
        values   = colour,
        na.value = "grey85",
        name     = fill_by
      )
  } else {
    p <- p +
      do.call(
        ggplot2::geom_point,
        c(list(data = nhood_umap, mapping = nhood_aes, fill = colour), nhood_fixed)
      )
  }

  if (!is.null(size_col)) {
    p <- p +
      ggplot2::scale_size_continuous(
        range = size_range,
        name  = gsub("_", " ", tools::toTitleCase(size_col))
      )
  }

  p
}
