#' Merge Neighbourhood Analysis Results
#'
#' Merges any number of neighbourhood analysis dataframes stored in the
#' metadata slot of a neighbourhood object. Automatically detects which
#' dataframes exist, validates neighbourhood consistency across all sources,
#' and handles redundant columns.
#'
#' @param object A \code{SingleCellExperiment} or \code{Milo} object containing
#'   neighbourhood analysis results in its metadata slot.
#' @param dataframe_names Character vector of metadata slot names to merge.
#'   Defaults to \code{c("nhood_diversity", "nhood_shm")}. Must contain at
#'   least two names.
#' @param merge_key Column to merge on (default: \code{"nhood"}).
#' @param size_col Column name containing neighbourhood size information, used
#'   to validate that all dataframes describe the same neighbourhoods. Can be a
#'   single string (same column name expected in all dataframes) or a named
#'   character vector mapping each dataframe name to its size column, e.g.
#'   \code{c(nhood_data = "nhood_size", nhood_shm = "n_cells")}. Dataframes
#'   not present in a named vector fall back to \code{"nhood_size"}.
#'   Set to \code{NULL} to skip the size consistency check entirely.
#' @param verbose Print merge details (default: \code{TRUE}).
#' @param result_name Name for the merged output in metadata
#'   (default: \code{"merged_nhood_analysis"}).
#'
#' @return Updated \code{object} with merged dataframe stored in
#'   \code{metadata(object)[[result_name]]}.
#'
#' @details
#' ## Dataframe availability:
#' \enumerate{
#'   \item \strong{Multiple dataframes found} — validates shared neighbourhoods
#'     then merges all sequentially on \code{merge_key} using inner joins,
#'     removing redundant columns at each step.
#'   \item \strong{Only one dataframe found} — halts with an error (a merge
#'   requires at least two inputs)
#'   \item \strong{No dataframes found} — halts with an error
#'   \item \strong{Result already exists} — overwrites with a warning.
#' }
#'
#' ## Neighbourhood validation (two-stage):
#' \enumerate{
#'   \item \emph{merge_key overlap} — checks that the \code{merge_key} entries
#'     in every dataframe overlap with those in the reference (first) dataframe.
#'     If no overlap exists for any pair the function halts. If only partial
#'     overlap exists a warning is issued and only shared neighbourhoods
#'     proceed.
#'   \item \emph{size consistency} — for neighbourhoods shared across
#'     dataframes, confirms that the size column contains identical values. Any
#'     mismatch halts the function as it indicates analyses were run on
#'     different data.
#' }
#'
#' ## Redundant column handling:
#' At each pairwise merge step, columns
#' shared between the two dataframes are compared value-by-value. Columns with
#' identical values are dropped from the right-hand dataframe. Columns with
#' different values are retained with \code{_1} / \code{_2} suffixes.
#'
#' @seealso
#' \code{\link{nhood_diversity}}, \code{\link{nhood_cdr3_convergence}},
#' \code{\link{nhood_shm}}, \code{\link{nhood_annotation}}
#'
#' @examples
#' \dontrun{
#' # Default: merge "nhood_diversity" and "nhood_shm"
#' object <- merge_nhood_results(object)
#'
#' # Merge three custom dataframes
#' object <- merge_nhood_results(
#'   object,
#'   dataframe_names = c("diversity_results", "nhood_shm", "spatial_features")
#' )
#'
#' # Dataframes use different column names for neighbourhood size
#' object <- merge_nhood_results(
#'   object,
#'   dataframe_names = c("diversity_results", "nhood_shm", "spatial_features"),
#'   size_col = c(diversity_results = "nhood_size",
#'                nhood_shm         = "n_cells",
#'                spatial_features  = "size")
#' )
#'
#' # Custom merge key and output name
#' object <- merge_nhood_results(
#'   object,
#'   dataframe_names = c("results_a", "results_b"),
#'   result_name     = "combined_results",
#'   merge_key       = "NhoodID")
#' }
#'
#' @export
#'
merge_nhood_results <- function(object,
                                dataframe_names = c("nhood_diversity", "nhood_shm"),
                                merge_key       = "nhood",
                                size_col        = "nhood_size",
                                verbose         = TRUE,
                                result_name     = "merged_nhood_analysis") {

  # ── 1. Object type ────────────────────────────────────────────────────────────
  is_milo     <- inherits(object, "Milo") && inherits(object, "SingleCellExperiment")
  is_sce_only <- !is_milo && inherits(object, "SingleCellExperiment")

  if (!is_milo && !is_sce_only) {
    stop("'object' must be a SingleCellExperiment (or Milo) object.")
  }

  # ── 2. Identify which requested dataframes are actually present ───────────────
  if (verbose) {
    message(strrep("=", 70))
    message("Merging Neighbourhood Analysis Results")
    message(strrep("=", 70))
    message("Checking available data (", length(dataframe_names), " requested):")
  }

  found <- vapply(dataframe_names, function(nm) {
    nm %in% names(S4Vectors::metadata(object))
  }, logical(1))

  if (verbose) {
    for (nm in dataframe_names) {
      status <- if (found[nm]) { "\u2713 Found" } else { "\u2717 Not found" }
      message("  ", nm, ": ", status)
    }
  }

  found_names   <- dataframe_names[found]
  missing_names <- dataframe_names[!found]

  if (length(found_names) == 0) {
    stop(
      "None of the requested dataframes were found in metadata(object).\n",
      "  Requested : ", paste(dataframe_names, collapse = ", "), "\n",
      "  Available : ", paste(names(S4Vectors::metadata(object)), collapse = ", ")
    )
  }

  if (length(missing_names) > 0 && verbose) {
    message("\u26a0 ", length(missing_names),
            " requested dataframe(s) not found and will be skipped: ",
            paste(missing_names, collapse = ", "))
  }

  # ── 3. Require at least two dataframes ───────────────────────────────────────
  if (length(found_names) == 1) {
    stop(
      "Only one dataframe was found ('", found_names, "').\n",
      "  A merge requires at least two dataframes.\n",
      "  Please check 'dataframe_names' or add further analyses to object."
    )
  }

  # ── 4. Retrieve all dataframes ────────────────────────────────────────────────
  df_list <- lapply(found_names, function(nm) S4Vectors::metadata(object)[[nm]])
  names(df_list) <- found_names

  if (verbose) {
    message("Input dimensions:")
    for (nm in found_names) {
      message(sprintf("  %-35s %d rows \u00d7 %d columns",
                      nm, nrow(df_list[[nm]]), ncol(df_list[[nm]])))
    }
  }

  # ── 5. Validate neighbourhood consistency across ALL dataframes ───────────────
  #    Stage 1 — merge_key overlap
  #    Stage 2 — size column consistency (on shared keys only)
  if (verbose) {
    message("Validating neighbourhood consistency...")
    message(strrep("-", 70))
  }

  ref_keys <- df_list[[1]][[merge_key]]
  ref_name <- found_names[1]

  for (nm in found_names[-1]) {
    cmp_keys <- df_list[[nm]][[merge_key]]

    shared   <- intersect(ref_keys, cmp_keys)
    n_shared <- length(shared)
    n_ref    <- length(ref_keys)
    n_cmp    <- length(cmp_keys)

    # ── Stage 1: merge_key overlap ─────────────────────────────────────────────
    if (n_shared == 0) {
      stop(
        "No shared '", merge_key, "' entries between '",
        ref_name, "' (", n_ref, " neighbourhoods) and '",
        nm, "' (", n_cmp, " neighbourhoods).\n",
        "  These dataframes appear to describe entirely different neighbourhoods.\n",
        "  Cannot merge \u2014 please check that analyses were run on the same data."
      )
    }

    if (n_shared < n_ref || n_shared < n_cmp) {
      if (verbose) {
        message(sprintf(
          "  \u26a0 '%s' vs '%s': partial overlap \u2014 %d/%d '%s' and %d/%d '%s' neighbourhoods shared",
          ref_name, nm,
          n_shared, n_ref, ref_name,
          n_shared, n_cmp, nm
        ))
      }
    } else {
      if (verbose) {
        message(sprintf("  \u2713 '%s' vs '%s': all %d neighbourhoods shared",
                        ref_name, nm, n_shared))
      }
    }

    # ── Stage 2: size consistency (shared keys only) ───────────────────────────
    ref_sc <- .get_size_col(ref_name, size_col)
    cmp_sc <- .get_size_col(nm, size_col)

    if (!is.null(ref_sc) && !is.null(cmp_sc)) {
      ref_has_sc <- ref_sc %in% colnames(df_list[[ref_name]])
      cmp_has_sc <- cmp_sc %in% colnames(df_list[[nm]])

      if (!ref_has_sc || !cmp_has_sc) {
        missing_info <- c(
          if (!ref_has_sc) { paste0("'", ref_name, "' (looking for '", ref_sc, "')") },
          if (!cmp_has_sc) { paste0("'", nm,       "' (looking for '", cmp_sc, "')") }
        )
        if (verbose) {
          message("  \u26a0 Size column not found in: ",
                  paste(missing_info, collapse = "; "),
                  " \u2014 skipping size check for this pair")
        }
      } else {
        ref_sub <- df_list[[ref_name]][
          df_list[[ref_name]][[merge_key]] %in% shared,
          c(merge_key, ref_sc), drop = FALSE
        ]
        cmp_sub <- df_list[[nm]][
          df_list[[nm]][[merge_key]] %in% shared,
          c(merge_key, cmp_sc), drop = FALSE
        ]
        colnames(ref_sub)[2] <- "size_ref"
        colnames(cmp_sub)[2] <- "size_cmp"

        size_cmp_df <- merge(ref_sub, cmp_sub, by = merge_key)
        diffs <- size_cmp_df$size_ref != size_cmp_df$size_cmp
        diffs[is.na(size_cmp_df$size_ref) & is.na(size_cmp_df$size_cmp)] <- FALSE
        n_diff <- sum(diffs, na.rm = TRUE)

        if (n_diff > 0) {
          stop(
            "Neighbourhood size differs between '", ref_name,
            "' (column '", ref_sc, "') and '", nm,
            "' (column '", cmp_sc, "')!\n",
            "  ", n_diff, " of ", nrow(size_cmp_df),
            " shared neighbourhoods have different cell counts.\n",
            "  This suggests the analyses were run on different datasets.\n",
            "  Example mismatches (first 5):\n",
            paste(utils::capture.output(print(utils::head(size_cmp_df[diffs, ], 5))),
                  collapse = "\n  "),
            "\n\n  Cannot merge \u2014 please re-run analyses on the same data."
          )
        } else {
          if (verbose) {
            sc_info <- if (ref_sc == cmp_sc) {
              paste0("'", ref_sc, "'")
            } else {
              paste0("'", ref_sc, "' / '", cmp_sc, "'")
            }
            message(sprintf(
              "  \u2713 '%s' vs '%s': size column (%s) IDENTICAL across all %d shared neighbourhoods",
              ref_name, nm, sc_info, n_shared
            ))
          }
        }
      }
    } else if (verbose) {
      message("  \u26a0 size_col = NULL \u2014 skipping size consistency check")
    }
  }

  if (verbose) {
    message(strrep("-", 70))
  }

  # ── 6. Sequential pairwise merge with redundancy checks ──────────────────────
  if (verbose) {
    message("Performing sequential merge on '", merge_key, "' column...")
  }

  accumulated      <- df_list[[1]]
  accumulated_name <- found_names[1]

  for (i in seq(2, length(found_names))) {
    right_name <- found_names[i]
    right_df   <- df_list[[i]]

    if (verbose) {
      message(sprintf("  Step %d: '%s'  +  '%s'", i - 1, accumulated_name, right_name))
    }

    right_df <- .check_and_clean_redundant(
      left_df    = accumulated,
      right_df   = right_df,
      left_name  = accumulated_name,
      right_name = right_name,
      merge_key  = merge_key,
      iter       = i - 1,
      verbose    = verbose
    )

    accumulated      <- merge(accumulated, right_df, by = merge_key, all = FALSE)
    accumulated_name <- paste0(accumulated_name, " + ", right_name)

    if (verbose) {
      message(sprintf("    \u2192 %d rows \u00d7 %d columns",
                      nrow(accumulated), ncol(accumulated)))
    }
  }

  # ── 7. Report rows lost to inner joining ──────────────────────────────────────
  if (verbose) {
    max_input_rows <- max(vapply(df_list, nrow, integer(1)))
    n_lost <- max_input_rows - nrow(accumulated)
    if (n_lost > 0) {
      message("\u26a0 Note: inner join(s) excluded ", n_lost,
              " neighbourhood(s) not present in all dataframes.")
    }
  }

  # ── 8. Store results in metadata and return ───────────────────────────────────
  if (result_name %in% names(S4Vectors::metadata(object))) {
    warning(sprintf(
      "metadata(object)[['%s']] already exists and will be overwritten.",
      result_name
    ))
  }

  S4Vectors::metadata(object)[[result_name]] <- accumulated

  message(sprintf(
    "Stored merged results (%d rows \u00d7 %d columns) in metadata(object)[['%s']].",
    nrow(accumulated), ncol(accumulated), result_name
  ))

  return(object)
}


# =============================================================================
# Internal helpers
# =============================================================================

#' Resolve the size-column name for a given dataframe
#'
#' Given the \code{size_col} argument passed to \code{merge_nhood_results},
#' returns the column name that should be used to look up neighbourhood sizes
#' for a specific dataframe.  Handles the three supported forms of
#' \code{size_col}: \code{NULL} (skip the check), an unnamed scalar/vector
#' (same column name for every dataframe), or a named character vector
#' (per-dataframe mapping with \code{"nhood_size"} as the fallback).
#'
#' @param df_name A single character string giving the name of the dataframe
#'   whose size column is being resolved.
#' @param size_col The \code{size_col} argument as received by
#'   \code{merge_nhood_results}: either \code{NULL}, an unnamed character
#'   vector (the first element is used for all dataframes), or a named
#'   character vector mapping dataframe names to column names.
#'
#' @return A single character string containing the column name to use, or
#'   \code{NULL} if \code{size_col} is \code{NULL} (i.e. the size check
#'   should be skipped entirely).
#'
#' @noRd
.get_size_col <- function(df_name, size_col) {
  if (is.null(size_col)) {
    return(NULL)
  }
  if (!is.null(names(size_col))) {
    if (df_name %in% names(size_col)) {
      size_col[[df_name]]
    } else {
      "nhood_size"
    }
  } else {
    size_col[[1]]
  }
}


#' Detect and remove redundant columns from the right-hand merge dataframe
#'
#' Before two dataframes are merged on \code{merge_key}, this helper
#' identifies every column (other than \code{merge_key}) that appears in
#' \emph{both} dataframes.  For each shared column it performs a row-wise
#' value comparison (after aligning on \code{merge_key}) and classifies it as
#' either identical or different:
#'
#' \itemize{
#'   \item \strong{Identical} columns are silently dropped from
#'     \code{right_df} to avoid duplicated information in the merged output.
#'   \item \strong{Different} columns are left in \code{right_df} unchanged;
#'     \code{base::merge} will append \code{"_1"} / \code{"_2"} suffixes
#'     (where the numbers reflect the current iteration) to disambiguate them.
#' }
#'
#' If there are no columns in common (excluding the merge key) the function
#' returns \code{right_df} unmodified.
#'
#' @param left_df A data frame — the accumulated result of all previous merge
#'   steps.
#' @param right_df A data frame — the next dataframe to be merged in.
#' @param left_name Character string; display name of \code{left_df} used in
#'   progress messages.
#' @param right_name Character string; display name of \code{right_df} used in
#'   progress messages.
#' @param merge_key Character string; name of the column that the two
#'   dataframes will be joined on.  This column is always excluded from the
#'   redundancy check.
#' @param iter Integer; the current merge-step index (1-based).  Used to
#'   construct informative suffix labels in the verbose message for differing
#'   columns.
#' @param verbose Logical; if \code{TRUE} a per-column summary is printed to
#'   the console.
#'
#' @return \code{right_df} with all value-identical shared columns removed.
#'   The \code{merge_key} column and any columns unique to \code{right_df} are
#'   always retained.
#'
#' @noRd
.check_and_clean_redundant <- function(left_df, right_df,
                                       left_name, right_name,
                                       merge_key, iter, verbose) {
  common_cols <- setdiff(intersect(colnames(left_df), colnames(right_df)),
                         merge_key)

  if (length(common_cols) == 0) {
    return(right_df)
  }

  if (verbose) {
    message("  Checking ", length(common_cols), " common column(s) between '",
            left_name, "' and '", right_name, "':")
  }

  identical_cols <- character(0)
  different_cols <- character(0)

  for (col in common_cols) {
    tmp <- merge(
      left_df[,  c(merge_key, col), drop = FALSE],
      right_df[, c(merge_key, col), drop = FALSE],
      by       = merge_key,
      suffixes = c("_l", "_r")
    )
    diffs <- tmp[[paste0(col, "_l")]] != tmp[[paste0(col, "_r")]]
    diffs[is.na(tmp[[paste0(col, "_l")]]) &
            is.na(tmp[[paste0(col, "_r")]])] <- FALSE
    n_diff <- sum(diffs, na.rm = TRUE)

    if (n_diff == 0) {
      identical_cols <- c(identical_cols, col)
      if (verbose) { message(sprintf("    \u2713 %-30s IDENTICAL", col)) }
    } else {
      different_cols <- c(different_cols, col)
      if (verbose) {
        message(sprintf("    \u2717 %-30s DIFFERENT (%d/%d rows)",
                        col, n_diff, nrow(tmp)))
      }
    }
  }

  if (length(identical_cols) > 0) {
    if (verbose) {
      message("  Removing ", length(identical_cols),
              " redundant column(s) from '", right_name, "': ",
              paste(identical_cols, collapse = ", "))
    }
    right_df <- right_df[,
                         !(colnames(right_df) %in% identical_cols), drop = FALSE]
  }

  if (length(different_cols) > 0 && verbose) {
    message("  \u26a0 ", length(different_cols),
            " column(s) with differing values will be kept with suffixes '_",
            iter, "' and '_", iter + 1, "': ",
            paste(different_cols, collapse = ", "))
  }

  right_df
}
