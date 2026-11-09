# =============================================================================
# Combine best occupancy probability plots
#
# This script:
#
#   1. Loads or extracts model outputs.
#   2. Selects the best models using WAIC:
#        - one best multiple-species model
#        - one best model per species for single-species models
#   3. Finds the corresponding occupancy probability PNG plots.
#   4. Verifies that each species has one plot from each model type.
#   5. Combines the two plots side-by-side for each species.
#   6. Combines all species plots into a single multi-page PDF.
#
# =============================================================================

# =============================================================================
# Function: combine_best_occupancy_plots
# =============================================================================

#' Select the best models by WAIC
#'
#' Lower WAIC indicates a better-fitting model. For the multiple-species
#' models the model with the lowest WAIC is selected; for the single-species
#' models the model with the lowest WAIC is selected separately for each
#' species. The selection is saved as `best_models.csv` when
#' `path_to_postProcessing` is given. This step does not need any predictions.
#'
#' @param outputs Model outputs with a `waic_vals` table.
#' @param path_to_postProcessing Folder for `best_models.csv` (or `NULL`).
#'
#' @return A data frame with the selected models.
#' @export
select_best_models <- function(outputs, path_to_postProcessing = NULL) {
  if (
    !is.data.frame(outputs$waic_vals) ||
      !all(c("WAIC", "single_species") %in% names(outputs$waic_vals))
  ) {
    # (also the case when no model files were found: the table is empty)
    stop(
      "`outputs$waic_vals` has no WAIC values. Extract the model outputs with ",
      "`modelChecks = TRUE` (the saved model_outputs.RData may come from a ",
      "run without model checks).",
      call. = FALSE
    )
  }

  best_models <- bind_rows(
    # Best multiple-species model
    outputs$waic_vals %>%
      filter(single_species == FALSE) %>%
      slice_min(
        WAIC,
        n = 1
      ),

    # Best model for each species
    outputs$waic_vals %>%
      filter(single_species == TRUE) %>%
      group_by(Species) %>%
      slice_min(
        WAIC,
        n = 1
      ) %>%
      ungroup()
  )

  if (!is.null(path_to_postProcessing)) {
    write.csv(
      best_models,
      file.path(
        path_to_postProcessing,
        "best_models.csv"
      ),
      row.names = FALSE
    )
  }

  best_models
}


combine_best_occupancy_plots <- function(
  outputs,
  path_to_postProcessing,
  output_dir_name = "combined_species_plots",
  final_filename = "all_species_occupancy_probability.pdf",
  strict = TRUE
) {
  # ===========================================================================
  # 1. Select the best models using WAIC and save them
  # ===========================================================================

  best_models <- select_best_models(
    outputs,
    path_to_postProcessing = path_to_postProcessing
  )

  # Result returned when the plots cannot be combined (strict = FALSE)
  best_models_only <- function(reason) {
    warning(
      reason,
      " The best models were saved, but the occupancy plots were not combined.",
      call. = FALSE
    )
    invisible(list(best_models = best_models))
  }

  # ===========================================================================
  # 2. Find occupancy probability PNG files
  # ===========================================================================
  #
  # Replace `modelFit` with `prediction_outputs` because the prediction plots
  # are stored under the prediction-output directory.
  #
  # Multiple-species model:
  #   Search recursively because species-specific plots are below the model
  #   directory.
  #
  # Single-species model:
  #   Search within:
  #
  #       <model path>/<Species>/plots/
  # ---------------------------------------------------------------------------

  # Predictions are written to "prediction" (see 02_predictions.R);
  # "prediction_outputs" is used as a fallback for older runs.
  resolve_prediction_path <- function(model_path) {
    candidates <- c(
      sub("modelFit", "prediction", model_path, fixed = TRUE),
      sub("modelFit", "prediction_outputs", model_path, fixed = TRUE)
    )

    existing <- candidates[dir.exists(candidates)]

    if (length(existing) > 0) existing[1] else candidates[1]
  }

  # The plots are named occupancy_probability.png and are stored in the
  # plots folder of each species (older runs: <Species>_occupancy_probability.png)
  plot_pattern <- "(^|_)occupancy_probability\\.png$"

  png_files <- unlist(
    lapply(
      seq_len(nrow(best_models)),
      function(i) {
        base_path <- resolve_prediction_path(best_models$path[i])

        # ---------------------------------------------------------------------
        # Multiple-species model
        # ---------------------------------------------------------------------

        if (is.na(best_models$Species[i])) {
          list.files(
            base_path,
            pattern = plot_pattern,
            full.names = TRUE,
            recursive = TRUE,
            ignore.case = TRUE
          )

          # ---------------------------------------------------------------------
          # Single-species model
          # ---------------------------------------------------------------------
        } else {
          plot_path <- file.path(
            base_path,
            best_models$Species[i],
            "plots"
          )

          list.files(
            plot_path,
            pattern = plot_pattern,
            full.names = TRUE,
            recursive = TRUE,
            ignore.case = TRUE
          )
        }
      }
    ),
    use.names = FALSE
  )

  # Remove duplicate paths
  png_files <- unique(png_files)

  # ===========================================================================
  # Check that PNG files were found
  # ===========================================================================

  if (length(png_files) == 0) {
    if (strict) {
      stop("No occupancy probability PNG files were found.", call. = FALSE)
    }

    return(
      best_models_only("No occupancy probability PNG files were found.")
    )
  }

  # ===========================================================================
  # 3. Create data frame describing the PNG files
  # ===========================================================================

  # The species is the folder above `plots`
  png_df <- data.frame(
    Species = basename(dirname(dirname(png_files))),

    png = png_files,

    model = ifelse(
      grepl(
        "multiple_species_model",
        png_files
      ),
      "multiple",
      "single"
    ),

    stringsAsFactors = FALSE
  )

  # ===========================================================================
  # 4. Order plots
  #
  #    Species
  #    Multiple-species model
  #    Single-species model
  # ===========================================================================

  png_df <- png_df %>%
    arrange(
      Species,
      factor(
        model,
        levels = c(
          "multiple",
          "single"
        )
      )
    )

  # ===========================================================================
  # 5. Check that every species has exactly two plots
  # ===========================================================================

  plot_counts <- png_df %>%
    count(
      Species,
      model
    )

  check_counts <- plot_counts %>%
    pivot_wider(
      names_from = model,
      values_from = n,
      values_fill = 0
    )

  # A model type without any plot has no column after pivoting
  for (col in c("multiple", "single")) {
    if (!col %in% names(check_counts)) {
      check_counts[[col]] <- 0L
    }
  }

  if (
    any(check_counts$multiple != 1) ||
      any(check_counts$single != 1)
  ) {
    print(plot_counts)

    msg <- paste(
      "Each species must have exactly one multiple-species",
      "and one single-species occupancy plot."
    )

    if (strict) {
      stop(msg, call. = FALSE)
    }

    return(best_models_only(msg))
  }

  # ===========================================================================
  # 6. Create output directory
  # ===========================================================================

  output_dir <- file.path(
    path_to_postProcessing,
    output_dir_name
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # ===========================================================================
  # 7. Create one combined PNG per species
  # ===========================================================================

  species_names <- unique(
    png_df$Species
  )

  combined_pages <- lapply(
    species_names,
    function(sp) {
      message(
        "Processing: ",
        sp
      )

      # -----------------------------------------------------------------------
      # Get the two PNG files for this species
      # -----------------------------------------------------------------------

      species_files <- png_df %>%
        filter(
          Species == sp
        ) %>%
        arrange(
          factor(
            model,
            levels = c(
              "multiple",
              "single"
            )
          )
        ) %>%
        pull(png)

      # -----------------------------------------------------------------------
      # Safety check
      # -----------------------------------------------------------------------

      if (length(species_files) != 2) {
        stop(
          "Expected 2 PNGs for ",
          sp,
          ", but found ",
          length(species_files)
        )
      }

      # -----------------------------------------------------------------------
      # Read the PNG files
      # -----------------------------------------------------------------------

      images <- lapply(
        species_files,
        magick::image_read
      )

      # -----------------------------------------------------------------------
      # Resize both plots to the same dimensions
      # -----------------------------------------------------------------------

      images <- lapply(
        images,
        function(img) {
          magick::image_resize(
            img,
            "1200x900!"
          )
        }
      )

      # -----------------------------------------------------------------------
      # Put the two plots side-by-side
      #
      # Left:
      #   Multiple-species model
      #
      # Right:
      #   Single-species model
      # -----------------------------------------------------------------------

      combined <- magick::image_append(
        c(
          images[[1]],
          images[[2]]
        ),
        stack = FALSE
      )

      # -----------------------------------------------------------------------
      # Save combined species plot
      # -----------------------------------------------------------------------

      page_file <- file.path(
        output_dir,
        paste0(
          sp,
          "_combined.png"
        )
      )

      magick::image_write(
        combined,
        path = page_file,
        format = "png"
      )

      message(
        "  Saved: ",
        page_file
      )

      page_file
    }
  )

  # ===========================================================================
  # 8. Combine all species plots into one PDF
  # ===========================================================================

  combined_pages <- unlist(
    combined_pages,
    use.names = FALSE
  )

  # Read all combined PNGs
  all_pages <- lapply(
    combined_pages,
    magick::image_read
  )

  # Join the images into a multi-page image
  all_pages <- magick::image_join(
    all_pages
  )

  # ---------------------------------------------------------------------------
  # Final PDF path
  # ---------------------------------------------------------------------------

  final_pdf <- file.path(
    output_dir,
    final_filename
  )

  # ---------------------------------------------------------------------------
  # Write the multi-page PDF
  # ---------------------------------------------------------------------------

  magick::image_write(
    all_pages,
    path = final_pdf,
    format = "pdf"
  )

  message(
    "\nFinal PDF saved to:\n",
    final_pdf
  )

  # ===========================================================================
  # Return results
  # ===========================================================================

  invisible(
    list(
      best_models = best_models,
      png_df = png_df,
      plot_counts = plot_counts,
      combined_pages = combined_pages,
      final_pdf = final_pdf
    )
  )
}
