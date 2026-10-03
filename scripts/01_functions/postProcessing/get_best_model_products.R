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

combine_best_occupancy_plots <- function(
  outputs,
  path_to_postProcessing,
  output_dir_name = "combined_species_plots",
  final_filename = "all_species_occupancy_probability.pdf"
) {
  # ===========================================================================
  # 1. Select the best models using WAIC
  # ===========================================================================
  #
  # Lower WAIC indicates a better-fitting model.
  #
  # Multiple-species models:
  #   Select the model with the lowest WAIC.
  #
  # Single-species models:
  #   Select the model with the lowest WAIC separately for each species.
  # ---------------------------------------------------------------------------

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

  # ===========================================================================
  # Save selected models
  # ===========================================================================

  write.csv(
    best_models,
    file.path(
      path_to_postProcessing,
      "best_models.csv"
    ),
    row.names = FALSE
  )

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

  png_files <- unlist(
    lapply(
      seq_len(nrow(best_models)),
      function(i) {
        base_path <- gsub(
          "modelFit",
          "prediction_outputs",
          best_models$path[i]
        )

        # ---------------------------------------------------------------------
        # Multiple-species model
        # ---------------------------------------------------------------------

        if (is.na(best_models$Species[i])) {
          list.files(
            base_path,
            pattern = "_occupancy_probability\\.png$",
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
            pattern = "_occupancy_probability\\.png$",
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
    stop(
      "No occupancy probability PNG files were found."
    )
  }

  # ===========================================================================
  # 3. Create data frame describing the PNG files
  # ===========================================================================

  png_df <- data.frame(
    Species = sub(
      "_occupancy_probability\\.png$",
      "",
      basename(png_files)
    ),

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

  if (
    any(check_counts$multiple != 1) ||
      any(check_counts$single != 1)
  ) {
    print(plot_counts)

    stop(
      paste(
        "Each species must have exactly one multiple-species",
        "and one single-species occupancy plot."
      )
    )
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
