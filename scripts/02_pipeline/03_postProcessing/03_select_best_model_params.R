# ============================================================
# Extract model outputs
# ============================================================

model_outputs_file <- file.path(path_to_postProcessing, "model_outputs.RData")

extract_outputs_for_selection <- function() {
  # The best models are selected with WAIC, which is only extracted when the
  # posterior predictive checks are included (modelChecks = TRUE)
  if (!isTRUE(modelChecks)) {
    stop(
      "Selecting the best models needs WAIC: set `modelChecks <- TRUE`.",
      call. = FALSE
    )
  }

  if (!compareMultipleOutputs) {
    extract_model_outputs_with_archive(path_to_results, modelChecks)
  } else {
    # If compare multiple outputs, load all the files for the ecosystem
    paths_df <- build_model_paths(
      output_root = output_root,
      version_labels = paste0("version_0", 1:2),
      ecosystemType = ecosystemType,
      existing_only = TRUE,
      include_single_species = c(TRUE, FALSE),
      include_spatial_model = c(TRUE, FALSE)
    )

    extract_model_outputs_with_archive(paths_df$path_to_results, modelChecks)
  }
}

if (file.exists(model_outputs_file)) {
  load(model_outputs_file)

  # Outputs saved from a run without model checks, or when no model files
  # were found, have no WAIC values
  if (!is.data.frame(outputs$waic_vals) || nrow(outputs$waic_vals) == 0) {
    message("Saved model outputs have no WAIC values; extracting them again.")
    outputs <- extract_outputs_for_selection()
    save(outputs, file = model_outputs_file)
  }
} else {
  outputs <- extract_outputs_for_selection()

  # Save extracted outputs
  save(outputs, file = model_outputs_file)
}

if (!is.data.frame(outputs$waic_vals) || nrow(outputs$waic_vals) == 0) {
  # Do not keep an empty cache: it would be reused by the next run
  unlink(model_outputs_file)

  stop(
    "No model results were found, so the best models cannot be selected. ",
    "Check that the model files (.rds) are in the modelFit folder and have ",
    "not been moved to the archive folder:\n",
    paste(
      if (compareMultipleOutputs) {
        file.path(output_root, "modelFit")
      } else {
        path_to_results
      },
      collapse = "\n"
    ),
    call. = FALSE
  )
}


# =============================================================================
# Configuration
# =============================================================================

# Path containing the extracted model outputs
# Modify this if needed.
# ------------------------------------------------------------------------------

# path_to_postProcessing <- ...

# Output directory
# ------------------------------------------------------------------------------

output_dir_name <- "combined_species_plots"

final_filename <- "all_species_occupancy_probability.pdf"


# =============================================================================
# Main execution
# =============================================================================

message(
  "\n============================================================\n",
  "Combining best occupancy probability plots\n",
  "============================================================"
)


combined_occupancy_plots <- combine_best_occupancy_plots(
  outputs = outputs,
  path_to_postProcessing = path_to_postProcessing,
  output_dir_name = output_dir_name,
  final_filename = final_filename,
  # The best models are saved even when the prediction plots do not exist yet
  # (02_predictions.R sources this script to get best_models.csv)
  strict = FALSE
)


# =============================================================================
# Print final output
# =============================================================================

message(
  "\n============================================================\n",
  "Completed\n",
  "============================================================"
)

if (!is.null(combined_occupancy_plots$final_pdf)) {
  message(
    "Final PDF:\n",
    combined_occupancy_plots$final_pdf
  )

  message(
    "\nNumber of species: ",
    length(
      unique(
        combined_occupancy_plots$png_df$Species
      )
    )
  )
}
