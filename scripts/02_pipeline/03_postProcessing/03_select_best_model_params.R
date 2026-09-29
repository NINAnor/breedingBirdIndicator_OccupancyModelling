# ============================================================
# Extract model outputs
# ============================================================

if (!file.exists(file.path(path_to_postProcessing, "model_outputs.RData"))) {
  if (!compareMultipleOutputs) {
    outputs <- extract_model_outputs(path_to_results, modelChecks)
  } else {
    paths_df <- build_model_paths(
      output_root = output_root,
      version_labels = paste0("version_0", 1:2),
      ecosystemType = ecosystemType,
      existing_only = TRUE,
      include_single_species = c(TRUE, FALSE),
      include_spatial_model = c(TRUE, FALSE)
    )

    outputs <- extract_model_outputs(paths_df$path_to_results, modelChecks)
  }

  # Save extracted outputs
  save(
    outputs,
    file = file.path(path_to_postProcessing, "model_outputs.RData")
  )
} else {
  load(file.path(path_to_postProcessing, "model_outputs.RData"))
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
  final_filename = final_filename
)


# =============================================================================
# Print final output
# =============================================================================

message(
  "\n============================================================\n",
  "Completed\n",
  "============================================================"
)

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
