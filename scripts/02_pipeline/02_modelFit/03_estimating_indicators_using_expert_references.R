################################################################################
# Reference values from experts
################################################################################

reference_file <- file.path(
  expertFolder,
  "referenceValues",
  paste0(tolower(ecosystemType), "_reference_vals.csv")
)

if (!file.exists(reference_file)) {
  stop(
    "Expert reference file does not exist:\n",
    reference_file,
    call. = FALSE
  )
}

species_reference_vals <- readr::read_csv2(
  reference_file,
  show_col_types = FALSE
) |>
  dplyr::mutate(
    Values = as.numeric(Values)
  )


################################################################################
# Nature regions
################################################################################

gdb <- file.path(
  geoSpatFolder,
  "AdministrativeUnits",
  "Norway_AdministrativeUnits",
  "Original",
  "Norway_County",
  "versjon2025",
  "Administrative enheter fylker FGDB-format",
  "Basisdata_0000_Norge_25833_Fylker_FGDB",
  "Basisdata_0000_Norge_25833_Fylker_FGDB.gdb"
)

bio_regions <- get_nature_regions(
  input_file = gdb,
  modelCrs = modelCrs,
  output_file = file.path(
    dataFolder,
    "processedData",
    "nature_regions",
    "nature_regions.parquet"
  ),
  overwrite = FALSE,
  plot = TRUE,
  plot_file = file.path(
    dataFolder,
    "processedData",
    "nature_regions",
    "figures"
  )
)


################################################################################
# Function to calculate scaled richness
################################################################################

calculate_indicator <- function(
  prediction_dir,
  output_dir,
  n.boot,
  n.Samps
) {
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  calculate_scaled_ecosystem_richness_with_expert_data(
    input_dir = prediction_dir,

    # Expert reference condition
    reference_data = species_reference_vals,

    # Nature regions for regional summaries
    categorisation_shape = bio_regions,

    ecosystem_type = ecosystemType,

    # Prediction column containing mean occupancy probability
    value_col = indicator_col,
    region_col = "nature_region",

    # Output
    output_dir = output_dir,
    save_outputs = TRUE,
    overwrite = TRUE,

    # Bootstrap
    rep_values_for_each_region = TRUE,
    n.boot = n.boot,
    boot.seed = 123,
    n.Samps = n.Samps
  )
}


################################################################################
# Calculate indicators
################################################################################

if (compareMultipleOutputs) {
  #---------------------------------------------------------------------------
  # Identify best models
  #---------------------------------------------------------------------------

  source(here::here(
    "scripts",
    "02_pipeline",
    "03_postProcessing",
    "03_select_best_model_params.R"
  ))

  best_model_details <- combined_occupancy_plots$best_models

  #---------------------------------------------------------------------------
  # Output directory
  #---------------------------------------------------------------------------

  ind_dir <- file.path(
    output_root,
    "indicators_outputs",
    "compareMultipleOutputs"
  )

  dir.create(
    ind_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  #---------------------------------------------------------------------------
  # Split models into single- and multiple-species models
  #---------------------------------------------------------------------------

  model_groups <- split(
    best_model_details,
    best_model_details$single_species
  )

  #---------------------------------------------------------------------------
  # Process each model type
  #---------------------------------------------------------------------------

  scaled_richness <- lapply(
    names(model_groups),
    function(sppInfo) {
      model_info <- model_groups[[sppInfo]]

      is_single_species <- identical(
        sppInfo,
        "TRUE"
      )

      # Prediction directory
      pp_for_preds <- if (is_single_species) {
        file.path(
          model_info$path,
          model_info$Species
        )
      } else {
        model_info$path
      }

      pp_for_preds <- gsub(
        "modelFit",
        "prediction_outputs",
        pp_for_preds,
        fixed = TRUE
      )

      # Output directory
      pp_for_indicators <- file.path(
        ind_dir,
        if (is_single_species) {
          "single_species"
        } else {
          "multiple_species"
        }
      )

      # Calculate indicator
      calculate_indicator(
        prediction_dir = pp_for_preds,
        output_dir = pp_for_indicators,

        # Keep your original comparison settings
        n.boot = 1000,
        n.Samps = 20000
      )
    }
  )

  # Give results meaningful names
  names(scaled_richness) <- names(model_groups)
} else {
  #---------------------------------------------------------------------------
  # Standard analysis
  #---------------------------------------------------------------------------

  scaled_richness <- calculate_indicator(
    prediction_dir = path_to_predictions,
    output_dir = path_to_indicators,

    # Original standard-analysis settings
    n.boot = 1000,
    n.Samps = 20000
  )
}
