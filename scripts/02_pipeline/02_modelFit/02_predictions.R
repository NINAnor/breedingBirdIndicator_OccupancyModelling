################################################################################
# Generate species occupancy predictions
#
# This script:
#   1. Creates a fine-resolution raster template (100 m).
#   2. Loads all fitted species models from the model output directory.
#   3. Predicts occupancy probabilities at the prediction locations.
#   4. Interpolates predictions onto a raster grid.
#   5. Saves prediction rasters and optional diagnostic plots.
#
# Outputs:
#   - Mean occupancy probability rasters for each species.
#   - Lower and upper credible interval rasters (2.5% and 97.5%).
#   - Species-specific occupancy maps (optional).
################################################################################

################################################################################
# Create raster template for prediction maps
################################################################################

# Predictions are generated on a finer spatial resolution (100 m)
# than that used for model fitting (typically 1 km). This raster serves
# as a common grid onto which species occupancy probabilities are projected.

if (ecosystemType %in% "Mountain") {
  input_file <- file.path(
    geoSpatFolder,
    "LandCover",
    "Norway_MountainAreas",
    "Original",
    "Fjell_norge_publisering.gdb"
  )

  pred_boundary <- get_mountain_layer(
    input_file = input_file,
    layer = "fjell_grid50",
    norway_boundary = norway_boundary,
    modelCrs = modelCrs,
    output_file = file.path(
      dataFolder,
      "processedData",
      "prediction_regions",
      "mountain.parquet"
    )
  )

  baseRasterForPreds <- create_template_raster(
    study_region = pred_boundary,
    resolution = pred_resolution, # 100-m prediction grid
    crs = modelCrs,
    output_file = file.path(
      dataFolder,
      "processedData",
      "prediction_grids",
      paste0("base_raster_for_predictions_", ecosystemType, "_100m.tif")
    ),
    overwrite = FALSE
  )
} else if (ecosystemType %in% "Forest") {
  if (
    !file.exists(file.path(
      dataFolder,
      "processedData",
      "prediction_grids",
      paste0("base_raster_for_predictions_", ecosystemType, "_100m.tif")
    ))
  ) {
    legend <- readxl::read_excel(
      paste0(dataFolder, "/Covariates/Landcover_F_M/legend.xlsx"),
      sheet = "NO_binnedClasses"
    )
    classes <- legend$`lc_GLAD.tif value`[
      legend$`For Forest mask - predictions` == "T"
    ]
    classes <- classes[!is.na(classes)]

    pred_boundary <- get_GLAD_landcover(
      input_file = file.path(
        dataFolder,
        "Covariates",
        "Landcover_F_M",
        "lc2020nor_mergedClasses.tif"
      ),
      model_crs = modelCrs,
      output_file = file.path(
        dataFolder,
        "processedData",
        "GLAD_lc",
        "LCglad.tif"
      )
    ) %>%
      resample_to_baseRaster(
        r = .,
        baseRaster = baseRaster,
        name = "glad_landCover"
      )

    base_template <- pred_boundary
    base_template[!base_template %in% classes] <- NA
    base_template[!is.na(base_template)] <- 1

    terra::writeRaster(
      base_template,
      filename = file.path(
        dataFolder,
        "processedData",
        "prediction_grids",
        paste0("base_raster_for_predictions_", ecosystemType, "_100m.tif")
      )
    )

    baseRasterForPreds <- base_template
  } else {
    baseRasterForPreds <- terra::rast(file.path(
      dataFolder,
      "processedData",
      "prediction_grids",
      paste0("base_raster_for_predictions_", ecosystemType, "_100m.tif")
    ))
  }
}

if (use_best_model_file) {
  # 03_select_best_model_params.R saves best_models.csv in
  # `path_to_postProcessing`: the multipleOutputs folder when models are
  # compared (compareMultipleOutputs = TRUE), otherwise the folder of the
  # current model configuration.
  best_model_file <- file.path(
    path_to_postProcessing,
    "best_models.csv"
  )

  if (!file.exists(best_model_file)) {
    source(here::here(
      "scripts",
      "02_pipeline",
      "03_postProcessing",
      "03_select_best_model_params.R"
    ))
  }

  best_models <- readr::read_csv(
    best_model_file,
    show_col_types = FALSE
  )

  path_to_results <- unique(best_models$path)
  path_to_predictions <- sub(
    "modelFit",
    "prediction_outputs",
    path_to_results,
    fixed = TRUE
  )

  # Order this so that the first is always the multi-species model
  path_to_predictions <- path_to_predictions[
    order(!grepl("multiple_species_model", path_to_predictions))
  ]

  # There are issues with the long name for the output folder.
  # So I created a temporary folder called model_tmp and save the results there.
  # Afterwards, I move it to the required folder.
  path_to_predictions <- sub(
    "^C:/Users/kwaku\\.adjei/OneDrive - NINA/112635 - Videreutvikling av indikatorer for hekkefugl - Documents/ninaProjects/breedingBirdsOutput",
    "C:/model_tmp",
    path_to_predictions
  )

  if (length(path_to_results) != length(path_to_predictions)) {
    stop("Mismatch between number of model results and predicted paths")
  }

  # Now lets save the model outputs we need and then delete the other files
  # to create space for the prediction outputs. Since we are using the best predictive
  # model for the indicators, we do not need to run the predictions for all the
  # iterations

  ## Get the model outputs
  # =============================================================================
  # 03. EXTRACT MODEL OUTPUTS
  # =============================================================================
  #
  # If model_outputs.RData already exists, use the previously extracted outputs.
  # Otherwise, extract them from the model results.
  #
  # When compareMultipleOutputs = TRUE, model paths are generated for:
  #   - version 01 and version 02
  #   - single- and multiple-species models
  #   - spatial and non-spatial models
  # -----------------------------------------------------------------------------

  model_outputs_file <- file.path(
    path_to_postProcessing,
    "model_outputs.RData"
  )

  # A cached file without any model (saved when no model files were found) is
  # not used
  if (file.exists(model_outputs_file)) {
    cached_outputs <- new.env()
    load(model_outputs_file, envir = cached_outputs)

    if (nrow(cached_outputs$outputs$meta) == 0) {
      unlink(model_outputs_file)
    }

    rm(cached_outputs)
  }

  if (!file.exists(model_outputs_file)) {
    if (!compareMultipleOutputs) {
      outputs <- extract_model_outputs_with_archive(
        path_to_results,
        modelChecks
      )
    } else {
      paths_df <- build_model_paths(
        output_root = output_root,
        version_labels = paste0("version_0", 1:2),
        ecosystemType = ecosystemType,
        existing_only = TRUE,
        include_single_species = c(TRUE, FALSE),
        include_spatial_model = c(TRUE, FALSE)
      )

      outputs <- extract_model_outputs_with_archive(
        paths_df$path_to_results,
        modelChecks
      )
    }

    # Save extracted outputs so that they do not need to be recalculated
    save(
      outputs,
      file = model_outputs_file
    )
  } else {
    # Load previously extracted model outputs
    load(model_outputs_file)
  }

  # Now let's delete the files in the best model
  all_files_in_results_folder <- list.files(
    path = file.path(output_root, "modelFit"),
    full.names = TRUE,
    recursive = TRUE
  )

  norm <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

  # Only the files of the current ecosystem are considered. The modelFit
  # folder also holds the models of the other ecosystems, which must not be
  # archived when the best models of this ecosystem are selected.
  all_files_in_results_folder <- all_files_in_results_folder[
    grepl(
      paste0("/", ecosystemType, "/"),
      norm(all_files_in_results_folder),
      fixed = TRUE
    )
  ]

  all_files_norm <- norm(all_files_in_results_folder)
  keep_dirs <- paste0(sub("/+$", "", norm(path_to_results)), "/")

  # Safety check: the folders of the best models must contain model files,
  # otherwise every model file would be treated as "not the best model"
  keep_has_models <- vapply(
    path_to_results,
    function(d) {
      length(list.files(d, pattern = "\\.rds$", recursive = TRUE)) > 0
    },
    logical(1)
  )

  if (!any(keep_has_models)) {
    stop(
      "None of the best-model folders contain model files; nothing was ",
      "archived:\n",
      paste(path_to_results, collapse = "\n"),
      call. = FALSE
    )
  }

  in_best <- vapply(
    all_files_norm,
    function(f) any(startsWith(f, keep_dirs)),
    logical(1)
  )

  files_not_in_results <- all_files_in_results_folder[!in_best]

  # Only archive if there is something left to archive. On a re-run the
  # non-best files have already been moved, so this block is skipped.
  if (length(files_not_in_results) == 0) {
    message("No files to archive; skipping archiving step.")
  } else {
    dir.create(
      file.path(archive_root, ecosystemType),
      recursive = TRUE,
      showWarnings = FALSE
    )

    # Identify the common root directory
    root_dir <- file.path(
      output_root,
      "modelFit"
    )

    # Generate relative paths
    relative_paths <- substring(
      files_not_in_results,
      nchar(root_dir) + 2
    )

    # Destination paths
    destination_files <- file.path(
      archive_root,
      relative_paths
    )

    # Create destination directories
    invisible(lapply(
      unique(dirname(destination_files)),
      dir.create,
      recursive = TRUE,
      showWarnings = FALSE
    ))

    # Copy only files not already archived (e.g. after an interrupted run)
    already_archived <- file.exists(destination_files)
    copied <- already_archived
    copied[!already_archived] <- file.copy(
      files_not_in_results[!already_archived],
      destination_files[!already_archived],
      overwrite = FALSE
    )

    # Only delete the original if the archived copy has the same size
    ok <- copied &
      file.exists(destination_files) &
      file.size(destination_files) == file.size(files_not_in_results)

    removed <- rep(FALSE, length(files_not_in_results))
    removed[ok] <- file.remove(files_not_in_results[ok])

    # Check results
    archive_summary <- data.frame(
      file = basename(files_not_in_results),
      copied = copied,
      verified = ok,
      removed = removed
    )
    print(archive_summary)
  }
}
################################################################################
# Generate occupancy predictions for all fitted species
################################################################################

# The function automatically:
#
#   - Reads all model objects (*.rds) in 'path_to_results'.
#   - Generates occupancy predictions for every species.
#   - Computes posterior summaries and credible intervals.
#   - Saves prediction rasters to 'path_to_predictions'.
#   - Optionally creates species-specific occupancy maps.
#
# Predictions are generated using the fixed effects only
# (ignore_RE = TRUE), which removes route-level random effects and
# provides smoother landscape-scale occupancy estimates.

# The function handles one model directory at a time, but `path_to_results`
# can hold several (one per selected best-model path), so loop over them.
# Results are already written to disk, so they are not kept in memory, and
# memory is released after each path.
for (i in seq_along(path_to_results)) {
  input_path <- path_to_results[i]
  output_path <- path_to_predictions[i]

  message(
    "Predicting path ",
    i,
    " of ",
    length(path_to_results),
    ": ",
    input_path
  )

  invisible(
    generate_occupancy_predictions_for_each_species_with_blocks(
      ##########################################################################
      # Input model directory
      ##########################################################################

      input_dir = input_path,

      ##########################################################################
      # Prediction data
      ##########################################################################

      prediction_points = predictionPoints_for_modelling,
      prediction_type = "both",

      det_covariates = det_covs,
      occ_covariates = occ_covs,

      ##########################################################################
      # Spatial information
      ##########################################################################

      base_raster = baseRasterForPreds,
      model_crs = modelCrs,
      ecosystem_type = ecosystemType,

      ##########################################################################
      # Output directory
      ##########################################################################

      output_dir = output_path,

      ##########################################################################
      # Output type
      ##########################################################################

      # Return predictions as raster layers
      return = "raster", #c("sf", "raster"), # options: sf or raster
      block_size = 1000,

      # Posterior credible intervals
      probs = c(0.025, 0.975),

      ##########################################################################
      # Random effects
      ##########################################################################

      # Route random effects were included during model fitting
      use_route_random_effect = use_route_random_effect,

      # Ignore route random effects during prediction
      # to obtain population-level occupancy estimates
      ignore_RE = TRUE,
      verbose = FALSE,

      ##########################################################################
      # Plotting options
      ##########################################################################

      save_plots = TRUE,

      plot_width = 8,
      plot_height = 6,
      plot_dpi = 300,
      n_omp_threads = n.omp.threads
    )
  )

  rm(input_path, output_path)
  gc(verbose = FALSE)
}

# prediction_results <- generate_occupancy_predictions_for_each_species(
#
#   ##########################################################################
#   # Input model directory
#   ##########################################################################
#
#   input_dir = path_to_results,
#
#
#   ##########################################################################
#   # Prediction data
#   ##########################################################################
#
#   prediction_points = predictionPoints_for_modelling[1:10000, ],
#
#   det_covariates = det_covs,
#   occ_covariates = occ_covs,
#
#
#   ##########################################################################
#   # Spatial information
#   ##########################################################################
#
#   base_raster = baseRasterForPreds,
#   model_crs   = modelCrs,
#   ecosystem_type = ecosystemType,
#
#
#   ##########################################################################
#   # Output directory
#   ##########################################################################
#
#   output_dir = path_to_predictions,
#
#
#   ##########################################################################
#   # Output type
#   ##########################################################################
#
#   # Return predictions as raster layers
#   return = "sf", # options: sf or raster
#
#   # Posterior credible intervals
#   probs = c(0.025, 0.975),
#
#
#   ##########################################################################
#   # Random effects
#   ##########################################################################
#
#   # Route random effects were included during model fitting
#   use_route_random_effect = TRUE,
#
#   # Ignore route random effects during prediction
#   # to obtain population-level occupancy estimates
#   ignore_RE = TRUE,
#
#
#   ##########################################################################
#   # Plotting options
#   ##########################################################################
#
#   save_plots = TRUE,
#
#   plot_width  = 8,
#   plot_height = 6,
#   plot_dpi    = 300
# )
