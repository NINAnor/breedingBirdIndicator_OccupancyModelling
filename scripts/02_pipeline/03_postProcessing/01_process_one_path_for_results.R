# =============================================================================
# POST-PROCESSING: BREEDING BIRD OCCUPANCY MODELS
#
# This script:
#   1. Extracts model outputs
#   2. Extracts environmental covariates at bird survey locations
#   3. Loads the selected/best models
#   4. Creates model diagnostic and parameter-effect plots
#   5. Validates expected richness against occupancy covariates
#   6. Reads and plots scaled indicator outputs
#
# =============================================================================

# =============================================================================
# 01. LOAD PACKAGES
# =============================================================================

library(sf)
library(dplyr)
library(ggplot2)
library(terra)
library(viridis)
library(readr)
library(sfarrow)


# =============================================================================
# 02. HELPER FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# Convert a Windows path to its short 8.3-format equivalent.
#
# This is useful when working with the long OneDrive paths used in this
# project, which can exceed the path-length limits of some R functions.
# -----------------------------------------------------------------------------

short_path <- function(path) {
  result <- system(
    sprintf(
      'cmd /c for %%I in ("%s") do @echo %%~sI',
      path
    ),
    intern = TRUE
  )

  result <- trimws(result)

  if (length(result) != 1 || !nzchar(result)) {
    stop(
      "Could not convert path to Windows short path:\n",
      path,
      call. = FALSE
    )
  }

  result
}


# -----------------------------------------------------------------------------
# Common plotting theme used for model-output figures.
# -----------------------------------------------------------------------------

theme_model <- function() {
  ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5
      ),
      axis.title = ggplot2::element_text(
        face = "bold"
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      legend.title = ggplot2::element_text(
        face = "bold"
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(
        linewidth = 0.2
      )
    )
}


# -----------------------------------------------------------------------------
# Save a ggplot in both PNG and PDF formats.
# -----------------------------------------------------------------------------

save_plot <- function(
  plot_obj,
  filename,
  width = 14,
  height = 10
) {
  ggplot2::ggsave(
    filename = file.path(
      path_to_postProcessing,
      paste0(filename, ".png")
    ),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )

  ggplot2::ggsave(
    filename = file.path(
      path_to_postProcessing,
      paste0(filename, ".pdf")
    ),
    plot = plot_obj,
    width = width,
    height = height
  )
}


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

if (!file.exists(model_outputs_file)) {
  if (!compareMultipleOutputs) {
    outputs <- extract_model_outputs(
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

    outputs <- extract_model_outputs(
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


# =============================================================================
# 04. EXTRACT COVARIATES AT BIRD SURVEY LOCATIONS
# =============================================================================
#
# Extract environmental covariates from the raster stack at each bird survey
# location. These data are used for model fitting and subsequent validation.
# -----------------------------------------------------------------------------

cov_stack_unscaled <- terra::rast(
  file.path(
    dataFolder,
    "processedData",
    "covariate_stack",
    paste0(
      "covariate_stack_unscaled_",
      ecosystemType,
      ".tif"
    )
  )
)


breeding_bird_data_for_modelling <-
  extract_model_covariates(
    x = ecosystemBreedingBird %>%
      dplyr::select(-DEM),
    modelCrs = modelCrs,
    all_covariates = all_covariates,
    lon_col = "Longitude",
    lat_col = "Latitude",
    input_crs = 4326,
    covariate_stack = cov_stack_unscaled,
    return = "data.frame"
  ) %>%

  # Add route and point identifiers
  dplyr::mutate(
    route_id = sub("_.*", "", Route_point_id),
    point_id = sub(".*_", "", Route_point_id)
  )


# =============================================================================
# 05. CONVERT CATEGORICAL COVARIATES TO FACTORS
# =============================================================================
#
# Only apply these transformations when the corresponding covariate is
# included in the model.
# -----------------------------------------------------------------------------

if ("dominantTree" %in% all_covariates) {
  breeding_bird_data_for_modelling <-
    breeding_bird_data_for_modelling %>%

    dplyr::mutate(
      dominantTree = dplyr::if_else(
        is.na(dominantTree),
        4,
        dominantTree
      ),
      dominantTree = factor(dominantTree)
    )
}


if ("glad_landcover" %in% all_covariates) {
  breeding_bird_data_for_modelling <-
    breeding_bird_data_for_modelling %>%

    dplyr::mutate(
      glad_landcover = factor(glad_landcover)
    )
}


if ("forestType" %in% all_covariates) {
  breeding_bird_data_for_modelling <-
    breeding_bird_data_for_modelling %>%

    dplyr::mutate(
      forestType = factor(forestType)
    )
}


# =============================================================================
# 06. LOAD BEST-MODEL INFORMATION
# =============================================================================
#
# best_models.csv contains the selected model for each species/model
# combination and is used to determine which model outputs should be plotted.
# -----------------------------------------------------------------------------

best_models <- readr::read_csv(
  file.path(
    output_root,
    "postProcessing_outputs",
    "multipleOutputs",
    ecosystemType,
    "best_models.csv"
  ),
  show_col_types = FALSE
)


# =============================================================================
# 07. MODEL EFFECT AND DIAGNOSTIC PLOTS
# =============================================================================

# -----------------------------------------------------------------------------
# 07.1 Species-level occupancy effects
# -----------------------------------------------------------------------------

p_species_occ <- plot(
  outputs,
  type = "error",
  components = c("occ_effects"),
  best_models = best_models
) +
  ggplot2::labs(
    title = "Occupancy effects",
    subtitle = ecosystemType
  ) +
  theme_model() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  theme(
    legend.position = "bottom"
  )


# -----------------------------------------------------------------------------
# 07.2 Species-level detection effects
# -----------------------------------------------------------------------------

p_species_det <- plot(
  outputs,
  type = "error",
  components = c("det_effects"),
  best_models = best_models
) +
  ggplot2::labs(
    title = "Detection effects",
    subtitle = ecosystemType
  ) +
  theme_model() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  theme(
    legend.position = "bottom"
  )


# -----------------------------------------------------------------------------
# 07.3 Bayesian posterior predictive checks
# -----------------------------------------------------------------------------

if (modelChecks) {
  p_bayes <- plot(
    outputs,
    type = "error",
    components = c("bayes_P"),
    best_models = best_models
  ) +
    ggplot2::labs(
      title = "Posterior predictive checks",
      subtitle = ecosystemType,
      x = "Bayesian p-value"
    ) +
    theme_model()
}


# -----------------------------------------------------------------------------
# 07.4 Community-level occupancy parameters
# -----------------------------------------------------------------------------

p_comm_occ <- plot(
  outputs,
  type = "error",
  components = c("beta_summary"),
  best_models = best_models
) +
  ggplot2::labs(
    title = "Community-level parameter summaries",
    subtitle = ecosystemType,
    x = "Posterior mean (95% credible interval)"
  ) +
  theme_model() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  theme(
    legend.position = "bottom"
  )


# -----------------------------------------------------------------------------
# 07.5 Community-level detection parameters
# -----------------------------------------------------------------------------

p_comm_det <- plot(
  outputs,
  type = "error",
  components = c("alpha_summary"),
  best_models = best_models
) +
  ggplot2::labs(
    title = "Community-level parameter summaries",
    subtitle = ecosystemType,
    x = "Posterior mean (95% credible interval)"
  ) +
  theme_model() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  theme(
    legend.position = "bottom"
  )


# =============================================================================
# 08. SAVE MODEL FIGURES
# =============================================================================

save_plot(
  p_species_occ,
  "occupancy_effects"
)

save_plot(
  p_species_det,
  "detection_effects"
)

if (modelChecks) {
  save_plot(
    p_bayes,
    "posterior_predictive_checks"
  )
}

save_plot(
  p_comm_occ,
  "occupancy_community_parameter_summaries"
)

save_plot(
  p_comm_det,
  "detection_community_parameter_summaries"
)


# =============================================================================
# 09. DISPLAY MODEL FIGURES
# =============================================================================

if (modelChecks) {
  print(p_bayes)
}

print(p_comm_occ)
print(p_comm_det)
print(p_species_occ)
print(p_species_det)


# =============================================================================
# 10. EXPECTED RICHNESS VALIDATION
# =============================================================================
#
# Compare expected richness from the indicator outputs with occupancy
# covariates at the bird survey/prediction locations.
#
# When compareMultipleOutputs = TRUE, expected-richness results are searched
# within the compareMultipleOutputs directory.
# -----------------------------------------------------------------------------

if (compareMultipleOutputs) {
  rds_files <- list.files(
    file.path(
      output_root,
      "indicators_outputs",
      "compareMultipleOutputs"
    ),
    pattern = "\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )
} else {
  rds_files <- list.files(
    file.path(
      output_root,
      "indicators_outputs"
    ),
    pattern = "\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )
}


# -----------------------------------------------------------------------------
# Keep only expected-richness result files
# -----------------------------------------------------------------------------

exp_rds_files <- rds_files[
  grepl(
    "[/\\\\]expected_richness[/\\\\]all_results\\.rds$",
    rds_files
  )
]


# -----------------------------------------------------------------------------
# Load expected-richness results
# -----------------------------------------------------------------------------

exp_richness_rds_files <- lapply(
  exp_rds_files,
  readRDS
)


# -----------------------------------------------------------------------------
# Create spatial prediction locations
#
# These locations are based on the bird modelling dataset and are transformed
# to the model CRS.
# -----------------------------------------------------------------------------

predCoords <- breeding_bird_data_for_modelling %>%
  sf::st_as_sf(
    coords = c("long", "lat"),
    crs = modelCrs
  ) %>%
  dplyr::select(
    dplyr::all_of(occ_covs)
  )


# -----------------------------------------------------------------------------
# Extract expected-richness median values at prediction locations
# -----------------------------------------------------------------------------

median_values <- lapply(
  exp_richness_rds_files,
  function(x) {
    expected <- x$expected_richness

    expected_xy <- sf::st_coordinates(
      expected
    )

    expected_vect <- terra::vect(
      expected
    )

    # Convert expected-richness prediction points to a raster.
    # A 2-km resolution is used for interpolation to the prediction locations.
    median_raster <- terra::rasterize(
      expected_vect,
      terra::rast(
        xmin = min(expected_xy[, 1]),
        xmax = max(expected_xy[, 1]),
        ymin = min(expected_xy[, 2]),
        ymax = max(expected_xy[, 2]),
        resolution = 2000,
        crs = sf::st_crs(expected)$wkt
      ),
      field = "median"
    )

    # Extract expected richness at bird prediction locations
    terra::extract(
      median_raster,
      terra::vect(predCoords),
      method = "bilinear"
    )[, 2]
  }
)


# -----------------------------------------------------------------------------
# Derive names for the expected-richness model outputs.
#
# Example:
#   .../multiple_species/expected_richness/all_results.rds
#       -> multiple_species
#
#   .../single_species/expected_richness/all_results.rds
#       -> single_species
# -----------------------------------------------------------------------------

model_names <- basename(
  dirname(
    dirname(exp_rds_files)
  )
)


# -----------------------------------------------------------------------------
# Add expected-richness predictions to prediction locations
# -----------------------------------------------------------------------------

predCoords <- dplyr::bind_cols(
  predCoords,
  setNames(
    median_values,
    model_names
  )
)


# -----------------------------------------------------------------------------
# Validate expected richness against occupancy covariates
# -----------------------------------------------------------------------------

validationPlots <- lapply(
  model_names,
  function(x) {
    ret <- plot_occupancy_covariates(
      data = predCoords,
      occ_covs = occ_covs,
      occupancy_col = x,
      ncol = 2,
      point_alpha = 0.2,
      point_size = 0.8,
      se = TRUE
    )

    # Save validation plot as PNG
    ggsave(
      filename = file.path(
        path_to_postProcessing,
        paste0(x, "_richCovRel.png")
      ),
      plot = ret,
      units = "in",
      width = 9,
      height = 9
    )

    # Save validation plot as PDF
    ggsave(
      filename = file.path(
        path_to_postProcessing,
        paste0(x, "_richCovRel.pdf")
      ),
      plot = ret,
      units = "in",
      width = 9,
      height = 9
    )

    ret
  }
)


# =============================================================================
# 11. PLOT FINAL INDICATOR VALUES
# =============================================================================
#
# Read the final indicator GeoParquet outputs and create maps of:
#   - Expected richness
#   - Scaled indicator
# -----------------------------------------------------------------------------

if (compareMultipleOutputs) {
  parquet_files <- list.files(
    file.path(
      output_root,
      "indicators_outputs",
      "compareMultipleOutputs"
    ),
    pattern = "\\.parquet$",
    full.names = TRUE,
    recursive = TRUE
  )
} else {
  parquet_files <- list.files(
    file.path(
      output_root,
      "indicators_outputs"
    ),
    pattern = "\\.parquet$",
    full.names = TRUE,
    recursive = TRUE
  )
}


# -----------------------------------------------------------------------------
# Keep indicator-value GeoParquet files
# -----------------------------------------------------------------------------

indicator_parquet_files <- parquet_files[
  grepl(
    "[/\\\\]indicator_vals\\.parquet$",
    parquet_files
  )
]


# -----------------------------------------------------------------------------
# Read indicator GeoParquet files
#
# short_path() is used because the project is stored in a long OneDrive path.
# -----------------------------------------------------------------------------

expected_richness_sf <- lapply(
  indicator_parquet_files,
  function(x) {
    sfarrow::st_read_parquet(
      short_path(x)
    )
  }
)


# -----------------------------------------------------------------------------
# Name the expected-richness datasets
# -----------------------------------------------------------------------------

names(expected_richness_sf) <- c(
  "multiple_species",
  "single_species"
)


# =============================================================================
# 12. SAVE INDICATOR MAPS
# =============================================================================

plot_dir <- file.path(
  path_to_postProcessing,
  "indicator_plots"
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# -----------------------------------------------------------------------------
# Plot expected richness and scaled indicator for each model type
# -----------------------------------------------------------------------------

model_names <- names(
  expected_richness_sf
)


for (i in seq_along(expected_richness_sf)) {
  model_name <- model_names[i]

  expected_sf <- expected_richness_sf[[i]]

  # ---------------------------------------------------------------------------
  # Expected richness
  # ---------------------------------------------------------------------------

  p_richness <- ggplot(
    expected_sf
  ) +

    geom_sf(
      aes(fill = indicator_vals),
      shape = 22,
      size = 1,
      stroke = 0
    ) +

    scale_fill_viridis_c(
      option = "viridis",
      limits = c(2, 13),
      breaks = seq(2, 12, by = 2),
      name = "Expected\nrichness",
      na.value = "transparent"
    ) +

    coord_sf(
      expand = FALSE
    ) +

    labs(
      title = model_name,
      x = "Easting",
      y = "Northing"
    ) +

    theme_minimal(
      base_size = 14
    ) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 16
      ),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  # Save expected-richness map
  ggsave(
    filename = file.path(
      plot_dir,
      paste0(
        model_name,
        "_expected_richness.png"
      )
    ),
    plot = p_richness,
    width = 9,
    height = 10,
    dpi = 300
  )

  # ---------------------------------------------------------------------------
  # Scaled indicator
  # ---------------------------------------------------------------------------

  p_scaled <- ggplot(
    expected_sf
  ) +

    geom_sf(
      aes(fill = scaled_indicator),
      shape = 22,
      size = 1,
      stroke = 0
    ) +

    scale_fill_viridis_c(
      option = "viridis",
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      name = "Scaled\nindicator",
      na.value = "transparent"
    ) +

    coord_sf(
      expand = FALSE
    ) +

    labs(
      title = model_name,
      x = "Easting",
      y = "Northing"
    ) +

    theme_minimal(
      base_size = 14
    ) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 16
      ),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  # Save scaled-indicator map
  ggsave(
    filename = file.path(
      plot_dir,
      paste0(
        model_name,
        "_scaled_indicator.png"
      )
    ),
    plot = p_scaled,
    width = 9,
    height = 10,
    dpi = 300
  )
}


# =============================================================================
# 13. OCCUPANCY COVARIATE EFFECTS
# =============================================================================
#
# Run the dedicated occupancy-covariate-effects post-processing script.
# -----------------------------------------------------------------------------

source(
  "scripts/02_pipeline/03_postProcessing/04_occupancy_covariate_effects.R"
)
