# =============================================================================
# Load packages
# =============================================================================

library(sf)
library(dplyr)
library(ggplot2)
library(viridis)
library(terra)


# =============================================================================
# Helper functions
# =============================================================================

# -----------------------------------------------------------------------------
# Convert a long Windows path to a short 8.3 path
#
# This is useful for readRDS() when Windows path-length limitations cause
# problems with long OneDrive/project paths.
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
# Find prediction RDS files for one model
#
# For:
#   single_species = FALSE:
#       search the complete prediction directory
#
#   single_species = TRUE:
#       search only the directory belonging to that species
# -----------------------------------------------------------------------------

get_prediction_files <- function(model_path, species = NA_character_) {
  search_path <- if (is.na(species)) {
    model_path
  } else {
    file.path(model_path, species)
  }

  if (!dir.exists(search_path)) {
    warning(
      "Prediction directory does not exist:\n",
      search_path,
      call. = FALSE
    )

    return(character(0))
  }

  list.files(
    search_path,
    pattern = "\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )
}


# -----------------------------------------------------------------------------
# Extract median predictions from one prediction RDS
#
# The prediction points are rasterized at 2 km resolution and then extracted
# at the bird prediction locations using bilinear interpolation.
# -----------------------------------------------------------------------------

extract_median_prediction <- function(
  rds_file,
  predCoords_vect
) {
  # Read prediction object
  dat <- readRDS(short_path(rds_file))

  # Keep only the median prediction. Prediction files with both occupancy and
  # detection use `occupancy_median`; older files use `median`.
  value_col <- intersect(c("occupancy_median", "median"), names(dat))[1]

  if (is.na(value_col)) {
    stop(
      "No median prediction column in: ",
      rds_file,
      call. = FALSE
    )
  }

  expected <- dat[value_col]
  names(expected)[names(expected) == value_col] <- "median"

  # Coordinates of prediction points
  expected_xy <- sf::st_coordinates(expected)

  # Convert to terra vector
  expected_vect <- terra::vect(expected)

  # ---------------------------------------------------------------------------
  # Create raster template
  # ---------------------------------------------------------------------------

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

  # ---------------------------------------------------------------------------
  # Extract values at prediction locations
  # ---------------------------------------------------------------------------

  terra::extract(
    median_raster,
    predCoords_vect,
    method = "bilinear"
  )[, 2]
}


# =============================================================================
# Prepare prediction paths
# =============================================================================

# Predictions are written to "prediction" (see 02_predictions.R);
# "prediction_outputs" is used as a fallback for older runs.
resolve_prediction_path <- function(model_path) {
  candidates <- c(
    gsub("modelFit", "prediction", model_path, fixed = TRUE),
    gsub("modelFit", "prediction_outputs", model_path, fixed = TRUE)
  )

  existing <- candidates[dir.exists(candidates)]

  if (length(existing) > 0) existing[1] else candidates[1]
}

best_models$prediction_path <- vapply(
  best_models$path,
  resolve_prediction_path,
  character(1),
  USE.NAMES = FALSE
)


# =============================================================================
# Find prediction RDS files
#
# This is done ONCE for each row of best_models.
#
# The resulting list-column has:
#
#   FALSE row -> all species RDS files for the multiple-species model
#
#   TRUE rows -> RDS files for the corresponding single species
# =============================================================================

best_models$prediction_files <- lapply(
  seq_len(nrow(best_models)),
  function(i) {
    get_prediction_files(
      model_path = best_models$prediction_path[i],
      species = best_models$Species[i]
    )
  }
)


# =============================================================================
# Prepare prediction coordinates once
# =============================================================================

# Convert predCoords to terra only once.
#
# Previously this conversion happened every time an RDS file was processed.

predCoords_vect <- terra::vect(predCoords)


# =============================================================================
# Process the two model types
#
# i = FALSE -> multiple-species model
# i = TRUE  -> single-species models
# =============================================================================

for (i in unique(best_models$single_species)) {
  cat(
    "\n==============================================\n",
    "Processing single_species = ",
    i,
    "\n",
    "==============================================\n",
    sep = ""
  )

  # ---------------------------------------------------------------------------
  # Select models belonging to this group
  # ---------------------------------------------------------------------------

  model_rows <- best_models |>
    filter(single_species == i)

  # ---------------------------------------------------------------------------
  # Define model type
  # ---------------------------------------------------------------------------

  modelType <- if (isTRUE(i)) {
    "single_species"
  } else {
    "multi_species"
  }

  # ---------------------------------------------------------------------------
  # Open one PDF for this model type
  #
  # FALSE:
  #   one multiple-species model -> multiple species pages
  #
  # TRUE:
  #   multiple single-species models -> one species page per model
  # ---------------------------------------------------------------------------

  pdf_file <- file.path(
    path_to_postProcessing,
    paste0(
      modelType,
      "_occ_cov_effects.pdf"
    )
  )

  pdf(
    file = pdf_file,
    width = 9,
    height = 9
  )

  # ===========================================================================
  # MULTIPLE-SPECIES MODEL
  # ===========================================================================

  if (!isTRUE(i)) {
    # -------------------------------------------------------------------------
    # There should be one multiple-species model.
    #
    # prediction_files[[1]] contains all RDS prediction files.
    # -------------------------------------------------------------------------

    rds_files <- model_rows$prediction_files[[1]]

    if (length(rds_files) == 0) {
      warning(
        "No prediction RDS files found for the multiple-species model."
      )

      dev.off()
      next
    }

    # -------------------------------------------------------------------------
    # Extract predictions for all species
    # -------------------------------------------------------------------------

    predictions <- lapply(
      rds_files,
      extract_median_prediction,
      predCoords_vect = predCoords_vect
    )

    # Species names are encoded in the parent directory
    species <- basename(dirname(rds_files))

    # -------------------------------------------------------------------------
    # Combine predictions
    # -------------------------------------------------------------------------

    prediction_df <- as.data.frame(
      setNames(
        predictions,
        species
      )
    )

    pred_vals <- bind_cols(
      predCoords,
      prediction_df
    )

    # -------------------------------------------------------------------------
    # Plot each species
    # -------------------------------------------------------------------------

    for (x in species) {
      ret <- plot_occupancy_covariates(
        data = pred_vals,
        occ_covs = occ_covs,
        occupancy_col = x,
        ncol = 2,
        point_alpha = 0.2,
        point_size = 0.8,
        se = TRUE,
        y_label = "Occupancy probability",
        covariate_names = covariateSelection
      ) +
        patchwork::plot_annotation(
          title = gsub("_", " ", x, fixed = TRUE)
        )

      print(ret)
    }

    # Free potentially large objects before continuing
    rm(
      predictions,
      prediction_df,
      pred_vals
    )

    gc()

    # ===========================================================================
    # SINGLE-SPECIES MODELS
    # ===========================================================================
  } else {
    # -------------------------------------------------------------------------
    # model_rows contains all single-species models.
    #
    # Each row corresponds to one species model.
    # -------------------------------------------------------------------------

    for (j in seq_len(nrow(model_rows))) {
      cat(
        "  Species model ",
        j,
        " of ",
        nrow(model_rows),
        ": ",
        model_rows$Species[j],
        "\n",
        sep = ""
      )

      # -----------------------------------------------------------------------
      # Get RDS files for this single-species model
      # -----------------------------------------------------------------------

      rds_files <- model_rows$prediction_files[[j]]

      if (length(rds_files) == 0) {
        warning(
          "No prediction RDS files found for: ",
          model_rows$Species[j]
        )

        next
      }

      # -----------------------------------------------------------------------
      # Extract predictions
      # -----------------------------------------------------------------------

      predictions <- lapply(
        rds_files,
        extract_median_prediction,
        predCoords_vect = predCoords_vect
      )

      # Species name
      species <- basename(dirname(rds_files))

      # -----------------------------------------------------------------------
      # Combine prediction with coordinates
      # -----------------------------------------------------------------------

      prediction_df <- as.data.frame(
        setNames(
          predictions,
          species
        )
      )

      pred_vals <- bind_cols(
        predCoords,
        prediction_df
      )

      # -----------------------------------------------------------------------
      # Plot species
      # -----------------------------------------------------------------------

      for (x in species) {
        ret <- plot_occupancy_covariates(
          data = pred_vals,
          occ_covs = occ_covs,
          occupancy_col = x,
          ncol = 2,
          point_alpha = 0.2,
          point_size = 0.8,
          se = TRUE,
          covariate_names = covariateSelection,
          y_label = "Occupancy probability"
        ) +
          patchwork::plot_annotation(
            title = gsub("_", " ", x, fixed = TRUE)
          )

        print(ret)
      }

      # -----------------------------------------------------------------------
      # Release large objects before processing the next species
      # -----------------------------------------------------------------------

      rm(
        predictions,
        prediction_df,
        pred_vals
      )

      gc()
    }
  }

  # ---------------------------------------------------------------------------
  # Close PDF
  # ---------------------------------------------------------------------------

  dev.off()

  cat(
    "\nSaved: ",
    pdf_file,
    "\n",
    sep = ""
  )
}
