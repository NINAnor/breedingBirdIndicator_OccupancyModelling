# =============================================================================
# Calculate regional and national median species predictions
# =============================================================================

# Required packages
library(sf)
library(dplyr)
library(tidyr)
library(terra)
library(readr)


# =============================================================================
# Helper functions
# =============================================================================

# -----------------------------------------------------------------------------
# Convert a long Windows path to a short 8.3 path
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

  if (length(result) != 1L || !nzchar(result)) {
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
# -----------------------------------------------------------------------------

get_prediction_files <- function(model_path, species = NA_character_) {
  search_path <- if (length(species) == 1L && is.na(species)) {
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
# Rasterize predictions at 2 km resolution, then extract at prediction
# locations using bilinear interpolation.
# -----------------------------------------------------------------------------

extract_median_prediction <- function(rds_file, predCoords_vect) {
  dat <- readRDS(short_path(rds_file))

  if (!inherits(dat, "sf")) {
    stop("Expected an sf object in: ", rds_file)
  }

  if (!"median" %in% names(dat)) {
    stop("Column 'median' not found in: ", rds_file)
  }

  if (nrow(dat) == 0L) {
    stop("Prediction file contains no rows: ", rds_file)
  }

  expected <- dat |>
    dplyr::select(median)

  expected_xy <- sf::st_coordinates(expected)

  expected_vect <- terra::vect(expected)

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

  extracted <- terra::extract(
    median_raster,
    predCoords_vect,
    method = "bilinear"
  )

  extracted[[2]]
}


# -----------------------------------------------------------------------------
# Calculate regional and national medians
# -----------------------------------------------------------------------------

calculate_species_medians <- function(pred_vals, species) {
  dat <- pred_vals |>
    sf::st_drop_geometry() |>
    dplyr::select(nature_region, dplyr::all_of(species))

  # Regional medians
  regional <- dat |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(species),
      names_to = "species",
      values_to = "value"
    ) |>
    dplyr::filter(!is.na(nature_region)) |>
    dplyr::group_by(species, nature_region) |>
    dplyr::summarise(
      median = if (all(is.na(value))) {
        NA_real_
      } else {
        median(value, na.rm = TRUE)
      },
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      names_from = nature_region,
      values_from = median
    )

  # National median across all prediction points, including points
  # without an assigned nature region
  national <- dat |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(species),
      names_to = "species",
      values_to = "value"
    ) |>
    dplyr::group_by(species) |>
    dplyr::summarise(
      National = if (all(is.na(value))) {
        NA_real_
      } else {
        median(value, na.rm = TRUE)
      },
      .groups = "drop"
    )

  # Combine regional and national medians
  regional |>
    dplyr::left_join(national, by = "species") |>
    dplyr::select(
      species,
      dplyr::any_of(c("Central", "East", "North", "South", "West")),
      National
    )
}


# =============================================================================
# Prepare prediction paths and locate prediction files
# =============================================================================

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

predCoords <- breeding_bird_data_for_modelling %>%
  sf::st_as_sf(
    coords = c("long", "lat"),
    crs = modelCrs
  ) %>%
  dplyr::select(
    dplyr::all_of(occ_covs)
  )


best_models$prediction_path <- gsub(
  "modelFit",
  "prediction_outputs",
  best_models$path,
  fixed = TRUE
)

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
# Load nature regions
# =============================================================================

gdb <- paste0(
  geoSpatFolder,
  "/AdministrativeUnits/",
  "Norway_AdministrativeUnits/Original/Norway_County/",
  "versjon2025/Administrative enheter fylker FGDB-format/",
  "Basisdata_0000_Norge_25833_Fylker_FGDB/",
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


# =============================================================================
# Prepare prediction coordinates and assign nature regions
# =============================================================================

predCoords_vect <- terra::vect(predCoords)

bio_regions_vect <- terra::vect(bio_regions)

if (!terra::same.crs(bio_regions_vect, predCoords_vect)) {
  bio_regions_vect <- terra::project(
    bio_regions_vect,
    terra::crs(predCoords_vect)
  )
}

region_id <- terra::extract(
  bio_regions_vect["nature_region"],
  predCoords_vect
)

if (nrow(region_id) != nrow(predCoords)) {
  stop("The number of extracted regions does not match prediction points.")
}

predCoords$nature_region <- region_id$nature_region


# =============================================================================
# Process multi-species and single-species models
# =============================================================================

median_results <- list()

for (i in unique(best_models$single_species)) {
  cat(
    "\n==============================================\n",
    "Processing single_species = ",
    i,
    "\n==============================================\n",
    sep = ""
  )

  model_rows <- best_models |>
    dplyr::filter(single_species == i)

  # ---------------------------------------------------------------------------
  # Multi-species model
  # ---------------------------------------------------------------------------

  if (!isTRUE(i)) {
    if (nrow(model_rows) == 0L) {
      next
    }

    # The original workflow expects one multi-species model
    rds_files <- model_rows$prediction_files[[1]]

    if (length(rds_files) == 0L) {
      warning("No prediction RDS files found for the multi-species model.")
      next
    }

    species <- basename(dirname(rds_files))

    if (anyDuplicated(species)) {
      stop(
        "Duplicate species names found among multi-species prediction files. ",
        "Resolve these before combining predictions."
      )
    }

    cat(
      "Extracting ",
      length(rds_files),
      " multi-species prediction files.\n",
      sep = ""
    )

    predictions <- lapply(
      rds_files,
      extract_median_prediction,
      predCoords_vect = predCoords_vect
    )

    prediction_df <- as.data.frame(
      stats::setNames(predictions, species),
      check.names = FALSE
    )

    if (nrow(prediction_df) != nrow(predCoords)) {
      stop("Prediction rows do not match the number of prediction coordinates.")
    }

    pred_vals <- dplyr::bind_cols(
      predCoords,
      prediction_df
    )

    medians <- calculate_species_medians(
      pred_vals = pred_vals,
      species = species
    )

    # Preserve model type so results from both model types are retained
    medians$model_type <- "multi_species"

    median_results[[length(median_results) + 1L]] <- medians

    rm(predictions, prediction_df, pred_vals, medians)
    gc()
  } else {
    # -------------------------------------------------------------------------
    # Single-species models
    # -------------------------------------------------------------------------

    for (j in seq_len(nrow(model_rows))) {
      current_species <- model_rows$Species[j]

      cat(
        "Processing single-species model ",
        j,
        " of ",
        nrow(model_rows),
        ": ",
        current_species,
        "\n",
        sep = ""
      )

      rds_files <- model_rows$prediction_files[[j]]

      if (length(rds_files) == 0L) {
        warning(
          "No prediction RDS files found for: ",
          current_species
        )
        next
      }

      # A single-species model should provide one species prediction
      predictions <- lapply(
        rds_files,
        extract_median_prediction,
        predCoords_vect = predCoords_vect
      )

      prediction_names <- basename(dirname(rds_files))

      if (length(unique(prediction_names)) != 1L) {
        stop(
          "Multiple species directories found for single-species model: ",
          current_species
        )
      }

      prediction_df <- as.data.frame(
        stats::setNames(predictions, prediction_names),
        check.names = FALSE
      )

      if (nrow(prediction_df) != nrow(predCoords)) {
        stop(
          "Prediction rows do not match prediction coordinates for: ",
          current_species
        )
      }

      pred_vals <- dplyr::bind_cols(
        predCoords,
        prediction_df
      )

      medians <- calculate_species_medians(
        pred_vals = pred_vals,
        species = prediction_names
      )

      medians$model_type <- "single_species"

      median_results[[length(median_results) + 1L]] <- medians

      rm(predictions, prediction_df, pred_vals, medians)
      gc()
    }
  }
}


# =============================================================================
# Combine all results and save one CSV file
# =============================================================================

if (length(median_results) == 0L) {
  stop("No median results were generated. Check the prediction files.")
}

all_species_medians <- dplyr::bind_rows(median_results) |>
  dplyr::select(
    model_type,
    species,
    dplyr::any_of(c("Central", "East", "North", "South", "West")),
    National
  ) |>
  dplyr::arrange(model_type, species)

# Check whether a species occurs in both model types
duplicate_species <- all_species_medians |>
  dplyr::count(species) |>
  dplyr::filter(n > 1L)

if (nrow(duplicate_species) > 0L) {
  warning(
    "Some species occur in both model types. They are retained as separate ",
    "rows and identified by the model_type column."
  )
}

# Output path
output_file <- file.path(
  output_root,
  "median_vals",
  paste0(ecosystemType, "_medians.xlsx")
)

# Create output directory if necessary
dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

# Save as Excel workbook
openxlsx::write.xlsx(
  all_species_medians,
  file = output_file,
  overwrite = TRUE,
  na.string = ""
)

cat(
  "\nMedian results saved to:\n",
  normalizePath(output_file, winslash = "/", mustWork = FALSE),
  "\nNumber of result rows: ",
  nrow(all_species_medians),
  "\n",
  sep = ""
)

print(all_species_medians)
