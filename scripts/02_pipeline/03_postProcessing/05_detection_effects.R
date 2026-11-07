# =============================================================================
# Detection probability against covariates and against occupancy probability
#
# For every species of the best models this script:
#   1. extracts the predicted median occupancy and detection probability at the
#      bird survey locations (2-km raster + bilinear interpolation);
#   2. plots detection probability against each detection covariate
#      (one PDF page per species);
#   3. plots detection probability against occupancy probability, with one
#      panel per species.
#
# It is sourced from 01_process_one_path_for_results.R after
# 04_occupancy_covariate_effects.R and reuses the objects and helper functions
# created there (best_models, short_path, get_prediction_files).
# =============================================================================

if (
  !exists("get_prediction_files") ||
    !exists("short_path") ||
    !"prediction_files" %in% names(best_models)
) {
  stop(
    "Run 04_occupancy_covariate_effects.R first: it creates the prediction ",
    "file list and helper functions used here.",
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# Predicted occupancy and detection at the survey locations for one file
#
# The file is read once and both values are extracted.
# -----------------------------------------------------------------------------

extract_occ_det_at_locations <- function(
  rds_file,
  locations_vect,
  occ_col = "occupancy_median",
  det_col = "detection_median"
) {
  dat <- readRDS(short_path(rds_file))

  missing_cols <- setdiff(c(occ_col, det_col), names(dat))
  if (length(missing_cols) > 0) {
    stop(
      "Column(s) missing in ",
      rds_file,
      ": ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  xy <- sf::st_coordinates(dat)
  template <- terra::rast(
    xmin = min(xy[, 1]),
    xmax = max(xy[, 1]),
    ymin = min(xy[, 2]),
    ymax = max(xy[, 2]),
    resolution = 2000,
    crs = sf::st_crs(dat)$wkt
  )

  pred_vect <- terra::vect(dat[c(occ_col, det_col)])

  extract_one <- function(field) {
    terra::extract(
      terra::rasterize(pred_vect, template, field = field),
      locations_vect,
      method = "bilinear"
    )[, 2]
  }

  list(
    occupancy = extract_one(occ_col),
    detection = extract_one(det_col)
  )
}


# -----------------------------------------------------------------------------
# Survey locations with the detection covariates
#
# Detection covariates that are stored per year (e.g. mean_temperature_2020)
# are averaged over the years.
# -----------------------------------------------------------------------------

det_location_data <- breeding_bird_data_for_modelling

for (cov in det_covs) {
  if (cov %in% names(det_location_data)) {
    next
  }

  year_cols <- grep(
    paste0("^", cov, "_?[0-9]{4}$"),
    names(det_location_data),
    value = TRUE
  )

  if (length(year_cols) == 0) {
    stop(
      "Detection covariate not found in the survey data: ",
      cov,
      call. = FALSE
    )
  }

  det_location_data[[cov]] <- rowMeans(
    det_location_data[, year_cols, drop = FALSE],
    na.rm = TRUE
  )
}

det_coords <- det_location_data |>
  sf::st_as_sf(
    coords = c("long", "lat"),
    crs = modelCrs
  ) |>
  dplyr::select(dplyr::all_of(det_covs))

det_coords_vect <- terra::vect(det_coords)


# -----------------------------------------------------------------------------
# Process the model types (multiple-species and single-species models)
# -----------------------------------------------------------------------------

for (i in unique(best_models$single_species)) {
  modelType <- if (isTRUE(i)) "single_species" else "multi_species"

  cat(
    "\n==============================================\n",
    "Detection effects: ",
    modelType,
    "\n==============================================\n",
    sep = ""
  )

  rds_files <- unlist(
    best_models$prediction_files[best_models$single_species == i]
  )

  if (length(rds_files) == 0) {
    warning(
      "No prediction RDS files found for ",
      modelType,
      "; skipping detection plots.",
      call. = FALSE
    )
    next
  }

  species <- basename(dirname(rds_files))

  # ---- Extract predicted occupancy and detection ----------------------------
  values <- lapply(
    rds_files,
    extract_occ_det_at_locations,
    locations_vect = det_coords_vect
  )
  names(values) <- species

  det_df <- as.data.frame(lapply(values, `[[`, "detection"))
  names(det_df) <- species

  occ_det_long <- dplyr::bind_rows(
    lapply(
      species,
      function(sp) {
        data.frame(
          species = sp,
          occupancy = values[[sp]]$occupancy,
          detection = values[[sp]]$detection
        )
      }
    )
  )

  # ---- Detection probability against each detection covariate ---------------
  det_pred_vals <- dplyr::bind_cols(det_coords, det_df)

  pdf_file <- file.path(
    path_to_postProcessing,
    paste0(modelType, "_detection_covariate_relationships.pdf")
  )

  pdf(file = pdf_file, width = 9, height = 5)

  for (sp in species) {
    ret <- plot_occupancy_covariates(
      data = det_pred_vals,
      occ_covs = det_covs,
      occupancy_col = sp,
      ncol = 2,
      point_alpha = 0.2,
      point_size = 0.8,
      se = TRUE,
      y_label = "Detection probability",
      covariate_names = covariateSelection
    ) +
      patchwork::plot_annotation(title = gsub("_", " ", sp, fixed = TRUE))

    print(ret)
  }

  dev.off()

  # ---- Detection probability against occupancy probability ------------------
  p_det_occ <- plot_detection_vs_occupancy(occ_det_long)

  n_rows <- ceiling(length(species) / 4)

  for (ext in c("png", "pdf")) {
    ggplot2::ggsave(
      filename = file.path(
        path_to_postProcessing,
        paste0(modelType, "_detection_vs_occupancy.", ext)
      ),
      plot = p_det_occ,
      width = 12,
      height = max(4, 3 * n_rows),
      dpi = 300
    )
  }

  rm(values, det_df, occ_det_long, det_pred_vals, p_det_occ)
  gc()
}
