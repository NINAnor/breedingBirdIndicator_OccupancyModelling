#' Calculate Scaled Ecosystem Richness Using Expert Reference Data
#'
#' Calculates multispecies expected richness from species prediction outputs,
#' scales richness against an expert-derived reference value, and optionally
#' estimates representative regional values using bootstrap beta regression.
#'
#' @param input_dir Character. Directory containing species prediction outputs.
#'   Each species should have a subdirectory containing a `.tif` or `.parquet`
#'   prediction file.
#' @param reference_data Data frame containing expert reference values.
#'   Must contain columns `Species` and `Values`.
#' @param categorisation_shape An `sf` object containing the regions used for
#'   regional summaries.
#' @param ecosystem_type Character. Ecosystem type passed to
#'   `calculate_ecosystem_richness()`.
#' @param value_col Character. Prediction layer used to calculate expected
#'   richness. Default is `"mean"`.
#' @param point_spacing Numeric. Sampling spacing in metres.
#' @param region_col Character. Region column name. Default is
#'   `"nature_region"`.
#' @param output_dir Character or NULL. Directory for outputs.
#' @param save_outputs Logical. Save results to disk?
#' @param overwrite Logical. Overwrite existing files?
#' @param rep_values_for_each_region Logical. Estimate representative values
#'   for each region using bootstrap beta regression?
#' @param n.boot Integer. Number of bootstrap iterations.
#' @param boot.seed Integer. Random seed for bootstrap.
#' @param n.Samps Integer. Maximum number of raster cells sampled per
#'   bootstrap replicate.
#'
#' @return A named list containing the species predictions, expected richness,
#'   scaled indicator, regional raster, and bootstrap summaries.
#'
#' @export
calculate_scaled_ecosystem_richness_with_expert_data <- function(
  input_dir,
  reference_data,
  categorisation_shape,
  ecosystem_type,
  value_col = "mean",
  point_spacing = 1000,
  region_col = "nature_region",
  output_dir = NULL,
  save_outputs = FALSE,
  overwrite = TRUE,
  rep_values_for_each_region = TRUE,
  n.boot = 5000,
  boot.seed = 123,
  n.Samps = 5000
) {
  # --------------------------------------------------------------------------
  # 1. Validate inputs
  # --------------------------------------------------------------------------

  if (!is.data.frame(reference_data)) {
    stop(
      "`reference_data` must be a data.frame.",
      call. = FALSE
    )
  }

  if (!inherits(categorisation_shape, "sf")) {
    stop(
      "`categorisation_shape` must be an sf object.",
      call. = FALSE
    )
  }

  if (!dir.exists(input_dir)) {
    stop(
      "`input_dir` does not exist: ",
      input_dir,
      call. = FALSE
    )
  }

  if (save_outputs && is.null(output_dir)) {
    stop(
      "`output_dir` must be supplied when `save_outputs = TRUE`.",
      call. = FALSE
    )
  }

  if (save_outputs) {
    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }

  # Required columns in expert reference data
  required_reference_cols <- c("Species", "Values")

  missing_reference_cols <- setdiff(
    required_reference_cols,
    names(reference_data)
  )

  if (length(missing_reference_cols) > 0) {
    stop(
      "Missing columns in `reference_data`: ",
      paste(missing_reference_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Required region column
  if (!region_col %in% names(categorisation_shape)) {
    stop(
      "`region_col` ('",
      region_col,
      "') is not present in `categorisation_shape`.",
      call. = FALSE
    )
  }

  # --------------------------------------------------------------------------
  # 2. Calculate multispecies expected richness
  # --------------------------------------------------------------------------

  cat("\nCalculating multispecies expected richness...\n")

  multispecies_richness <- calculate_ecosystem_richness(
    input_dir = input_dir,
    ecosystem_type = ecosystem_type,
    value_col = value_col,
    na.rm = TRUE,
    output_dir = if (save_outputs) {
      file.path(output_dir, "expected_richness")
    } else {
      NULL
    },
    save_outputs = save_outputs,
    overwrite = overwrite
  )

  # Extract species results
  results <- multispecies_richness$results

  if (length(results) == 0) {
    stop(
      "No valid species prediction outputs were found.",
      call. = FALSE
    )
  }

  # --------------------------------------------------------------------------
  # 3. Match expert species to modelled species
  # --------------------------------------------------------------------------

  cat("\nCalculating expert reference value...\n")

  species_cols <- safe_name(names(results))

  reference_data_filtered <- reference_data |>
    dplyr::mutate(
      Species_col = gsub(
        " ",
        "_",
        Species
      )
    ) |>
    dplyr::filter(
      Species_col %in% species_cols
    )

  if (nrow(reference_data_filtered) == 0) {
    stop(
      "None of the species in `reference_data` matched the modelled species.",
      call. = FALSE
    )
  }

  # Check for species with missing expert values
  if (any(!is.finite(reference_data_filtered$Values))) {
    stop(
      "Expert reference values contain NA, NaN, Inf, or -Inf.",
      call. = FALSE
    )
  }

  # Expert-derived reference richness
  reference_predictions <- sum(
    reference_data_filtered$Values,
    na.rm = TRUE
  )

  if (
    !is.finite(reference_predictions) ||
      reference_predictions <= 0
  ) {
    stop(
      "The expert reference richness must be greater than zero.",
      call. = FALSE
    )
  }

  cat(
    "Expert reference richness:",
    round(reference_predictions, 4),
    "\n"
  )

  # --------------------------------------------------------------------------
  # 4. Extract expected richness
  # --------------------------------------------------------------------------

  expected_richness <- multispecies_richness$expected_richness

  if (!inherits(expected_richness, "SpatRaster")) {
    stop(
      "`multispecies_richness$expected_richness` must be a SpatRaster.",
      call. = FALSE
    )
  }

  if (!value_col %in% names(expected_richness)) {
    stop(
      "Layer `",
      value_col,
      "` was not found in `expected_richness`.",
      call. = FALSE
    )
  }

  # --------------------------------------------------------------------------
  # 5. Scale richness against expert reference
  # --------------------------------------------------------------------------

  cat("\nScaling richness against expert reference...\n")

  scaled_richness <- (expected_richness[[value_col]] /
    reference_predictions)

  # Restrict to [0, 1]
  scaled_richness <- terra::clamp(
    scaled_richness,
    lower = 0,
    upper = 1,
    values = TRUE
  )

  # Give layers clear names
  indicator_vals <- c(
    expected_richness[[value_col]],
    scaled_richness
  )

  names(indicator_vals) <- c(
    "indicator_vals",
    "scaled_indicator"
  )

  # --------------------------------------------------------------------------
  # 6. Add biological/nature regions to raster
  # --------------------------------------------------------------------------

  cat("\nRasterising regions...\n")

  regions_vect <- terra::vect(categorisation_shape)

  # Numeric ID for each region
  regions_vect$region_id <- seq_len(
    nrow(regions_vect)
  )

  # Rasterise using same grid as indicator
  region_raster <- terra::rasterize(
    regions_vect,
    indicator_vals[[1]],
    field = "region_id"
  )

  # Add region names as factor levels
  region_levels <- data.frame(
    ID = seq_len(nrow(regions_vect)),
    nature_region = as.character(
      regions_vect[[region_col]]
    )
  )

  levels(region_raster) <- region_levels

  names(region_raster) <- "bio_region"

  # Add region raster
  indicator_vals <- c(
    indicator_vals,
    region_raster
  )

  # --------------------------------------------------------------------------
  # 7. Bootstrap representative values for regions
  # --------------------------------------------------------------------------

  representative_indicator_vals <- NULL

  if (rep_values_for_each_region) {
    cat("\nEstimating representative values for each region...\n")

    # Helper for extracting one region
    extract_region_values <- function(
      indicator_vals,
      regions,
      region,
      region_col
    ) {
      region_shape <- regions[
        as.character(regions[[region_col]]) == region,
        ,
        drop = FALSE
      ]

      if (nrow(region_shape) == 0) {
        stop(
          "Region not found: ",
          region,
          call. = FALSE
        )
      }

      terra::extract(
        indicator_vals,
        terra::vect(region_shape)
      )
    }

    # Unique region names
    region_names <- unique(
      as.character(categorisation_shape[[region_col]])
    )

    region_names <- region_names[
      !is.na(region_names)
    ]

    regions_for_study <- c(
      "National",
      region_names
    )

    # ------------------------------------------------------------------------
    # Bootstrap each region
    # ------------------------------------------------------------------------

    representative_indicator_vals <- lapply(
      regions_for_study,
      function(region) {
        cat(
          "\n----------------------------------------\n"
        )

        cat(
          "Running bootstrap model for region: ",
          region,
          "\n",
          sep = ""
        )

        # --------------------------------------------------------------
        # Extract raster values
        # --------------------------------------------------------------

        if (region == "National") {
          selected_region_vals <- terra::extract(
            indicator_vals,
            terra::vect(categorisation_shape)
          )
        } else {
          selected_region_vals <- extract_region_values(
            indicator_vals = indicator_vals,
            regions = categorisation_shape,
            region = region,
            region_col = region_col
          )
        }

        # Remove NA values
        selected_region_vals <- selected_region_vals |>
          dplyr::filter(
            is.finite(scaled_indicator)
          )

        # Check sample size
        if (nrow(selected_region_vals) == 0) {
          warning(
            "No valid scaled indicator values for region: ",
            region
          )

          return(
            data.frame(
              mean = NA_real_,
              median = NA_real_,
              lower = NA_real_,
              upper = NA_real_,
              sd = NA_real_,
              n = 0,
              region = region
            )
          )
        }

        cat(
          "Number of available cells: ",
          nrow(selected_region_vals),
          "\n",
          sep = ""
        )

        # --------------------------------------------------------------
        # Bootstrap beta regression
        # --------------------------------------------------------------

        summarised_vals <- bootstrap_beta_intercept(
          data = selected_region_vals,
          value_col = "scaled_indicator",
          n_boot = n.boot,
          seed = boot.seed,
          nSamps = n.Samps
        )

        # Add sample size and region
        data.frame(
          summarised_vals$summary
        ) |>
          dplyr::mutate(
            n = summarised_vals$n,
            region = region
          )
      }
    ) |>
      dplyr::bind_rows()
  }

  # --------------------------------------------------------------------------
  # 8. Create output object
  # --------------------------------------------------------------------------

  output <- list(
    results = results,

    expected_richness = expected_richness,

    indicator_vals = indicator_vals,

    reference_data = reference_data_filtered,

    reference_predictions = reference_predictions,

    representative_indicator_vals = representative_indicator_vals
  )

  # --------------------------------------------------------------------------
  # 9. Save outputs
  # --------------------------------------------------------------------------

  if (save_outputs) {
    cat("\nSaving outputs...\n")

    save_object <- function(
      x,
      filename_base,
      dir = output_dir
    ) {
      if (inherits(x, "sf")) {
        sfarrow::st_write_parquet(
          obj = x,
          dsn = file.path(
            dir,
            paste0(filename_base, ".parquet")
          ),
          compression = "snappy"
        )
      } else if (inherits(x, "SpatRaster")) {
        terra::writeRaster(
          x,
          filename = file.path(
            dir,
            paste0(filename_base, ".tif")
          ),
          overwrite = overwrite
        )
      } else if (inherits(x, "data.frame")) {
        utils::write.csv(
          x,
          file = file.path(
            dir,
            paste0(filename_base, ".csv")
          ),
          row.names = FALSE
        )
      } else {
        saveRDS(
          x,
          file = file.path(
            dir,
            paste0(filename_base, ".rds")
          )
        )
      }
    }

    save_object(
      expected_richness,
      "expected_richness"
    )

    save_object(
      indicator_vals,
      "indicator_vals"
    )

    save_object(
      reference_data_filtered,
      "reference_data"
    )

    if (!is.null(representative_indicator_vals)) {
      save_object(
        representative_indicator_vals,
        "representative_indicator_vals"
      )
    }

    # Save individual species results
    species_output_dir <- file.path(
      output_dir,
      "species_results"
    )

    # Save complete object
    saveRDS(
      output,
      file = file.path(
        output_dir,
        "all_results.rds"
      )
    )
  }

  # --------------------------------------------------------------------------
  # 10. Return
  # --------------------------------------------------------------------------

  output
}
