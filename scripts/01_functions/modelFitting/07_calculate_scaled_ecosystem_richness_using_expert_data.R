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
  n.Samps = 5000,
  method = "rds"
) {
  # ============================================================
  # 1. Validate inputs
  # ============================================================

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

  if (!method %in% c("raster", "sf", "rds")) {
    stop(
      "`method` must be one of: 'raster', 'sf', or 'rds'.",
      call. = FALSE
    )
  }

  if (length(input_dir) == 0) {
    stop(
      "`input_dir` contains no paths.",
      call. = FALSE
    )
  }

  missing_dirs <- input_dir[!dir.exists(input_dir)]

  if (length(missing_dirs) > 0) {
    stop(
      "The following input directories do not exist:\n",
      paste(missing_dirs, collapse = "\n"),
      call. = FALSE
    )
  }

  if (!region_col %in% names(categorisation_shape)) {
    stop(
      "`region_col` ('",
      region_col,
      "') is not present in `categorisation_shape`.",
      call. = FALSE
    )
  }

  required_reference_cols <- c(
    "Species",
    "Values"
  )

  missing_reference_cols <- setdiff(
    required_reference_cols,
    names(reference_data)
  )

  if (length(missing_reference_cols) > 0) {
    stop(
      "Missing columns in `reference_data`: ",
      paste(
        missing_reference_cols,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  if (save_outputs && is.null(output_dir)) {
    stop(
      "`output_dir` must be supplied when ",
      "`save_outputs = TRUE`.",
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

  # ============================================================
  # 2. Calculate multispecies expected richness
  # ============================================================

  cat(
    "\nCalculating multispecies expected richness...\n"
  )

  multispecies_richness <-
    calculate_ecosystem_richness(
      input_dir = input_dir,
      ecosystem_type = ecosystem_type,
      value_col = value_col,
      na.rm = TRUE,
      output_dir = if (save_outputs) {
        file.path(
          output_dir,
          "expected_richness"
        )
      } else {
        NULL
      },
      save_outputs = save_outputs,
      overwrite = overwrite,
      method = method
    )

  results <- multispecies_richness$results

  if (length(results) == 0) {
    stop(
      "No valid species prediction outputs were found.",
      call. = FALSE
    )
  }

  # ============================================================
  # 3. Match expert reference species
  # ============================================================

  cat(
    "\nCalculating expert reference value...\n"
  )

  model_species <- names(results)

  reference_data_filtered <- reference_data |>
    dplyr::mutate(
      Species_col = safe_name(
        as.character(Species)
      )
    ) |>
    dplyr::filter(
      Species_col %in% model_species
    )

  if (nrow(reference_data_filtered) == 0) {
    stop(
      "None of the species in `reference_data` ",
      "matched the modelled species.",
      call. = FALSE
    )
  }

  # Check values
  if (
    any(
      !is.finite(
        reference_data_filtered$Values
      )
    )
  ) {
    stop(
      "Expert reference values contain NA, NaN, ",
      "Inf, or -Inf.",
      call. = FALSE
    )
  }

  # Optional diagnostic for duplicate species
  duplicate_species <- reference_data_filtered |>
    dplyr::count(Species_col) |>
    dplyr::filter(n > 1)

  if (nrow(duplicate_species) > 0) {
    warning(
      "Multiple expert reference values were found for: ",
      paste(
        duplicate_species$Species_col,
        collapse = ", "
      ),
      ". All matching values will be included in the sum.",
      call. = FALSE
    )
  }

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
    "Expert reference richness: ",
    round(reference_predictions, 4),
    "\n",
    sep = ""
  )

  # ============================================================
  # 4. Extract expected richness
  # ============================================================

  expected_richness <-
    multispecies_richness$expected_richness

  if (!value_col %in% names(expected_richness)) {
    stop(
      "Layer/column `",
      value_col,
      "` was not found in expected richness.",
      call. = FALSE
    )
  }

  # ============================================================
  # 5. Scale expected richness
  # ============================================================

  cat(
    "\nScaling richness against expert reference...\n"
  )

  if (inherits(expected_richness, "SpatRaster")) {
    scaled_richness <-
      expected_richness[[value_col]] /
      reference_predictions

    scaled_richness <- terra::clamp(
      scaled_richness,
      lower = 0,
      upper = 1,
      values = TRUE
    )

    indicator_vals <- c(
      expected_richness[[value_col]],
      scaled_richness
    )

    names(indicator_vals) <- c(
      "indicator_vals",
      "scaled_indicator"
    )
  } else if (inherits(expected_richness, "sf")) {
    expected_richness$scaled_indicator <-
      expected_richness[[value_col]] /
      reference_predictions

    expected_richness$scaled_indicator <-
      pmin(
        pmax(
          expected_richness$scaled_indicator,
          0
        ),
        1
      )

    indicator_vals <- expected_richness[
      c(
        value_col,
        "scaled_indicator",
        attr(
          expected_richness,
          "sf_column"
        )
      )
    ]

    names(indicator_vals)[
      names(indicator_vals) == value_col
    ] <- "indicator_vals"
  } else {
    stop(
      "`expected_richness` must be either a ",
      "SpatRaster or sf object.",
      call. = FALSE
    )
  }

  # ============================================================
  # 6. Assign biological/nature regions
  # ============================================================

  region_names <- unique(
    as.character(
      categorisation_shape[[region_col]]
    )
  )

  region_names <- region_names[
    !is.na(region_names)
  ]

  if (inherits(indicator_vals, "SpatRaster")) {
    cat(
      "\nRasterising regions...\n"
    )

    regions_vect <- terra::vect(
      categorisation_shape
    )

    regions_vect$region_id <-
      seq_len(nrow(regions_vect))

    region_raster <- terra::rasterize(
      regions_vect,
      indicator_vals[[1]],
      field = "region_id"
    )

    region_levels <- data.frame(
      ID = seq_len(
        nrow(regions_vect)
      ),
      bio_region = as.character(
        regions_vect[[region_col]]
      )
    )

    levels(region_raster) <-
      region_levels

    names(region_raster) <-
      "bio_region"

    indicator_vals <- c(
      indicator_vals,
      region_raster
    )
  } else {
    cat(
      "\nAssigning regions to prediction points...\n"
    )

    regions_sf <- categorisation_shape[
      c(
        region_col,
        attr(
          categorisation_shape,
          "sf_column"
        )
      )
    ]

    names(regions_sf)[
      names(regions_sf) == region_col
    ] <- "bio_region"

    indicator_vals <- sf::st_join(
      indicator_vals,
      regions_sf,
      join = sf::st_within,
      left = TRUE
    )
  }

  # ============================================================
  # 7. Estimate representative regional values
  # ============================================================

  representative_indicator_vals <- NULL

  if (rep_values_for_each_region) {
    cat(
      "\nEstimating representative values for each region...\n"
    )

    regions_for_study <- c(
      "National",
      region_names
    )

    # ----------------------------------------------------------
    # Helper for one region
    # ----------------------------------------------------------

    bootstrap_region <- function(
      values,
      region
    ) {
      values <- values[
        is.finite(values)
      ]

      if (length(values) == 0) {
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
        "  ",
        region,
        ": ",
        length(values),
        " valid values\n",
        sep = ""
      )

      bootstrap_data <- data.frame(
        scaled_indicator = values
      )

      summarised_vals <-
        bootstrap_beta_intercept(
          data = bootstrap_data,
          value_col = "scaled_indicator",
          n_boot = n.boot,
          seed = boot.seed,
          nSamps = n.Samps
        )

      data.frame(
        summarised_vals$summary
      ) |>
        dplyr::mutate(
          n = summarised_vals$n,
          region = region
        )
    }

    # ==========================================================
    # Raster
    # ==========================================================

    if (inherits(indicator_vals, "SpatRaster")) {
      representative_indicator_vals <-
        lapply(
          regions_for_study,
          function(region) {
            cat(
              "\n----------------------------------------\n",
              "Region: ",
              region,
              "\n",
              sep = ""
            )

            if (region == "National") {
              selected <- terra::values(
                indicator_vals[["scaled_indicator"]],
                mat = FALSE
              )
            } else {
              region_shape <-
                categorisation_shape[
                  as.character(
                    categorisation_shape[[region_col]]
                  ) ==
                    region,
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

              selected <- terra::extract(
                indicator_vals[["scaled_indicator"]],
                terra::vect(region_shape)
              )[[2]]
            }

            bootstrap_region(
              selected,
              region
            )
          }
        ) |>
        dplyr::bind_rows()

      # ==========================================================
      # sf
      # ==========================================================
    } else {
      # --------------------------------------------------------
      # IMPORTANT:
      # Do not repeatedly subset the full sf object.
      #
      # Extract only the two columns needed for bootstrap.
      # This avoids copying geometry for every region.
      # --------------------------------------------------------

      regional_values <- indicator_vals |>
        sf::st_drop_geometry() |>
        dplyr::select(
          scaled_indicator,
          bio_region
        )

      # Remove invalid values once
      regional_values <- regional_values[
        is.finite(
          regional_values$scaled_indicator
        ),
        ,
        drop = FALSE
      ]

      # Split once by region
      regional_values_split <- split(
        regional_values$scaled_indicator,
        regional_values$bio_region
      )

      representative_indicator_vals <-
        lapply(
          regions_for_study,
          function(region) {
            cat(
              "\n----------------------------------------\n",
              "Region: ",
              region,
              "\n",
              sep = ""
            )

            if (region == "National") {
              values <-
                regional_values$scaled_indicator
            } else {
              values <-
                regional_values_split[[region]]

              if (is.null(values)) {
                values <- numeric(0)
              }
            }

            bootstrap_region(
              values,
              region
            )
          }
        ) |>
        dplyr::bind_rows()
    }
  }

  # ============================================================
  # 8. Create output
  # ============================================================

  output <- list(
    results = results,

    expected_richness = expected_richness,

    indicator_vals = indicator_vals,

    reference_data = reference_data_filtered,

    reference_predictions = reference_predictions,

    representative_indicator_vals = representative_indicator_vals
  )

  # ============================================================
  # 9. Save outputs
  # ============================================================

  if (save_outputs) {
    cat(
      "\nSaving outputs...\n"
    )

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
            paste0(
              filename_base,
              ".parquet"
            )
          ),
          compression = "snappy"
        )
      } else if (inherits(x, "SpatRaster")) {
        terra::writeRaster(
          x,
          filename = file.path(
            dir,
            paste0(
              filename_base,
              ".tif"
            )
          ),
          overwrite = overwrite
        )
      } else if (inherits(x, "data.frame")) {
        utils::write.csv(
          x,
          file = file.path(
            dir,
            paste0(
              filename_base,
              ".csv"
            )
          ),
          row.names = FALSE
        )
      } else {
        saveRDS(
          x,
          file = file.path(
            dir,
            paste0(
              filename_base,
              ".rds"
            )
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

    if (
      !is.null(
        representative_indicator_vals
      )
    ) {
      save_object(
        representative_indicator_vals,
        "representative_indicator_vals"
      )
    }

    saveRDS(
      output,
      file = file.path(
        output_dir,
        "all_results.rds"
      )
    )
  }

  # ============================================================
  # 10. Return
  # ============================================================

  output
}
