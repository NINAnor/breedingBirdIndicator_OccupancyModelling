#' Generate Block-wise Occupancy Predictions for Each Species
#'
#' Runs occupancy predictions in blocks to avoid holding the complete
#' posterior prediction object in memory. Each block is predicted and
#' summarized immediately. The summarized block results are then merged
#' and saved as one final prediction file per species.
#'
#' @param input_dir Character. Directory containing fitted model RDS files.
#' @param prediction_points Data frame containing prediction locations and
#'   occupancy covariates.
#' @param det_covariates Detection covariates used by the fitted models.
#' @param occ_covariates Occupancy covariates used by the fitted models.
#' @param base_raster Character or terra::SpatRaster raster template.
#' @param model_crs Character or numeric CRS used by the model.
#' @param ecosystem_type Character ecosystem identifier for plot labels.
#' @param output_dir Character output directory.
#' @param return Character. Either `"sf"` or `"raster"`.
#' @param probs Numeric vector of lower and upper posterior probabilities.
#' @param use_route_random_effect Logical.
#' @param ignore_RE Logical. Passed to `.run_predict()`.
#' @param save_plots Logical.
#' @param verbose Logical.
#' @param plot_width Numeric.
#' @param plot_height Numeric.
#' @param plot_dpi Numeric.
#' @param block_size Integer. Number of prediction locations per block.
#'
#' @return Invisibly returns a named list containing the merged
#'   species-level predictions.
#'
#' @export
generate_occupancy_predictions_for_each_species_with_blocks <- function(
  input_dir,
  prediction_points,
  det_covariates,
  occ_covariates,
  base_raster = NULL,
  model_crs,
  ecosystem_type,
  output_dir,
  return = c("sf", "raster"),
  probs = c(0.025, 0.975),
  use_route_random_effect = FALSE,
  ignore_RE = TRUE,
  save_plots = TRUE,
  verbose = FALSE,
  plot_width = 8,
  plot_height = 6,
  plot_dpi = 300,
  block_size = 5000
) {
  # ========================================================================
  # 1. Validate inputs
  # ========================================================================

  return <- match.arg(return)

  if (!dir.exists(input_dir)) {
    stop(
      "`input_dir` does not exist or is not a directory: ",
      input_dir,
      call. = FALSE
    )
  }

  if (return == "raster" && is.null(base_raster)) {
    stop(
      "`base_raster` must be supplied when ",
      "`return = 'raster'`.",
      call. = FALSE
    )
  }

  if (
    length(block_size) != 1L ||
      !is.numeric(block_size) ||
      is.na(block_size) ||
      block_size < 1
  ) {
    stop(
      "`block_size` must be a positive integer.",
      call. = FALSE
    )
  }

  block_size <- as.integer(block_size)

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  output_dir <- normalizePath(
    output_dir,
    winslash = "/",
    mustWork = FALSE
  )

  # ========================================================================
  # 2. Find fitted model files
  # ========================================================================

  rds_files <- list.files(
    input_dir,
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(rds_files) == 0) {
    stop(
      "No .rds files found in `input_dir`.",
      call. = FALSE
    )
  }

  # ========================================================================
  # 3. Prepare prediction data
  # ========================================================================

  message(
    "Preparing prediction data."
  )

  prediction_points <- prediction_points %>%
    dplyr::filter(
      stats::complete.cases(.)
    )

  # ------------------------------------------------------------------------
  # Route random effect
  # ------------------------------------------------------------------------

  if (use_route_random_effect) {
    prediction_points <- prediction_points %>%
      dplyr::mutate(
        route_id = 1
      )
  }

  # ------------------------------------------------------------------------
  # Temporal means
  # ------------------------------------------------------------------------

  prediction_points <- prediction_points %>%
    dplyr::mutate(
      mean_temperature = rowMeans(
        dplyr::across(
          dplyr::matches(
            "^mean_temperature_[0-9]{4}$"
          )
        ),
        na.rm = TRUE
      ),

      mean_precipitation = rowMeans(
        dplyr::across(
          dplyr::matches(
            "^mean_precipitation_[0-9]{4}$"
          )
        ),
        na.rm = TRUE
      ),

      snowmelt = rowMeans(
        dplyr::across(
          dplyr::matches(
            "^snowmelt_[0-9]{4}$"
          )
        ),
        na.rm = TRUE
      )
    )

  # ========================================================================
  # 4. Occupancy formula
  # ========================================================================

  occ_rhs <- make_rhs(
    occ_covariates,
    use_route_random_effect
  )

  occ_formula <- stats::as.formula(
    paste(
      "~",
      occ_rhs
    )
  )

  # ========================================================================
  # 5. Prediction blocks
  # ========================================================================

  n_prediction <- nrow(
    prediction_points
  )

  block_starts <- seq(
    1L,
    n_prediction,
    by = block_size
  )

  n_blocks <- length(
    block_starts
  )

  message(
    "\nPrediction locations: ",
    format(
      n_prediction,
      big.mark = ","
    )
  )

  message(
    "Block size: ",
    format(
      block_size,
      big.mark = ","
    )
  )

  message(
    "Number of prediction blocks: ",
    n_blocks
  )

  # ========================================================================
  # 6. Final results
  # ========================================================================
  #
  # Only summarized predictions are retained.
  #
  # pred itself is deleted after every block.
  #
  # ========================================================================

  all_results <- list()

  # =========================================================================
  # 7. Loop over models
  # =========================================================================

  for (rds_file in rds_files) {
    cat(
      "\n",
      paste(
        rep("=", 72),
        collapse = ""
      ),
      "\n",
      sep = ""
    )

    cat(
      "Reading: ",
      rds_file,
      "\n",
      sep = ""
    )

    cat(
      paste(
        rep("=", 72),
        collapse = ""
      ),
      "\n",
      sep = ""
    )

    # ----------------------------------------------------------------------
    # Read model
    # ----------------------------------------------------------------------

    output <- tryCatch(
      readRDS(
        rds_file
      ),

      error = function(e) {
        stop(
          "Failed to read:\n",
          rds_file,
          "\n\n",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    if (is.null(output$model_products$model)) {
      warning(
        "Skipping file without ",
        "`output$model_products$model`: ",
        rds_file
      )

      rm(output)

      next
    }

    model <- output$model_products$model

    # ----------------------------------------------------------------------
    # Model metadata
    # ----------------------------------------------------------------------

    meta <- tryCatch(
      .resolve_model_meta(
        output,
        rds_file
      ),

      error = function(e) {
        warning(
          conditionMessage(e),
          "\nSkipping: ",
          rds_file
        )

        NULL
      }
    )

    if (is.null(meta)) {
      rm(
        output,
        model
      )

      gc(
        verbose = FALSE
      )

      next
    }

    is_single <- meta$single_species

    sp_names <- meta$species_names

    x_names <- meta$x_names

    # ----------------------------------------------------------------------
    # Spatial model
    # ----------------------------------------------------------------------

    is_spatial <- if (
      !is.null(
        output$fit_spatial_model
      )
    ) {
      isTRUE(
        output$fit_spatial_model
      )
    } else {
      inherits(
        model,
        c(
          "spPGOcc",
          "spMsPGOcc"
        )
      )
    }

    cat(
      sprintf(
        "Model: %s | %s | %d species\n",
        if (is_single) {
          "single-species"
        } else {
          "multi-species"
        },
        if (is_spatial) {
          "spatial"
        } else {
          "non-spatial"
        },
        length(sp_names)
      )
    )

    cat(
      "Species: ",
      paste(
        sp_names,
        collapse = ", "
      ),
      "\n",
      sep = ""
    )

    # ======================================================================
    # 8. Check raw occupancy variables
    # ======================================================================

    required_vars <- all.vars(
      occ_formula
    )

    missing_vars <- setdiff(
      required_vars,
      colnames(prediction_points)
    )

    if (length(missing_vars) > 0) {
      stop(
        "Prediction data are missing variables required by the ",
        "occupancy formula for ",
        basename(rds_file),
        ":\n",
        paste(
          missing_vars,
          collapse = ", "
        ),
        call. = FALSE
      )
    }

    # ======================================================================
    # 9. Preserve factor structure
    # ======================================================================

    factor_columns <- names(
      prediction_points
    )[
      vapply(
        prediction_points,
        is.factor,
        logical(1)
      )
    ]

    if (length(factor_columns) > 0) {
      for (v in factor_columns) {
        prediction_points[[v]] <- factor(
          prediction_points[[v]],
          levels = levels(
            prediction_points[[v]]
          )
        )
      }
    }

    # ======================================================================
    # 10. Block summary accumulator
    # ======================================================================

    species_blocks <- setNames(
      vector(
        "list",
        length(sp_names)
      ),
      sp_names
    )

    # =========================================================================
    # 11. Prediction blocks
    # =========================================================================

    for (b in seq_len(n_blocks)) {
      block_start <- block_starts[b]

      block_end <- min(
        block_start + block_size - 1L,
        n_prediction
      )

      block_idx <- block_start:block_end

      cat(
        sprintf(
          "\n  Block %d/%d: rows %s-%s (%s locations)\n",
          b,
          n_blocks,
          format(
            block_start,
            big.mark = ","
          ),
          format(
            block_end,
            big.mark = ","
          ),
          format(
            length(block_idx),
            big.mark = ","
          )
        )
      )

      # =====================================================================
      # 11.1 Extract block
      # =====================================================================

      prediction_block <- prediction_points[
        block_idx,
        ,
        drop = FALSE
      ]

      # =====================================================================
      # 11.2 Construct model matrix
      # =====================================================================

      X_0 <- tryCatch(
        stats::model.matrix(
          occ_formula,
          data = prediction_block
        ),

        error = function(e) {
          stop(
            "Could not construct model.matrix() for block ",
            b,
            " of ",
            basename(rds_file),
            ":\n",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )

      cat(
        "    Model matrix: ",
        nrow(X_0),
        " x ",
        ncol(X_0),
        "\n",
        sep = ""
      )

      # =====================================================================
      # 11.3 Add model columns absent from this block
      # =====================================================================
      #
      # This handles factor levels that do not occur in an individual
      # prediction block.
      #
      # If the full prediction data contain a level but the current block
      # does not, model.matrix() can omit that dummy column.
      #
      # The fitted model still expects the column, so add it as zero.
      #
      # =====================================================================

      missing_cols <- setdiff(
        x_names,
        colnames(X_0)
      )

      if (length(missing_cols) > 0) {
        cat(
          "    Adding missing model columns: ",
          paste(
            missing_cols,
            collapse = ", "
          ),
          "\n",
          sep = ""
        )

        for (col in missing_cols) {
          X_0[[col]] <- 0
        }
      }

      # =====================================================================
      # 11.4 Remove/reorder columns to exactly match fitted model
      # =====================================================================

      X_0 <- X_0[,
        x_names,
        drop = FALSE
      ]

      if (
        !identical(
          colnames(X_0),
          x_names
        )
      ) {
        stop(
          "Prediction model matrix does not match fitted model ",
          "after alignment in block ",
          b,
          ".",
          call. = FALSE
        )
      }

      # =====================================================================
      # 11.5 Coordinates
      # =====================================================================

      coords_0_mat <- as.matrix(
        prediction_block[,
          c(
            "long",
            "lat"
          ),
          drop = FALSE
        ]
      )

      coords_0_df <- prediction_block[,
        c(
          "long",
          "lat"
        ),
        drop = FALSE
      ]

      # =====================================================================
      # 11.6 Run prediction
      # =====================================================================

      cat(
        "    Running .run_predict()...\n"
      )

      pred <- tryCatch(
        .run_predict(
          model = model,
          X_0 = X_0,
          coords_0_mat = coords_0_mat,
          is_single = is_single,
          is_spatial = is_spatial,
          ignore_RE = ignore_RE,
          verbose = verbose
        ),

        error = function(e) {
          warning(
            "Prediction failed for ",
            basename(rds_file),
            ", block ",
            b,
            ":\n",
            conditionMessage(e)
          )

          NULL
        }
      )

      if (is.null(pred)) {
        rm(
          prediction_block,
          X_0,
          coords_0_mat,
          coords_0_df
        )

        gc(
          verbose = FALSE
        )

        next
      }

      # =====================================================================
      # 11.7 Temporary output directory
      # =====================================================================
      #
      # summarize_occupancy_predictions() requires a valid path because
      # it currently calls normalizePath(output_dir), even when the user
      # does not intend to save from that function.
      #
      # We therefore give each block its own temporary directory.
      #
      # =====================================================================

      model_safe <- safe_name(
        tools::file_path_sans_ext(
          basename(rds_file)
        )
      )

      block_output_dir <- file.path(
        output_dir,
        "_block_summaries",
        model_safe,
        paste0(
          "block_",
          b
        )
      )

      dir.create(
        block_output_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )

      # =====================================================================
      # 11.8 Summarize this block
      # =====================================================================

      cat(
        "    Summarizing posterior predictions...\n"
      )

      block_results <- tryCatch(
        summarize_occupancy_predictions(
          pred = pred,
          coords = coords_0_df,
          species_names = sp_names,
          template = NULL,
          return = "sf",
          crs = model_crs,
          probs = probs,
          output_dir = block_output_dir
        ),

        error = function(e) {
          warning(
            "Could not summarize block ",
            b,
            " for ",
            basename(rds_file),
            ":\n",
            conditionMessage(e)
          )

          NULL
        }
      )

      # =====================================================================
      # 11.9 Keep summarized block results
      # =====================================================================

      if (!is.null(block_results)) {
        for (sp in sp_names) {
          sp_block <- dplyr::filter(
            block_results,
            species == sp
          )

          if (nrow(sp_block) > 0) {
            species_blocks[[sp]][[
              length(
                species_blocks[[sp]]
              ) +
                1L
            ]] <- sp_block
          }
        }
      }

      # =====================================================================
      # 11.10 Remove temporary block output
      # =====================================================================

      if (dir.exists(block_output_dir)) {
        unlink(
          block_output_dir,
          recursive = TRUE,
          force = TRUE
        )
      }

      # =====================================================================
      # 11.11 Release large objects immediately
      # =====================================================================

      rm(
        pred,
        block_results,
        prediction_block,
        X_0,
        coords_0_mat,
        coords_0_df
      )

      gc(
        verbose = FALSE
      )
    }

    # =========================================================================
    # 12. MERGE BLOCK PREDICTIONS AND SAVE
    # =========================================================================

    cat(
      "\nMerging block predictions for model:\n",
      basename(rds_file),
      "\n"
    )

    for (sp in sp_names) {
      sp_safe <- safe_name(
        sp
      )

      # ----------------------------------------------------------------------
      # Check for successful blocks
      # ----------------------------------------------------------------------

      if (
        length(
          species_blocks[[sp]]
        ) ==
          0
      ) {
        warning(
          "No successful prediction blocks for species: ",
          sp
        )

        next
      }

      # ----------------------------------------------------------------------
      # Merge all summarized blocks
      # ----------------------------------------------------------------------

      cat(
        "  Merging ",
        length(
          species_blocks[[sp]]
        ),
        " blocks for ",
        sp,
        "\n",
        sep = ""
      )

      sp_result <- do.call(
        rbind,
        species_blocks[[sp]]
      )

      # ----------------------------------------------------------------------
      # Rebuild geometry after rbind if necessary
      # ----------------------------------------------------------------------

      if (
        !inherits(
          sp_result,
          "sf"
        )
      ) {
        sp_result <- sf::st_as_sf(
          sp_result,
          coords = c(
            "x",
            "y"
          ),
          crs = model_crs,
          remove = FALSE
        )
      }

      # ----------------------------------------------------------------------
      # Order by location
      # ----------------------------------------------------------------------

      sp_result <- sp_result[
        order(
          sp_result$location_id
        ),
        ,
        drop = FALSE
      ]

      # ----------------------------------------------------------------------
      # Final species directory
      # ----------------------------------------------------------------------

      species_output_dir <- file.path(
        output_dir,
        sp_safe
      )

      dir.create(
        species_output_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )

      # =====================================================================
      # 12.1 Save merged SF / GeoParquet
      # =====================================================================

      merged_parquet <- file.path(
        species_output_dir,
        paste0(
          sp_safe,
          "_occupancy_predictions.parquet"
        )
      )

      cat(
        "  Saving merged predictions:\n",
        "  ",
        merged_parquet,
        "\n",
        sep = ""
      )

      sfarrow::st_write_parquet(
        sp_result,
        merged_parquet,
        compression = "snappy"
      )

      # =====================================================================
      # 12.2 Save merged RDS
      # =====================================================================

      merged_rds <- file.path(
        species_output_dir,
        paste0(
          sp_safe,
          "_",
          return,
          "_occupancy_predictions.rds"
        )
      )

      saveRDS(
        sp_result,
        merged_rds
      )

      cat(
        "  Merged predictions: ",
        format(
          nrow(sp_result),
          big.mark = ","
        ),
        " rows\n",
        sep = ""
      )

      # =====================================================================
      # 12.3 Plot merged predictions
      # =====================================================================

      if (save_plots) {
        plot_dir <- file.path(
          species_output_dir,
          "plots"
        )

        dir.create(
          plot_dir,
          recursive = TRUE,
          showWarnings = FALSE
        )

        p <- ggplot2::ggplot(
          sp_result
        ) +

          ggplot2::geom_sf(
            ggplot2::aes(
              colour = mean
            ),
            size = 0.15,
            alpha = 0.8
          ) +

          ggplot2::scale_color_viridis_c(
            name = "Occupancy\nprobability",
            limits = c(
              0,
              1
            ),
            na.value = "grey80"
          ) +

          ggplot2::coord_sf() +

          ggplot2::theme_bw() +

          ggplot2::labs(
            title = sp,
            subtitle = paste0(
              ecosystem_type,
              " | ",
              if (is_single) {
                "single-species"
              } else {
                "multi-species"
              },
              " | ",
              if (is_spatial) {
                "spatial"
              } else {
                "non-spatial"
              }
            )
          )

        plot_base <- file.path(
          plot_dir,
          paste0(
            sp_safe,
            "_occupancy_probability"
          )
        )

        ggplot2::ggsave(
          filename = paste0(
            plot_base,
            ".png"
          ),
          plot = p,
          width = plot_width,
          height = plot_height,
          dpi = plot_dpi
        )

        ggplot2::ggsave(
          filename = paste0(
            plot_base,
            ".pdf"
          ),
          plot = p,
          width = plot_width,
          height = plot_height
        )
      }

      # ----------------------------------------------------------------------
      # Store final merged result
      # ----------------------------------------------------------------------

      all_results[[sp]] <- sp_result

      # ----------------------------------------------------------------------
      # Release species block objects
      # ----------------------------------------------------------------------

      species_blocks[[sp]] <- NULL

      rm(
        sp_result
      )

      gc(
        verbose = FALSE
      )
    }

    # =========================================================================
    # 13. Clean model
    # =========================================================================

    rm(
      output,
      model,
      species_blocks
    )

    gc(
      verbose = FALSE
    )
  }

  # =========================================================================
  # 14. Return merged SF results
  # =========================================================================

  if (return == "sf") {
    return(
      invisible(
        all_results
      )
    )
  }

  # =========================================================================
  # 15. Convert merged predictions to raster
  # =========================================================================

  if (
    !inherits(
      base_raster,
      "terra::SpatRaster"
    )
  ) {
    base_raster <- terra::rast(
      base_raster
    )
  }

  raster_results <- list()

  for (sp in names(all_results)) {
    sp_safe <- safe_name(
      sp
    )

    sp_sf <- all_results[[sp]]

    cat(
      "\nRasterizing merged predictions for: ",
      sp,
      "\n",
      sep = ""
    )

    # ----------------------------------------------------------------------
    # Rasterize summaries
    # ----------------------------------------------------------------------

    mean_raster <- terra::rasterize(
      terra::vect(sp_sf),
      base_raster,
      field = "mean"
    )

    median_raster <- terra::rasterize(
      terra::vect(sp_sf),
      base_raster,
      field = "median"
    )

    sd_raster <- terra::rasterize(
      terra::vect(sp_sf),
      base_raster[[1]],
      field = "sd"
    )

    lower_raster <- terra::rasterize(
      terra::vect(sp_sf),
      base_raster[[1]],
      field = "lower"
    )

    upper_raster <- terra::rasterize(
      terra::vect(sp_sf),
      base_raster[[1]],
      field = "upper"
    )

    names(mean_raster) <- "mean"
    names(median_raster) <- "median"
    names(sd_raster) <- "sd"
    names(lower_raster) <- "lower"
    names(upper_raster) <- "upper"

    sp_raster <- c(
      mean_raster,
      median_raster,
      sd_raster,
      lower_raster,
      upper_raster
    )

    # ----------------------------------------------------------------------
    # Save raster
    # ----------------------------------------------------------------------

    species_output_dir <- file.path(
      output_dir,
      sp_safe
    )

    raster_file <- file.path(
      species_output_dir,
      paste0(
        sp_safe,
        "_occupancy_predictions.tif"
      )
    )

    terra::writeRaster(
      sp_raster,
      raster_file,
      overwrite = TRUE,
      wopt = list(
        datatype = "FLT4S",
        gdal = c(
          "COMPRESS=DEFLATE",
          "PREDICTOR=2"
        )
      )
    )

    raster_results[[sp]] <- sp_raster

    rm(
      sp_sf,
      mean_raster,
      median_raster,
      sd_raster,
      lower_raster,
      upper_raster,
      sp_raster
    )

    gc(
      verbose = FALSE
    )
  }

  # =========================================================================
  # 16. Return raster results
  # =========================================================================

  invisible(
    raster_results
  )
}
