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
#' Generate block-wise occupancy and/or detection predictions
#'
#' Requires the existing package helpers `make_rhs()`, `safe_name()`, and
#' `.resolve_model_meta()`. The metadata helper must return:
#' `single_species`, `species_names`, `x_names` (occupancy design columns) and
#' `x_p_names` (detection design columns, needed for detection predictions).
#'
#' @param prediction_type One of "occupancy", "detection", or "both".
#' @param det_covariates Detection-model covariates (optional; if the fitted
#'   detection columns are year-specific, e.g. `mean_temperature2020`, the
#'   design matrix is built directly from the fitted column names).
#' @param occ_covariates Occupancy-model covariates.
#' @param return Either "sf" or "raster".
#' @param n_omp_threads Threads passed to spatial `predict()` methods.
#' @importFrom rlang .data
#' @export
generate_occupancy_predictions_for_each_species_with_blocks <- function(
  input_dir,
  prediction_points,
  det_covariates = NULL,
  occ_covariates = NULL,
  base_raster = NULL,
  model_crs,
  ecosystem_type,
  output_dir,
  prediction_type = c("occupancy", "detection", "both"),
  return = c("sf", "raster"),
  probs = c(0.025, 0.975),
  use_route_random_effect = FALSE,
  ignore_RE = TRUE,
  save_plots = TRUE,
  verbose = FALSE,
  plot_width = 8,
  plot_height = 6,
  plot_dpi = 300,
  block_size = 5000,
  n_omp_threads = 1L
) {
  prediction_type <- match.arg(prediction_type)
  return <- match.arg(return)

  need_occ <- prediction_type %in% c("occupancy", "both")
  need_det <- prediction_type %in% c("detection", "both")

  # ---- Argument checks ------------------------------------------------------
  if (!dir.exists(input_dir)) {
    stop("`input_dir` does not exist: ", input_dir, call. = FALSE)
  }
  if (return == "raster" && is.null(base_raster)) {
    stop("Supply `base_raster` when `return = 'raster'`.", call. = FALSE)
  }
  if (need_occ && is.null(occ_covariates)) {
    stop(
      "`occ_covariates` is required for occupancy predictions.",
      call. = FALSE
    )
  }
  if (
    length(block_size) != 1L ||
      !is.numeric(block_size) ||
      is.na(block_size) ||
      block_size < 1 ||
      block_size != as.integer(block_size)
  ) {
    stop("`block_size` must be a positive integer.", call. = FALSE)
  }
  if (
    length(probs) != 2L ||
      anyNA(probs) ||
      any(probs < 0 | probs > 1) ||
      probs[1] >= probs[2]
  ) {
    stop(
      "`probs` must contain two increasing probabilities in [0, 1].",
      call. = FALSE
    )
  }

  block_size <- as.integer(block_size)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

  rds_files <- list.files(
    input_dir,
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (!length(rds_files)) {
    stop("No .rds files found in `input_dir`.", call. = FALSE)
  }

  prediction_points <- as.data.frame(prediction_points)
  if (!all(c("long", "lat") %in% names(prediction_points))) {
    stop(
      "`prediction_points` must contain coordinate columns `long` and `lat`.",
      call. = FALSE
    )
  }

  if (use_route_random_effect) {
    prediction_points$route_id <- 1
  }

  # ---- Aggregate yearly covariates (used by occupancy formulas) -------------
  temporal_groups <- list(
    mean_temperature = "^mean_temperature_[0-9]{4}$",
    mean_precipitation = "^mean_precipitation_[0-9]{4}$",
    snowmelt = "^snowmelt_[0-9]{4}$"
  )
  for (nm in names(temporal_groups)) {
    cols <- grep(temporal_groups[[nm]], names(prediction_points), value = TRUE)
    if (length(cols)) {
      prediction_points[[nm]] <- rowMeans(
        prediction_points[, cols, drop = FALSE],
        na.rm = TRUE
      )
    }
  }

  occ_formula <- if (need_occ) {
    stats::as.formula(
      paste("~", make_rhs(occ_covariates, use_route_random_effect))
    )
  }
  det_formula <- if (need_det && !is.null(det_covariates)) {
    stats::as.formula(
      paste("~", make_rhs(det_covariates, use_route_random_effect))
    )
  }

  # Only drop rows missing variables needed by the requested prediction(s).
  required_vars <- unique(c(
    if (need_occ) all.vars(occ_formula),
    if (need_det && !is.null(det_formula)) all.vars(det_formula),
    "long",
    "lat"
  ))
  missing_vars <- setdiff(required_vars, names(prediction_points))
  if (length(missing_vars)) {
    stop(
      "Prediction data are missing required variables: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  keep <- stats::complete.cases(prediction_points[,
    required_vars,
    drop = FALSE
  ])
  if (any(!keep)) {
    warning(
      sum(!keep),
      " prediction rows removed because required variables are missing."
    )
    prediction_points <- prediction_points[keep, , drop = FALSE]
  }
  if (!nrow(prediction_points)) {
    stop("No complete prediction locations remain.", call. = FALSE)
  }

  n_prediction <- nrow(prediction_points)
  block_starts <- seq.int(1L, n_prediction, by = block_size)
  n_blocks <- length(block_starts)
  message(
    "Prediction locations: ",
    format(n_prediction, big.mark = ","),
    "; blocks: ",
    n_blocks,
    "; block size: ",
    block_size
  )

  requested_types <- if (prediction_type == "both") {
    c("occupancy", "detection")
  } else {
    prediction_type
  }

  if (return == "raster" && !inherits(base_raster, "SpatRaster")) {
    base_raster <- terra::rast(base_raster)
  }

  all_results <- list() # sf results (return = "sf") or rasters (return = "raster")

  # ---- Loop over fitted models ----------------------------------------------
  for (rds_file in rds_files) {
    message(
      "\n",
      strrep("=", 72),
      "\nReading: ",
      rds_file,
      "\n",
      strrep("=", 72)
    )

    output <- tryCatch(readRDS(rds_file), error = function(e) {
      stop(
        "Failed to read ",
        rds_file,
        ": ",
        conditionMessage(e),
        call. = FALSE
      )
    })

    if (is.null(output$model_products$model)) {
      warning("Skipping file without `output$model_products$model`: ", rds_file)
      next
    }
    model <- output$model_products$model

    meta <- tryCatch(
      .resolve_model_meta(output, rds_file),
      error = function(e) {
        warning(
          "Could not resolve model metadata for ",
          basename(rds_file),
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
    if (is.null(meta)) {
      next
    }

    is_single <- isTRUE(meta$single_species)
    species_names <- meta$species_names
    occ_x_names <- meta$x_names
    det_x_names <- meta$x_p_names

    if (need_det && (is.null(det_x_names) || !length(det_x_names))) {
      warning(
        "No detection design column names for ",
        basename(rds_file),
        "; skipping."
      )
      next
    }

    is_spatial <- if (!is.null(output$fit_spatial_model)) {
      isTRUE(output$fit_spatial_model)
    } else {
      inherits(model, c("spPGOcc", "spMsPGOcc"))
    }

    message(
      "Model: ",
      if (is_single) "single-species" else "multi-species",
      " | ",
      if (is_spatial) "spatial" else "non-spatial",
      " | ",
      length(species_names),
      " species"
    )

    # Plain data.frames per species; sf geometry is built once at the end.
    species_blocks <- stats::setNames(
      vector("list", length(species_names)),
      species_names
    )

    for (b in seq_along(block_starts)) {
      block_start <- block_starts[b]
      block_end <- min(block_start + block_size - 1L, n_prediction)
      block_idx <- block_start:block_end
      prediction_block <- prediction_points[block_idx, , drop = FALSE]
      message(
        "  Block ",
        b,
        "/",
        n_blocks,
        " (rows ",
        block_start,
        "-",
        block_end,
        ")"
      )

      coords_df <- prediction_block[, c("long", "lat"), drop = FALSE]
      coords_mat <- as.matrix(coords_df)

      block_summary <- NULL # named list: species -> data.frame
      block_failed <- FALSE

      # One type at a time so both posterior objects are never in memory together.
      for (ptype in requested_types) {
        if (ptype == "occupancy") {
          frm <- occ_formula
          expected <- occ_x_names
          sample_field <- "psi.0.samples"
        } else {
          expected <- det_x_names
          # Year-specific fitted columns: build directly from fitted names.
          frm <- if (any(grepl(.year_pattern, expected, perl = TRUE))) {
            NULL
          } else {
            det_formula
          }
          sample_field <- "p.0.samples"
        }

        X <- tryCatch(
          .make_prediction_matrix(
            formula = frm,
            data = prediction_block,
            expected_names = expected,
            type = ptype,
            label = paste(ptype, "prediction")
          ),
          error = function(e) {
            stop(
              "Failed to build ",
              ptype,
              " matrix for ",
              basename(rds_file),
              ", block ",
              b,
              ": ",
              conditionMessage(e),
              call. = FALSE
            )
          }
        )

        pred <- tryCatch(
          .run_predict(
            model = model,
            X_0 = X,
            coords_0_mat = coords_mat,
            is_single = is_single,
            is_spatial = is_spatial,
            ignore_RE = ignore_RE,
            verbose = verbose,
            type = ptype,
            n_omp_threads = n_omp_threads
          ),
          error = function(e) {
            warning(
              "Prediction failed (",
              ptype,
              ") for ",
              basename(rds_file),
              ", block ",
              b,
              ": ",
              conditionMessage(e)
            )
            NULL
          }
        )
        draws <- if (!is.null(pred)) pred[[sample_field]]
        rm(X, pred)

        if (is.null(draws)) {
          warning(
            "No `$",
            sample_field,
            "` for ",
            basename(rds_file),
            ", block ",
            b,
            "; block skipped."
          )
          block_failed <- TRUE
          break
        }

        smry <- .summarise_draws(
          draws = draws,
          species_names = species_names,
          coords = coords_df,
          location_ids = block_idx,
          probs = probs,
          is_single = is_single,
          prefix = if (prediction_type == "both") paste0(ptype, "_") else ""
        )
        rm(draws)

        if (is.null(block_summary)) {
          block_summary <- smry
        } else {
          # Same species and same row order, so columns can be bound directly.
          id_cols <- c("species", "location_id", "x", "y")
          for (sp in names(smry)) {
            block_summary[[sp]] <- cbind(
              block_summary[[sp]],
              smry[[sp]][, setdiff(names(smry[[sp]]), id_cols), drop = FALSE]
            )
          }
        }
        rm(smry)
      }

      # A block is kept only if every requested type succeeded, so columns
      # are consistent across blocks and rbind() cannot fail later.
      if (!block_failed && !is.null(block_summary)) {
        for (sp in names(block_summary)) {
          species_blocks[[sp]][[length(species_blocks[[sp]]) + 1L]] <-
            block_summary[[sp]]
        }
      }
      rm(block_summary, prediction_block, coords_df, coords_mat)
      gc(verbose = FALSE)
    }

    # ---- Per-species output -------------------------------------------------
    for (sp in species_names) {
      if (!length(species_blocks[[sp]])) {
        warning("No successful prediction blocks for species: ", sp)
        next
      }

      sp_df <- do.call(rbind, species_blocks[[sp]])
      species_blocks[[sp]] <- NULL
      sp_df <- sp_df[order(sp_df$location_id), , drop = FALSE]
      sp_result <- sf::st_as_sf(
        sp_df,
        coords = c("x", "y"),
        crs = model_crs,
        remove = FALSE
      )
      rm(sp_df)

      sp_safe <- safe_name(sp)
      species_dir <- file.path(output_dir, sp_safe)
      dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)

      ok <- dir.create(species_dir, recursive = TRUE, showWarnings = FALSE) ||
        dir.exists(species_dir)
      if (!ok) {
        stop("Could not create directory: ", species_dir, call. = FALSE)
      }

      parquet_file <- file.path(
        species_dir,
        paste0(sp_safe, "_", type_suffix, ".parquet")
      )
      if (.Platform$OS.type == "windows" && nchar(parquet_file) >= 250) {
        stop(
          "Output path is ",
          nchar(parquet_file),
          " characters, above the Windows 260 limit. Use a shorter `output_dir`:\n",
          parquet_file,
          call. = FALSE
        )
      }

      type_suffix <- paste0(prediction_type, "_predictions")
      sfarrow::st_write_parquet(
        sp_result,
        file.path(species_dir, paste0(type_suffix, ".parquet")),
        compression = "snappy"
      )

      saveRDS(
        sp_result,
        file.path(
          species_dir,
          paste0(sp_safe, "_", return, "_", type_suffix, ".rds")
        )
      )

      if (save_plots) {
        plot_dir <- file.path(species_dir, "plots")
        dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
        for (ptype in requested_types) {
          value_col <- if (prediction_type == "both") {
            paste0(ptype, "_mean")
          } else {
            "mean"
          }
          p <- ggplot2::ggplot(sp_result) +
            ggplot2::geom_sf(
              ggplot2::aes(colour = .data[[value_col]]),
              size = 0.15,
              alpha = 0.8
            ) +
            ggplot2::scale_colour_viridis_c(
              name = paste0(tools::toTitleCase(ptype), "\nprobability"),
              limits = c(0, 1),
              na.value = "grey80"
            ) +
            ggplot2::coord_sf() +
            ggplot2::theme_bw() +
            ggplot2::labs(
              title = sp,
              subtitle = paste(
                ecosystem_type,
                if (is_single) "single-species" else "multi-species",
                if (is_spatial) "spatial" else "non-spatial",
                sep = " | "
              )
            )
          plot_base <- file.path(
            plot_dir,
            paste0(sp_safe, "_", ptype, "_probability")
          )
          ggplot2::ggsave(
            paste0(plot_base, ".png"),
            p,
            width = plot_width,
            height = plot_height,
            dpi = plot_dpi
          )
          ggplot2::ggsave(
            paste0(plot_base, ".pdf"),
            p,
            width = plot_width,
            height = plot_height
          )
          rm(p)
        }
      }

      key <- paste(basename(rds_file), sp, sep = "::")
      if (return == "sf") {
        all_results[[key]] <- sp_result
      } else {
        # Rasterise immediately so full sf results are not held in memory.
        all_results[[key]] <- .write_species_raster(
          sp_result,
          base_raster,
          prediction_type,
          file.path(species_dir, paste0(sp_safe, "_", type_suffix, ".tif"))
        )
      }
      rm(sp_result)
      gc(verbose = FALSE)
    }
    rm(output, model, species_blocks)
    gc(verbose = FALSE)
  }

  invisible(all_results)
}
