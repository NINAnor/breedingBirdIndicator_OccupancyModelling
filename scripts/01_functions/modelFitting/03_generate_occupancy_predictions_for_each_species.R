# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

make_rhs <- function(covariates, use_route_random_effect = FALSE) {
  if (length(covariates) == 0) {
    rhs <- "1"
  } else {
    rhs <- paste(covariates, collapse = " + ")
  }

  if (use_route_random_effect) {
    rhs <- paste(rhs, "+ (1 | route_id)")
  }

  rhs
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]+", "_", x)
  gsub("^_+|_+$", "", x)
}


# -------------------------------------------------------------------------
# Internal: resolve species identity and model type from a loaded RDS output
#
# Returns a list with:
#   $single_species  logical
#   $species_names   character vector (length 1 for single-species)
#   $x_names         character vector of occupancy predictor column names
# -------------------------------------------------------------------------

.resolve_model_meta <- function(output, rds_file) {

  model <- output$model_products$model

  if (is.null(model)) {
    stop(
      "No fitted model found in `output$model_products$model` for RDS: ",
      rds_file,
      call. = FALSE
    )
  }

  # ---- Determine single vs multi-species ----------------------------------

  is_single <- if (!is.null(output$single_species)) {
    isTRUE(output$single_species)
  } else {
    is.null(model$sp.names)
  }

  # ---- Species name(s) ----------------------------------------------------

  if (is_single) {

    sp_names <- output$species

    if (is.null(sp_names) || length(sp_names) != 1L || is.na(sp_names)) {

      base <- tools::file_path_sans_ext(basename(rds_file))

      sp_names <- sub(
        "^(spPGOcc|PGOcc)_[^_]+_",
        "",
        base
      )

      sp_names <- gsub("_", " ", sp_names)

      if (is.null(sp_names) || length(sp_names) != 1L || !nzchar(sp_names)) {
        stop(
          "Cannot determine species name from single-species RDS: ",
          rds_file,
          "\nExpected `output$species` to be a length-1 character vector, ",
          "or filename like `spPGOcc_Ecosystem_Genus_species.rds`.",
          call. = FALSE
        )
      }
    }

  } else {

    sp_names <- if (!is.null(model$sp.names)) {
      model$sp.names
    } else {
      output$species
    }

    if (is.null(sp_names) || length(sp_names) == 0L) {
      stop(
        "Cannot determine species names from multi-species RDS: ",
        rds_file,
        "\nExpected `model$sp.names` or `output$species` to be non-NULL.",
        call. = FALSE
      )
    }
  }

  # ---- Occupancy predictor matrix column names ----------------------------

  if (is.null(output$formulas$occ_formula)) {
    stop(
      "No occupancy formula found at `output$formulas$occ_formula` for RDS: ",
      rds_file,
      call. = FALSE
    )
  }

  if (is.null(output$data_list$occ.covs)) {
    stop(
      "No occupancy covariates found at `output$data_list$occ.covs` for RDS: ",
      rds_file,
      call. = FALSE
    )
  }

  occ_formula_chr <- paste(deparse(output$formulas$occ_formula), collapse = " ")

  occ_formula_chr <- gsub(
    "\\+\\s*\\(1\\s*\\|\\s*route_id\\)",
    "",
    occ_formula_chr
  )

  occ_formula_chr <- gsub(
    "\\(1\\s*\\|\\s*route_id\\)",
    "",
    occ_formula_chr
  )

  occ_formula_fixed <- stats::as.formula(occ_formula_chr)

  x_names <- colnames(
    stats::model.matrix(
      occ_formula_fixed,
      data = output$data_list$occ.covs
    )
  )

  if (is.null(x_names) || length(x_names) == 0L) {
    stop(
      "Cannot determine model-matrix predictor columns for RDS: ",
      rds_file,
      call. = FALSE
    )
  }

  list(
    single_species = is_single,
    species_names  = sp_names,
    x_names        = x_names,
    occ_formula_fixed = occ_formula_fixed
  )
}


# -------------------------------------------------------------------------
# Internal: call predict() with the right signature for the model class
#
# Key differences between single and multi-species predict methods:
#   - spMsPGOcc / msPGOcc : accept `ignore.RE`
#   - spPGOcc  / PGOcc    : do NOT accept `ignore.RE`
#   - spMsPGOcc / spPGOcc : require `coords.0` as a plain matrix (s.x, s.y)
#   - msPGOcc  / PGOcc    : do NOT accept `coords.0`
# -------------------------------------------------------------------------

.run_predict <- function(model, X_0, coords_0_mat, is_single, is_spatial,
                         ignore_RE = TRUE, verbose = TRUE) {
  if (is_single) {
    if (is_spatial) {
      # spPGOcc::predict — needs coords.0, no ignore.RE
      predict(
        model,
        X.0      = X_0,
        coords.0 = coords_0_mat,
        ignore.RE = ignore_RE,
        verbose  = verbose
      )
    } else {
      # PGOcc::predict — no coords.0, no ignore.RE
      predict(
        model,
        X.0     = X_0,
        ignore.RE = ignore_RE,
        verbose = verbose
      )
    }
  } else {
    if (is_spatial) {
      # spMsPGOcc::predict — needs coords.0 and accepts ignore.RE
      predict(
        model,
        X.0      = X_0,
        coords.0 = coords_0_mat,
        ignore.RE = ignore_RE,
        verbose  = verbose
      )
    } else {
      # msPGOcc::predict — no coords.0, accepts ignore.RE
      predict(
        model,
        X.0       = X_0,
        ignore.RE = ignore_RE,
        verbose   = verbose
      )
    }
  }
}


# -------------------------------------------------------------------------
# Main batch prediction function
# -------------------------------------------------------------------------

generate_occupancy_predictions_for_each_species <- function(
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
    ignore_RE = TRUE,          # passed to multi-species predict only
    save_plots = TRUE,
    verbose = FALSE,
    plot_width  = 8,
    plot_height = 6,
    plot_dpi    = 300
) {
  return <- match.arg(return)

  if (!dir.exists(input_dir)) {
    stop("`input_dir` does not exist or is not a directory: ", input_dir,
         call. = FALSE)
  }

  if (return == "raster" && is.null(base_raster)) {
    stop("`base_raster` must be supplied when `return = 'raster'`.",
         call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  rds_files <- list.files(
    input_dir,
    pattern   = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(rds_files) == 0) {
    stop("No .rds files found in `input_dir`.", call. = FALSE)
  }

  # Build the occupancy design matrix once — columns are aligned per-model
  # below using model$x.names, so this only needs to contain a superset of
  # all possible predictors.
  occ_rhs     <- make_rhs(occ_covariates, use_route_random_effect)
  occ_formula <- as.formula(paste("~", occ_rhs))

  if(use_route_random_effect) prediction_points <- prediction_points %>%
    dplyr::mutate(route_id = 1)

  X_0_full <- model.matrix(occ_formula, data = prediction_points)

  # Coordinates as a plain numeric matrix — used by spatial predict methods.
  # Column names must be present but the spatial predict generics accept any
  # two-column matrix; we name them to match the convention used during fitting.
  coords_0_mat <- as.matrix(
    prediction_points[, c("longitude", "latitude")]
  )

  # sf / terra coordinate data frame used by summarize_occupancy_predictions()
  coords_0_df <- prediction_points[, c("longitude", "latitude")]

  all_results <- list()

  # -------------------------------------------------------------------------
  # Loop over RDS files
  # -------------------------------------------------------------------------

  for (rds_file in rds_files) {

    cat("\nReading:", rds_file, "\n")

    output <- tryCatch(
      readRDS(rds_file),
      error = function(e) {
        stop("Failed to read:\n", rds_file, "\n\n", conditionMessage(e))
      }
    )

    if (is.null(output$model_products$model)) {
      warning("Skipping file without `output$model_products$model`: ", rds_file)
      next
    }

    model <- output$model_products$model

    # ------------------------------------------------------------------
    # Resolve metadata (single vs multi-species, species names, x.names)
    # ------------------------------------------------------------------

    meta <- tryCatch(
      .resolve_model_meta(output, rds_file),
      error = function(e) {
        warning(conditionMessage(e), "\nSkipping: ", rds_file)
        NULL
      }
    )

    if (is.null(meta)) next

    is_single  <- meta$single_species
    sp_names   <- meta$species_names
    x_names    <- meta$x_names

    # Detect whether this was a spatial model from the saved flag, or fall
    # back to inspecting the model class.
    is_spatial <- if (!is.null(output$fit_spatial_model)) {
      isTRUE(output$fit_spatial_model)
    } else {
      inherits(model, c("spPGOcc", "spMsPGOcc"))
    }

    # ------------------------------------------------------------------
    # Align prediction matrix columns to what the model was fitted with
    # ------------------------------------------------------------------

    missing_cols <- setdiff(x_names, colnames(X_0_full))

    if (length(missing_cols) > 0) {
      stop(
        "Prediction matrix is missing required model columns for ",
        basename(rds_file), ": ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }


    X_0 <- X_0_full[, x_names, drop = FALSE]


    # ------------------------------------------------------------------
    # Run prediction
    # ------------------------------------------------------------------

    cat(sprintf(
      "  Predicting [%s | %s | %d species]: %s\n",
      if (is_single) "single-species" else "multi-species",
      if (is_spatial) "spatial" else "non-spatial",
      length(sp_names),
      paste(sp_names, collapse = ", ")
    ))

    pred <- tryCatch(
      .run_predict(
        model        = model,
        X_0          = X_0,
        coords_0_mat = coords_0_mat,
        is_single    = is_single,
        is_spatial   = is_spatial,
        ignore_RE    = ignore_RE,
        verbose      = verbose
      ),
      error = function(e) {
        warning(
          "predict() failed for ", basename(rds_file),
          ": ", conditionMessage(e),
          "\nSkipping."
        )
        NULL
      }
    )

    if (is.null(pred)) next

    # ------------------------------------------------------------------
    # Summarize predictions
    # summarize_occupancy_predictions() now handles both 2-D (single) and
    # 3-D (multi-species) psi.0.samples automatically.
    # ------------------------------------------------------------------

    model_results <- summarize_occupancy_predictions(
      pred          = pred,
      coords        = coords_0_df,
      species_names = sp_names,
      template      = base_raster,
      return        = return,
      crs           = model_crs,
      probs         = probs,
      output_dir    = output_dir #normalizePath(output_dir, winslash = "/", mustWork = FALSE)
    )

    # ------------------------------------------------------------------
    # Per-species output: save RDS + optional plots
    # ------------------------------------------------------------------

    for (sp in sp_names) {
      sp_safe <- safe_name(sp)

      cat(sprintf(
        "  [%s] [%d/%d] %s\n",
        basename(rds_file),
        match(sp, sp_names),
        length(sp_names),
        sp
      ))

      species_output_dir <- normalizePath(
        file.path(output_dir, ecosystem_type, sp_safe),
        winslash = "/",
        mustWork = FALSE
      )
      dir.create(species_output_dir, recursive = TRUE, showWarnings = FALSE)

      # Extract this species' result from the combined output
      if (return == "sf") {
        sp_result <- dplyr::filter(model_results, species == sp)
      } else {
        sp_result <- model_results[[sp]]
      }

      all_results[[sp]] <- sp_result

      # Save RDS
      rds_outfile <- file.path(
        species_output_dir,
        paste0(sp_safe, "_", return, "_occupancy_predictions.rds")
      )

      cat("  Saving to:", rds_outfile, "\n")
      saveRDS(sp_result, file = rds_outfile)

      # ------------------------------------------------------------------
      # Plots
      # ------------------------------------------------------------------

      if (save_plots) {
        plot_dir <- file.path(species_output_dir, "plots")
        dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

        if (return == "sf") {
          p <- ggplot2::ggplot(sp_result) +
            ggplot2::geom_sf(
              ggplot2::aes(color = mean),
              size  = 0.15,
              alpha = 0.8
            ) +
            ggplot2::scale_color_viridis_c(
              name   = "Occupancy\nprobability",
              limits = c(0, 1)
            ) +
            ggplot2::coord_sf() +
            ggplot2::theme_bw() +
            ggplot2::labs(
              title    = sp,
              subtitle = paste0(
                ecosystem_type, " | ",
                if (is_single) "single-species" else "multi-species", " | ",
                if (is_spatial) "spatial" else "non-spatial"
              )
            )

        } else {
          r_mean <- sp_result[["mean"]]
          r_df   <- terra::as.data.frame(r_mean, xy = TRUE, na.rm = TRUE)
          names(r_df) <- c("x", "y", "mean")

          p <- ggplot2::ggplot(r_df) +
            ggplot2::geom_raster(
              ggplot2::aes(x = x, y = y, fill = mean)
            ) +
            ggplot2::scale_fill_viridis_c(
              name   = "Occupancy\nprobability",
              limits = c(0, 1)
            ) +
            ggplot2::coord_equal() +
            ggplot2::theme_bw() +
            ggplot2::labs(
              title    = sp,
              subtitle = paste0(
                ecosystem_type, " | ",
                if (is_single) "single-species" else "multi-species", " | ",
                if (is_spatial) "spatial" else "non-spatial"
              ),
              x = NULL,
              y = NULL
            )
        }

        plot_base <- file.path(plot_dir, paste0(sp_safe, "_occupancy_probability"))

        ggplot2::ggsave(
          filename = paste0(plot_base, ".png"),
          plot     = p,
          width    = plot_width,
          height   = plot_height,
          dpi      = plot_dpi
        )

        ggplot2::ggsave(
          filename = paste0(plot_base, ".pdf"),
          plot     = p,
          width    = plot_width,
          height   = plot_height
        )
      }
    }
  }

  invisible(all_results)
}
