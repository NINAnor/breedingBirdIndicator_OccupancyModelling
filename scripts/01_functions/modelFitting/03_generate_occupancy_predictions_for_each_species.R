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
    ignore_RE = TRUE,
    save_plots = TRUE,
    verbose = FALSE,
    plot_width = 8,
    plot_height = 6,
    plot_dpi = 300
) {
  return <- match.arg(return)

  if (!dir.exists(input_dir)) {
    stop("`input_dir` does not exist or is not a directory: ", input_dir, call. = FALSE)
  }

  if (return == "raster" && is.null(base_raster)) {
    stop("`base_raster` must be supplied when `return = 'raster'`.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  rds_files <- list.files(
    input_dir,
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(rds_files) == 0) {
    stop("No .rds files found in `input_dir`.", call. = FALSE)
  }

  all_results <- list()

  coords_0 <- prediction_points |>
    dplyr::select(longitude, latitude)

  occ_rhs <- make_rhs(
    covariates = occ_covariates,
    use_route_random_effect = use_route_random_effect
  )

  occ_formula <- as.formula(paste("~", occ_rhs))

  for (rds_file in rds_files) {

    cat("\nReading:", rds_file, "\n")

    output <- tryCatch(
      readRDS(rds_file),
      error = function(e) {
        stop(
          "Failed to read:\n",
          rds_file,
          "\n\n",
          conditionMessage(e)
        )
      }
    )

    if (is.null(output$model_products$model)) {
      warning("Skipping file without `output$model_products$model`: ", rds_file)
      next
    }

    model <- output$model_products$model

    if (is.null(model$sp.names) || is.null(model$x.names)) {
      warning("Skipping model missing `sp.names` or `x.names`: ", rds_file)
      next
    }

    X_0_full <- model.matrix(
      occ_formula,
      data = prediction_points
    )

    missing_cols <- setdiff(model$x.names, colnames(X_0_full))

    if (length(missing_cols) > 0) {
      stop(
        "Prediction matrix is missing required model columns for ",
        basename(rds_file),
        ": ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }

    X_0 <- X_0_full[, model$x.names, drop = FALSE]

    pred <- predict(
      model,
      X.0 = X_0,
      coords.0 = coords_0,
      ignore.RE = ignore_RE,
      verbose = verbose
    )

    model_results <- summarize_occupancy_predictions(
      pred = pred,
      coords = coords_0,
      species_names = model$sp.names,
      template = base_raster,
      return = return,
      crs = model_crs,
      probs = probs,
      output_dir = file.path(output_dir, ecosystem_type)
    )

    for (sp in model$sp.names) {
      sp_safe <- safe_name(sp)

      cat(
        sprintf(
          "[Model %s] [%d/%d] %s\n",
          basename(rds_file),
          match(sp, model$sp.names),
          length(model$sp.names),
          sp
        )
      )

      species_output_dir <- file.path(
        output_dir,
        ecosystem_type,
        sp_safe
      )

      dir.create(species_output_dir, recursive = TRUE, showWarnings = FALSE)

      if (return == "sf") {
        sp_result <- dplyr::filter(model_results, species == sp)
      } else {
        sp_result <- model_results[[sp]]
      }

      all_results[[sp]] <- sp_result

      species_output_dir <- normalizePath(
        species_output_dir,
        winslash = "/",
        mustWork = TRUE
      )

      rds_outfile <- file.path(
        species_output_dir,
        paste0(sp_safe, "_", return, "_occupancy_predictions.rds")
      )

      cat("Saving to:", rds_outfile, "\n")
      cat("Directory exists:", dir.exists(dirname(rds_outfile)), "\n")

      saveRDS(
        sp_result,
        file = rds_outfile
      )

      if (save_plots) {
        plot_dir <- file.path(species_output_dir, "plots")
        dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

        if (return == "sf") {
          p <- ggplot2::ggplot(sp_result) +
            ggplot2::geom_sf(
              ggplot2::aes(color = mean),
              size = 0.15,
              alpha = 0.8
            ) +
            ggplot2::scale_color_viridis_c(
              name = "Occupancy\nprobability",
              limits = c(0, 1)
            ) +
            ggplot2::coord_sf() +
            ggplot2::theme_bw() +
            ggplot2::labs(
              title = sp,
              subtitle = ecosystem_type
            )
        } else {
          r_mean <- sp_result[["mean"]]

          r_df <- terra::as.data.frame(
            r_mean,
            xy = TRUE,
            na.rm = TRUE
          )

          names(r_df) <- c("x", "y", "mean")

          p <- ggplot2::ggplot(r_df) +
            ggplot2::geom_raster(
              ggplot2::aes(x = x, y = y, fill = mean)
            ) +
            ggplot2::scale_fill_viridis_c(
              name = "Occupancy\nprobability",
              limits = c(0, 1)
            ) +
            ggplot2::coord_equal() +
            ggplot2::theme_bw() +
            ggplot2::labs(
              title = sp,
              subtitle = ecosystem_type,
              x = NULL,
              y = NULL
            )
        }

        ggplot2::ggsave(
          filename = file.path(plot_dir, paste0(sp_safe, "_occupancy_probability.png")),
          plot = p,
          width = plot_width,
          height = plot_height,
          dpi = plot_dpi
        )

        ggplot2::ggsave(
          filename = file.path(plot_dir, paste0(sp_safe, "_occupancy_probability.pdf")),
          plot = p,
          width = plot_width,
          height = plot_height
        )
      }
    }
  }

  invisible(all_results)
}
