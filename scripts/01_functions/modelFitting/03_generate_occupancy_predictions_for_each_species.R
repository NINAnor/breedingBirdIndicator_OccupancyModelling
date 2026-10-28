# -------------------------------------------------------------------------
# Main batch prediction function
# -------------------------------------------------------------------------
#' Generate Occupancy Predictions for Each Species from Fitted Models
#'
#' Loads fitted occupancy models from RDS files and generates predictions for
#' each species at specified locations. Supports single- and multi-species models,
#' spatial and non-spatial variants, and produces outputs as either point-based
#' \pkg{sf} objects or raster layers.
#'
#' This function supports FAIR principles:
#' \itemize{
#'   \item \strong{Findable}: Outputs are organized by ecosystem and species.
#'   \item \strong{Accessible}: Results saved as RDS, plots (PNG/PDF), and spatial formats.
#'   \item \strong{Interoperable}: Uses standard formats (\pkg{sf}, GeoTIFF).
#'   \item \strong{Reusable}: Includes metadata-rich outputs and consistent naming.
#' }
#'
#' @param input_dir Character. Directory containing `.rds` model files.
#' @param prediction_points Data frame with coordinates and covariates for prediction.
#' @param det_covariates Character vector of detection covariates (not currently used).
#' @param occ_covariates Character vector of occupancy covariates.
#' @param base_raster Optional \pkg{terra} raster object used when `return = "raster"`.
#' @param model_crs Coordinate reference system for output data.
#' @param ecosystem_type Character. Ecosystem label used in output paths.
#' @param output_dir Character. Directory where outputs should be saved.
#' @param return Character. Output type: `"sf"` or `"raster"`.
#' @param probs Numeric vector of quantiles for uncertainty intervals.
#' @param use_route_random_effect Logical. Include random intercept for route.
#' @param ignore_RE Logical. Whether to ignore random effects during prediction
#'   (used in multi-species models).
#' @param save_plots Logical. If `TRUE`, saves prediction plots.
#' @param verbose Logical. Print detailed progress messages.
#' @param plot_width,plot_height Numeric. Plot dimensions in inches.
#' @param plot_dpi Numeric. Plot resolution.
#'
#' @return Invisibly returns a named list of per-species prediction outputs.
#'
#' @details
#' The function performs:
#' \enumerate{
#'   \item Loads model objects from `.rds` files
#'   \item Resolves species identity and model structure
#'   \item Builds a prediction design matrix
#'   \item Aligns predictor variables with model requirements
#'   \item Runs predictions using appropriate model methods
#'   \item Summarizes predictions into mean and uncertainty metrics
#'   \item Outputs results as `sf` or raster formats
#'   \item Saves per-species outputs and optional plots
#' }
#'
#' Predictions include occupancy probability estimates summarized across:
#' \itemize{
#'   \item mean occupancy probability
#'   \item uncertainty intervals (defined by \code{probs})
#' }
#'
#' Output structure:
#' \itemize{
#'   \item One folder per species
#'   \item Serialized results (`.rds`)
#'   \item Optional visualization plots (`.png`, `.pdf`)
#' }
#'
#' @examples
#' \dontrun{
#' results <- generate_occupancy_predictions_for_each_species(
#'   input_dir = "models/",
#'   prediction_points = df_points,
#'   occ_covariates = c("elevation", "forest_cover"),
#'   model_crs = 25833,
#'   ecosystem_type = "Forest",
#'   output_dir = "outputs/",
#'   return = "sf"
#' )
#' }
#'
#' @seealso
#' [predict()], [sf::st_as_sf()], [terra::rast()]
#'
#' @export

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
  ignore_RE = TRUE, # passed to multi-species predict only
  save_plots = TRUE,
  verbose = FALSE,
  plot_width = 8,
  plot_height = 6,
  plot_dpi = 300
) {
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
      "`base_raster` must be supplied when `return = 'raster'`.",
      call. = FALSE
    )
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

  # Drop any NAs in the prediction data. It causes problems with aliggnment of prediction data

  prediction_points <- prediction_points %>%
    na.omit()

  # Build the occupancy design matrix once — columns are aligned per-model
  # below using model$x.names, so this only needs to contain a superset of
  # all possible predictors.
  occ_rhs <- make_rhs(occ_covariates, use_route_random_effect)
  occ_formula <- as.formula(paste("~", occ_rhs))

  if (use_route_random_effect) {
    prediction_points <- prediction_points %>%
      dplyr::mutate(route_id = 1)
  }

  prediction_points <- prediction_points %>%
    dplyr::mutate(
      mean_temperature = rowMeans(
        dplyr::across(
          dplyr::matches("^mean_temperature_[0-9]{4}$")
        ),
        na.rm = TRUE
      ),

      mean_precipitation = rowMeans(
        dplyr::across(
          dplyr::matches("^mean_precipitation_[0-9]{4}$")
        ),
        na.rm = TRUE
      )
    )

  X_0_full <- model.matrix(occ_formula, data = prediction_points)

  # Coordinates as a plain numeric matrix — used by spatial predict methods.
  # Column names must be present but the spatial predict generics accept any
  # two-column matrix; we name them to match the convention used during fitting.
  coords_0_mat <- as.matrix(
    prediction_points[, c("long", "lat")]
  )

  # sf / terra coordinate data frame used by summarize_occupancy_predictions()
  coords_0_df <- prediction_points[, c("long", "lat")]

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

    if (is.null(meta)) {
      next
    }

    is_single <- meta$single_species
    sp_names <- meta$species_names
    x_names <- meta$x_names

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
        basename(rds_file),
        ": ",
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
          "predict() failed for ",
          basename(rds_file),
          ": ",
          conditionMessage(e),
          "\nSkipping."
        )
        NULL
      }
    )

    if (is.null(pred)) {
      next
    }

    # ------------------------------------------------------------------
    # Summarize predictions
    # summarize_occupancy_predictions() now handles both 2-D (single) and
    # 3-D (multi-species) psi.0.samples automatically.
    # ------------------------------------------------------------------

    model_results <- summarize_occupancy_predictions(
      pred = pred,
      coords = coords_0_df,
      species_names = sp_names,
      template = base_raster,
      return = return,
      crs = model_crs,
      probs = probs,
      output_dir = output_dir #normalizePath(output_dir, winslash = "/", mustWork = FALSE)
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
        file.path(output_dir, sp_safe),
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
              subtitle = paste0(
                ecosystem_type,
                " | ",
                if (is_single) "single-species" else "multi-species",
                " | ",
                if (is_spatial) "spatial" else "non-spatial"
              )
            )
        } else {
          r_mean <- sp_result[["mean"]]
          r_df <- terra::as.data.frame(r_mean, xy = TRUE, na.rm = TRUE)
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
              subtitle = paste0(
                ecosystem_type,
                " | ",
                if (is_single) "single-species" else "multi-species",
                " | ",
                if (is_spatial) "spatial" else "non-spatial"
              ),
              x = NULL,
              y = NULL
            )
        }

        plot_base <- file.path(
          plot_dir,
          paste0(sp_safe, "_occupancy_probability")
        )

        ggplot2::ggsave(
          filename = paste0(plot_base, ".png"),
          plot = p,
          width = plot_width,
          height = plot_height,
          dpi = plot_dpi
        )

        ggplot2::ggsave(
          filename = paste0(plot_base, ".pdf"),
          plot = p,
          width = plot_width,
          height = plot_height
        )
      }
    }
  }

  invisible(all_results)
}
