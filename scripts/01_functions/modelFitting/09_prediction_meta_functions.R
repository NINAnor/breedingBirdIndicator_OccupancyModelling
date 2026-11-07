# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------
#' Construct Right-Hand Side of a Model Formula
#'
#' Builds the right-hand side (RHS) of a regression or occupancy model formula
#' from a set of covariates, optionally adding a random effect term for route.
#'
#' @param covariates Character vector of covariate names.
#' @param use_route_random_effect Logical. If `TRUE`, appends a random intercept
#'   term \code{(1 | route_id)} to the formula.
#'
#' @return Character string representing the RHS of a model formula.
#'
#' @details
#' If no covariates are provided, the function returns `"1"` (intercept-only model).
#' If \code{use_route_random_effect = TRUE}, the random intercept is added.
#'
#' @examples
#' make_rhs(c("elevation", "forest_cover"))
#' make_rhs(c("elevation"), use_route_random_effect = TRUE)
#'
#' @keywords internal

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


#' Sanitize Strings for Safe File or Object Names
#'
#' Cleans a character string by replacing non-alphanumeric characters with
#' underscores and trimming leading/trailing underscores.
#'
#' @param x Character vector.
#'
#' @return Character vector with safe names suitable for filenames or column names.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Replaces all non-alphanumeric characters (except `-` and `_`) with `_`
#'   \item Removes leading and trailing underscores
#' }
#'
#' @examples
#' safe_name("Haliaeetus albicilla")
#' safe_name("species (test)-1")
#'
#' @keywords internal

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

#' Resolve species names and design-matrix column names from a fitted-model RDS
#'
#' Design-matrix column names are taken from the fitted model's coefficient
#' samples (the authoritative source). If those are unavailable, they are
#' rebuilt from the stored formula and covariates.
#'
#' Note: detection covariates are stored as site x occasion matrices. A single
#' coefficient is estimated per covariate, so the detection design has one
#' column per covariate, not one per covariate-year.
#'
#' @param output List. Loaded RDS object containing model outputs.
#' @param rds_file Character. Path to the RDS file (used for error messages).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{single_species}: Logical indicating model type
#'   \item \code{species_names}: Character vector of species names
#'   \item \code{x_names}: Predictor column names for model matrix
#'   \item \code{occ_formula_fixed}: Processed occupancy formula
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Determines whether the model is single- or multi-species
#'   \item Extracts species names from model object or filename
#'   \item Removes random effects from the occupancy formula
#'   \item Reconstructs the model matrix to identify predictor columns
#' }
#'
#' It ensures compatibility between stored model structure and prediction inputs.
#'
#' @keywords internal
.resolve_model_meta <- function(output, rds_file) {
  model <- output$model_products$model

  if (is.null(model)) {
    stop(
      "No fitted model found in `output$model_products$model` for RDS: ",
      rds_file,
      call. = FALSE
    )
  }

  # ---- Single vs multi-species ------------------------------------------------

  is_single <- if (!is.null(output$single_species)) {
    isTRUE(output$single_species)
  } else {
    is.null(model$sp.names)
  }

  # ---- Species name(s) --------------------------------------------------------

  if (is_single) {
    sp_names <- output$species

    if (is.null(sp_names) || length(sp_names) != 1L || is.na(sp_names)) {
      base <- tools::file_path_sans_ext(basename(rds_file))
      sp_names <- sub("^(spPGOcc|PGOcc)_[^_]+_", "", base)
      sp_names <- gsub("_", " ", sp_names)

      if (length(sp_names) != 1L || !nzchar(sp_names)) {
        stop(
          "Cannot determine species name from single-species RDS: ",
          rds_file,
          "\nExpected `output$species` to be a length-1 character vector, ",
          "or a filename like `spPGOcc_Ecosystem_Genus_species.rds`.",
          call. = FALSE
        )
      }
    }
  } else {
    sp_names <- if (!is.null(model$sp.names)) model$sp.names else output$species

    if (is.null(sp_names) || length(sp_names) == 0L) {
      stop(
        "Cannot determine species names from multi-species RDS: ",
        rds_file,
        "\nExpected `model$sp.names` or `output$species` to be non-NULL.",
        call. = FALSE
      )
    }
  }

  # ---- Formulas (random-effect terms removed) ---------------------------------

  occ_formula_fixed <- .fixed_effect_formula(
    output$formulas$occ_formula,
    "occupancy",
    rds_file
  )
  det_formula_fixed <- .fixed_effect_formula(
    output$formulas$det_formula,
    "detection",
    rds_file
  )

  # ---- Design-matrix column names ---------------------------------------------

  # Fitted coefficient names: preferred source.
  occ_samples <- if (is_single) model$beta.samples else model$beta.comm.samples
  det_samples <- if (is_single) {
    model$alpha.samples
  } else {
    model$alpha.comm.samples
  }

  x_names <- colnames(occ_samples)
  x_p_names <- colnames(det_samples)

  # Fallback: rebuild from stored covariates.
  if (is.null(x_names) || !length(x_names)) {
    if (is.null(output$data_list$occ.covs)) {
      stop(
        "No occupancy covariates found at `output$data_list$occ.covs` for RDS: ",
        rds_file,
        call. = FALSE
      )
    }
    x_names <- colnames(
      stats::model.matrix(occ_formula_fixed, data = output$data_list$occ.covs)
    )
  }

  if (is.null(x_p_names) || !length(x_p_names)) {
    if (is.null(output$data_list$det.covs)) {
      stop(
        "No detection covariates found at `output$data_list$det.covs` for RDS: ",
        rds_file,
        call. = FALSE
      )
    }
    # Flatten site x occasion matrices to vectors, as spOccupancy does.
    det_df <- as.data.frame(lapply(output$data_list$det.covs, as.vector))
    chr <- vapply(det_df, is.character, logical(1))
    det_df[chr] <- lapply(det_df[chr], as.factor)
    x_p_names <- colnames(stats::model.matrix(det_formula_fixed, data = det_df))
  }

  if (!length(x_names) || !length(x_p_names)) {
    stop(
      "Cannot determine model-matrix predictor columns for RDS: ",
      rds_file,
      call. = FALSE
    )
  }

  list(
    single_species = is_single,
    species_names = sp_names,
    x_names = x_names,
    x_p_names = x_p_names,
    occ_formula_fixed = occ_formula_fixed,
    det_formula_fixed = det_formula_fixed
  )
}


# Drop random-effect terms such as (1 | route_id), keeping the intercept state.
.fixed_effect_formula <- function(formula, label, rds_file) {
  if (is.null(formula)) {
    stop(
      "No ",
      label,
      " formula found at `output$formulas` for RDS: ",
      rds_file,
      call. = FALSE
    )
  }

  tt <- stats::terms(formula)
  labs <- attr(tt, "term.labels")
  labs <- labs[!grepl("|", labs, fixed = TRUE)]
  if (!length(labs)) {
    labs <- "1"
  }

  out <- stats::reformulate(labs, intercept = attr(tt, "intercept") == 1L)
  environment(out) <- environment(formula)
  out
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

#' Internal Wrapper for Occupancy Model Prediction
#'
#' Calls the appropriate \code{predict()} method depending on whether the model
#' is single- or multi-species and spatial or non-spatial.
#'
#' @param model Fitted occupancy model object.
#' @param X_0 Numeric matrix of covariates for prediction.
#' @param coords_0_mat Numeric matrix of coordinates (longitude, latitude).
#' @param is_single Logical. Whether the model is single-species.
#' @param is_spatial Logical. Whether the model includes spatial structure.
#' @param ignore_RE Logical. Whether to ignore random effects during prediction.
#' @param verbose Logical. Whether to print progress messages.
#' @param n_omp_threads Integer. OpenMP threads for spatial models.
#'
#' @return Prediction object returned by \code{predict()}.
#'
#' @details
#' Handles differences between model classes:
#' \itemize{
#'   \item \code{spPGOcc}, \code{PGOcc} (single-species)
#'   \item \code{spMsPGOcc}, \code{msPGOcc} (multi-species)
#' }
#'
#' Ensures correct handling of:
#' \itemize{
#'   \item \code{coords.0} argument for spatial models
#'   \item \code{ignore.RE} argument for supported models
#' }
#'
#' @keywords internal
#'

.run_predict <- function(
  model,
  X_0,
  coords_0_mat = NULL,
  is_single,
  is_spatial,
  ignore_RE = TRUE,
  verbose = TRUE,
  type = c("occupancy", "detection"),
  n_omp_threads = 1L
) {
  type <- match.arg(type)

  # Common prediction arguments
  predict_args <- list(
    X.0 = X_0,
    type = type,
    ignore.RE = ignore_RE,
    verbose = verbose
  )

  # Spatial models require prediction coordinates
  if (is_spatial) {
    if (is.null(coords_0_mat)) {
      stop(
        "`coords_0_mat` must be supplied for spatial models.",
        call. = FALSE
      )
    }

    predict_args$coords.0 <- coords_0_mat
    # Only the spatial predict() methods accept n.omp.threads
    predict_args$n.omp.threads <- as.integer(n_omp_threads)
  }

  # Dispatch to the predict method for the fitted spOccupancy model
  do.call(
    stats::predict,
    c(list(model), predict_args)
  )
}


# ---- Helpers ------------------------------------------------------------------

# Year-suffixed covariate names, with or without an underscore before the year.
.year_pattern <- "^(.*?)_?([0-9]{4})$"

#' Build a fixed-effect design matrix aligned to the fitted columns.
#'
#' If `formula` is NULL the formula is derived from `expected_names`
#' (used for year-specific detection covariates).
.make_prediction_matrix <- function(
  formula = NULL,
  data,
  expected_names,
  type = c("occupancy", "detection"),
  label = type
) {
  type <- match.arg(type)

  if (is.null(expected_names) || !length(expected_names)) {
    stop(
      "`expected_names` must contain the fitted model matrix columns.",
      call. = FALSE
    )
  }
  expected_names <- as.character(expected_names)

  # Detection: create year-specific columns the fitted model expects
  # (mean_temperature2020 <- mean_temperature_2020, or vice versa).
  if (type == "detection") {
    year_names <- grep(.year_pattern, expected_names, value = TRUE, perl = TRUE)
    for (nm in year_names) {
      if (nm %in% names(data)) {
        next
      }
      parts <- regmatches(nm, regexec(.year_pattern, nm, perl = TRUE))[[1]]
      if (length(parts) != 3L) {
        next
      }
      candidates <- c(
        paste0(parts[2], "_", parts[3]), # underscore form preferred
        paste0(parts[2], parts[3])
      )
      available <- candidates[candidates %in% names(data)]
      if (!length(available)) {
        stop(
          label,
          ": cannot find a prediction variable for expected column '",
          nm,
          "'. Tried: ",
          paste(candidates, collapse = ", "),
          call. = FALSE
        )
      }
      data[[nm]] <- data[[available[1L]]]
    }
  }

  # Fixed-effect formula (random-effect terms such as (1 | route_id) removed).
  if (is.null(formula)) {
    term_chr <- setdiff(expected_names, "(Intercept)")
    rhs <- if (length(term_chr)) {
      paste0("`", term_chr, "`", collapse = " + ")
    } else {
      "1"
    }
    if (!"(Intercept)" %in% expected_names) {
      rhs <- paste(rhs, "- 1")
    }
    fixed_formula <- stats::as.formula(paste("~", rhs))
  } else {
    tt <- stats::terms(formula)
    labs <- attr(tt, "term.labels")
    labs <- labs[!grepl("|", labs, fixed = TRUE)]
    if (!length(labs)) {
      labs <- "1"
    }
    fixed_formula <- stats::reformulate(
      labs,
      intercept = attr(tt, "intercept") == 1L
    )
    environment(fixed_formula) <- environment(formula)
  }

  missing_vars <- setdiff(all.vars(fixed_formula), names(data))
  if (length(missing_vars)) {
    stop(
      label,
      ": formula variables missing from prediction data: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  # na.pass keeps row alignment with the coordinates (no silent row dropping).
  mf <- stats::model.frame(
    fixed_formula,
    data = data,
    na.action = stats::na.pass
  )
  X <- stats::model.matrix(fixed_formula, data = mf)

  if (nrow(X) != nrow(data)) {
    stop(
      label,
      ": design matrix has ",
      nrow(X),
      " rows but data has ",
      nrow(data),
      ".",
      call. = FALSE
    )
  }
  if (anyNA(X)) {
    stop(label, ": design matrix contains NA values.", call. = FALSE)
  }

  missing_cols <- setdiff(expected_names, colnames(X))
  if (length(missing_cols)) {
    stop(
      label,
      ": expected columns could not be constructed: ",
      paste(missing_cols, collapse = ", "),
      "\nGenerated columns: ",
      paste(colnames(X), collapse = ", "),
      call. = FALSE
    )
  }
  extra_cols <- setdiff(colnames(X), expected_names)
  if (length(extra_cols)) {
    warning(
      label,
      ": dropping unexpected columns: ",
      paste(extra_cols, collapse = ", "),
      call. = FALSE
    )
  }

  X[, expected_names, drop = FALSE]
}


# Summarise posterior probability draws for one block.
# Single-species draws: samples x locations.
# Multi-species draws: samples x species x locations.
# Returns a named list (by species) of plain data.frames; no sf work here.
.summarise_draws <- function(
  draws,
  species_names,
  coords,
  location_ids,
  probs,
  is_single,
  prefix = ""
) {
  if (inherits(draws, "mcmc") || inherits(draws, "mcmc.list")) {
    draws <- as.matrix(draws)
  }
  d <- dim(draws)

  if (is.null(d) || length(d) == 2L) {
    if (!is_single) {
      stop("Expected a 3-D posterior array for a multi-species model.")
    }
    mat <- as.matrix(draws)
    if (ncol(mat) != nrow(coords)) {
      stop(
        "Single-species draw matrix has ",
        ncol(mat),
        " locations but coordinates have ",
        nrow(coords),
        " rows."
      )
    }
    by_species <- stats::setNames(list(mat), species_names[[1L]])
  } else if (length(d) == 3L) {
    if (d[2] != length(species_names) || d[3] != nrow(coords)) {
      stop("Posterior array dimensions do not match species and coordinates.")
    }
    by_species <- lapply(seq_along(species_names), function(i) {
      matrix(draws[, i, ], nrow = d[1L], ncol = d[3L]) # safe for 1 location
    })
    names(by_species) <- species_names
  } else {
    stop("Unsupported posterior draw dimensions: ", paste(d, collapse = " x "))
  }

  use_ms <- requireNamespace("matrixStats", quietly = TRUE)
  value_names <- paste0(prefix, c("mean", "median", "sd", "lower", "upper"))

  lapply(names(by_species), function(sp) {
    mat <- by_species[[sp]]
    if (use_ms) {
      q <- matrixStats::colQuantiles(
        mat,
        probs = probs,
        na.rm = TRUE,
        useNames = FALSE
      )
      vals <- list(
        matrixStats::colMeans2(mat, na.rm = TRUE),
        matrixStats::colMedians(mat, na.rm = TRUE),
        matrixStats::colSds(mat, na.rm = TRUE),
        q[, 1L],
        q[, 2L]
      )
    } else {
      # One quantile pass for lower, median and upper together.
      q <- apply(
        mat,
        2,
        stats::quantile,
        probs = c(probs[1], 0.5, probs[2]),
        na.rm = TRUE,
        names = FALSE
      )
      vals <- list(
        colMeans(mat, na.rm = TRUE),
        q[2L, ],
        apply(mat, 2, stats::sd, na.rm = TRUE),
        q[1L, ],
        q[3L, ]
      )
    }
    names(vals) <- value_names
    data.frame(
      species = sp,
      location_id = location_ids,
      x = coords$long,
      y = coords$lat,
      vals,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }) |>
    stats::setNames(names(by_species))
}


# Rasterise one species' sf results and write a compressed GeoTIFF.
.write_species_raster <- function(
  sp_sf,
  base_raster,
  prediction_type,
  raster_file
) {
  type_cols <- if (prediction_type == "both") {
    as.vector(outer(
      c("occupancy_", "detection_"),
      c("mean", "median", "sd", "lower", "upper"),
      paste0
    ))
  } else {
    c("mean", "median", "sd", "lower", "upper")
  }
  type_cols <- intersect(type_cols, names(sp_sf))
  if (!length(type_cols)) {
    return(NULL)
  }

  v <- terra::vect(sp_sf)
  rasters <- lapply(type_cols, function(nm) {
    terra::rasterize(v, base_raster, field = nm, fun = "mean")
  })
  sp_raster <- terra::rast(rasters)
  names(sp_raster) <- type_cols

  terra::writeRaster(
    sp_raster,
    raster_file,
    overwrite = TRUE,
    wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
  )
  sp_raster
}
