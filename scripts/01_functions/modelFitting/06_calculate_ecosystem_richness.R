#' Calculate Multispecies Ecosystem Richness from Species Occupancy Predictions
#'
#' Aggregates species-specific occupancy prediction outputs to calculate
#' multispecies expected ecosystem richness. For each species, the function
#' extracts the mean, median, lower, and upper occupancy prediction layers
#' and sums them cell-by-cell across all species.
#'
#' The resulting raster contains four layers:
#' \itemize{
#'   \item \strong{mean}: Sum of species-specific mean occupancy predictions.
#'   \item \strong{median}: Sum of species-specific median occupancy predictions.
#'   \item \strong{lower}: Sum of species-specific lower occupancy predictions.
#'   \item \strong{upper}: Sum of species-specific upper occupancy predictions.
#' }
#'
#' The resulting mean layer represents the expected number of species occurring
#' at each prediction location, based on the species-specific occupancy
#' probabilities. The lower and upper layers represent the corresponding
#' cell-wise sums of the species-specific uncertainty bounds.
#'
#' The function supports species prediction outputs stored as GeoTIFF rasters
#' or GeoParquet files and preserves the spatial reference and resolution of
#' the input prediction surfaces.
#'
#' This function supports FAIR principles:
#' \itemize{
#'   \item \strong{Findable}: Outputs use descriptive filenames and structured
#'         metadata.
#'   \item \strong{Accessible}: Raster outputs are written as GeoTIFF and
#'         intermediate results can be serialized as RDS.
#'   \item \strong{Interoperable}: Spatial outputs preserve CRS information
#'         and use standard spatial formats.
#'   \item \strong{Reusable}: Processing parameters, species names, input
#'         locations, software versions, and output formats are recorded in
#'         JSON metadata.
#' }
#'
#' @param input_dir Character. Directory containing species prediction outputs.
#'   Each species should have a subdirectory containing a GeoTIFF (`.tif`) or
#'   GeoParquet (`.parquet`) prediction file.
#' @param reference_shape An `sf` object defining reference locations.
#'   Retained for compatibility with the broader workflow but not used in the
#'   current multispecies raster aggregation.
#' @param categorisation_shape An `sf` object containing region categories.
#'   Retained for compatibility with the broader workflow but not used in the
#'   current multispecies raster aggregation.
#' @param ecosystem_type Character. Ecosystem type associated with the
#'   prediction outputs.
#' @param value_col Character. Name of the prediction layer. Default is
#'   `"mean"`. Retained for compatibility with the broader workflow; the
#'   current implementation aggregates the `mean`, `median`, `lower`, and
#'   `upper` layers.
#' @param point_spacing Numeric. Sampling spacing in metres. Retained for
#'   compatibility with the broader workflow.
#' @param region_col Character. Column name defining regions. Retained for
#'   compatibility with the broader workflow.
#' @param set_lower_bound Numeric or `NULL`. Retained for compatibility with
#'   the broader workflow.
#' @param forest_only_strict Logical. Retained for compatibility with the
#'   broader workflow.
#' @param na.rm Logical. Logical indicating whether missing values should be
#'   removed during aggregation. Default is `TRUE`.
#' @param output_dir Character or `NULL`. Directory where outputs and metadata
#'   are saved when `save_outputs = TRUE`.
#' @param save_outputs Logical. If `TRUE`, saves the aggregated richness
#'   raster, metadata, and complete results object to disk.
#' @param overwrite Logical. If `TRUE`, existing output files are overwritten.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{results}: Named list of species-specific prediction objects.
#'   \item \code{expected_richness}: A `SpatRaster` containing four layers:
#'         `mean`, `median`, `lower`, and `upper`, representing the
#'         cell-wise sum of the corresponding prediction layers across all
#'         species.
#' }
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{results}: List of species prediction objects.
#'   \item \code{reference_shape_with_categories}: Reference features with
#'         assigned region categories.
#'   \item \code{reference_predictions}: Extracted species predictions at reference locations.
#'   \item \code{reference_predictions_expected_richness}: Regional summary of
#'         expected richness for reference conditions.
#'   \item \code{prediction_values_scaled_richness}: Per-point scaled richness values.
#'   \item \code{scaled_richness_summary}: Regional summary statistics.
#'   \item \code{lower_bound}: Minimum expected richness used in scaling.
#'   \item \code{regional_reference}: Maximum expected richness per region.
#' }
#'
#' @details
#' The function performs the following workflow:
#' \enumerate{
#'   \item Validates the input directory.
#'   \item Identifies species-specific subdirectories.
#'   \item Loads the first available GeoTIFF or GeoParquet prediction file
#'         from each species directory.
#'   \item Extracts the `mean`, `median`, `lower`, and `upper` prediction
#'         layers from each species.
#'   \item Sums each prediction layer across all species on a cell-by-cell
#'         basis.
#'   \item Combines the four summed layers into a single `SpatRaster`.
#'   \item Names the output layers `mean`, `median`, `lower`, and `upper`.
#'   \item Optionally writes the aggregated richness raster as a GeoTIFF.
#'   \item Optionally writes processing metadata as `metadata.json` and the
#'         complete output object as `all_results.rds`.
#' }
#'
#' The multispecies expected richness is calculated as:
#' \deqn{
#'   R_{mean}(x) = \sum_{s=1}^{S} p_{s,mean}(x)
#' }
#' where \eqn{p_{s,mean}(x)} is the mean predicted occupancy probability
#' for species \eqn{s} at location \eqn{x}, and \eqn{S} is the number of
#' species.
#'
#' The median, lower, and upper layers are calculated analogously by summing
#' the corresponding species-specific prediction layers.
#'
#' Note that the summed lower and upper layers are the cell-wise sums of the
#' species-specific marginal uncertainty bounds. They should therefore not be
#' interpreted as the exact posterior credible interval for the multispecies
#' sum unless the underlying posterior samples are aggregated before
#' calculating the interval.
#'
#' Output files (if \code{save_outputs = TRUE}) include:
#' \itemize{
#'   \item Reference and prediction datasets (GeoParquet / GeoTIFF)
#'   \item Regional summaries (CSV)
#'   \item Species-specific outputs
#'   \item Metadata file (`metadata.json`)
#'   \item Serialized results (`all_results.rds`)
#' }
#'
#' @examples
#' \dontrun{
#' result <- calculate_scaled_ecosystem_richness(
#'   input_dir = "data/species_predictions",
#'   reference_shape = ref_sf,
#'   categorisation_shape = region_sf,
#'   ecosystem_type = "Forest",
#'   output_dir = "outputs/richness",
#'   save_outputs = TRUE
#' )
#'
#' result$scaled_richness_summary
#' }
#'
#' @seealso
#' [terra::rast()], [sf::st_intersection()], [sfarrow::st_read_parquet()],
#' [dplyr::summarise()]
#'
#' @export

calculate_ecosystem_richness <- function(
  input_dir,
  ecosystem_type,
  value_col = "mean",
  na.rm = TRUE,
  output_dir = NULL,
  save_outputs = FALSE,
  overwrite = TRUE,
  method = "rds"
) {
  # ============================================================
  # 1. Validate arguments
  # ============================================================

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

  # ============================================================
  # 2. Windows short-path helper
  #
  # Required because the OneDrive paths are too long for some
  # R/Arrow/Windows file operations.
  # ============================================================

  short_path <- function(path) {
    if (length(path) != 1 || is.na(path)) {
      stop(
        "Invalid file path supplied to `short_path()`.",
        call. = FALSE
      )
    }

    result <- system(
      sprintf(
        'cmd /c for %%I in ("%s") do @echo %%~sI',
        path
      ),
      intern = TRUE
    )

    result <- trimws(result)

    if (length(result) != 1 || !nzchar(result)) {
      stop(
        "Could not convert path to Windows short path:\n",
        path,
        call. = FALSE
      )
    }

    result
  }

  # ============================================================
  # 3. Identify species directories
  # ============================================================

  if (length(input_dir) == 1) {
    species_dirs <- list.dirs(
      input_dir,
      recursive = FALSE,
      full.names = TRUE
    )
  } else {
    # input_dir already contains species directories
    species_dirs <- input_dir
  }

  # Remove internal/output directories such as _block_summaries
  species_dirs <- species_dirs[
    !grepl("^_", basename(species_dirs))
  ]

  if (length(species_dirs) == 0) {
    stop(
      "No species directories were found.",
      call. = FALSE
    )
  }

  # ============================================================
  # 4. Function for locating the prediction file
  # ============================================================

  find_prediction_file <- function(sp_dir) {
    files <- list.files(
      sp_dir,
      full.names = TRUE,
      recursive = FALSE
    )

    pattern <- switch(
      method,

      raster = "\\.tif$",

      sf = "\\.parquet$",

      rds = "\\.rds$"
    )

    prediction_files <- files[
      grepl(
        pattern,
        files,
        ignore.case = TRUE
      )
    ]

    if (length(prediction_files) == 0) {
      stop(
        "No ",
        method,
        " prediction file found in:\n",
        sp_dir,
        call. = FALSE
      )
    }

    # Use the first matching prediction file
    prediction_files[1]
  }

  # ============================================================
  # 5. Read one species prediction
  # ============================================================

  read_prediction <- function(sp_dir) {
    prediction_file <- find_prediction_file(sp_dir)

    # Convert long OneDrive path to Windows short path
    file_short <- short_path(prediction_file)

    switch(
      method,

      raster = {
        terra::rast(file_short)
      },

      sf = {
        sfarrow::st_read_parquet(file_short)
      },

      rds = {
        readRDS(file_short)
      }
    )
  }

  # ============================================================
  # 6. Read species predictions
  # ============================================================

  results <- setNames(
    lapply(
      species_dirs,
      read_prediction
    ),
    basename(species_dirs)
  )

  if (length(results) == 0) {
    stop(
      "No species predictions were successfully loaded.",
      call. = FALSE
    )
  }

  # ============================================================
  # 7. Validate prediction objects
  # ============================================================

  required_layers <- c(
    "mean",
    "median",
    "lower",
    "upper"
  )

  if (method == "raster") {
    invalid <- names(results)[
      !vapply(
        results,
        function(x) {
          all(required_layers %in% names(x))
        },
        logical(1)
      )
    ]

    if (length(invalid) > 0) {
      stop(
        "The following species are missing one or more required ",
        "raster layers:\n",
        paste(invalid, collapse = "\n"),
        call. = FALSE
      )
    }
  }

  # ============================================================
  # 8. Aggregate species predictions
  # ============================================================

  if (method == "raster") {
    # ----------------------------------------------------------
    # Raster aggregation
    # ----------------------------------------------------------

    sum_occupancy <- results[[1]][
      c("mean", "median", "lower", "upper")
    ]

    names(sum_occupancy) <- required_layers

    # Add remaining species
    if (length(results) > 1) {
      for (i in seq.int(2, length(results))) {
        x <- results[[i]][required_layers]

        # Make sure all rasters have matching geometry
        if (
          !terra::compareGeom(
            sum_occupancy,
            x,
            stopOnError = FALSE
          )
        ) {
          stop(
            "Spatial geometry does not match for species: ",
            names(results)[i],
            call. = FALSE
          )
        }

        if (na.rm) {
          for (j in seq_along(required_layers)) {
            sum_occupancy[[j]] <- terra::ifel(
              is.na(sum_occupancy[[j]]),
              x[[j]],
              sum_occupancy[[j]]
            )

            sum_occupancy[[j]] <- terra::ifel(
              !is.na(x[[j]]),
              sum_occupancy[[j]] + x[[j]],
              sum_occupancy[[j]]
            )
          }
        } else {
          sum_occupancy <- sum_occupancy + x
        }
      }
    }
  } else {
    # ----------------------------------------------------------
    # sf / RDS aggregation
    # ----------------------------------------------------------

    sum_occupancy <- results[[1]]

    # Check that required columns exist
    missing_columns <- setdiff(
      required_layers,
      names(sum_occupancy)
    )

    if (length(missing_columns) > 0) {
      stop(
        "The first prediction object is missing: ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }

    # Check that all species have the same number of prediction
    # locations and the same geometry
    if (length(results) > 1) {
      for (i in seq.int(2, length(results))) {
        x <- results[[i]]

        if (nrow(x) != nrow(sum_occupancy)) {
          stop(
            "Number of prediction locations differs for species: ",
            names(results)[i],
            call. = FALSE
          )
        }

        # if (
        #   inherits(sum_occupancy, "sf") &&
        #     inherits(x, "sf")
        # ) {
        #   if (
        #     !identical(
        #       sf::st_crs(sum_occupancy),
        #       sf::st_crs(x)
        #     )
        #   ) {
        #     stop(
        #       "CRS differs for species: ",
        #       names(results)[i],
        #       call. = FALSE
        #     )
        #   }

        #  if (
        #    !isTRUE(
        #      all(
        #        sf::st_equals(
        #          sf::st_geometry(sum_occupancy),
        #          sf::st_geometry(x),
        #          sparse = FALSE
        #        )
        #      )
        #    )
        #  ) {
        #    stop(
        #      "Prediction geometries differ for species: ",
        #      names(results)[i],
        #      call. = FALSE
        #    )
        #  }
        #}

        # Add each prediction layer
        for (layer in required_layers) {
          if (na.rm) {
            sum_occupancy[[layer]] <- rowSums(
              cbind(
                sum_occupancy[[layer]],
                x[[layer]]
              ),
              na.rm = TRUE
            )
          } else {
            sum_occupancy[[layer]] <-
              sum_occupancy[[layer]] +
              x[[layer]]
          }
        }
      }
    }

    # Keep only the four prediction layers and geometry
    if (inherits(sum_occupancy, "sf")) {
      sum_occupancy <- sum_occupancy[
        c(required_layers, attr(sum_occupancy, "sf_column"))
      ]
    } else {
      sum_occupancy <- sum_occupancy[
        required_layers
      ]
    }
  }

  # ============================================================
  # 9. Create output object
  # ============================================================

  output <- list(
    results = results,
    expected_richness = sum_occupancy
  )

  # ============================================================
  # 10. Save outputs
  # ============================================================

  if (save_outputs) {
    # ----------------------------------------------------------
    # Determine output filename/format
    # ----------------------------------------------------------

    if (inherits(sum_occupancy, "SpatRaster")) {
      output_file <- file.path(
        output_dir,
        "expected_richness.tif"
      )

      terra::writeRaster(
        sum_occupancy,
        filename = output_file,
        overwrite = overwrite
      )

      output_format <- "GeoTIFF"
    } else if (inherits(sum_occupancy, "sf")) {
      output_file <- file.path(
        output_dir,
        "expected_richness.parquet"
      )

      sfarrow::st_write_parquet(
        obj = sum_occupancy,
        dsn = output_file,
        compression = "snappy"
      )

      output_format <- "GeoParquet"
    } else {
      output_file <- file.path(
        output_dir,
        "expected_richness.rds"
      )

      saveRDS(
        sum_occupancy,
        output_file
      )

      output_format <- "RDS"
    }

    # ----------------------------------------------------------
    # Metadata
    # ----------------------------------------------------------

    metadata <- list(
      title = "Multispecies ecosystem richness from occupancy predictions",

      description = paste(
        "Species-specific occupancy predictions were aggregated",
        "across species by summing mean, median, lower, and upper",
        "prediction values at each prediction location."
      ),

      created = as.character(Sys.time()),

      ecosystem_type = ecosystem_type,

      method = method,

      input_dir = normalizePath(
        input_dir,
        mustWork = FALSE
      ),

      output_dir = normalizePath(
        output_dir,
        mustWork = FALSE
      ),

      n_species = length(results),

      species = names(results),

      prediction_layers = required_layers,

      aggregation = list(
        method = "sum across species",
        mean = "sum of species-specific mean predictions",
        median = "sum of species-specific median predictions",
        lower = "sum of species-specific lower predictions",
        upper = "sum of species-specific upper predictions"
      ),

      output = list(
        filename = basename(output_file),
        format = output_format,
        layers = required_layers
      ),

      uncertainty_note = paste(
        "The lower and upper layers are sums of species-specific",
        "marginal uncertainty bounds. They are not equivalent to",
        "a posterior credible interval for the multispecies sum",
        "unless posterior samples are aggregated before calculating",
        "the interval."
      ),

      software_versions = list(
        R = R.version.string,
        sf = if (requireNamespace("sf", quietly = TRUE)) {
          as.character(utils::packageVersion("sf"))
        },
        terra = if (requireNamespace("terra", quietly = TRUE)) {
          as.character(utils::packageVersion("terra"))
        },
        sfarrow = if (requireNamespace("sfarrow", quietly = TRUE)) {
          as.character(utils::packageVersion("sfarrow"))
        },
        jsonlite = if (requireNamespace("jsonlite", quietly = TRUE)) {
          as.character(utils::packageVersion("jsonlite"))
        }
      )
    )

    jsonlite::write_json(
      metadata,
      path = file.path(
        output_dir,
        "metadata.json"
      ),
      pretty = TRUE,
      auto_unbox = TRUE
    )

    saveRDS(
      output,
      file = file.path(
        output_dir,
        "all_results.rds"
      )
    )
  }

  # ============================================================
  # 11. Return
  # ============================================================

  output
}
