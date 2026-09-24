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
  overwrite = TRUE
) {
  if (!dir.exists(input_dir)) {
    stop("`input_dir` does not exist: ", input_dir, call. = FALSE)
  }

  if (save_outputs && is.null(output_dir)) {
    stop(
      "`output_dir` must be supplied when `save_outputs = TRUE`.",
      call. = FALSE
    )
  }

  if (save_outputs) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  species_dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)

  results <- setNames(
    lapply(species_dirs, function(sp_dir) {
      files <- list.files(sp_dir, full.names = TRUE)

      tif_files <- files[grepl("\\.tif$", files)]
      parquet_files <- files[grepl("\\.parquet$", files)]

      if (length(tif_files) > 0) {
        terra::rast(tif_files[1])
      } else if (length(parquet_files) > 0) {
        sfarrow::st_read_parquet(parquet_files[1])
      } else {
        NULL
      }
    }),
    basename(species_dirs)
  )

  results <- Filter(Negate(is.null), results)

  if (length(results) == 0) {
    stop("No valid .tif or .parquet species outputs were found.", call. = FALSE)
  }

  species_cols <- safe_name(names(results))

  # Extract mean, lower, and upper rasters for all species
  mean_rasters <- lapply(results, function(x) x[["mean"]])
  median_rasters <- lapply(results, function(x) x[["median"]])
  lower_rasters <- lapply(results, function(x) x[["lower"]])
  upper_rasters <- lapply(results, function(x) x[["upper"]])

  # Sum across all species
  sum_mean <- Reduce(`+`, mean_rasters)
  sum_median <- Reduce(`+`, median_rasters)
  sum_lower <- Reduce(`+`, lower_rasters)
  sum_upper <- Reduce(`+`, upper_rasters)

  # Combine into one SpatRaster
  sum_occupancy <- c(
    sum_mean,
    sum_median,
    sum_lower,
    sum_upper
  )

  # Give meaningful names
  names(sum_occupancy) <- c(
    "mean",
    "median",
    "lower",
    "upper"
  )

  output <- list(
    results = results,
    expected_richness = sum_occupancy
  )

  if (save_outputs) {
    save_object <- function(x, filename_base, dir = output_dir) {
      if (inherits(x, "sf")) {
        sfarrow::st_write_parquet(
          obj = x,
          dsn = file.path(dir, paste0(filename_base, ".parquet")),
          compression = "snappy"
        )
      } else if (inherits(x, "SpatRaster")) {
        terra::writeRaster(
          x,
          filename = file.path(dir, paste0(filename_base, ".tif")),
          overwrite = overwrite
        )
      } else if (inherits(x, "data.frame")) {
        utils::write.csv(
          x,
          file = file.path(dir, paste0(filename_base, ".csv")),
          row.names = FALSE
        )
      } else {
        saveRDS(
          x,
          file = file.path(dir, paste0(filename_base, ".rds"))
        )
      }
    }

    save_object(
      sum_occupancy,
      "expected_richness"
    )

    metadata <- list(
      title = "Multispecies ecosystem richness from occupancy predictions",
      description = paste(
        "Cell-wise aggregation of species-specific occupancy predictions.",
        "For each species, mean, median, lower, and upper occupancy",
        "prediction layers are summed across species to produce a",
        "multispecies expected richness surface."
      ),
      created = as.character(Sys.time()),

      ecosystem_type = ecosystem_type,

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
      species_columns = safe_name(names(results)),

      prediction_layers = c(
        "mean",
        "median",
        "lower",
        "upper"
      ),

      aggregation = list(
        method = "cell-wise sum across species",
        mean = "sum of species-specific mean occupancy predictions",
        median = "sum of species-specific median occupancy predictions",
        lower = "sum of species-specific lower occupancy predictions",
        upper = "sum of species-specific upper occupancy predictions"
      ),

      output = list(
        filename = "expected_richness.tif",
        format = "GeoTIFF",
        layers = c(
          "mean",
          "median",
          "lower",
          "upper"
        ),
        interpretation = list(
          mean = "Expected multispecies richness",
          median = "Sum of species-specific median occupancy predictions",
          lower = "Sum of species-specific lower occupancy predictions",
          upper = "Sum of species-specific upper occupancy predictions"
        )
      ),

      uncertainty_note = paste(
        "Lower and upper layers are sums of species-specific marginal",
        "uncertainty bounds and are not equivalent to posterior credible",
        "intervals for the multispecies sum unless posterior samples are",
        "aggregated before calculating quantiles."
      ),

      file_format_rules = list(
        SpatRaster = "GeoTIFF",
        other = "RDS"
      ),

      FAIR_notes = list(
        findable = paste(
          "Outputs use descriptive filenames and metadata containing",
          "species names and processing information."
        ),
        accessible = paste(
          "Raster outputs are provided as GeoTIFF and complete results",
          "are serialized as RDS."
        ),
        interoperable = paste(
          "Raster outputs preserve the coordinate reference system,",
          "resolution, extent, and standard GeoTIFF structure."
        ),
        reusable = paste(
          "Metadata records input locations, species, prediction layers,",
          "aggregation method, parameters, software versions, and",
          "uncertainty interpretation."
        )
      ),

      software_versions = list(
        R = R.version.string,
        sf = as.character(utils::packageVersion("sf")),
        terra = as.character(utils::packageVersion("terra")),
        dplyr = as.character(utils::packageVersion("dplyr")),
        sfarrow = as.character(utils::packageVersion("sfarrow")),
        jsonlite = as.character(utils::packageVersion("jsonlite"))
      )
    )

    jsonlite::write_json(
      metadata,
      path = file.path(output_dir, "metadata.json"),
      pretty = TRUE,
      auto_unbox = TRUE
    )

    saveRDS(
      output,
      file = file.path(output_dir, "all_results.rds")
    )
  }

  output
}
