#' Create a Human-Readable Label from a Model Output Path
#'
#' Extracts model configuration information from a model output directory path
#' and converts it into a concise label describing the model type.
#'
#' The function assumes the directory structure follows the convention:
#'
#' \preformatted{
#' .../<species_model>/<spatial_model>
#' }
#'
#' where:
#' \itemize{
#'   \item \code{single_species_model} or \code{multiple_species_model}
#'   \item \code{spatialModel} or \code{non_spatialModel}
#' }
#'
#' These components are converted to a compact label suitable for plot legends
#' and summaries.
#'
#' Examples:
#'
#' \preformatted{
#' single-species;spatial
#' single-species;non-spatial
#' multi-species;spatial
#' multi-species;non-spatial
#' }
#'
#' @param path Character string containing the full path to a model output
#'   directory.
#'
#' @return A character string describing the model configuration.
#'
#' @examples
#' extract_run_label(
#'   "modelFit/version_01/Forest/single_species_model/spatialModel"
#' )
#'
#' extract_run_label(
#'   "modelFit/version_01/Forest/multiple_species_model/non_spatialModel"
#' )
#'
#' @export
extract_run_label <- function(path) {
  parts <- strsplit(normalizePath(path, winslash = "/"), "/")[[1]]

  species <- tail(parts, 2)[1]
  spatial <- tail(parts, 1)
  version_no <- tail(parts, 4)[1]

  species_short <- if (species == "multiple_species_model") {
    "multi-speices"
  } else {
    "single-species"
  }

  spatial_short <- if (spatial == "spatialModel") {
    "spatial"
  } else {
    "non-spatial"
  }

  paste(species_short, spatial_short, version_no, sep = ";")
}


#' Build Model Output Paths for One or More Model Configurations
#'
#' Constructs a table of expected model output directories based on combinations
#' of species-model and spatial-model settings.
#'
#' This function is useful when comparing multiple model configurations or
#' programmatically locating model output folders for downstream extraction,
#' summarisation, and plotting.
#'
#' The generated paths follow the directory structure:
#'
#' \preformatted{
#' <output_root>/modelFit/<version_label>/<ecosystemType>/
#'   <species_model>/<spatial_model>
#' }
#'
#' where:
#' \itemize{
#'   \item \code{single_species_model} or \code{multiple_species_model}
#'   \item \code{spatialModel} or \code{non_spatialModel}
#' }
#'
#' Optionally, non-existing directories can be removed from the output.
#'
#' @param output_root Character. Root project output directory.
#' @param version_label Character. Model version identifier.
#' @param ecosystemType Character. Ecosystem or analysis type used in the model.
#' @param include_single_species Logical vector indicating whether to include
#'   single-species model configurations. Defaults to both \code{TRUE} and
#'   \code{FALSE}.
#' @param include_spatial_model Logical vector indicating whether to include
#'   spatial model configurations. Defaults to both \code{TRUE} and
#'   \code{FALSE}.
#' @param existing_only Logical. If \code{TRUE}, only paths that currently
#'   exist on disk are returned.
#'
#' @return A tibble containing:
#' \describe{
#'   \item{fit_single_species}{Logical indicating model type.}
#'   \item{fit_spatial_model}{Logical indicating spatial configuration.}
#'   \item{species_dir}{Species-model directory name.}
#'   \item{spatial_dir}{Spatial-model directory name.}
#'   \item{path_to_results}{Full path to the model output directory.}
#' }
#'
#' @examples
#' paths_df <- build_model_paths(
#'   output_root = "outputs",
#'   version_label = "version_01",
#'   ecosystemType = "Forest"
#' )
#'
#' paths_df
#'
#' existing_paths <- build_model_paths(
#'   output_root = "outputs",
#'   version_label = "version_01",
#'   ecosystemType = "Forest",
#'   existing_only = TRUE
#' )
#'
#' @importFrom dplyr mutate filter if_else
#' @importFrom tidyr crossing
#'
#' @export
build_model_paths <- function(
  output_root,
  version_labels,
  ecosystemType,
  include_single_species = c(TRUE, FALSE),
  include_spatial_model = c(TRUE, FALSE),
  existing_only = FALSE
) {
  df <- tidyr::crossing(
    fit_single_species = include_single_species,
    fit_spatial_model = include_spatial_model,
    version_label = version_labels
  ) %>%
    dplyr::mutate(
      species_dir = dplyr::if_else(
        fit_single_species,
        "single_species_model",
        "multiple_species_model"
      ),
      spatial_dir = dplyr::if_else(
        fit_spatial_model,
        "spatialModel",
        "non_spatialModel"
      ),
      path_to_results = file.path(
        output_root,
        "modelFit",
        version_label,
        ecosystemType,
        species_dir,
        spatial_dir
      )
    )

  df <- df %>%
    dplyr::filter(!is.na(path_to_results), nzchar(path_to_results))

  if (existing_only) {
    df <- df %>% dplyr::filter(file.exists(path_to_results))
  }

  df
}
