#' Read species occupancy predictions (value column and geometry only)
#'
#' Reads one prediction file per species folder, for the species in
#' `reference_data`. Only `value_col` and the geometry are kept, which keeps
#' memory use low.
#'
#' @param input_dir Character. Folder with one sub-folder per species, or the
#'   species folders themselves (more than one path).
#' @param reference_data Data frame with a `Species` column. Only these
#'   species are read.
#' @param value_col Character. Occupancy column to keep.
#' @param method One of `"raster"` (`.tif`), `"sf"` (`.parquet`), `"rds"`.
#'
#' @return Named list (names are the species folder names) of `sf` objects
#'   or `SpatRaster` layers.
#' @export
read_species_predictions <- function(
  input_dir,
  reference_data,
  value_col = "occupancy_mean",
  method = c("raster", "sf", "rds")
) {
  method <- match.arg(method)

  if (length(input_dir) == 0) {
    stop("`input_dir` contains no paths.", call. = FALSE)
  }

  missing_dirs <- input_dir[!dir.exists(input_dir)]
  if (length(missing_dirs) > 0) {
    warning(
      "The following input directories do not exist:\n",
      paste(missing_dirs, collapse = "\n"),
      call. = FALSE
    )
  }
  input_dir <- input_dir[dir.exists(input_dir)]
  if (length(input_dir) == 0) {
    stop("None of the `input_dir` paths exist.", call. = FALSE)
  }

  # OneDrive paths are too long for some R/Arrow/Windows file operations
  short_path <- function(path) {
    if (.Platform$OS.type != "windows") {
      return(path)
    }
    result <- trimws(system(
      sprintf('cmd /c for %%I in ("%s") do @echo %%~sI', path),
      intern = TRUE
    ))
    if (length(result) != 1 || !nzchar(result)) {
      stop("Could not convert path to Windows short path:\n", path, call. = FALSE)
    }
    result
  }

  # Species folders (only species with a reference)
  species_dirs <- if (length(input_dir) == 1) {
    list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
  } else {
    input_dir
  }
  species_dirs <- species_dirs[!grepl("^_", basename(species_dirs))]

  dir_species <- gsub("_", " ", basename(species_dirs), fixed = TRUE)

  no_dir <- setdiff(reference_data$Species, dir_species)
  if (length(no_dir)) {
    warning(
      "No prediction folder for reference species: ",
      paste(no_dir, collapse = ", "),
      call. = FALSE
    )
  }

  species_dirs <- species_dirs[dir_species %in% reference_data$Species]
  if (length(species_dirs) == 0) {
    stop("No species folders match `reference_data`.", call. = FALSE)
  }

  pattern <- switch(method, raster = "\\.tif$", sf = "\\.parquet$", rds = "\\.rds$")

  read_prediction <- function(sp_dir) {
    files <- list.files(sp_dir, full.names = TRUE, recursive = FALSE)
    files <- files[grepl(pattern, files, ignore.case = TRUE)]
    if (length(files) == 0) {
      stop("No ", method, " prediction file found in:\n", sp_dir, call. = FALSE)
    }
    file_short <- short_path(files[1])

    x <- switch(
      method,
      raster = terra::rast(file_short),
      sf = sfarrow::st_read_parquet(file_short),
      rds = readRDS(file_short)
    )

    if (!value_col %in% names(x)) {
      stop("`", value_col, "` not found for ", basename(sp_dir), call. = FALSE)
    }
    if (inherits(x, "SpatRaster")) x[[value_col]] else x[value_col]
  }

  stats::setNames(
    lapply(species_dirs, read_prediction),
    basename(species_dirs)
  )
}


#' Scale species occupancy by region-specific reference values
#'
#' Reads the species-level occupancy predictions for the species in
#' `reference_data` (or uses `results`), scales them with
#' [scale_species_occupancy_by_region()], and averages the scaled values
#' across species.
#'
#' @param input_dir Character. Folder with one sub-folder per species, or the
#'   species folders themselves. Not needed when `results` is supplied.
#' @param reference_data Data frame with `Species` and `reference_quantile`
#'   (see [get_reference_quantiles()]). Only these species are used.
#' @param categorisation_shape `sf` polygons defining the regions.
#' @param region_col Character. Region name column in `categorisation_shape`.
#' @param ecosystem_type Character. Ecosystem type (informational).
#' @param value_col Character. Occupancy column to scale.
#' @param truncation Logical. Passed to [ecTools::ec_normalise()].
#' @param na.rm Logical. Ignore missing species values when computing the mean
#'   across species.
#' @param method One of `"raster"` (`.tif`), `"sf"` (`.parquet`), `"rds"`.
#' @param save_outputs Logical. Save results to `output_dir`.
#' @param output_dir Character. Required when `save_outputs = TRUE`.
#' @param overwrite Logical. Overwrite existing outputs.
#' @param results Optional list from [read_species_predictions()]. Supplying
#'   it avoids reading the files again, e.g. when scaling the same predictions
#'   with several reference quantiles.
#' @param ... Further arguments for [ecTools::ec_normalise()].
#'
#' @return A list with
#'   \describe{
#'     \item{scaled}{`sf` with a scaled column per species, `mean_scaled`
#'       (mean across species), `bio_region` and geometry.}
#'     \item{reference_values}{Reference value (`x100`) per region and
#'       species.}
#'     \item{mean_raw}{Unscaled mean occupancy across species.}
#'   }
#' @export
species_specific_scaling <- function(
  input_dir = NULL,
  reference_data,
  categorisation_shape,
  region_col,
  ecosystem_type = NA_character_,
  value_col = "occupancy_mean",
  truncation = TRUE,
  na.rm = TRUE,
  method = c("raster", "sf", "rds"),
  save_outputs = FALSE,
  output_dir = NULL,
  overwrite = FALSE,
  results = NULL,
  ...
) {
  method <- match.arg(method)

  if (save_outputs && is.null(output_dir)) {
    stop("`output_dir` must be supplied when `save_outputs = TRUE`.", call. = FALSE)
  }

  # ---- Read predictions (unless supplied) -------------------------------------
  if (is.null(results)) {
    results <- read_species_predictions(
      input_dir = input_dir,
      reference_data = reference_data,
      value_col = value_col,
      method = method
    )
  }

  # ---- Scale by region-specific reference -------------------------------------
  output <- scale_species_occupancy_by_region(
    results = results,
    reference_data = reference_data,
    categorisation_shape = categorisation_shape,
    region_col = region_col,
    value_col = value_col,
    truncation = truncation,
    ...
  )
  rm(results)

  # ---- Mean across species ----------------------------------------------------
  species_cols <- setdiff(names(output$scaled), c("bio_region", "geometry"))

  output$scaled$mean_scaled <- rowMeans(
    sf::st_drop_geometry(output$scaled)[species_cols],
    na.rm = na.rm
  )
  # rowMeans() gives NaN when every species is missing
  output$scaled$mean_scaled[is.nan(output$scaled$mean_scaled)] <- NA_real_

  output$scaled <- output$scaled[
    c(species_cols, "mean_scaled", "bio_region", "geometry")
  ]

  # ---- Save -------------------------------------------------------------------
  if (save_outputs) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    scaled_file <- file.path(output_dir, "scaled_occupancy.parquet")
    ref_file <- file.path(output_dir, "reference_values.csv")
    if (!overwrite && any(file.exists(c(scaled_file, ref_file)))) {
      stop(
        "Outputs already exist in `output_dir`; use `overwrite = TRUE`.",
        call. = FALSE
      )
    }

    sfarrow::st_write_parquet(output$scaled, scaled_file, compression = "snappy")
    utils::write.csv(output$reference_values, ref_file, row.names = FALSE)
    saveRDS(output, file.path(output_dir, "all_results.rds"))
  }

  output
}
