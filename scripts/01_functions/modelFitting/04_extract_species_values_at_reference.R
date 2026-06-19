#' Extract species prediction values at reference polygons
#'
#' @description
#' Extracts species-specific occupancy prediction summaries from a named list
#' of prediction outputs and appends them as columns to a reference `sf`
#' object. Prediction outputs may be provided either as `sf` objects or as
#' `terra::SpatRaster` objects. One column is created per species, using
#' sanitized species names as column names.
#'
#' @details
#' This function supports FAIR (Findable, Accessible, Interoperable, Reusable)
#' data principles.
#'
#' **Findable**
#'
#' Species predictions are stored in consistently named columns derived from
#' species names. A metadata attribute containing extraction information and
#' software versions is attached to the returned object.
#'
#' **Accessible**
#'
#' The returned object is an `sf` object and can be written to open formats
#' such as GeoPackage or GeoParquet.
#'
#' **Interoperable**
#'
#' Coordinate reference systems are checked and transformed where necessary
#' before spatial extraction. The function supports both vector (`sf`) and
#' raster (`terra`) prediction outputs.
#'
#' **Reusable**
#'
#' Metadata describing extraction settings, species names, summary functions,
#' and software versions are stored as an attribute named `"FAIR_metadata"`.
#'
#' For `sf` prediction objects, values are transferred using
#' `sf::st_join()` and summarized over intersecting features.
#'
#' For raster predictions, values are extracted using `terra::extract()`
#' and summarized across the reference polygons.
#'
#' @param allResults Named list of species prediction objects. Each element
#'   must be either an `sf` object or a `terra::SpatRaster`. List names are
#'   interpreted as species names.
#' @param reference_shape_with_regions An `sf` object representing the
#'   reference polygons to which species prediction values will be appended.
#' @param value_col Character string giving the prediction layer or column
#'   to extract. Defaults to `"mean"`.
#' @param summarise_fun Function used to summarize multiple values intersecting
#'   the same reference polygon. Defaults to `mean`.
#' @param na.rm Logical. Should missing values be removed during
#'   summarization? Defaults to `TRUE`.
#' @param reference_id_col Character string identifying the column used
#'   internally to uniquely identify reference polygons. If the column is
#'   absent, a unique identifier is created automatically. Defaults to
#'   `".reference_id"`.
#' @param keep_fair_metadata Logical. Should FAIR metadata describing the
#'   extraction process be attached to the returned object? Defaults to `TRUE`.
#'
#' @return
#' An `sf` object containing all original columns and geometry from
#' `reference_shape_with_regions`, with one additional numeric column per
#' species.
#'
#' If `keep_fair_metadata = TRUE`, a `"FAIR_metadata"` attribute is attached
#' containing:
#'
#' - extraction date and time;
#' - extracted prediction summary;
#' - summary function used;
#' - species names;
#' - extraction methods;
#' - software versions;
#' - FAIR documentation notes.
#'
#' @section Supported prediction types:
#'
#' - `sf`: values are joined spatially using `sf::st_intersects()`.
#' - `terra::SpatRaster`: polygon values are extracted using
#'   `terra::extract()`.
#'
#' @section FAIR recommendations:
#'
#' For improved reproducibility and reuse, it is recommended to export the
#' returned object to GeoPackage or GeoParquet and preserve the associated
#' `"FAIR_metadata"` attribute in a companion metadata file.
#'
#' Metadata files should document:
#'
#' - prediction source and version;
#' - species names;
#' - coordinate reference systems;
#' - extracted summary statistic;
#' - software versions;
#' - date of extraction;
#' - citations or identifiers for source data.
#'
#' @examples
#' \dontrun{
#' species_values <- extract_species_values_at_reference(
#'   allResults = allResults,
#'   reference_shape_with_regions = protected_regions,
#'   value_col = "mean"
#' )
#'
#' attr(species_values, "FAIR_metadata")
#'
#' sf::write_sf(
#'   species_values,
#'   "species_predictions.gpkg"
#' )
#' }
#'
#' @export extract_species_values_at_reference

extract_species_values_at_reference <- function(
    allResults,
    reference_shape_with_regions,
    value_col = "mean",
    point_spacing,
    region_col = "bioGeographic_region",
    summarise_fun = mean,
    na.rm = TRUE,
    keep_fair_metadata = TRUE
) {
  stopifnot(inherits(reference_shape_with_regions, "sf"))

  if (!is.list(allResults) || is.null(names(allResults))) {
    stop("`allResults` must be a named list of sf or terra::SpatRaster objects.", call. = FALSE)
  }

  if (missing(point_spacing)) {
    stop("Provide `point_spacing`, in the map units of `reference_shape_with_regions`.", call. = FALSE)
  }

  if (!region_col %in% names(reference_shape_with_regions)) {
    stop("`region_col` not found in `reference_shape_with_regions`: ", region_col, call. = FALSE)
  }

  reference_shape_with_regions <- sf::st_make_valid(reference_shape_with_regions)

  reference_shape_with_regions$.reference_id <- seq_len(nrow(reference_shape_with_regions))

  polygons <- reference_shape_with_regions[
    sf::st_geometry_type(reference_shape_with_regions) %in% c("POLYGON", "MULTIPOLYGON"),
  ]

  if (nrow(polygons) == 0) {
    stop("`reference_shape_with_regions` must contain POLYGON or MULTIPOLYGON geometries.", call. = FALSE)
  }

  cat("Generating regular points inside reference polygons...\n")

  point_geom <- sf::st_make_grid(
    polygons,
    cellsize = point_spacing,
    what = "centers",
    square = TRUE
  )

  points_sf <- sf::st_sf(
    geometry = point_geom,
    crs = sf::st_crs(polygons)
  )

  points_sf <- sf::st_join(
    points_sf,
    polygons[, c(".reference_id", region_col), drop = FALSE],
    join = sf::st_within,
    left = FALSE
  )

  points_sf$.point_id <- seq_len(nrow(points_sf))

  extraction_log <- list()

  for (sp in names(allResults)) {

    cat("Extracting:", sp, "\n")

    x <- allResults[[sp]]
    sp_safe <- safe_name(sp)

    extraction_log[[sp_safe]] <- list(
      species_original_name = sp,
      species_column_name = sp_safe,
      source_object_class = paste(class(x), collapse = ", "),
      extracted_value = value_col,
      status = "started",
      extraction_method = NA_character_
    )

    if (inherits(x, "sf")) {

      if (!value_col %in% names(x)) {
        warning("Skipping ", sp, ": column not found: ", value_col)
        extraction_log[[sp_safe]]$status <- "skipped_column_not_found"
        next
      }

      if (sf::st_crs(x) != sf::st_crs(points_sf)) {
        x <- sf::st_transform(x, sf::st_crs(points_sf))
      }

      nearest_idx <- sf::st_nearest_feature(points_sf, x)

      points_sf[[sp_safe]] <- x[[value_col]][nearest_idx]

      extraction_log[[sp_safe]]$status <- "completed"
      extraction_log[[sp_safe]]$extraction_method <- "sf::st_nearest_feature from regular reference points"
    }

    else if (inherits(x, "SpatRaster")) {

      if (!value_col %in% names(x)) {
        warning("Skipping ", sp, ": raster layer not found: ", value_col)
        extraction_log[[sp_safe]]$status <- "skipped_layer_not_found"
        next
      }

      points_vect <- terra::vect(points_sf)

      if (!terra::same.crs(x, points_vect)) {
        points_vect <- terra::project(points_vect, terra::crs(x))
      }

      vals <- terra::extract(
        x[[value_col]],
        points_vect
      )

      points_sf[[sp_safe]] <- vals[[value_col]]

      extraction_log[[sp_safe]]$status <- "completed"
      extraction_log[[sp_safe]]$extraction_method <- "terra::extract at regular reference points"
    }

    else {
      warning("Skipping ", sp, ": unsupported result type.")
      extraction_log[[sp_safe]]$status <- "skipped_unsupported_type"
    }
  }

  if (keep_fair_metadata) {
    attr(points_sf, "FAIR_metadata") <- list(
      title = "Species occupancy prediction values extracted at regular points within reference polygons",
      created = as.character(Sys.time()),
      value_column_extracted = value_col,
      point_spacing = point_spacing,
      region_column = region_col,
      species_columns = safe_name(names(allResults)),
      extraction_log = extraction_log,
      software_versions = list(
        R = R.version.string,
        sf = as.character(utils::packageVersion("sf")),
        terra = as.character(utils::packageVersion("terra")),
        dplyr = as.character(utils::packageVersion("dplyr"))
      ),
      FAIR_notes = list(
        findable = "Each species is stored in a clearly named column using a sanitized species name.",
        accessible = "The returned object is an sf point object and can be written to GeoPackage or GeoParquet.",
        interoperable = "CRS is preserved and transformed before extraction where required.",
        reusable = "Metadata records point spacing, extraction value, software versions, region column, and extraction log."
      )
    )
  }

  points_sf
}
