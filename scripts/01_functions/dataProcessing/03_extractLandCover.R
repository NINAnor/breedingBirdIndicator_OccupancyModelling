#' Extract CORINE land-cover labels for breeding-bird route points
#'
#' This function converts breeding-bird observations with longitude/latitude
#' coordinates into an `sf` object, projects them to the model CRS, extracts
#' CORINE land-cover raster values at each point, joins those raster codes to
#' the CORINE land-cover lookup table, and returns selected identifiers and
#' land-cover labels.
#'
#' FAIR notes:
#' - Findable: input and covariate file paths are constructed explicitly via
#'   `path_finder("dataFolder")`.
#' - Accessible: required external files are documented in the function body.
#' - Interoperable: spatial data are handled using standard `sf` and `terra`
#'   classes with explicit CRS definitions.
#' - Reusable: output columns are stable and limited to documented identifiers
#'   and land-cover label fields.
#'
#' @param ecosystemBreedingBird A data frame containing at least
#'   `Longitude`, `Latitude`, `Route_point_id`, and `Year`.
#' @param modelCrs Target coordinate reference system as an EPSG code,
#'   for example `25833`.
#'
#' @return A data frame with columns:
#'   `Route_point_id`, `Year`, `LABEL1`, `LABEL2`, and `LABEL`.
#'
#' @examples
#' landcover_labels <- get_landcover_labels(
#'   ecosystemBreedingBird = ecosystemBreedingBird,
#'   modelCrs = 25833
#' )
#'
get_landcover_labels <- function(ecosystemBreedingBird,
                                 modelCrs,
                                 cov_folder,
                                 output = "data.frame") {

  # Keep only records with valid coordinates.
  # FAIR: improves data quality and makes spatial conversion reproducible.
  ecosystemBreedingBird_sf <- ecosystemBreedingBird |>
    dplyr::filter(!is.na(Longitude), !is.na(Latitude)) |>
    dplyr::select(Longitude, Latitude, Year, Route_point_id) |>
    dplyr::group_by(Route_point_id, Year) |>
    dplyr::distinct() |>
    sf::st_as_sf(
      coords = c("Longitude", "Latitude"),
      crs = 4326,        # WGS84 input coordinates
      remove = FALSE
    ) |>
    sf::st_transform(modelCrs)


  # Load and project CORINE land-cover raster.
  # FAIR: CRS is explicitly harmonised with the model CRS.
  landCoverRaster <- terra::rast(
    file.path(cov_folder, "CorineLandCover2018_1km.tif")
  ) |>
    terra::project(
      paste0("EPSG:", modelCrs),
      method = "near"    # nearest neighbour for categorical raster data
    )

  # Extract raster values at point locations.
  vals <- terra::extract(
    landCoverRaster,
    terra::vect(ecosystemBreedingBird_sf)
  )

  # Load CORINE code lookup table.
  # FAIR: maps numeric raster codes to interpretable land-cover labels.
  Corine_Land_Cover_codes <- readxl::read_excel(
    file.path(cov_folder, "Corine_Land_Cover_codes.xlsx")
  )

  # Join extracted raster codes to land-cover labels and return stable columns.
  if(output %in% "data.frame"){
  result <- cbind(ecosystemBreedingBird_sf, vals) |>
    dplyr::left_join(
      Corine_Land_Cover_codes,
      by = c("clc18_kode" = "CLC_CODE")
    ) |>
    sf::st_drop_geometry() |>
    dplyr::select(
      Route_point_id,
      Year,
      LABEL1,
      LABEL2,
      LABEL3
    )
  } else if(output %in% "raster"){
    result <- landCoverRaster |>
      tidyterra::mutate(LABEL1 = Corine_Land_Cover_codes$LABEL1[
        match(clc18_kode, Corine_Land_Cover_codes$CLC_CODE)
      ],
      LABEL2 = Corine_Land_Cover_codes$LABEL2[
        match(clc18_kode, Corine_Land_Cover_codes$CLC_CODE)
      ],
      LABEL3 = Corine_Land_Cover_codes$LABEL3[
        match(clc18_kode, Corine_Land_Cover_codes$CLC_CODE)
      ])|>
      tidyterra::select(
        clc18_kode,
        LABEL1,
        LABEL2,
        LABEL3
      )
  }

  return(result)
}

