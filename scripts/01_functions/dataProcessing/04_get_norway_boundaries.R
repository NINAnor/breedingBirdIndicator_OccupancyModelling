#' Retrieve and Cache Norway Administrative Boundaries
#'
#' Downloads Norway administrative boundaries from GADM using the
#' `geodata` package and caches them locally as GeoParquet files.
#'
#' If cached files already exist, they are read from disk unless
#' `overwrite = TRUE`.
#'
#' The function follows FAIR principles:
#' \itemize{
#'   \item \strong{Findable}: data and metadata are stored with descriptive names.
#'   \item \strong{Accessible}: boundaries are stored in open GeoParquet format.
#'   \item \strong{Interoperable}: boundaries are returned as `sf` objects with CRS.
#'   \item \strong{Reusable}: provenance, CRS, source, and retrieval date are
#'         recorded in a YAML metadata file.
#' }
#'
#' @param out_dir Character. Directory for storing GeoParquet files and metadata.
#'   Defaults to `"norway_boundaries"`.
#' @param temp_dir Character or NULL. Temporary download directory. If NULL,
#'   a temporary directory is created automatically.
#' @param source Character. Source label. Defaults to `"GADM"`.
#' @param overwrite Logical. If TRUE, forces re-download and overwrites cached
#'   files. Defaults to FALSE.
#'
#' @details
#' Administrative levels:
#' \describe{
#'   \item{Level 0}{National boundary of Norway}
#'   \item{Level 1}{County (fylke) boundaries}
#' }
#'
#' Metadata are stored in YAML format to improve machine readability and
#' reproducibility.
#'
#' @return A named list containing:
#' \describe{
#'   \item{national_boundary}{An sf object of Norway's national boundary.}
#'   \item{county_boundaries}{An sf object of county boundaries.}
#'   \item{files}{Paths to the GeoParquet and metadata files.}
#' }
#'
#' @references
#' Hijmans, R.J. (2025). geodata: Download Geographic Data.
#' https://rspatial.org/geodata/
#'
#' GADM Database of Global Administrative Areas.
#' https://gadm.org/
#'
#' @examples
#' \dontrun{
#' norway <- get_norway_boundaries()
#'
#' plot(sf::st_geometry(norway$national_boundary))
#' plot(sf::st_geometry(norway$county_boundaries))
#'
#' norway <- get_norway_boundaries(overwrite = TRUE)
#' }
#'
#' @export
get_norway_boundaries <- function(
    out_dir = "norway_boundaries",
    temp_dir = NULL,
    source = "GADM",
    overwrite = FALSE
) {

  required_pkgs <- c(
    "sf",
    "geodata",
    "terra",
    "sfarrow",
    "yaml"
  )

  missing_pkgs <- required_pkgs[
    !vapply(
      required_pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing_pkgs) > 0) {
    stop(
      "Install missing packages first: ",
      paste(missing_pkgs, collapse = ", "),
      call. = FALSE
    )
  }

  # Create output directory
  out_dir <- normalizePath(
    out_dir,
    winslash = "/",
    mustWork = FALSE
  )

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Create temporary directory if needed
  if (is.null(temp_dir)) {
    temp_dir <- tempfile("norway_boundaries_")
  }

  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

  # FAIR-compliant file names
  national_file <- normalizePath(
    file.path(out_dir, "norway_boundary_level0.parquet"),
    winslash = "/",
    mustWork = FALSE
  )

  county_file <- normalizePath(
    file.path(out_dir, "norway_boundary_level1_counties.parquet"),
    winslash = "/",
    mustWork = FALSE
  )

  metadata_file <- normalizePath(
    file.path(out_dir, "norway_boundaries_metadata.yml"),
    winslash = "/",
    mustWork = FALSE
  )

  files_exist <- all(vapply(
    c(national_file, county_file, metadata_file),
    file.exists,
    logical(1)
  ))

  if (files_exist && !overwrite) {

    message("Reading cached GeoParquet files.")

    national <- sfarrow::st_read_parquet(
      national_file
    )

    counties <- sfarrow::st_read_parquet(
      county_file
    )

  } else {

    message("Downloading Norway boundaries from ", source)

    # Download data
    national_raw <- geodata::gadm(
      country = "NOR",
      level = 0,
      path = temp_dir
    )

    counties_raw <- geodata::gadm(
      country = "NOR",
      level = 1,
      path = temp_dir
    )

    # Convert to sf
    national <- sf::st_make_valid(
      sf::st_as_sf(national_raw)
    )

    counties <- sf::st_make_valid(
      sf::st_as_sf(counties_raw)
    )

    # Save as GeoParquet
    sfarrow::st_write_parquet(
      national,
      national_file
    )

    sfarrow::st_write_parquet(
      counties,
      county_file
    )

    # FAIR metadata
    metadata <- list(
      title = "Norway administrative boundaries",
      description = "National and county boundaries for Norway",
      country = "Norway",
      country_iso3 = "NOR",
      source = source,
      source_url = "https://gadm.org/",
      retrieval_date = as.character(Sys.time()),
      format = "GeoParquet",
      administrative_levels = list(
        level0 = "National boundary",
        level1 = "County boundaries"
      ),
      coordinate_reference_system = list(
        level0 = sf::st_crs(national)$input,
        level1 = sf::st_crs(counties)$input
      ),
      files = list(
        national_boundary = basename(national_file),
        county_boundaries = basename(county_file)
      ),
      license = "See GADM license",
      FAIR = list(
        Findable = TRUE,
        Accessible = TRUE,
        Interoperable = TRUE,
        Reusable = TRUE
      )
    )

    metadata_file <- normalizePath(
      file.path(out_dir, "norway_boundaries_metadata.yml"),
      mustWork = FALSE
    )

    dir.create(
      dirname(metadata_file),
      recursive = TRUE,
      showWarnings = FALSE
    )

    if (!dir.exists(dirname(metadata_file))) {
      stop("Metadata directory does not exist: ", dirname(metadata_file))
    }

    if (file.access(dirname(metadata_file), mode = 2) != 0) {
      stop("No write permission for metadata directory: ", dirname(metadata_file))
    }

    yaml::write_yaml(
      metadata,
      file = metadata_file,
      fileEncoding = "UTF-8"
    )

    # Clean temporary files
    unlink(
      list.files(
        temp_dir,
        full.names = TRUE,
        all.files = TRUE,
        no.. = TRUE
      ),
      recursive = TRUE,
      force = TRUE
    )
  }

  return(
    list(
      national_boundary = national,
      county_boundaries = counties,
      files = list(
        national = national_file,
        counties = county_file,
        metadata = metadata_file
      )
    )
  )
}
