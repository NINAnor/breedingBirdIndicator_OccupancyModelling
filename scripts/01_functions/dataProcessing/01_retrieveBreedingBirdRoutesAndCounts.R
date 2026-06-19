###############################################################################
# Script: retrieve_breedingBirds_routesAndCounts.R
#
# Purpose:
#   Retrieve, clean, document, and export breeding bird point-count data,
#   route information, and species metadata for FIAR/TOV-E analyses.
#
# FAIR compliance:
#   Findable:
#     - Uses descriptive output filenames.
#     - Creates machine-readable metadata.
#     - Preserves route-point and EURING identifiers.
#
#   Accessible:
#     - Saves outputs to a user-defined data folder.
#     - Exports data as CSV, Parquet, and GeoPackage.
#
#   Interoperable:
#     - Uses tidy tabular data.
#     - Uses EURING species codes.
#     - Stores spatial data with explicit CRS information.
#
#   Reusable:
#     - Includes field-level metadata.
#     - Records provenance.
#     - Performs quality-control checks.
#     - Documents flock-to-pair conversion rules.
###############################################################################


retrieve_breedingBirds_routesAndCounts <- function(
    con = NULL,
    dataFolder,
    overwrite_cache = FALSE,
    save_outputs = TRUE,
    crs
) {

  ###########################################################################
  # 1. Install and load required packages
  ###########################################################################

  required_pkgs <- c(
    "dplyr",
    "tidyr",
    "readr",
    "stringr",
    "sf",
    "DBI",
    "dbplyr",
    "arrow",
    "tibble"
  )

  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs, dependencies = TRUE)
  }

  invisible(
    lapply(required_pkgs, library, character.only = TRUE)
  )


  ###########################################################################
  # 2. Validate user inputs
  ###########################################################################

  if (missing(dataFolder) || is.null(dataFolder)) {
    stop("'dataFolder' must be supplied.")
  }

  if (!dir.exists(dataFolder)) {
    dir.create(dataFolder, recursive = TRUE)
  }

  if (!is.logical(overwrite_cache)) {
    stop("'overwrite_cache' must be TRUE or FALSE.")
  }

  if (!is.logical(save_outputs)) {
    stop("'save_outputs' must be TRUE or FALSE.")
  }


  ###########################################################################
  # 3. Helper function: retrieve data from database or cached processed CSV
  #
  # Important:
  #   Cached CSV files are assumed to already be processed.
  #   Therefore, renaming and Route_point_id creation are only applied when
  #   data are retrieved directly from the database.
  ###########################################################################

  get_or_read_table <- function(table_name, file_name) {

    csv_path <- file.path(dataFolder, file_name)

    if (!is.null(con) && overwrite_cache) {

      message("Retrieving raw table from database: ", table_name)

      dat <- dplyr::tbl(con, table_name) |>
        dplyr::collect()

      readr::write_csv(dat, csv_path)

      return(dat)
    }

    if (file.exists(csv_path)) {

      message("Reading cached processed file: ", csv_path)

      return(
        readr::read_csv(csv_path, show_col_types = FALSE)
      )
    }

    if (!is.null(con)) {

      message("Cached file not found. Retrieving raw table from database: ", table_name)

      dat <- dplyr::tbl(con, table_name) |>
        dplyr::collect()

      readr::write_csv(dat, csv_path)

      return(dat)
    }

    stop(
      "Cached file not found and no database connection supplied: ",
      csv_path
    )
  }


  ###########################################################################
  # 4. Retrieve sampling-point information
  #
  # Source table:
  #   TakseringspunktInfo
  #
  # Scientific role:
  #   Contains route-point information and site-level variables such as
  #   coordinates, elevation, and other sampling-point attributes.
  ###########################################################################

  Sampling_point_info <- get_or_read_table(
    table_name = "TakseringspunktInfo",
    file_name  = "TakseringspunktInfo.csv"
  )

  if (!is.null(con) && overwrite_cache) {

    Sampling_point_info <- Sampling_point_info |>
      dplyr::rename(
        PointID = PktID,
        RouteID = RuteID
      ) |>
      tidyr::unite(
        col = "Route_point_id",
        RouteID,
        PointID,
        remove = FALSE
      )
  }


  ###########################################################################
  # 5. Check sampling-point columns
  ###########################################################################

  required_sampling_cols <- c(
    "Route_point_id",
    "RouteID",
    "PointID",
    "Latitude",
    "Longitude",
    "Hoh",
    "Hoh_DEM"
  )

  missing_sampling_cols <- setdiff(
    required_sampling_cols,
    names(Sampling_point_info)
  )

  if (length(missing_sampling_cols) > 0) {
    stop(
      "Sampling-point data are missing required columns: ",
      paste(missing_sampling_cols, collapse = ", ")
    )
  }


  ###########################################################################
  # 6. Retrieve verified breeding bird point-count data
  #
  # Source table:
  #   Punkttaksering_Verifisert
  #
  # Scientific role:
  #   Contains verified bird observations from route-based point counts.
  ###########################################################################

  Bird_data_points <- get_or_read_table(
    table_name = "Punkttaksering_Verifisert",
    file_name  = "Punkttaksering_Verifisert.csv"
  )

    Bird_data_points <- Bird_data_points |>
      tidyr::unite(
        col = "Route_point_id",
        RouteID,
        PointID,
        remove = FALSE
      )



  ###########################################################################
  # 7. Check bird-count columns
  ###########################################################################

  required_bird_cols <- c(
    "Route_point_id",
    "RouteID",
    "PointID",
    "YEAR",
    "EURINGCode",
    "ScientificName",
    "Obs<50m",
    "Obs>50m",
    "ObsFlokk"
  )

  missing_bird_cols <- setdiff(
    required_bird_cols,
    names(Bird_data_points)
  )

  if (length(missing_bird_cols) > 0) {
    stop(
      "Bird point-count data are missing required columns: ",
      paste(missing_bird_cols, collapse = ", ")
    )
  }


  ###########################################################################
  # 8. Join sampling-point information to bird observations
  #
  # FAIR relevance:
  #   Links each bird observation to a stable route-point identifier.
  #   This improves traceability and reuse in later spatial or ecological models.
  ###########################################################################

  Bird_data_points <- Bird_data_points |>
    dplyr::left_join(
      Sampling_point_info,
      by = "Route_point_id",
      suffix = c("", "_site")
    )


  ###########################################################################
  # 9. Retrieve species metadata
  #
  # Source table:
  #   Art
  #
  # Scientific role:
  #   Contains species names, EURING codes, and flock-conversion information.
  ###########################################################################

  SppID <- get_or_read_table(
    table_name = "Art",
    file_name  = "Art.csv"
  )

  if (!all(c("Species", "EURINGCode", "Artsnavn_Lat", "FK_Kode_Flokk") %in% names(SppID))) {

    SppID <- SppID |>
      dplyr::mutate(
        Species = as.numeric(ArtsID)
      ) |>
      dplyr::select(
        Species,
        EURINGCode,
        Artsnavn_Lat,
        FK_Kode_Flokk
      ) |>
      dplyr::filter(!is.na(Species), Species > 0) |>
      dplyr::arrange(Species)
  }


  ###########################################################################
  # 10. Check species metadata columns
  ###########################################################################

  required_species_cols <- c(
    "Species",
    "EURINGCode",
    "Artsnavn_Lat",
    "FK_Kode_Flokk"
  )

  missing_species_cols <- setdiff(
    required_species_cols,
    names(SppID)
  )

  if (length(missing_species_cols) > 0) {
    stop(
      "Species metadata are missing required columns: ",
      paste(missing_species_cols, collapse = ", ")
    )
  }


  ###########################################################################
  # 11. Join species metadata and clean species names
  #
  # FAIR relevance:
  #   EURING codes improve taxonomic interoperability.
  #   Scientific names are trimmed to avoid duplicate values caused by
  #   trailing spaces.
  ###########################################################################

  Bird_data <- Bird_data_points |>
    dplyr::left_join(
      SppID,
      by = "EURINGCode"
    ) |>
    dplyr::mutate(
      Species_name = stringr::str_trim(Artsnavn_Lat)
    )


  ###########################################################################
  # 12. Convert flock observations to estimated breeding pairs
  #
  # Scientific rule:
  #   Observers may record birds as flocks. Because the analytical unit is
  #   breeding pairs, flock counts are converted using species-specific codes:
  #
  #     FK_Kode_Flokk == 1: ObsFlokk / 5
  #     FK_Kode_Flokk == 2: ObsFlokk / 2
  #     Otherwise:          ObsFlokk unchanged
  ###########################################################################

  Bird_data <- Bird_data |>
    dplyr::mutate(
      dplyr::across(
        c(`Obs<50m`, `Obs>50m`, ObsFlokk),
        ~ suppressWarnings(as.numeric(.x))
      ),
      FK_Kode_Flokk = suppressWarnings(as.numeric(FK_Kode_Flokk)),
      Pairs_in_flock = dplyr::case_when(
        FK_Kode_Flokk == 1 ~ ObsFlokk / 5,
        FK_Kode_Flokk == 2 ~ ObsFlokk / 2,
        TRUE               ~ ObsFlokk
      ),
      TotalCount = `Obs<50m` + `Obs>50m` + Pairs_in_flock
    )


  ###########################################################################
  # 13. Quality-control checks
  #
  # FAIR relevance:
  #   Reusable data should include basic validation and transparent warnings.
  ###########################################################################

  if (any(Bird_data$TotalCount < 0, na.rm = TRUE)) {
    warning("Negative total counts detected.")
  }

  if (any(is.na(Bird_data$Species_name))) {
    warning("Some records are missing scientific species names.")
  }

  if (any(is.na(Bird_data$EURINGCode))) {
    warning("Some records are missing EURING codes.")
  }

  n_missing_coords <- Bird_data |>
    dplyr::filter(is.na(Longitude) | is.na(Latitude)) |>
    nrow()

  if (n_missing_coords > 0) {
    warning(n_missing_coords, " records have missing coordinates.")
  }

  n_missing_total <- sum(is.na(Bird_data$TotalCount))

  if (n_missing_total > 0) {
    warning(n_missing_total, " records have missing total counts.")
  }


  ###########################################################################
  # 14. Create spatial sf object
  #
  # CRS:
  #   Input coordinates are assumed to be longitude/latitude in WGS84.
  #   Output is transformed to EPSG:25833, ETRS89 / UTM zone 33N.
  ###########################################################################

  Bird_data_sf <- Bird_data |>
    dplyr::filter(!is.na(Longitude), !is.na(Latitude)) |>
    sf::st_as_sf(
      coords = c("Longitude", "Latitude"),
      crs = 4326,
      remove = FALSE
    ) |>
    sf::st_transform(crs)


  ###########################################################################
  # 15. Create route summary object
  ###########################################################################

  routes <- Bird_data_sf |>
    dplyr::select(RouteID, Route_point_id) |>
    dplyr::distinct()


  ###########################################################################
  # 16. Create final analytical dataset
  #
  # This is the main cleaned tabular dataset for modelling.
  ###########################################################################

  bird.dat <- Bird_data |>
    dplyr::select(
      Route_point_id,
     # RouteID,
      Year = YEAR,
      Species = Species_name,
      Euring = EURINGCode,
      `Obs<50m`,
      `Obs>50m`,
      Pairs_in_flock,
      TotalCount,
      Latitude,
      Longitude,
      DEM = Hoh_DEM
    )


  ###########################################################################
  # 17. Create field-level metadata
  #
  # FAIR relevance:
  #   Metadata make the dataset easier to find, interpret, and reuse.
  ###########################################################################

  field_metadata <- tibble::tibble(
    field = names(bird.dat),
    description = c(
      "Unique identifier combining survey route and point.",
      "Unique identifier for the ruote.",
      "Year in which the bird survey was conducted.",
      "Scientific species name.",
      "EURING species code used for taxonomic interoperability.",
      "Observed count within 50 metres of the sampling point.",
      "Observed count beyond 50 metres of the sampling point.",
      "Flock observations converted to estimated breeding pairs.",
      "Total estimated count: Obs<50m + Obs>50m + Pairs_in_flock.",
      "Latitude in WGS84 geographic coordinates.",
      "Longitude in WGS84 geographic coordinates.",
      "Elevation from digital elevation model."
    ),
    unit = c(
      "identifier",
      "year",
      "scientific name",
      "EURING code",
      "count",
      "count",
      "estimated pairs",
      "estimated pairs",
      "decimal degrees",
      "decimal degrees",
      "metres"
    )
  )


  ###########################################################################
  # 18. Create dataset-level provenance metadata
  #
  # FAIR relevance:
  #   Provenance records where the data came from, when it was processed,
  #   and what transformations were applied.
  ###########################################################################

  dataset_metadata <- tibble::tibble(
    metadata_field = c(
      "dataset_name",
      "processing_date",
      "processed_by",
      "source_tables",
      "cache_folder",
      "coordinate_reference_system_input",
      "coordinate_reference_system_output",
      "flock_conversion_rule",
      "output_formats",
      "script_function"
    ),
    value = c(
      "Cleaned breeding bird point-count data",
      as.character(Sys.time()),
      Sys.info()[["user"]],
      "TakseringspunktInfo; Punkttaksering_Verifisert; Art",
      normalizePath(dataFolder, winslash = "/", mustWork = FALSE),
      "EPSG:4326",
      paste0("EPSG:", crs),
      "FK_Kode_Flokk 1 = ObsFlokk/5; FK_Kode_Flokk 2 = ObsFlokk/2; otherwise unchanged",
      "CSV; Parquet; GeoPackage",
      "retrieve_breedingBirds_routesAndCounts"
    )
  )


  ###########################################################################
  # 19. Save FAIR-compliant outputs
  ###########################################################################

  if (save_outputs) {

    readr::write_csv(
      bird.dat,
      file.path(dataFolder, "breeding_bird_counts_clean.csv")
    )

    arrow::write_parquet(
      bird.dat,
      file.path(dataFolder, "breeding_bird_counts_clean.parquet")
    )

    sf::st_write(
      Bird_data_sf,
      file.path(dataFolder, "breeding_bird_counts_spatial.gpkg"),
      delete_dsn = TRUE,
      quiet = TRUE
    )

    readr::write_csv(
      field_metadata,
      file.path(dataFolder, "breeding_bird_counts_field_metadata.csv")
    )

    readr::write_csv(
      dataset_metadata,
      file.path(dataFolder, "breeding_bird_counts_dataset_metadata.csv")
    )
  }


  ###########################################################################
  # 20. Return objects to R
  ###########################################################################

  return(
    list(
      breedingBirdData = bird.dat,
      breedingBirdData_sf = Bird_data_sf,
      routes = routes,
      field_metadata = field_metadata,
      dataset_metadata = dataset_metadata
    )
  )
}
