#' Create Latitude Raster
#'
#' Creates a raster containing geographic latitude in decimal degrees for
#' each raster cell centre. Cell coordinates are transformed from the
#' model CRS to EPSG:4326 in memory-efficient blocks, avoiding creation
#' of a coordinate object for the entire raster.
#'
#' The output is written incrementally to disk, making the function suitable
#' for very large rasters containing hundreds of millions of cells.
#'
#' @param boundary An `sf` polygon containing the spatial extent for the
#'   latitude raster.
#' @param model_crs Character or numeric. CRS used by the model,
#'   e.g. `"EPSG:25833"` or `25833`.
#' @param template Character or `terra::SpatRaster`. Raster template defining
#'   the extent, resolution, dimensions, and CRS of the output.
#' @param output_file Character. Path to the output GeoTIFF.
#' @param metadata_file Character or `NULL`. Path to YAML metadata file.
#' @param overwrite Logical. If `TRUE`, recreate the raster.
#' @param block_rows Integer. Number of raster rows processed at a time.
#'   Defaults to `1000`.
#'
#' @return A `terra::SpatRaster` containing latitude in decimal degrees.
#'
#' @export
get_latitude_block <- function(
    boundary,
    model_crs,
    template,
    output_file,
    metadata_file = NULL,
    overwrite = FALSE,
    block_rows = 1000
) {

  # --------------------------------------------------------
  # Required packages
  # --------------------------------------------------------

  required_pkgs <- c(
    "terra",
    "sf",
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
      paste(
        missing_pkgs,
        collapse = ", "
      ),
      call. = FALSE
    )
  }


  # --------------------------------------------------------
  # Validate block size
  # --------------------------------------------------------

  if (
    length(block_rows) != 1 ||
    !is.numeric(block_rows) ||
    is.na(block_rows) ||
    block_rows < 1
  ) {

    stop(
      "`block_rows` must be a positive integer.",
      call. = FALSE
    )
  }

  block_rows <- as.integer(
    block_rows
  )


  # --------------------------------------------------------
  # Validate model CRS
  # --------------------------------------------------------

  if (is.numeric(model_crs)) {

    model_crs <- paste0(
      "EPSG:",
      model_crs
    )
  }

  model_crs <- as.character(
    model_crs
  )

  model_crs_obj <- tryCatch(

    sf::st_crs(model_crs),

    error = function(e) {

      stop(
        "Invalid model CRS: ",
        model_crs,
        "\n",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (is.na(model_crs_obj)) {

    stop(
      "Invalid model CRS: ",
      model_crs,
      call. = FALSE
    )
  }


  # --------------------------------------------------------
  # Validate boundary
  # --------------------------------------------------------

  if (!inherits(boundary, "sf")) {

    stop(
      "`boundary` must be an sf object.",
      call. = FALSE
    )
  }

  if (is.na(sf::st_crs(boundary))) {

    stop(
      "`boundary` does not have a valid CRS.",
      call. = FALSE
    )
  }


  # --------------------------------------------------------
  # Transform boundary to model CRS
  # --------------------------------------------------------

  if (
    !sf::st_crs(boundary) == model_crs_obj
  ) {

    message(
      "Transforming boundary to model CRS: ",
      model_crs
    )

    boundary <- sf::st_transform(
      boundary,
      model_crs_obj
    )
  }

  boundary <- sf::st_make_valid(
    boundary
  )


  # --------------------------------------------------------
  # Normalise output path
  # --------------------------------------------------------

  output_file <- normalizePath(
    output_file,
    winslash = "/",
    mustWork = FALSE
  )

  dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
  )


  # --------------------------------------------------------
  # Read cached raster
  # --------------------------------------------------------

  if (
    file.exists(output_file) &&
    !overwrite
  ) {

    message(
      "Reading cached latitude raster: ",
      output_file
    )

    cached_latitude <- terra::rast(
      output_file
    )

    if (
      terra::same.crs(
        cached_latitude,
        model_crs
      )
    ) {

      message(
        "Cached latitude raster already uses model CRS."
      )

      return(
        cached_latitude
      )
    }

    message(
      "Cached latitude raster has an unexpected CRS."
    )

    message(
      "Recreating raster."
    )
  }


  # --------------------------------------------------------
  # Read template
  # --------------------------------------------------------

  if (inherits(template, "SpatRaster")) {

    template_raster <- template

  } else if (is.character(template)) {

    template <- normalizePath(
      template,
      winslash = "/",
      mustWork = TRUE
    )

    template_raster <- terra::rast(
      template
    )

  } else {

    stop(
      "`template` must be either a terra::SpatRaster ",
      "or a path to a raster file.",
      call. = FALSE
    )
  }


  # --------------------------------------------------------
  # Validate template CRS
  # --------------------------------------------------------

  if (
    !terra::same.crs(
      template_raster,
      model_crs
    )
  ) {

    stop(
      "The supplied raster template does not use model_crs.\n",
      "Template CRS: ",
      terra::crs(template_raster),
      "\n",
      "Model CRS: ",
      model_crs,
      call. = FALSE
    )
  }


  # --------------------------------------------------------
  # Create output raster
  # --------------------------------------------------------

  latitude_raster <- template_raster[[1]]

  names(latitude_raster) <- "latitude"


  # --------------------------------------------------------
  # Crop template to boundary extent
  #
  # This can substantially reduce the number of cells if the
  # template extends beyond the study area.
  # --------------------------------------------------------

  boundary_vect <- terra::vect(
    boundary
  )

  latitude_raster <- terra::crop(
    latitude_raster,
    boundary_vect
  )

  names(latitude_raster) <- "latitude"


  # --------------------------------------------------------
  # Raster dimensions
  # --------------------------------------------------------

  n_rows <- terra::nrow(
    latitude_raster
  )

  n_cols <- terra::ncol(
    latitude_raster
  )

  n_cells <- terra::ncell(
    latitude_raster
  )

  message(
    "Latitude raster dimensions: ",
    n_rows,
    " rows x ",
    n_cols,
    " columns"
  )

  message(
    "Cells to process: ",
    format(
      n_cells,
      big.mark = ","
    )
  )


  # --------------------------------------------------------
  # Temporary output raster
  #
  # FLT4S is sufficient because latitude is rounded to four
  # decimal places.
  # --------------------------------------------------------

  latitude_raster <- terra::rast(
    latitude_raster
  )

  names(latitude_raster) <- "latitude"


  # --------------------------------------------------------
  # Start block-wise writing
  # --------------------------------------------------------

  latitude_raster <- terra::writeStart(
    latitude_raster,
    filename = output_file,
    overwrite = TRUE,
    datatype = "FLT4S",
    gdal = c(
      "COMPRESS=DEFLATE",
      "PREDICTOR=2",
      "ZLEVEL=6"
    )
  )


  # --------------------------------------------------------
  # Number of blocks
  # --------------------------------------------------------

  n_blocks <- ceiling(
    n_rows / block_rows
  )


  message(
    "Processing ",
    n_blocks,
    " blocks of approximately ",
    block_rows,
    " rows."
  )


  # --------------------------------------------------------
  # Process raster block-by-block
  # --------------------------------------------------------

  for (block_id in seq_len(n_blocks)) {

    row_start <- (
      (block_id - 1L) * block_rows
    ) + 1L

    row_end <- min(
      block_id * block_rows,
      n_rows
    )

    rows <- row_start:row_end


    # ------------------------------------------------------
    # Cell indices for current block
    # ------------------------------------------------------

    cells <- terra::cellFromRowCol(
      latitude_raster,
      row = rep(
        rows,
        each = n_cols
      ),
      col = rep(
        seq_len(n_cols),
        times = length(rows)
      )
    )


    # ------------------------------------------------------
    # Get projected coordinates
    # ------------------------------------------------------

    xy <- terra::xyFromCell(
      latitude_raster,
      cells
    )


    # ------------------------------------------------------
    # Transform coordinates directly
    #
    # terra::project() transforms the coordinate matrix
    # without creating an sf object containing millions
    # of point geometries.
    # ------------------------------------------------------

    xy_lonlat <- terra::project(
      xy,
      from = model_crs,
      to = "EPSG:4326"
    )


    # ------------------------------------------------------
    # Extract latitude
    # ------------------------------------------------------

    latitude_values <- round(
      xy_lonlat[, 2],
      4
    )


    # ------------------------------------------------------
    # Mask cells outside boundary
    # ------------------------------------------------------

    block_raster <- terra::rast(
      latitude_raster,
      nrows = length(rows),
      ncols = n_cols,
      extent = terra::ext(
        latitude_raster
      )
    )

    # Correct block extent
    block_raster <- terra::crop(
      latitude_raster,
      terra::ext(
        latitude_raster,
        row = rows
      )
    )

    terra::values(
      block_raster
    ) <- latitude_values


    block_mask <- terra::rasterize(
      boundary_vect,
      block_raster,
      field = 1,
      background = NA
    )

    latitude_values[
      is.na(
        terra::values(block_mask)
      )
    ] <- NA_real_


    # ------------------------------------------------------
    # Write block
    # ------------------------------------------------------

    latitude_raster <- terra::writeValues(
      latitude_raster,
      latitude_values,
      start = row_start
    )


    # ------------------------------------------------------
    # Progress
    # ------------------------------------------------------

    if (
      block_id == 1L ||
      block_id %% 10L == 0L ||
      block_id == n_blocks
    ) {

      message(
        sprintf(
          "Block %d/%d (rows %d-%d)",
          block_id,
          n_blocks,
          row_start,
          row_end
        )
      )
    }
  }


  # --------------------------------------------------------
  # Finish writing
  # --------------------------------------------------------

  latitude_raster <- terra::writeStop(
    latitude_raster
  )


  # --------------------------------------------------------
  # Metadata path
  # --------------------------------------------------------

  if (is.null(metadata_file)) {

    metadata_file <- file.path(
      dirname(output_file),
      paste0(
        tools::file_path_sans_ext(
          basename(output_file)
        ),
        "_metadata.yml"
      )
    )
  }

  metadata_file <- normalizePath(
    metadata_file,
    winslash = "/",
    mustWork = FALSE
  )

  dir.create(
    dirname(metadata_file),
    recursive = TRUE,
    showWarnings = FALSE
  )


  # --------------------------------------------------------
  # Metadata
  # --------------------------------------------------------

  metadata <- list(

    title = "Latitude raster",

    description = paste(
      "Raster containing geographic latitude at each raster",
      "cell centre, calculated from the model CRS."
    ),

    created = as.character(
      Sys.time()
    ),

    input = list(

      boundary_crs =
        sf::st_crs(boundary)$input,

      model_crs =
        model_crs,

      geographic_crs =
        "EPSG:4326",

      template =
        if (inherits(template, "SpatRaster")) {
          "terra::SpatRaster"
        } else {
          template
        }
    ),

    output = list(

      file =
        output_file,

      file_format =
        "GeoTIFF",

      variable =
        "latitude",

      units =
        "decimal degrees",

      decimal_places =
        4,

      datatype =
        "FLT4S"
    ),

    processing_method = list(

      coordinate_transformation =
        "terra::project()",

      cell_coordinates =
        "terra::xyFromCell()",

      processing =
        "Block-wise raster processing",

      geographic_crs =
        "EPSG:4326",

      rounding =
        "Latitude values rounded to 4 decimal places.",

      masking =
        "Raster masked to the supplied boundary."
    ),

    raster_properties = list(

      extent =
        as.character(
          terra::ext(
            latitude_raster
          )
        ),

      resolution =
        as.numeric(
          terra::res(
            latitude_raster
          )
        ),

      dimensions = c(

        nrow =
          terra::nrow(
            latitude_raster
          ),

        ncol =
          terra::ncol(
            latitude_raster
          ),

        nlyr =
          terra::nlyr(
            latitude_raster
          )
      ),

      names =
        names(
          latitude_raster
        )
    )
  )


  # --------------------------------------------------------
  # Write metadata
  # --------------------------------------------------------

  tryCatch(

    yaml::write_yaml(
      metadata,
      metadata_file
    ),

    error = function(e) {

      warning(
        "Latitude raster was saved, but metadata could not ",
        "be written: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )


  message(
    "Latitude raster saved to: ",
    output_file
  )

  message(
    "Metadata saved to: ",
    metadata_file
  )


  # --------------------------------------------------------
  # Return
  # --------------------------------------------------------

  latitude_raster
}
