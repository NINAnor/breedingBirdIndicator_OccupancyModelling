#' Build and Align Model Covariate Raster Stack
#'
#' Loads requested model covariates, aligns them to a common template raster,
#' and returns a single `terra::SpatRaster` stack.
#'
#' This function supports FAIR principles by:
#' \itemize{
#'   \item \strong{Findable}: covariates are loaded from explicit file paths.
#'   \item \strong{Accessible}: outputs are standard `terra::SpatRaster` objects.
#'   \item \strong{Interoperable}: all rasters are aligned to one CRS, extent,
#'         resolution, and origin.
#'   \item \strong{Reusable}: the function records deterministic processing steps
#'         through explicit arguments and consistent layer names.
#' }
#'
#' @param x Object passed to `get_landcover_labels()`, usually the ecosystem or
#'   study object used to derive land-cover covariates.
#' @param modelCrs CRS used by the model and covariate rasters.
#' @param all_covariates Character vector of requested covariate names, e.g.
#'   `c("DEM", "LABEL1")`.
#' @param template A `terra::SpatRaster` used as the base raster for CRS,
#'   extent, resolution, and origin.
#' @param cov_folder Character. Folder containing land-cover covariate files.
#' @param data_folder Character. Base data folder. If `NULL`, uses
#'   `path_finder("dataFolder")`.
#'
#' @return A `terra::SpatRaster` containing aligned covariate layers.
#'
#' @export
put_covariates_together <- function(
    x,
    modelCrs,
    all_covariates,
    template,
    output_file,
    cov_folder = file.path(path_finder("dataFolder"), "Covariates"),
    data_folder = path_finder("dataFolder"),
    overwrite = FALSE,
    scale_vars = TRUE,
    norway_boundary,
    modelYear
) {
  required_pkgs <- c("sf", "terra", "dplyr", "tidyterra")

  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_pkgs) > 0) {
    stop(
      "Install missing packages first: ",
      paste(missing_pkgs, collapse = ", "),
      call. = FALSE
    )
  }

  if (!inherits(template, "SpatRaster")) {
    stop("`template` must be a terra::SpatRaster.", call. = FALSE)
  }

  if (!is.character(all_covariates)) {
    stop("`all_covariates` must be a character vector.", call. = FALSE)
  }

  years <- seq(modelYear[1], modelYear[2])

  output_file <- normalizePath(output_file, winslash = "/", mustWork = FALSE)

  if (file.exists(output_file) && !overwrite) {
    message("Found cached covariate stack: ", output_file)

    cached_stack <- terra::rast(output_file)

    same_as_template <- terra::compareGeom(
      cached_stack,
      template,
      stopOnError = FALSE
    )

    raster_names <- names(cached_stack)

    matched_cov_names <- vapply(
      all_covariates,
      function(cov_name) {
        any(startsWith(raster_names, cov_name))
      },
      logical(1)
    )

    all_matched <- all(matched_cov_names)

    if (same_as_template && all_matched) {

      message(
        "Cached covariate stack matches template and all requested covariates."
      )

      matched_rasters <- unlist(
        lapply(
          all_covariates,
          function(cov_name) {
            raster_names[
              startsWith(raster_names, cov_name)
            ]
          }
        ),
        use.names = FALSE
      )

      return(
        cached_stack[[matched_rasters]]
      )
    }

    if (!same_as_template) {
      message(
        "Cached covariate stack does not match template."
      )
    }

    if (!all_matched) {

      missing_covariates <- all_covariates[
        !matched_cov_names
      ]

      message(
        "Cached covariate stack is missing: ",
        paste(
          missing_covariates,
          collapse = ", "
        )
      )
    }

    message(
      "Rebuilding covariate stack."
    )

    # matched_cov_names <- all(all_covariates %in% names(cached_stack))
    #
    # if (same_as_template && matched_cov_names) {
    #   message("Cached covariate stack matches template and requested covariates.")
    #   return(cached_stack[[all_covariates]])
    # }
    #
    # if (!same_as_template) {
    #   message("Cached covariate stack does not match template.")
    # }
    #
    # if (!matched_cov_names) {
    #   missing_covariates <- setdiff(all_covariates, names(cached_stack))
    #   message(
    #     "Cached covariate stack is missing: ",
    #     paste(missing_covariates, collapse = ", ")
    #   )
    # }
    #
    # message("Rebuilding covariate stack.")
  }

  covariate_stack <- NULL

  ## =====================================
  ## DEM
  ## ====================================

  if ("DEM" %in% all_covariates) {
    dem <- get_DEM(
      output_file = file.path(
        data_folder,
        "processedData",
        "DEM",
        "DEM.tif"
      )
    )

    dem <- resample_to_baseRaster(
      r = dem,
      baseRaster = template,
      name = "DEM",
      scale_continuous = scale_vars
    )

    covariate_stack <- dem
  }

  ## =====================================
  ## Slope
  ## ====================================

  if("slope" %in% all_covariates){
    if(!exists("dem")){
      dem <- get_DEM(
        output_file = file.path(
          data_folder,
          "processedData",
          "DEM",
          "DEM.tif"
        )
      )
    }

    slope <- get_slope(dem,
                       output_file = file.path(
                         data_folder,
                         "processedData",
                         "DEM",
                         "slope.tif"
                       )
                       )

    slope <- resample_to_baseRaster(
      r = slope,
      baseRaster = template,
      name = "slope",
      scale_continuous = scale_vars
    )

    covariate_stack <- if (is.null(covariate_stack)) {
      slope
    } else {
      c(covariate_stack, slope)
    }
  }


  ## =====================================
  ## Aspect
  ## ====================================

  if("aspect" %in% all_covariates){
    if(!exists("dem")){
      dem <- get_DEM(
        output_file = file.path(
          data_folder,
          "processedData",
          "DEM",
          "DEM.tif"
        )
      )
    }

    aspect <- get_aspect(dem,
                       output_file = file.path(
                         data_folder,
                         "processedData",
                         "DEM",
                         "aspect.tif"
                       )
    )

    aspect <- resample_to_baseRaster(
      r = aspect,
      baseRaster = template,
      name = "aspect",
      scale_continuous = scale_vars
    )

    covariate_stack <- if (is.null(covariate_stack)) {
      aspect
    } else {
      c(covariate_stack, aspect)
    }
  }




  ## =====================================
  ## CORINE Landcover
  ## ====================================

  if ("corine_landCover" %in% all_covariates) {
    landCoverVars <- get_landcover_labels(
      x,
      modelCrs = modelCrs,
      cov_folder = cov_folder,
      output = "raster"
    ) |>
      tidyterra::select(LABEL1)

    landCoverVars <- resample_to_baseRaster(
      r = landCoverVars,
      baseRaster = template,
      name = "corine_landCover"
    )

    covariate_stack <- if (is.null(covariate_stack)) {
      landCoverVars
    } else {
      c(covariate_stack, landCoverVars)
    }
  }


  ## =====================================
  ## GLAD Landcover
  ## ====================================
  if("glad_landcover"  %in% all_covariates) {
    GLADRaster <- get_GLAD_landcover(
      input_file = file.path(
        dataFolder,
        "Covariates",
        "Landcover_F_M",
        "lc2020nor_mergedClasses.tif"
      ),
      model_crs = modelCrs,
      output_file = file.path(
        dataFolder,
        "processedData",
        "GLAD_lc",
        "LCglad.tif"
      )
    )

    GLADRaster <- resample_to_baseRaster(
      r = GLADRaster,
      baseRaster = template,
      name = "glad_landCover"
    )

    covariate_stack <- if (is.null(covariate_stack)) {
      GLADRaster
    } else {
      c(covariate_stack, GLADRaster)
    }
  }

  ## =====================================
  ## Snow melt
  ## ====================================
  if("snowmelt" %in% all_covariates){
    snow_dir <- file.path(
      path_finder("ninaGeospatialFiles"),
      "Meteorology",
      "Norway_SeNorge_snow",
      "Original",
      "sd"
    )

    output_dir <- file.path(
      dataFolder,
      "processedData",
      "Snowmelt"
    )

    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    snowmelt <- lapply(
      years,
      function(y) {

       get_snowmelt(
          snow_dir = snow_dir,
          output_file = file.path(
            output_dir,
            paste0("snowmelt_", y, ".tif")
          ),
          year = y,
          snow_threshold = 0,
          run_threshold = 14,
          target_crs = modelCrs,
          overwrite = FALSE
        ) %>%
          resample_to_baseRaster(
          r = .,
          baseRaster = template,
          name = paste0("snowmelt_", y),
          scale_continuous = scale_vars,
          max_unique = 10
        )
      }
    ) |>
      rast()

    covariate_stack <- if (is.null(covariate_stack)) {
      snowmelt
    } else {
      c(covariate_stack, snowmelt)
    }

  }

  ## =====================================
  ## Dominant treeline
  ## ====================================
  if("dominantTree" %in% all_covariates){

    sr16_dir <- file.path(
      dataFolder,
      "Covariates",
      "SR16_forest_F"
    )

    # output_dir <- file.path(
    #   dataFolder,
    #   "processedData",
    #   "TreeSpp"
    # )
    #
    # output_file <- file.path(
    #   output_dir,
    #   "TreeSpp.tif"
    # )

    TreeSppRaster <- get_TreeSpp(
      input_file = file.path(
        sr16_dir,
        "dominantTreeSpecies.tif"
      ),
      output_file = file.path(
        dataFolder,
        "processedData",
        "TreeSpp",
        "TreeSpp.tif"
      ),
      target_crs = paste0("EPSG:", modelCrs),
      overwrite = FALSE
    ) |>
      tidyterra::select(dominantTree)

    TreeSppRaster <- resample_to_baseRaster(
        r = TreeSppRaster,
        baseRaster = template,
        name = "dominantTree",
        scale_continuous = scale_vars,
        max_unique = 3
      )

    covariate_stack <- if (is.null(covariate_stack)) {
      TreeSppRaster
    } else {
      c(covariate_stack, TreeSppRaster)
    }
  }


  ## =====================================
  ## Mean temperature and precipitation
  ## ====================================

  if(any(c("mean_precipitation", "mean_temperature") %in% all_covariates)){
    prec_temp <- lapply(
      years,
      function(yr) {

        get_prec_temp(
          nc_dir = paste0(geoSpatFolder,
                          "/Meteorology/Norway_SeNorge2018_v22.09/Original/rr_tm"),
          output_dir = file.path(
            dataFolder,
            "processedData",
            "PrecTemp"
          ),
          year = yr,
          nc_pattern = "seNorge2018_%d.nc",
          tg_sub = "tg",
          rr_sub = "rr",
          breeding_season = 135:166,
          raster_crs = "EPSG:25833",
          overwrite = FALSE
        )
      }
    )

    if("mean_temperature" %in% all_covariates){
      temp <- lapply(prec_temp, function(dat){
        dat$temperature %>%
          resample_to_baseRaster(
            r = .,
            baseRaster = template,
            #name = "dominantTree",
            scale_continuous = scale_vars,
            max_unique = 10
          )
      })|>
        rast()

      covariate_stack <- if (is.null(covariate_stack)) {
        temp
      } else {
        c(covariate_stack, temp )
      }
    }

    if("mean_precipitation" %in% all_covariates){
      prec <- lapply(prec_temp, function(dat){
        dat$precipitation %>%
          resample_to_baseRaster(
            r = .,
            baseRaster = template,
            #name = "dominantTree",
            scale_continuous = scale_vars,
            max_unique = 10
          )
      })|>
        rast()

      covariate_stack <- if (is.null(covariate_stack)) {
        prec
      } else {
        c(covariate_stack, prec )
      }
    }

  }

  if("latitude" %in% all_covariates){
    latitude_raster <- get_latitude(
      boundary = norway_boundary,
      model_crs = modelCrs,
      template = template,
      output_file = file.path(
        dataFolder,
        "processedData",
        "Latitude",
        "latitude.tif"
      ),
      overwrite = FALSE
    ) %>%
      resample_to_baseRaster(
        r = .,
        baseRaster = template,
        #name = "dominantTree",
        scale_continuous = scale_vars,
        max_unique = 10
      )

    covariate_stack <- if (is.null(covariate_stack)) {
      latitude_raster
    } else {
      c(covariate_stack, latitude_raster )
    }
  }


  if (is.null(covariate_stack)) {
    stop(
      "No supported covariates found in `all_covariates`. Supported values currently include: DEM, LABEL1.",
      call. = FALSE
    )
  }

  #covariate_stack <- terra::rast(covariate_stack)
  names(covariate_stack) <- make.unique(names(covariate_stack))

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  terra::writeRaster(
    covariate_stack,
    output_file,
    overwrite = TRUE
  )

  covariate_stack
}
