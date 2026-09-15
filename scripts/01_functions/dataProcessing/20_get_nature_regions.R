#' Define and Cache Norwegian Nature Regions
#'
#' Reads Norwegian county polygons from a geodatabase, assigns counties
#' to predefined nature regions, dissolves county boundaries within each
#' region, transforms the result to the model CRS, and caches the result
#' as a GeoParquet file.
#'
#' @param input_file Character. Path to the Norwegian county geodatabase.
#' @param layer Character. Name of the county layer. Defaults to `"fylke"`.
#' @param modelCrs Coordinate reference system for the output.
#' @param output_file Character. Path to the cached GeoParquet file.
#' @param overwrite Logical. If `TRUE`, overwrite an existing cache.
#' @param plot Logical. If `TRUE`, save PDF and PNG maps of the regions.
#' @param plot_file Character. Output filename without extension for maps.
#'
#' @return An `sf` object containing the dissolved nature regions.
#'
#' @export
get_nature_regions <- function(
    input_file,
    layer = "fylke",
    modelCrs,
    output_file,
    overwrite = FALSE,
    plot = FALSE,
    plot_file = NULL
) {

  # ------------------------------------------------------------------
  # Required packages
  # ------------------------------------------------------------------

  required_pkgs <- c(
    "sf",
    "dplyr",
    "sfarrow",
    "ggplot2"
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
      "Missing packages: ",
      paste(missing_pkgs, collapse = ", "),
      call. = FALSE
    )
  }


  # ------------------------------------------------------------------
  # Validate input
  # ------------------------------------------------------------------

  if (!file.exists(input_file)) {
    stop(
      "Input geodatabase does not exist: ",
      input_file,
      call. = FALSE
    )
  }


  # ------------------------------------------------------------------
  # Prepare output directory
  # ------------------------------------------------------------------

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


  # ------------------------------------------------------------------
  # Read cached regions if available
  # ------------------------------------------------------------------

  if (file.exists(output_file) && !overwrite) {

    message("Reading cached nature regions.")

    nature_regions <- sfarrow::st_read_parquet(
      output_file
    )

  } else {

    message("Processing Norwegian nature regions.")


    # --------------------------------------------------------------
    # Read county polygons
    # --------------------------------------------------------------

    counties <- sf::st_read(
      input_file,
      layer = layer,
      quiet = TRUE
    )


    # --------------------------------------------------------------
    # Assign counties to nature regions
    # --------------------------------------------------------------

    nature_regions <- counties %>%

      dplyr::mutate(
        fylkesnummer = as.character(fylkesnummer),

        nature_region = dplyr::case_when(

          # South
          fylkesnummer %in% c(
            "40",
            "42"
          ) ~ "South",

          # East
          fylkesnummer %in% c(
            "31",
            "32",
            "03",
            "39",
            "33",
            "34"
          ) ~ "East",

          # West
          fylkesnummer %in% c(
            "46",
            "11"
          ) ~ "West",

          # Central
          fylkesnummer %in% c(
            "50",
            "15"
          ) ~ "Central",

          # North
          fylkesnummer %in% c(
            "18",
            "55",
            "56"
          ) ~ "North",

          TRUE ~ NA_character_
        )
      ) %>%

      dplyr::filter(
        !is.na(nature_region)
      )


    # --------------------------------------------------------------
    # Transform to model CRS
    # --------------------------------------------------------------

    nature_regions <- nature_regions %>%
      sf::st_transform(modelCrs) %>%
      sf::st_make_valid()


    # --------------------------------------------------------------
    # Dissolve counties within each nature region
    # --------------------------------------------------------------

    nature_regions <- nature_regions %>%

      dplyr::group_by(
        nature_region
      ) %>%

      dplyr::summarise(
        .groups = "drop"
      )


    # --------------------------------------------------------------
    # Write GeoParquet cache
    # --------------------------------------------------------------

    sfarrow::st_write_parquet(
      nature_regions,
      output_file
    )

    message(
      "Nature regions saved to: ",
      output_file
    )
  }


  # ------------------------------------------------------------------
  # Optional plotting
  # ------------------------------------------------------------------

  if (plot) {

    if (is.null(plot_file)) {
      plot_file <- tools::file_path_sans_ext(
        output_file
      )
    }

    plot_file <- normalizePath(
      plot_file,
      winslash = "/",
      mustWork = FALSE
    )

    dir.create(
      dirname(plot_file),
      recursive = TRUE,
      showWarnings = FALSE
    )


    p <- ggplot2::ggplot() +

      ggplot2::geom_sf(
        data = nature_regions,
        ggplot2::aes(
          fill = nature_region
        ),
        colour = "black",
        linewidth = 0.4
      ) +

      ggplot2::labs(
        fill = "Nature region"
      ) +

      ggplot2::theme_bw() +

      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        legend.position = "right"
      )


    # PDF
    ggplot2::ggsave(
      filename = paste0(
        plot_file,
        ".pdf"
      ),
      plot = p,
      width = 8,
      height = 8,
      units = "in"
    )


    # PNG
    ggplot2::ggsave(
      filename = paste0(
        plot_file,
        ".png"
      ),
      plot = p,
      width = 8,
      height = 8,
      units = "in",
      dpi = 300
    )

    message(
      "Nature region plots saved to: ",
      plot_file
    )
  }


  # ------------------------------------------------------------------
  # Return
  # ------------------------------------------------------------------

  return(nature_regions)
}
