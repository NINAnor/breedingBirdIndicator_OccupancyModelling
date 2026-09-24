plot_occupancy_covariates <- function(
  data,
  occ_covs,
  occupancy_col = "mean",
  ncol = 2,
  point_alpha = 0.2,
  point_size = 0.8,
  se = TRUE
) {
  # ---------------------------------------------------------------------------
  # Check inputs
  # ---------------------------------------------------------------------------

  if (!inherits(data, "sf")) {
    stop("`data` must be an sf object.", call. = FALSE)
  }

  missing_covs <- setdiff(occ_covs, names(data))

  if (length(missing_covs) > 0) {
    stop(
      "The following covariates are missing from `data`: ",
      paste(missing_covs, collapse = ", "),
      call. = FALSE
    )
  }

  if (!occupancy_col %in% names(data)) {
    stop(
      "`occupancy_col` ('",
      occupancy_col,
      "') is not present in `data`.",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Determine whether each covariate is binary or continuous
  # ---------------------------------------------------------------------------

  covariate_type <- purrr::map_chr(
    occ_covs,
    function(x) {
      values <- data[[x]]

      # Remove NA
      values <- values[!is.na(values)]

      # Exactly two unique values = binary
      if (length(unique(values)) == 2) {
        "binary"
      } else {
        "continuous"
      }
    }
  )

  covariate_types <- tibble::tibble(
    covariate = occ_covs,
    type = covariate_type
  )

  print(covariate_types)

  # ---------------------------------------------------------------------------
  # Create individual plots
  # ---------------------------------------------------------------------------

  plots <- purrr::map2(
    occ_covs,
    covariate_type,
    function(cov, type) {
      plot_data <- data |>
        sf::st_drop_geometry() |>
        dplyr::select(
          covariate = dplyr::all_of(cov),
          occupancy = dplyr::all_of(occupancy_col)
        ) |>
        dplyr::filter(
          !is.na(covariate),
          !is.na(occupancy)
        )

      # -----------------------------------------------------------------------
      # Binary covariate
      # -----------------------------------------------------------------------

      if (type == "binary") {
        ggplot2::ggplot(
          plot_data,
          ggplot2::aes(
            x = factor(covariate),
            y = occupancy
          )
        ) +

          ggplot2::geom_boxplot(
            outlier.shape = NA
          ) +

          ggplot2::geom_jitter(
            width = 0.15,
            alpha = point_alpha,
            size = point_size
          ) +

          ggplot2::labs(
            x = cov,
            y = "Sum of mean occupancy"
          ) +

          ggplot2::theme_bw() +

          ggplot2::theme(
            strip.text = ggplot2::element_text(face = "bold"),
            panel.grid.minor = ggplot2::element_blank()
          )
      } else {
        # ---------------------------------------------------------------------
        # Continuous covariate
        # ---------------------------------------------------------------------

        ggplot2::ggplot(
          plot_data,
          ggplot2::aes(
            x = covariate,
            y = occupancy
          )
        ) +

          ggplot2::geom_point(
            alpha = point_alpha,
            size = point_size
          ) +

          ggplot2::geom_smooth(
            method = "gam",
            formula = y ~ s(x, bs = "cs"),
            se = se
          ) +

          ggplot2::labs(
            x = cov,
            y = "Sum of mean occupancy"
          ) +

          ggplot2::theme_bw() +

          ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank()
          )
      }
    }
  )

  # ---------------------------------------------------------------------------
  # Combine plots
  # ---------------------------------------------------------------------------

  patchwork::wrap_plots(
    plots,
    ncol = ncol
  )
}
