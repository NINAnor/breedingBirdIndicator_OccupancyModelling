plot_occupancy_covariates <- function(
  data,
  occ_covs,
  occupancy_col = "mean",
  ncol = 2,
  point_alpha = 0.2,
  point_size = 0.8,
  se = TRUE,
  y_label = "Occupancy probability",
  covariate_names = NULL,
  smooth_method = c("bam", "loess"),
  n_sample = 5000,
  seed = 2871
) {
  smooth_method <- match.arg(smooth_method)

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
  # Use a random sample of locations (without replacement) for the points and
  # the smoother. This keeps the plots fast and the PDFs small.
  # ---------------------------------------------------------------------------

  if (!is.null(n_sample)) {
    valid_rows <- which(!is.na(data[[occupancy_col]]))

    if (length(valid_rows) > n_sample) {
      if (!is.null(seed)) {
        # Restore the caller's random-number state afterwards
        old_seed <- if (exists(".Random.seed", envir = globalenv())) {
          get(".Random.seed", envir = globalenv())
        }
        on.exit(
          if (!is.null(old_seed)) {
            assign(".Random.seed", old_seed, envir = globalenv())
          },
          add = TRUE
        )
        set.seed(seed)
      }

      data <- data[
        sample(valid_rows, size = n_sample, replace = FALSE),
      ]
    }
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
      if (length(unique(values)) <= 7) {
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

  # ---------------------------------------------------------------------------
  # Smoother: fast GAM (mgcv::bam) by default, loess optionally
  # ---------------------------------------------------------------------------

  smoother_layer <- if (smooth_method == "bam") {
    ggplot2::geom_smooth(
      method = mgcv::bam,
      formula = y ~ s(x, bs = "cs"),
      method.args = list(method = "fREML", discrete = TRUE),
      se = se
    )
  } else {
    ggplot2::geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = se
    )
  }

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

          # Wrap long category names so that they do not overlap
          ggplot2::scale_x_discrete(
            labels = scales::label_wrap(10)
          ) +

          ggplot2::labs(
            x = if (is.null(covariate_names)) {
              cov
            } else {
              covariate_names$`Covariate definition`[
                covariate_names$covariate == cov
              ]
            },
            y = y_label
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

          smoother_layer +

          ggplot2::labs(
            x = if (is.null(covariate_names)) {
              cov
            } else {
              covariate_names$`Covariate definition`[
                covariate_names$covariate == cov
              ]
            },
            y = y_label
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
