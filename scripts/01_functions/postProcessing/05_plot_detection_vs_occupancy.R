#' Plot detection probability against occupancy probability for all species
#'
#' One panel per species, showing a random sample of locations and a smoother.
#'
#' @param data Data frame with columns `species`, `occupancy` and `detection`.
#' @param n_sample Integer. Maximum number of distinct locations per species
#'   (sampled without replacement). `NULL` uses all of them.
#' @param smooth_method `"bam"` (default, fast GAM) or `"loess"`.
#' @param se Logical. Show the confidence band of the smoother.
#' @param ncol Number of panel columns.
#' @param point_alpha,point_size Point appearance.
#' @param seed Integer. Seed for the sample (the caller's random-number state is
#'   restored afterwards).
#'
#' @return A ggplot object.
#' @export
plot_detection_vs_occupancy <- function(
  data,
  n_sample = 2000,
  smooth_method = c("bam", "loess"),
  se = TRUE,
  ncol = 4,
  point_alpha = 0.3,
  point_size = 0.6,
  seed = 2871
) {
  smooth_method <- match.arg(smooth_method)

  missing_cols <- setdiff(c("species", "occupancy", "detection"), names(data))
  if (length(missing_cols) > 0) {
    stop(
      "`data` is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  data <- data[stats::complete.cases(data[c("occupancy", "detection")]), ]

  # Repeated surveys share the same prediction, so keep distinct locations
  data <- unique(data[c("species", "occupancy", "detection")])

  if (!is.null(n_sample)) {
    if (!is.null(seed)) {
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

    keep <- unlist(
      lapply(
        split(seq_len(nrow(data)), data$species),
        function(idx) {
          if (length(idx) > n_sample) {
            idx[sample.int(length(idx), n_sample)]
          } else {
            idx
          }
        }
      ),
      use.names = FALSE
    )
    data <- data[sort(keep), ]
  }

  smoother_layer <- if (smooth_method == "bam") {
    ggplot2::geom_smooth(
      method = mgcv::bam,
      formula = y ~ s(x, bs = "cs", k = 5),
      method.args = list(method = "fREML", discrete = TRUE),
      se = se
    )
  } else {
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = se)
  }

  data$species <- gsub("_", " ", data$species, fixed = TRUE)

  ggplot2::ggplot(
    data,
    ggplot2::aes(x = occupancy, y = detection)
  ) +
    ggplot2::geom_point(alpha = point_alpha, size = point_size) +
    smoother_layer +
    ggplot2::facet_wrap(~species, ncol = ncol, scales = "free") +
    ggplot2::labs(
      x = "Occupancy probability",
      y = "Detection probability"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}
