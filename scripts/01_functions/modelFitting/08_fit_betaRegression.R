#' Squeeze values from [0, 1] into (0, 1) for beta regression
#'
#' Smithson and Verkuilen (2006): `(y * (n - 1) + 0.5) / n`. Applied only when
#' at least one value is exactly 0 or 1.
#'
#' @param y Numeric vector of finite values in [0, 1].
#' @param n Sample size used in the transformation.
#' @keywords internal
.squeeze_unit_interval <- function(y, n = length(y)) {
  if (any(y < 0 | y > 1)) {
    stop("Values must lie between 0 and 1.", call. = FALSE)
  }
  if (any(y <= 0 | y >= 1)) {
    y <- (y * (n - 1) + 0.5) / n
  }
  y
}

bootstrap_beta_intercept <- function(
  data,
  value_col = "scaled_indicator",
  n_boot = 5000,
  seed = 123,
  nSamps = 5000
) {
  set.seed(seed)

  # Extract response
  y <- data[[value_col]]

  # Remove non-finite values
  y <- y[is.finite(y)]

  if (length(y) == 0) {
    stop(
      "No finite values available for bootstrap.",
      call. = FALSE
    )
  }

  # Beta regression requires 0 < y < 1; values of exactly 0 or 1 are
  # squeezed into the open interval
  y <- .squeeze_unit_interval(y)

  # Number of observations per bootstrap replicate
  n <- min(
    nSamps,
    length(y)
  )

  # Storage
  bootstrap_estimates <- numeric(n_boot)

  # Progress bar
  pb <- utils::txtProgressBar(
    min = 0,
    max = n_boot,
    style = 3
  )

  on.exit(
    close(pb),
    add = TRUE
  )

  # Bootstrap
  for (i in seq_len(n_boot)) {
    # Conventional bootstrap:
    # sampling WITH replacement
    y_boot <- sample(
      y,
      size = n,
      replace = TRUE
    )

    # Fit intercept-only beta regression
    model <- betareg::betareg(
      y_boot ~ 1
    )

    # Transform logit intercept to response scale
    bootstrap_estimates[i] <- plogis(
      coef(model)[["(Intercept)"]]
    )

    # Update progress bar
    utils::setTxtProgressBar(
      pb,
      i
    )
  }

  # Bootstrap summary
  summary <- data.frame(
    mean = mean(
      bootstrap_estimates
    ),

    median = median(
      bootstrap_estimates
    ),

    lower = unname(
      quantile(
        bootstrap_estimates,
        probs = 0.025
      )
    ),

    upper = unname(
      quantile(
        bootstrap_estimates,
        probs = 0.975
      )
    ),

    sd = sd(
      bootstrap_estimates
    )
  )

  list(
    summary = summary,
    bootstrap_estimates = bootstrap_estimates,
    n = n,
    n_boot = n_boot
  )
}


#' Bootstrap the national intercept of a beta mixed model
#'
#' Pools the values of all regions and fits an intercept-only beta mixed model
#' with region as a random effect, `y ~ 1 + (1 | region)`, to each bootstrap
#' replicate. Each replicate resamples within regions (stratified bootstrap),
#' so every region keeps its own weight regardless of its size.
#'
#' @param data Data frame with the response and region columns.
#' @param value_col Character. Response column (values in [0, 1]).
#' @param region_col Character. Region column.
#' @param n_boot Integer. Number of bootstrap replicates.
#' @param seed Integer. Random seed.
#' @param nSamps Integer. Maximum observations drawn per region and replicate.
#'
#' @return A list like [bootstrap_beta_intercept()] with `summary` (mean,
#'   median, 2.5% and 97.5% quantiles, sd), `bootstrap_estimates`, `n` (total
#'   observations per replicate), `n_boot`, `n_failed` (replicates whose model
#'   did not fit) and `n_regions`.
#'
#' @details The estimate is the back-transformed fixed intercept, i.e. the
#'   typical region (random effect = 0) on the response scale. Values of
#'   exactly 0 or 1 are squeezed into (0, 1) first.
#' @export
bootstrap_beta_mixed_intercept <- function(
  data,
  value_col = "scaled_indicator",
  region_col = "bio_region",
  n_boot = 5000,
  seed = 123,
  nSamps = 5000
) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("Package `glmmTMB` is required.", call. = FALSE)
  }

  set.seed(seed)

  y <- data[[value_col]]
  region <- as.character(data[[region_col]])

  keep <- is.finite(y) & !is.na(region)
  y <- y[keep]
  region <- region[keep]

  if (length(y) == 0) {
    stop("No finite values available for bootstrap.", call. = FALSE)
  }
  if (length(unique(region)) < 2) {
    stop(
      "At least two regions are needed to fit a random effect for region.",
      call. = FALSE
    )
  }

  y <- .squeeze_unit_interval(y)

  idx_by_region <- split(seq_along(y), region)
  n_per_region <- pmin(nSamps, lengths(idx_by_region))

  bootstrap_estimates <- rep(NA_real_, n_boot)

  pb <- utils::txtProgressBar(min = 0, max = n_boot, style = 3)
  on.exit(close(pb), add = TRUE)

  for (i in seq_len(n_boot)) {
    sampled <- unlist(
      Map(
        function(idx, m) idx[sample.int(length(idx), size = m, replace = TRUE)],
        idx_by_region,
        n_per_region
      ),
      use.names = FALSE
    )

    boot_data <- data.frame(y = y[sampled], region = factor(region[sampled]))

    fit <- tryCatch(
      suppressWarnings(
        glmmTMB::glmmTMB(
          y ~ 1 + (1 | region),
          data = boot_data,
          family = glmmTMB::beta_family()
        )
      ),
      error = function(e) NULL
    )

    if (!is.null(fit)) {
      est <- tryCatch(
        stats::plogis(glmmTMB::fixef(fit)$cond[["(Intercept)"]]),
        error = function(e) NA_real_
      )
      bootstrap_estimates[i] <- est
    }

    utils::setTxtProgressBar(pb, i)
  }

  n_failed <- sum(is.na(bootstrap_estimates))
  if (n_failed == n_boot) {
    stop("All bootstrap replicates failed to fit.", call. = FALSE)
  }
  if (n_failed > 0) {
    warning(
      n_failed, " of ", n_boot, " bootstrap replicates failed to fit.",
      call. = FALSE
    )
  }

  ok <- bootstrap_estimates[!is.na(bootstrap_estimates)]

  list(
    summary = data.frame(
      mean = mean(ok),
      median = stats::median(ok),
      lower = unname(stats::quantile(ok, probs = 0.025)),
      upper = unname(stats::quantile(ok, probs = 0.975)),
      sd = stats::sd(ok)
    ),
    bootstrap_estimates = bootstrap_estimates,
    n = sum(n_per_region),
    n_boot = n_boot,
    n_failed = n_failed,
    n_regions = length(idx_by_region)
  )
}
