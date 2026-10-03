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

  # Beta regression requires 0 < y < 1
  if (any(y <= 0 | y >= 1)) {
    stop(
      paste0(
        "Beta regression requires all values to be strictly ",
        "between 0 and 1. Found values outside this range."
      ),
      call. = FALSE
    )
  }

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
