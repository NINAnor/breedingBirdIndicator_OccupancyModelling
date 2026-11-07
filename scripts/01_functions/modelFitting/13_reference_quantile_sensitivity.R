#' Sensitivity of the indicator to the reference quantile of a category
#'
#' Runs the expert-reference indicator for several reference quantiles of one
#' category (default `samme`) and collects the `x100` reference values and the
#' bootstrapped intercept estimates. The species predictions are read once and
#' the same bootstrap seed is used for every value, so differences between
#' values reflect the quantile and not bootstrap noise.
#'
#' @param quantile_values Numeric vector of quantiles to test for `category`.
#' @param category Character. Category to vary (`"samme"` or `"h\u00F8yere"`);
#'   the other category keeps its default quantile.
#' @param input_dir Species prediction folders.
#' @param expert_folder,ecosystem_type Passed to [get_reference_quantiles()].
#' @param categorisation_shape,region_col Regions.
#' @param value_col,method Prediction column and file type.
#' @param n.boot,n.Samps,boot.seed Bootstrap settings.
#' @param output_dir Optional folder; if given, the combined tables are saved
#'   there as csv.
#'
#' @return A list of data frames, each with a `quantile` column:
#'   `x100_values`, `representative_indicator_vals` and `species_summary`.
#' @export
reference_quantile_sensitivity <- function(
  quantile_values = c(0.5, 0.6, 0.7, 0.8, 0.9),
  category = "samme",
  input_dir,
  expert_folder,
  ecosystem_type,
  categorisation_shape,
  region_col = "nature_region",
  value_col = "occupancy_mean",
  method = "rds",
  n.boot = 100,
  n.Samps = 5000,
  boot.seed = 123,
  output_dir = NULL
) {
  build_reference <- function(q) {
    overrides <- stats::setNames(list(q), category)
    get_reference_quantiles(
      expert_folder = expert_folder,
      ecosystem_type = ecosystem_type,
      ref_quantiles = do.call(reference_quantile_lookup, overrides)
    )
  }

  # Read the predictions once (species are the same for every quantile)
  cat("Reading species predictions once...\n")
  results <- read_species_predictions(
    input_dir = input_dir,
    reference_data = build_reference(quantile_values[1]),
    value_col = value_col,
    method = method
  )

  runs <- lapply(quantile_values, function(q) {
    cat("\n==== ", category, " quantile = ", q, " ====\n", sep = "")

    out <- calculate_scaled_ecosystem_richness_with_expert_data(
      reference_data = build_reference(q),
      categorisation_shape = categorisation_shape,
      ecosystem_type = ecosystem_type,
      value_col = value_col,
      region_col = region_col,
      save_outputs = FALSE,
      rep_values_for_each_region = TRUE,
      n.boot = n.boot,
      boot.seed = boot.seed,
      n.Samps = n.Samps,
      method = method,
      results = results
    )

    add_q <- function(x) {
      x <- as.data.frame(x)
      x$quantile <- q
      x
    }

    # Keep only the small tables (the per-location sf objects are dropped)
    list(
      x100_values = add_q(out$x100_values),
      representative_indicator_vals = add_q(out$representative_indicator_vals),
      species_summary = add_q(out$species_summary)
    )
  })

  combined <- lapply(
    stats::setNames(nm = names(runs[[1]])),
    function(el) dplyr::bind_rows(lapply(runs, `[[`, el))
  )

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    for (nm in names(combined)) {
      utils::write.csv(
        combined[[nm]],
        file.path(output_dir, paste0("sensitivity_", nm, ".csv")),
        row.names = FALSE
      )
    }
  }

  combined
}


#' Plot the reference quantile sensitivity results
#'
#' @param sensitivity Output of [reference_quantile_sensitivity()].
#' @param category Character. Category that was varied (used in axis labels).
#' @param region_order Character vector giving the order of the panels.
#'
#' @return A list of two ggplot objects: `estimates` (bootstrap mean intercept
#'   and 95% interval against the quantile, by region) and `x100` (reference
#'   value against the quantile for each species of the varied category, by
#'   region).
#' @export
plot_reference_quantile_sensitivity <- function(
  sensitivity,
  category = "samme",
  region_order = c("National", "Central", "East", "North", "South", "West")
) {
  x_label <- paste0("Reference quantile for '", category, "'")

  estimates <- sensitivity$representative_indicator_vals
  estimates$region <- factor(
    estimates$region,
    levels = union(region_order, unique(estimates$region))
  )

  p_estimates <- ggplot2::ggplot(
    estimates,
    ggplot2::aes(quantile, mean)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.3
    ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~region) +
    ggplot2::labs(
      x = x_label,
      y = "Bootstrap mean intercept (scaled indicator)"
    )

  # Species whose reference quantile is the one being varied
  x100 <- sensitivity$x100_values
  x100 <- x100[x100$reference_quantile == x100$quantile, , drop = FALSE]

  p_x100 <- ggplot2::ggplot(
    x100,
    ggplot2::aes(quantile, reference_value, colour = species)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1) +
    ggplot2::facet_wrap(~bio_region) +
    ggplot2::labs(
      x = x_label,
      y = "x100 (reference value)",
      colour = "Species"
    )

  list(estimates = p_estimates, x100 = p_x100)
}
