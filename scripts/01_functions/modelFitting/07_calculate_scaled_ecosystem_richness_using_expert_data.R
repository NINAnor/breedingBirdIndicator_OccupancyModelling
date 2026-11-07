#' Scaled ecosystem indicator using expert reference quantiles
#'
#' Scales each species' occupancy against a region-specific reference value
#' (the species' reference quantile of occupancy within the region, used as
#' `x100` in [ecTools::ec_normalise()]), averages the scaled values across
#' species at every location, and summarises the result by region and
#' nationally with bootstrapped beta models.
#'
#' @param input_dir Species prediction folders (see
#'   [species_specific_scaling()]).
#' @param reference_data Data frame with `Species` and `reference_quantile`
#'   (see [get_reference_quantiles()]).
#' @param categorisation_shape `sf` polygons defining the regions.
#' @param ecosystem_type Character. Ecosystem type.
#' @param value_col Character. Occupancy column in the prediction files.
#' @param region_col Character. Region column in `categorisation_shape`.
#' @param output_dir Output folder (needed when `save_outputs = TRUE`).
#' @param save_outputs,overwrite Logical. Save results; allow overwriting.
#' @param rep_values_for_each_region Logical. Run the bootstrap summaries.
#' @param n.boot Integer. Bootstrap replicates.
#' @param boot.seed Integer. Random seed.
#' @param n.Samps Integer. Maximum observations drawn per replicate (per
#'   region for the national model).
#' @param method One of `"raster"`, `"sf"`, `"rds"`: the prediction file type.
#' @param results Optional list from [read_species_predictions()], used
#'   instead of reading `input_dir` again.
#'
#' @return A list with
#'   \describe{
#'     \item{indicator_vals}{`sf` per location: `indicator_vals` (unscaled
#'       mean occupancy across species), `scaled_indicator` (mean of the
#'       scaled species values), `bio_region`, geometry.}
#'     \item{scaled_species}{`sf` per location with the scaled value of every
#'       species, `mean_scaled`, `bio_region` and geometry.}
#'     \item{x100_values}{Reference value (`x100`) per region and species.}
#'     \item{species_summary}{Mean, median, 2.5% and 97.5% quantiles of the
#'       scaled values per species and region (and `National`).}
#'     \item{representative_indicator_vals}{Bootstrap mean, median, 2.5% and
#'       97.5% quantiles of the intercept per region (beta regression) and
#'       `National` (beta mixed model with region as random effect).}
#'     \item{reference_data}{The reference data used.}
#'   }
#' @export
calculate_scaled_ecosystem_richness_with_expert_data <- function(
  input_dir = NULL,
  reference_data,
  categorisation_shape,
  ecosystem_type,
  value_col = "occupancy_mean",
  region_col = "nature_region",
  output_dir = NULL,
  save_outputs = FALSE,
  overwrite = TRUE,
  rep_values_for_each_region = TRUE,
  n.boot = 5000,
  boot.seed = 123,
  n.Samps = 5000,
  method = "rds",
  results = NULL
) {
  # ============================================================
  # 1. Validate inputs
  # ============================================================

  if (!is.data.frame(reference_data)) {
    stop("`reference_data` must be a data.frame.", call. = FALSE)
  }

  if (!inherits(categorisation_shape, "sf")) {
    stop("`categorisation_shape` must be an sf object.", call. = FALSE)
  }

  if (!method %in% c("raster", "sf", "rds")) {
    stop("`method` must be one of: 'raster', 'sf', or 'rds'.", call. = FALSE)
  }

  if (is.null(results) && length(input_dir) == 0) {
    stop("Supply `input_dir` or `results`.", call. = FALSE)
  }

  if (!region_col %in% names(categorisation_shape)) {
    stop(
      "`region_col` ('",
      region_col,
      "') is not present in `categorisation_shape`.",
      call. = FALSE
    )
  }

  missing_reference_cols <- setdiff(
    c("Species", "reference_quantile"),
    names(reference_data)
  )
  if (length(missing_reference_cols) > 0) {
    stop(
      "Missing columns in `reference_data`: ",
      paste(missing_reference_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (save_outputs && is.null(output_dir)) {
    stop("`output_dir` must be supplied when `save_outputs = TRUE`.", call. = FALSE)
  }

  if (save_outputs) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # ============================================================
  # 2. Scale each species by its regional reference value (x100)
  # ============================================================

  cat("\nScaling each species by its regional reference quantile...\n")

  scaling <- species_specific_scaling(
    input_dir = input_dir,
    reference_data = reference_data,
    categorisation_shape = categorisation_shape,
    region_col = region_col,
    ecosystem_type = ecosystem_type,
    value_col = value_col,
    method = method,
    save_outputs = FALSE,
    results = results
  )

  scaled_species <- scaling$scaled
  x100_values <- scaling$reference_values

  species_cols <- setdiff(
    names(scaled_species),
    c("mean_scaled", "bio_region", "geometry")
  )

  cat(
    "Species scaled: ",
    length(species_cols),
    "\n",
    sep = ""
  )

  # Per-location table for maps (compatible with the post-processing scripts)
  indicator_vals <- sf::st_sf(
    indicator_vals = scaling$mean_raw,
    scaled_indicator = scaled_species$mean_scaled,
    bio_region = scaled_species$bio_region,
    geometry = sf::st_geometry(scaled_species)
  )

  # Drop geometry once; bootstrap and summaries only need the values
  scaled_values <- sf::st_drop_geometry(scaled_species)

  region_names <- sort(unique(stats::na.omit(scaled_values$bio_region)))
  regions_for_study <- c("National", region_names)

  # ============================================================
  # 3. Species-specific summaries by region and nationally
  # ============================================================

  cat("\nSummarising species-specific scaled values...\n")

  summarise_values <- function(v) {
    v <- v[is.finite(v)]
    if (length(v) == 0) {
      return(
        data.frame(
          mean = NA_real_,
          median = NA_real_,
          lower = NA_real_,
          upper = NA_real_,
          n = 0L
        )
      )
    }
    q <- stats::quantile(v, probs = c(0.025, 0.5, 0.975), names = FALSE)
    data.frame(
      mean = mean(v),
      median = q[2],
      lower = q[1],
      upper = q[3],
      n = length(v)
    )
  }

  region_index <- c(
    list(National = seq_len(nrow(scaled_values))),
    split(seq_len(nrow(scaled_values)), scaled_values$bio_region)
  )

  species_summary <- do.call(
    rbind,
    lapply(species_cols, function(sp) {
      do.call(
        rbind,
        lapply(names(region_index), function(r) {
          cbind(
            species = gsub("_", " ", sp, fixed = TRUE),
            region = r,
            summarise_values(scaled_values[[sp]][region_index[[r]]]),
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )
  rownames(species_summary) <- NULL

  # ============================================================
  # 4. Representative values: region (beta regression) and
  #    national (beta mixed model, region as random effect)
  # ============================================================

  representative_indicator_vals <- NULL

  if (rep_values_for_each_region) {
    cat("\nEstimating representative values (bootstrap)...\n")

    indicator_df <- data.frame(
      scaled_indicator = scaled_values$mean_scaled,
      bio_region = scaled_values$bio_region
    )
    indicator_df <- indicator_df[
      is.finite(indicator_df$scaled_indicator) & !is.na(indicator_df$bio_region),
      ,
      drop = FALSE
    ]

    empty_row <- function(region, model) {
      data.frame(
        mean = NA_real_,
        median = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        sd = NA_real_,
        n = 0,
        region = region,
        model = model
      )
    }

    # ---- Regions: intercept-only beta regression ---------------------------
    regional_rows <- lapply(region_names, function(region) {
      cat(
        "\n----------------------------------------\n",
        "Region: ",
        region,
        "\n",
        sep = ""
      )

      values <- indicator_df$scaled_indicator[indicator_df$bio_region == region]

      if (length(values) == 0) {
        warning("No valid scaled indicator values for region: ", region)
        return(empty_row(region, "beta intercept"))
      }

      cat("  ", length(values), " valid values\n", sep = "")

      boot <- bootstrap_beta_intercept(
        data = data.frame(scaled_indicator = values),
        value_col = "scaled_indicator",
        n_boot = n.boot,
        seed = boot.seed,
        nSamps = n.Samps
      )

      data.frame(boot$summary) |>
        dplyr::mutate(n = boot$n, region = region, model = "beta intercept")
    })

    # ---- National: beta mixed model with region random effect ---------------
    cat(
      "\n----------------------------------------\n",
      "Region: National\n",
      sep = ""
    )

    national_row <- tryCatch(
      {
        boot <- bootstrap_beta_mixed_intercept(
          data = indicator_df,
          value_col = "scaled_indicator",
          region_col = "bio_region",
          n_boot = n.boot,
          seed = boot.seed,
          nSamps = n.Samps
        )
        data.frame(boot$summary) |>
          dplyr::mutate(
            n = boot$n,
            region = "National",
            model = "beta mixed intercept (region random effect)"
          )
      },
      error = function(e) {
        warning("National estimate failed: ", conditionMessage(e), call. = FALSE)
        empty_row("National", "beta mixed intercept (region random effect)")
      }
    )

    representative_indicator_vals <- dplyr::bind_rows(
      national_row,
      regional_rows
    )
  }

  # ============================================================
  # 5. Create output
  # ============================================================

  reference_data_used <- reference_data[
    reference_data$Species %in% gsub("_", " ", species_cols, fixed = TRUE),
    ,
    drop = FALSE
  ]

  output <- list(
    indicator_vals = indicator_vals,
    scaled_species = scaled_species,
    x100_values = x100_values,
    species_summary = species_summary,
    representative_indicator_vals = representative_indicator_vals,
    reference_data = reference_data_used
  )

  # ============================================================
  # 6. Save outputs
  # ============================================================

  if (save_outputs) {
    cat("\nSaving outputs...\n")

    save_object <- function(x, filename_base, dir = output_dir) {
      if (inherits(x, "sf")) {
        sfarrow::st_write_parquet(
          obj = x,
          dsn = file.path(dir, paste0(filename_base, ".parquet")),
          compression = "snappy"
        )
      } else {
        utils::write.csv(
          x,
          file = file.path(dir, paste0(filename_base, ".csv")),
          row.names = FALSE
        )
      }
    }

    save_object(indicator_vals, "indicator_vals")
    save_object(scaled_species, "scaled_species")
    save_object(x100_values, "x100_values")
    save_object(species_summary, "species_summary")
    save_object(reference_data_used, "reference_data")

    if (!is.null(representative_indicator_vals)) {
      save_object(representative_indicator_vals, "representative_indicator_vals")
    }

    saveRDS(output, file = file.path(output_dir, "all_results.rds"))
  }

  # ============================================================
  # 7. Return
  # ============================================================

  output
}
