################################################################################
# Estimate ecosystem indicators using expert reference values
#
# This script:
#   1. Gets the expert reference quantile for each species.
#   2. Loads the nature regions.
#   3. For each set of best models, scales every species against its regional
#      reference value (x100), averages across species at each location, and
#      bootstraps regional (beta regression) and national (beta mixed model
#      with region as random effect) intercepts.
#   4. Saves the results (parquet for spatial data, csv for tables) and a
#      combined summary across model groups.
#
# Outputs per model group (see calculate_scaled_ecosystem_richness_with_expert_data):
#   indicator_vals.parquet, scaled_species.parquet, x100_values.csv,
#   species_summary.csv, representative_indicator_vals.csv, reference_data.csv,
#   all_results.rds
################################################################################

################################################################################
# Settings
################################################################################

# Prediction column holding mean occupancy probability
indicator_value_col <- "occupancy_mean"

# Sensitivity of the indicator to the reference quantile of "samme".
# Set runSensitivity <- FALSE before sourcing this script to skip it.
run_sensitivity <- if (exists("runSensitivity")) {
  isTRUE(runSensitivity)
} else {
  TRUE
}

################################################################################
# Reference values from experts
################################################################################

# Defaults: samme = 0.9, høyere = 0.99. Override with e.g. samme = 0.6.
reference_data <- get_reference_quantiles(
  expert_folder = expertFolder,
  ecosystem_type = ecosystemType,
  ref_quantiles = reference_quantile_lookup()
)

################################################################################
# Nature regions
################################################################################

gdb <- file.path(
  geoSpatFolder,
  "AdministrativeUnits",
  "Norway_AdministrativeUnits",
  "Original",
  "Norway_County",
  "versjon2025",
  "Administrative enheter fylker FGDB-format",
  "Basisdata_0000_Norge_25833_Fylker_FGDB",
  "Basisdata_0000_Norge_25833_Fylker_FGDB.gdb"
)

bio_regions <- get_nature_regions(
  input_file = gdb,
  modelCrs = modelCrs,
  output_file = file.path(
    dataFolder,
    "processedData",
    "nature_regions",
    "nature_regions.parquet"
  ),
  overwrite = FALSE,
  plot = TRUE,
  plot_file = file.path(
    dataFolder,
    "processedData",
    "nature_regions",
    "figures"
  )
)

################################################################################
# Helper functions
################################################################################

norm_path <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

# Prediction folder for a model-fit folder. Predictions are written to
# "prediction" (see 02_predictions.R); "prediction_outputs" is used as a
# fallback for older runs. The first folder that exists is returned.
resolve_prediction_dir <- function(model_path) {
  fit_root <- norm_path(file.path(output_root, "modelFit"))

  vapply(
    norm_path(model_path),
    function(p) {
      candidates <- c(
        sub(
          fit_root,
          norm_path(file.path(output_root, "prediction")),
          p,
          fixed = TRUE
        ),
        sub(
          fit_root,
          norm_path(file.path(output_root, "prediction_outputs")),
          p,
          fixed = TRUE
        )
      )
      existing <- candidates[dir.exists(candidates)]
      if (length(existing)) existing[1] else candidates[1]
    },
    character(1),
    USE.NAMES = FALSE
  )
}

calculate_indicator <- function(
  prediction_dir,
  output_dir,
  n.boot,
  n.Samps
) {
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  calculate_scaled_ecosystem_richness_with_expert_data(
    input_dir = prediction_dir,

    # Expert reference quantiles per species
    reference_data = reference_data,

    # Nature regions for regional summaries
    categorisation_shape = bio_regions,

    ecosystem_type = ecosystemType,

    value_col = indicator_value_col,
    region_col = "nature_region",

    # Output
    output_dir = output_dir,
    save_outputs = TRUE,
    overwrite = TRUE,

    # Bootstrap
    rep_values_for_each_region = TRUE,
    n.boot = n.boot,
    boot.seed = 123,
    n.Samps = n.Samps
  )
}

# Stack a table from each result and save it as one csv
combine_and_save <- function(results, element, group_name, file) {
  tables <- lapply(names(results), function(nm) {
    x <- results[[nm]][[element]]
    if (is.null(x)) {
      return(NULL)
    }
    x <- as.data.frame(x)
    x[[group_name]] <- nm
    x
  })

  combined <- dplyr::bind_rows(tables)
  if (nrow(combined) > 0) {
    utils::write.csv(combined, file, row.names = FALSE)
  }
  invisible(combined)
}

################################################################################
# Calculate indicators
################################################################################

if (compareMultipleOutputs) {
  #---------------------------------------------------------------------------
  # Best models (read the file, or create it if it does not exist)
  #---------------------------------------------------------------------------

  best_model_file <- file.path(
    output_root,
    "postProcessing_outputs",
    "multipleOutputs",
    ecosystemType,
    "best_models.csv"
  )

  if (!file.exists(best_model_file)) {
    source(here::here(
      "scripts",
      "02_pipeline",
      "03_postProcessing",
      "03_select_best_model_params.R"
    ))
  }

  best_model_details <- readr::read_csv(
    best_model_file,
    show_col_types = FALSE
  )

  #---------------------------------------------------------------------------
  # Output directory
  #---------------------------------------------------------------------------

  ind_dir <- file.path(
    output_root,
    "indicators_outputs",
    "compareMultipleOutputs",
    ecosystemType
  )

  dir.create(
    ind_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  #---------------------------------------------------------------------------
  # Split models into single- and multiple-species models
  #---------------------------------------------------------------------------

  model_groups <- split(
    best_model_details,
    ifelse(
      best_model_details$single_species,
      "single_species",
      "multiple_species"
    )
  )

  #---------------------------------------------------------------------------
  # Process each model type
  #---------------------------------------------------------------------------

  scaled_richness <- lapply(
    names(model_groups),
    function(group) {
      model_info <- model_groups[[group]]

      # Prediction directory: one folder per species for single-species
      # models, the model folder itself for the multi-species model
      pp_for_preds <- if (group == "single_species") {
        file.path(
          resolve_prediction_dir(model_info$path),
          model_info$Species
        )
      } else {
        unique(resolve_prediction_dir(model_info$path))
      }

      missing_preds <- pp_for_preds[!dir.exists(pp_for_preds)]
      if (length(missing_preds) > 0) {
        warning(
          "Prediction folders not found for ",
          group,
          ":\n",
          paste(missing_preds, collapse = "\n"),
          call. = FALSE
        )
      }
      pp_for_preds <- pp_for_preds[dir.exists(pp_for_preds)]

      if (length(pp_for_preds) == 0) {
        warning("No predictions to process for ", group, call. = FALSE)
        return(NULL)
      }

      calculate_indicator(
        prediction_dir = pp_for_preds,
        output_dir = file.path(ind_dir, group),
        n.boot = boot_settings_compare$n.boot,
        n.Samps = boot_settings_compare$n.Samps
      )
    }
  )

  names(scaled_richness) <- names(model_groups)
  scaled_richness <- Filter(Negate(is.null), scaled_richness)

  #---------------------------------------------------------------------------
  # Combined summaries across model groups
  #---------------------------------------------------------------------------

  combine_and_save(
    scaled_richness,
    "representative_indicator_vals",
    "model_group",
    file.path(ind_dir, "indicator_summary.csv")
  )
  combine_and_save(
    scaled_richness,
    "x100_values",
    "model_group",
    file.path(ind_dir, "x100_values_all_groups.csv")
  )
  combine_and_save(
    scaled_richness,
    "species_summary",
    "model_group",
    file.path(ind_dir, "species_summary_all_groups.csv")
  )

  # The sensitivity analysis uses the best multi-species model
  sens_pred_dir <- if ("multiple_species" %in% names(model_groups)) {
    unique(resolve_prediction_dir(model_groups$multiple_species$path))
  } else {
    NULL
  }
} else {
  #---------------------------------------------------------------------------
  # Standard analysis: the model of the current configuration
  #---------------------------------------------------------------------------

  scaled_richness <- calculate_indicator(
    prediction_dir = path_to_predictions,
    output_dir = path_to_indicators,
    n.boot = boot_settings_standard$n.boot,
    n.Samps = boot_settings_standard$n.Samps
  )

  sens_pred_dir <- path_to_predictions
}

################################################################################
# Sensitivity of the indicator to the reference quantile of "samme"
################################################################################

if (run_sensitivity) {
  if (is.null(sens_pred_dir) || !all(dir.exists(sens_pred_dir))) {
    warning(
      "Prediction folder for the sensitivity analysis not found; skipping.",
      call. = FALSE
    )
  } else {
    sens_dir <- file.path(
      output_root,
      "indicators_outputs",
      "sensitivity",
      ecosystemType,
      "samme_quantile"
    )
    sens_fig_dir <- file.path(sens_dir, "figures")
    dir.create(sens_fig_dir, recursive = TRUE, showWarnings = FALSE)

    sensitivity_results <- reference_quantile_sensitivity(
      quantile_values = sensitivity_quantiles,
      category = "samme",
      input_dir = sens_pred_dir,
      expert_folder = expertFolder,
      ecosystem_type = ecosystemType,
      categorisation_shape = bio_regions,
      region_col = "nature_region",
      value_col = indicator_value_col,
      n.boot = boot_settings_sensitivity$n.boot,
      n.Samps = boot_settings_sensitivity$n.Samps,
      boot.seed = 123,
      # Saves the combined tables as csv
      output_dir = sens_dir
    )

    # -----------------------------------------------------------------------------
    # Common plotting theme used for model-output figures.
    # -----------------------------------------------------------------------------

    theme_model <- function() {
      ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            face = "bold",
            hjust = 0.5
          ),
          plot.subtitle = ggplot2::element_text(
            hjust = 0.5
          ),
          axis.title = ggplot2::element_text(
            face = "bold"
          ),
          strip.text = ggplot2::element_text(
            face = "bold"
          ),
          legend.title = ggplot2::element_text(
            face = "bold"
          ),
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_line(
            linewidth = 0.2
          )
        )
    }

    sensitivity_plots <- plot_reference_quantile_sensitivity(
      sensitivity_results,
      category = "samme"
    )

    ggplot2::ggsave(
      file.path(sens_fig_dir, "indicator_estimate_vs_samme_quantile.png"),
      sensitivity_plots$estimates + theme_model(),
      width = 9,
      height = 6,
      dpi = 300
    )

    ggplot2::ggsave(
      file.path(sens_fig_dir, "x100_vs_samme_quantile_by_species.png"),
      sensitivity_plots$x100 + theme_model(),
      width = 9,
      height = 6,
      dpi = 300
    )
  }
}
