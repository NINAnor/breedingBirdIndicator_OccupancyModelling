# ============================================================
# Extract model outputs
# ============================================================

if (!compareMultipleOutputs) {
  outputs <- extract_model_outputs(path_to_results, modelChecks)
} else {
  paths_df <- build_model_paths(
    output_root = output_root,
    version_labels = paste0("version_0", 1:2),
    ecosystemType = ecosystemType,
    existing_only = TRUE,
    include_single_species = c(TRUE, FALSE),
    include_spatial_model = c(TRUE, FALSE)
  )

  outputs <- extract_model_outputs(paths_df$path_to_results, modelChecks)
}

# ============================================================
# Save extracted outputs
# ============================================================

save(
  outputs,
  file = file.path(path_to_postProcessing, "model_outputs.RData")
)

print(outputs)

# ============================================================
# Common plotting theme
# ============================================================

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
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.2)
    )
}

# ============================================================
# Create plots
# ============================================================

p_species_occ <- plot(
  outputs,
  type = "error",
  components = c("occ_effects")
) +
  ggplot2::labs(
    title = "Occupancy effects",
    subtitle = ecosystemType
  ) +
  theme_model() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "bottom")


p_species_det <- plot(
  outputs,
  type = "error",
  components = c("det_effects")
) +
  ggplot2::labs(
    title = "Detection effects",
    subtitle = ecosystemType
  ) +
  theme_model() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "bottom")

if (modelChecks) {
  p_bayes <- plot(
    outputs,
    type = "error",
    components = c("bayes_P")
  ) +
    ggplot2::labs(
      title = "Posterior predictive checks",
      subtitle = ecosystemType,
      x = "Bayesian p-value"
    ) +
    theme_model()
}

p_comm_occ <- plot(
  outputs,
  type = "error",
  components = c("beta_summary") #, "alpha_summary")
) +
  ggplot2::labs(
    title = "Community-level parameter summaries",
    subtitle = ecosystemType,
    x = "Posterior mean (95% credible interval)"
  ) +
  theme_model() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "bottom")


p_comm_det <- plot(
  outputs,
  type = "error",
  components = c("alpha_summary")
) +
  ggplot2::labs(
    title = "Community-level parameter summaries",
    subtitle = ecosystemType,
    x = "Posterior mean (95% credible interval)"
  ) +
  theme_model() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "bottom")

# ============================================================
# Save helper
# ============================================================

save_plot <- function(plot_obj, filename, width = 14, height = 10) {
  ggplot2::ggsave(
    filename = file.path(
      path_to_postProcessing,
      paste0(filename, ".png")
    ),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )

  ggplot2::ggsave(
    filename = file.path(
      path_to_postProcessing,
      paste0(filename, ".pdf")
    ),
    plot = plot_obj,
    width = width,
    height = height
  )
}

# ============================================================
# Save figures
# ============================================================

save_plot(
  p_species_occ,
  "occupancy_effects"
)

save_plot(
  p_species_det,
  "detection_effects"
)

if (modelChecks) {
  save_plot(
    p_bayes,
    "posterior_predictive_checks"
  )
}

save_plot(
  p_comm_occ,
  "occupancy_community_parameter_summaries"
)

save_plot(
  p_comm_det,
  "detection_community_parameter_summaries"
)

# ============================================================
# Display plots in current session
# ============================================================

if (modelChecks) {
  print(p_bayes)
}

print(p_comm_occ)
print(p_comm_det)
print(p_species_occ)
print(p_species_det)

## =======================================
## Expected richness vrs expected richness
## ======================================

library(sf)
library(dplyr)
library(ggplot2)
library(viridis)

tif_files <- list.files(
  file.path(output_root, "indicators_outputs"),
  pattern = "\\.tif$",
  full.names = TRUE,
  recursive = TRUE
)

exp_richness_tif_files <- list.files(
  file.path(output_root, "indicators_outputs"),
  pattern = "expected_richness\\.tif$",
  full.names = TRUE,
  recursive = TRUE
)

exp_richness_tif_files <- exp_richness_tif_files[
  grepl(
    "[/\\\\]expected_richness[/\\\\]expected_richness\\.tif$",
    exp_richness_tif_files
  )
]

exp_richness_tif_rasters <- lapply(exp_richness_tif_files, terra::rast)

predCoords <- breeding_bird_data_for_modelling %>%
  sf::st_as_sf(coords = c("long", "lat"), crs = modelCrs) %>%
  dplyr::select(dplyr::all_of(occ_covs))

median_values <- lapply(
  exp_richness_tif_rasters,
  function(r) {
    terra::extract(r[["median"]], predCoords)[, 2]
  }
)

predCoords <- dplyr::bind_cols(
  predCoords,
  setNames(
    median_values,
    c(
      "richness_v1_multi_median",
      "richness_v1_single_median",
      "richness_v2_multi_median",
      "richness_v2_single_median"
    )
  )
)

validationPlots <- lapply(
  c(
    "richness_v1_multi_median",
    "richness_v1_single_median",
    "richness_v2_multi_median",
    "richness_v2_single_median"
  ),
  function(x) {
    ret <- plot_occupancy_covariates(
      data = predCoords,
      occ_covs = occ_covs,
      occupancy_col = x,
      ncol = 2,
      point_alpha = 0.2,
      point_size = 0.8,
      se = TRUE
    )

    ggsave(
      filename = file.path(path_to_postProcessing, paste0(x, ".png")),
      plot = ret,
      units = "in",
      width = 9,
      height = 9
    )

    ggsave(
      filename = file.path(path_to_postProcessing, paste0(x, ".pdf")),
      plot = ret,
      units = "in",
      width = 9,
      height = 9
    )
  }
)


#######################################
# Richness plot
#########################################

library(terra)
library(ggplot2)
# Output directory
plot_dir <- file.path(
  path_to_postProcessing,
  "indicator_plots"
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

raster_names <- c(
  "version_01_multiple_species",
  "version_01_single_species",
  "version_02_multiple_species",
  "version_02_single_species"
)

for (i in seq_along(exp_richness_tif_rasters)) {
  # --------------------------------------------------------------------------
  # Extract mean layer
  # --------------------------------------------------------------------------

  mean_raster <- exp_richness_tif_rasters[[i]][["mean"]]

  # Remove outer rows/columns containing only NA values
  mean_raster <- terra::trim(mean_raster)

  # --------------------------------------------------------------------------
  # Convert raster to data frame
  # --------------------------------------------------------------------------

  mean_df <- terra::as.data.frame(
    mean_raster,
    xy = TRUE,
    na.rm = TRUE
  )

  names(mean_df)[3] <- "expected_richness"

  # --------------------------------------------------------------------------
  # Plot
  # --------------------------------------------------------------------------

  p <- ggplot(
    mean_df,
    aes(
      x = x,
      y = y,
      fill = expected_richness
    )
  ) +
    geom_tile(
      width = 500,
      height = 500
    ) +
    coord_equal(
      expand = FALSE
    ) +
    scale_fill_viridis_c(
      option = "viridis",
      limits = c(2, 13),
      breaks = seq(2, 12, by = 2),
      name = "Expected\nrichness",
      na.value = "transparent"
    ) +
    labs(
      title = raster_names[i],
      x = "Easting",
      y = "Northing"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 16
      ),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  # --------------------------------------------------------------------------
  # Save
  # --------------------------------------------------------------------------

  ggsave(
    filename = file.path(
      plot_dir,
      paste0(
        raster_names[i],
        "_expected_richness_mean.png"
      )
    ),
    plot = p,
    width = 9,
    height = 10,
    dpi = 300
  )
}

###################################
# Plot indicator values
###################################

indicator_richness_tif_files <- list.files(
  file.path(output_root, "indicators_outputs"),
  pattern = "indicator_vals\\.tif$",
  full.names = TRUE,
  recursive = TRUE
)

indicator_richness_tif_rasters <- lapply(
  indicator_richness_tif_files,
  terra::rast
)


raster_names <- c(
  "version_01_multiple_species",
  "version_01_single_species",
  "version_02_multiple_species",
  "version_02_single_species"
)

for (i in seq_along(indicator_richness_tif_rasters)) {
  # --------------------------------------------------------------------------
  # Extract mean layer
  # --------------------------------------------------------------------------

  mean_raster <- indicator_richness_tif_rasters[[i]][["scaled_indicator"]]

  # Remove outer rows/columns containing only NA values
  mean_raster <- terra::trim(mean_raster)

  # --------------------------------------------------------------------------
  # Convert raster to data frame
  # --------------------------------------------------------------------------

  mean_df <- terra::as.data.frame(
    mean_raster,
    xy = TRUE,
    na.rm = TRUE
  )

  names(mean_df)[3] <- "expected_richness"

  # --------------------------------------------------------------------------
  # Plot
  # --------------------------------------------------------------------------

  p <- ggplot(
    mean_df,
    aes(
      x = x,
      y = y,
      fill = expected_richness
    )
  ) +
    geom_tile(
      width = 500,
      height = 500
    ) +
    coord_equal(
      expand = FALSE
    ) +
    scale_fill_viridis_c(
      option = "viridis",
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      name = "Indicator",
      na.value = "transparent"
    ) +
    labs(
      title = raster_names[i],
      x = "Easting",
      y = "Northing"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 16
      ),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  # --------------------------------------------------------------------------
  # Save
  # --------------------------------------------------------------------------

  ggsave(
    filename = file.path(
      plot_dir,
      paste0(
        raster_names[i],
        "_indicator_mean.png"
      )
    ),
    plot = p,
    width = 9,
    height = 10,
    dpi = 300
  )
}
