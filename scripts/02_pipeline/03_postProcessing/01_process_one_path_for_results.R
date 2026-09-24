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
##
## ======================================

library(sf)
library(dplyr)
library(ggplot2)
library(viridis)

# ==========================================================
# Prepare data
# ==========================================================

richness_sf <- scaled_richness$prediction_values_scaled_richness %>%
  dplyr::select(
    nature_region,
    scaled_richness,
    geometry
  ) %>%
  sf::st_as_sf()


# Check regions
unique(
  richness_sf$nature_region
)


# ==========================================================
# Common plot theme
# ==========================================================

richness_theme <- theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    strip.text = element_text(
      face = "bold",
      size = 11
    )
  )


# ==========================================================
# 1. ALL BIOGEOGRAPHIC REGIONS TOGETHER
# ==========================================================

p_all_regions <- ggplot(
  richness_sf
) +

  geom_sf(
    aes(
      colour = scaled_richness
    ),
    size = 0.7,
    alpha = 0.8
  ) +

  scale_colour_viridis_c(
    name = "Scaled richness",
    limits = c(0, 1),
    na.value = "grey80"
  ) +

  labs(
    title = "Scaled bird richness",
    subtitle = "All biogeographic regions"
  ) +

  richness_theme

ggsave(
  filename = file.path(
    path_to_postProcessing,
    paste0(
      "scaled_richness_",
      "all",
      ".png"
    )
  ),
  plot = p_all_regions,
  width = 8,
  height = 8,
  dpi = 300
)

ggsave(
  filename = file.path(
    path_to_postProcessing,
    paste0(
      "scaled_richness_",
      "all",
      ".pdf"
    )
  ),
  plot = p_all_regions,
  width = 8,
  height = 8,
  dpi = 300
)


# ==========================================================
# 2. ONE PANEL PER BIOGEOGRAPHIC REGION
# ==========================================================

# p_regions <- ggplot(
#   richness_sf
# ) +
#
#   geom_sf(
#     aes(
#       colour = scaled_richness
#     ),
#     size = 0.6,
#     alpha = 0.8
#   ) +
#
#   scale_colour_viridis_c(
#     name = "Scaled richness",
#     limits = c(0, 1),
#     na.value = "grey80"
#   ) +
#
#   facet_wrap(
#     ~ bioGeographic_region
#   ) +
#
#   labs(
#     title = "Scaled bird richness by biogeographic region"
#   ) +
#
#   richness_theme
#
#
# p_regions

# ==========================================================
# Individual map for each region
# ==========================================================

regions <- unique(
  richness_sf$nature_region
)

for (reg in regions) {
  dat_reg <- richness_sf %>%
    dplyr::filter(
      nature_region == reg
    )

  p <- ggplot() +

    # All nature regions
    # geom_sf(
    #   data = bio_regions,
    #   fill = "grey95",
    #   colour = "grey40",
    #   linewidth = 0.4
    # ) +

    # Selected region boundary
    geom_sf(
      data = bio_regions %>%
        dplyr::filter(
          nature_region == reg
        ),
      fill = NA,
      colour = "black",
      linewidth = 0.8
    ) +

    # Richness observations
    geom_sf(
      data = dat_reg,
      aes(
        colour = scaled_richness
      ),
      size = 0.7,
      alpha = 0.8
    ) +

    scale_colour_viridis_c(
      name = "Scaled richness",
      limits = c(0, 1),
      na.value = "grey80"
    ) +

    labs(
      title = "Scaled bird richness",
      subtitle = reg
    ) +

    richness_theme

  reg_file <- gsub(
    "[^A-Za-z0-9_-]",
    "_",
    reg
  )

  ggsave(
    filename = file.path(
      path_to_postProcessing,
      paste0(
        "scaled_richness_",
        reg_file,
        ".png"
      )
    ),
    plot = p,
    width = 8,
    height = 8,
    dpi = 300
  )

  ggsave(
    filename = file.path(
      path_to_postProcessing,
      paste0(
        "scaled_richness_",
        reg_file,
        ".pdf"
      )
    ),
    plot = p,
    width = 8,
    height = 8,
    dpi = 300
  )
}
