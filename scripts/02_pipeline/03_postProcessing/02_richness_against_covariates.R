predCoords <- predictionPoints_for_modelling %>%
  sf::st_as_sf(coords = c("long", "lat"), crs = modelCrs) %>%
  select(all_of(occ_covs)) %>%
  bind_cols(., terra::extract(multispecies_richness$expected_richness, .))

plot_occupancy_covariates(
  data = predCoords,
  occ_covs = occ_covs
)


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
