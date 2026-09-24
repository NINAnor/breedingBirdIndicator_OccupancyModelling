predCoords <- predictionPoints_for_modelling %>%
  sf::st_as_sf(coords = c("long", "lat"), crs = modelCrs) %>%
  select(all_of(occ_covs)) %>%
  bind_cols(., terra::extract(multispecies_richness$expected_richness, .))

plot_occupancy_covariates(
  data = predCoords,
  occ_covs = occ_covs
)
