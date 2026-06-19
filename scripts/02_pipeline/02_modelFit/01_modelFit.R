
# ---------------------------------------------------------------------------
# Load species grouping file
# ---------------------------------------------------------------------------
species_grouping_dir = file.path(path_finder("expertData"),
                                 "processedData",
                                 "speciesGroupings")

grouping_file <- switch(
  ecosystemType,
  "Forest"   = file.path(species_grouping_dir, "forestBirdSpecies.xlsx"),
  "Mountain" = file.path(species_grouping_dir, "mountainBirdSpecies.xlsx"),
  stop("ecosystemType must be either 'Forest' or 'Mountain'")
)

species_groupings <- readxl::read_excel(grouping_file)

if (!all(c("Species", "broad_guild") %in% names(species_groupings))) {
  stop("Species grouping file must contain columns: Species and broad_guild")
}

res_forest <- fit_ecosystem_spoccupancy(
  ecosystemBreedingBird = breeding_bird_data_for_modelling,
  ecosystemType = ecosystemType,

  occ_covariates = occ_covs,
  det_covariates = det_covs,

  species_groupings = NULL,
  maxSpeciesInGrps = 3,

  site_id_col  = "Route_point_id",
  species_col  = "Species",
  year_col     = "Year",
  count_col    = "TotalCount",
  lon_col      = "longitude",
  lat_col      = "latitude",

  route_id_col = NULL,
  point_id_col = NULL,

  use_route_random_effect = TRUE,

  # --- Model type switches ---
  single_species    = fit_single_species,   # TRUE -> PGOcc / spPGOcc per species
  fit_spatial_model = fit_spatial_model,    # TRUE -> spatial NNGP variant

  covariates_to_scale = NULL,
  remove_routes_with_missing_coords = isTRUE(fit_spatial_model),

  n.batch       = 5,
  batch.length  = 800,
  n.burn        = 2000,
  n.thin        = 10,
  n.chains      = 3,
  n.neighbors   = 5,
  n.omp.threads = 5,
  cov.model     = "exponential",

  fit_models   = TRUE,
  save_results = TRUE,
  verbose      = TRUE,
  output_dir   = path_to_results
)
