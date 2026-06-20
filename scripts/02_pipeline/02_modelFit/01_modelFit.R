################################################################################
# Load species grouping information
#
# Species groupings can be used to fit groups of ecologically similar species
# together. Currently, species are fitted independently (species_groupings = NULL),
# but grouping files are loaded here for future analyses.
################################################################################

# Directory containing ecosystem-specific species grouping files
species_grouping_dir <- file.path(
  expertFolder,
  "processedData",
  "speciesGroupings"
)

# Select the appropriate grouping file for the ecosystem of interest
grouping_file <- switch(
  ecosystemType,
  "Forest"   = file.path(species_grouping_dir, "forestBirdSpecies.xlsx"),
  "Mountain" = file.path(species_grouping_dir, "mountainBirdSpecies.xlsx"),
  stop("ecosystemType must be either 'Forest' or 'Mountain'")
)

# Read species grouping table
species_groupings <- readxl::read_excel(grouping_file)

# Ensure required columns are present
required_cols <- c("Species", "broad_guild")

if (!all(required_cols %in% names(species_groupings))) {

  stop(
    paste(
      "Species grouping file must contain the columns:",
      paste(required_cols, collapse = ", ")
    ),
    call. = FALSE
  )

}


################################################################################
# Fit occupancy models
#
# Depending on the options selected above, the function fits either:
#
#  - Single-species models (PGOcc or spPGOcc)
#  - Spatial or non-spatial models
#
# Results are automatically saved to disk.
################################################################################

res_forest <- fit_ecosystem_spoccupancy(

  ##########################################################################
  # Input data
  ##########################################################################

  ecosystemBreedingBird = breeding_bird_data_for_modelling,
  ecosystemType         = ecosystemType,

  occ_covariates = occ_covs,
  det_covariates = det_covs,


  ##########################################################################
  # Species grouping options
  ##########################################################################

  # Set species_groupings = NULL to fit each species independently.
  species_groupings = NULL,

  # Maximum number of species per ecological group.
  maxSpeciesInGrps = maxSpeciesInGrps,


  ##########################################################################
  # Column names
  ##########################################################################

  site_id_col = "Route_point_id",
  species_col = "Species",
  year_col    = "Year",
  count_col   = "TotalCount",

  lon_col = "longitude",
  lat_col = "latitude",

  route_id_col = NULL,
  point_id_col = NULL,


  ##########################################################################
  # Random effects
  ##########################################################################

  # Include a random effect for survey route
  use_route_random_effect = TRUE,


  ##########################################################################
  # Model type
  ##########################################################################

  # TRUE: fit one model per species (PGOcc or spPGOcc)
  # FALSE: fit a multi-species occupancy model
  single_species = fit_single_species,

  # TRUE: fit spatial NNGP models
  # FALSE: fit non-spatial models
  fit_spatial_model = fit_spatial_model,


  ##########################################################################
  # Data preprocessing
  ##########################################################################

  covariates_to_scale = NULL,
  fill_detection_vals = TRUE,

  # Spatial models require coordinates for all routes
  remove_routes_with_missing_coords = isTRUE(fit_spatial_model),


  ##########################################################################
  # MCMC settings
  ##########################################################################

  n.batch      = n.batch,
  batch.length = batch.length,
  n.burn       = n.burn,
  n.thin       = n.thin,
  n.chains     = n.chains,


  ##########################################################################
  # Spatial settings
  ##########################################################################

  n.neighbors = n.neighbors,
  cov.model   = "exponential",


  ##########################################################################
  # Parallel processing
  ##########################################################################

  n.omp.threads = n.omp.threads,


  ##########################################################################
  # Output options
  ##########################################################################

  fit_models   = TRUE,
  save_results = TRUE,
  verbose      = TRUE,

  output_dir = path_to_results
)
