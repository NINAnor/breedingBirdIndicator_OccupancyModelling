################################################################################
# Load spatial reference layers and calculate ecosystem indicators
#
# This script:
#   1. Loads and processes biogeographic regions.
#   2. Loads and processes protected areas.
#   3. Calculates scaled ecosystem richness from species occupancy predictions.
#
# Assumptions:
#   - Occupancy predictions represent landscape-level mean occupancy.
#   - Protected areas approximate reference ecological conditions.
#   - Indicators are summarized within biogeographic regions.
#   - Richness estimates are calculated on a 1-km grid.
################################################################################


################################################################################
# Load biogeographic regions
################################################################################

# R:\GeoSpatialData\BiogeographicalRegions\Norway_Biogeography\Original
# bio_regions <- get_biogeographic_region(
#   input_file = file.path(
#     geoSpatFolder,
#     "BiogeographicalRegions",
#     "Norway_Biogeography",
#     "Original",
#     "biogeografiske_regioner_omrade.shp"
#   ),
#   norway_boundary = norway_boundary,
#   modelCrs = modelCrs,
#   output_file = file.path(
#     dataFolder,
#     "processedData",
#     "biogeographic_regions",
#     "biogeographic_regions.parquet"
#   ),
#   overwrite = FALSE
# )


gdb <- paste0(geoSpatFolder,
  "/AdministrativeUnits/",
  "Norway_AdministrativeUnits/Original/Norway_County/",
  "versjon2025/Administrative enheter fylker FGDB-format/",
  "Basisdata_0000_Norge_25833_Fylker_FGDB/",
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
# Load protected areas
################################################################################
# The file path is from NINA's R drive
# R:\GeoSpatialData\ProtectedSites\Norway_ProtectedAreas\Original\Natur_Vern_norge_med_svalbard_25833.gdb
protected_regions <- get_protected_areas(
  input_file = file.path(
    geoSpatFolder,
    "ProtectedSites",
    "Norway_ProtectedAreas",
    "Original",
    "Natur_Vern_norge_med_svalbard_25833.gdb"
  ),
  norway_boundary = norway_boundary,
  modelCrs = modelCrs,
  output_file = file.path(
    dataFolder,
    "processedData",
    "protected_regions",
    "protected_regions.parquet"
  ),
  layer = "Naturvernområde",
  overwrite = FALSE
)


################################################################################
# Calculate scaled ecosystem richness
################################################################################

scaled_richness <- calculate_scaled_ecosystem_richness(
  input_dir = path_to_predictions,

  # Protected areas are used as the reference condition.
  reference_shape = protected_regions,

  # Biogeographic regions are used for regional summaries.
  categorisation_shape = bio_regions,

  ecosystem_type = ecosystemType,

  # Column in the prediction outputs containing mean occupancy probability.
  value_col = indicator_col,
  region_col = "nature_region",

  # Can set 0 as the lower bound if 0 refers to the worst condition
  set_lower_bound = 0,

  # Aggregation grid spacing in metres.
  point_spacing = aggregation_spacing,

  # Restrict calculations to forest pixels when modelling forest ecosystems.
  forest_only_strict = ecosystemType == "Forest",

  na.rm = TRUE,

  output_dir = path_to_indicators,
  save_outputs = TRUE,
  overwrite = FALSE
)
