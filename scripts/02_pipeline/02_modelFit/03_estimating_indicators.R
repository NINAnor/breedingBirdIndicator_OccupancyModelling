# Lets the load the biogeopraphic region
# The file path is from NINA's R drive
# R:\GeoSpatialData\BiogeographicalRegions\Norway_Biogeography\Original
# We load and crop to the boundary of Norway

bio_regions <- get_biogeographic_region(
  input_file = file.path(path_finder("ninaGeospatialFiles"),
                         "BiogeographicalRegions",
                         "Norway_Biogeography",
                         "Original",
                         "biogeografiske_regioner_omrade.shp"),
  norway_boundary,
  modelCrs,
  output_file = file.path(path_finder("dataFolder"),
                          "processedData",
                          "biogeopragic_regions",
                          "biogeographic_region.parquet"),
  overwrite = FALSE
)


# Lets the load the protection regions
# The file path is from NINA's R drive
# R:\GeoSpatialData\ProtectedSites\Norway_ProtectedAreas\Original\Natur_Vern_norge_med_svalbard_25833.gdb
# We load and crop to the boundary of Norway

protectedRegions <- get_protected_areas(input_file = file.path(path_finder("ninaGeospatialFiles"),
                                                               "ProtectedSites",
                                                               "Norway_ProtectedAreas",
                                                               "Original",
                                                               "Natur_Vern_norge_med_svalbard_25833.gdb"),
                                        norway_boundary,
                                        modelCrs,
                                        output_file = file.path(path_finder("dataFolder"),
                                                                "processedData",
                                                                "protected_regions",
                                                                "protected_region.parquet"),
                                        layer = "Naturvernområde",
                                        overwrite = FALSE)


# Estimate indicators
forest_scaled_richness <- calculate_scaled_ecosystem_richness(
  input_dir = file.path(path_finder("results"),
                        "prediction_outputs",
                        ecosystemType),
  reference_shape = protectedRegions,
  categorisation_shape = bio_regions,
  ecosystem_type = ecosystemType,
  value_col = "mean",
  point_spacing = 1000,
  forest_only_strict = TRUE,
  na.rm = TRUE,
  output_dir = file.path(path_finder("results"),
                         "indicator_outputs",
                         ecosystemType),
  save_outputs = TRUE,
  overwrite = FALSE
)
