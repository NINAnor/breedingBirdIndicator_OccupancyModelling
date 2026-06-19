# Lets the load the biogeopraphic region
# The file path is from NINA's R drive
# R:\GeoSpatialData\BiogeographicalRegions\Norway_Biogeography\Original
# We load and crop to the boundary of Norway
input_file <- file.path(path_finder("ninaGeospatialFiles"),
                        "BiogeographicalRegions",
                        "Norway_Biogeography",
                        "Original",
                        "biogeografiske_regioner_omrade.shp")

bio_regions <- get_biogeographic_region(
  input_file,
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

input_file <- file.path(path_finder("ninaGeospatialFiles"),
                        "ProtectedSites",
                        "Norway_ProtectedAreas",
                        "Original",
                        "Natur_Vern_norge_med_svalbard_25833.gdb")

protectedRegions <- get_protected_areas(input_file,
                                        norway_boundary,
                                        modelCrs,
                                        output_file = file.path(path_finder("dataFolder"),
                                                                "processedData",
                                                                "protected_regions",
                                                                "protected_region.parquet"),
                                        layer = "Naturverngrense",
                                        overwrite = FALSE)
