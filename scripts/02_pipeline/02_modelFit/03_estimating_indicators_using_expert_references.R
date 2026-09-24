# Obtain the reference values from experts

if (ecosystemType %in% "Mountain") {
  species_reference_vals <- readr::read_csv2(
    file.path(
      expertFolder,
      "referenceValues",
      "mountain_reference_vals.csv"
    )
  ) %>%
    dplyr::mutate(Values = as.numeric(Values))
} else if (ecosystemType %in% "Forest") {
  species_reference_vals <- readr::read_csv2(
    file.path(
      expertFolder,
      "referenceValues",
      "forest_reference_vals.csv"
    )
  ) %>%
    dplyr::mutate(Values = as.numeric(Values))
}

# Estimate scaled indicators
gdb <- paste0(
  geoSpatFolder,
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


scaled_richness <- calculate_scaled_ecosystem_richness_with_expert_data(
  input_dir = path_to_predictions,

  # Protected areas are used as the reference condition.
  reference_data = species_reference_vals,

  # Biogeographic regions are used for regional summaries.
  categorisation_shape = bio_regions,

  ecosystem_type = ecosystemType,

  # Column in the prediction outputs containing mean occupancy probability.
  value_col = indicator_col,
  region_col = "nature_region",

  output_dir = path_to_indicators,
  save_outputs = TRUE,
  overwrite = TRUE,

  # Bootstrap results
  rep_values_for_each_region = TRUE,
  n.boot = 1000,
  boot.seed = 123,
  n.Samps = 20000
)
