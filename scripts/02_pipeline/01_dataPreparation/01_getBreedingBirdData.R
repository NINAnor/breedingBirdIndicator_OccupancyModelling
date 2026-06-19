#### Connecting to the DB ####
sort(unique(odbcListDrivers()[[1]]))

con <- tryCatch(

  DBI::dbConnect(
    odbc::odbc(),
    Driver   = "SQL Server",
    Server   = "ninsql07.nina.no",
    Database = "TOVTaksering",
   # UID      = username,
  #  PWD      = password,
    Trusted_Connection = "True"
  ),

  error = function(e) {

    warning(
      "Unable to connect to TOVTaksering database.\n",
      "Proceeding with cached local files.\n",
      "Original error: ", e$message
    )

    NULL
  }

)

# Load the processed breeding bird data
if(!file.exists(file.path(path_finder("dataFolder"), "pointCounts", "breeding_bird_counts_clean.csv"))){
  # This dataset contains the breedingBirdData, breedingbirdData in sf format, the routes (sf format) and the
  # meta data associated with the field and dataset.

breedingBirdData <- retrieve_breedingBirds_routesAndCounts(
    con = con,
    dataFolder = file.path(path_finder("dataFolder"), "pointCounts"),
    overwrite_cache = FALSE,
    save_outputs = TRUE,
    crs = modelCrs
)
} else {
  breedingBirdData <- readr::read_csv(file.path(path_finder("dataFolder"),
                                                "pointCounts",
                                                "breeding_bird_counts_clean.csv"))
}

# Get the expert species names for the ecosystem type
expertSpecies <- retrieve_expertSpeciesList(path_to_file = path_finder("expertData"),
                                             ecosystem = ecosystemType)


expertSpecies <- expertSpecies$expertSpeciesList$Species

cat("Number of species for ", ecosystemType, "=", length(expertSpecies))

# Filter the breeding bird data for these species
ecosystemBreedingBird <- breedingBirdData %>%
  dplyr::filter(Species %in% expertSpecies) %>%
  dplyr::filter(Year >= modelYear[1] & Year <= modelYear[2])

cat("Number of species for ", ecosystemType, "for selected years=", length(expertSpecies))




