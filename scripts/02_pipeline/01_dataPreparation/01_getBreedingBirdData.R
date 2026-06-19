#### Load breeding bird data for selected ecosystem ####

# This script:
# 1. Connects to the NINA SQL Server database if available.
# 2. Falls back to cached local files if the database is unavailable.
# 3. Loads the expert-defined species list for the selected ecosystem.
# 4. Filters breeding bird count data to the selected species and model years.


#### 1. Connect to database ####

# Optional: inspect available ODBC drivers on the current machine.
# Useful for debugging database connection problems.
sort(unique(odbc::odbcListDrivers()[[1]]))

con <- tryCatch(
  {
    DBI::dbConnect(
      odbc::odbc(),
      Driver   = "SQL Server",
      Server   = "ninsql07.nina.no",
      Database = "TOVTaksering",
      Trusted_Connection = "True"
    )
  },
  error = function(e) {
    warning(
      paste(
        "Unable to connect to the TOVTaksering database.",
        "Proceeding with cached local files instead.",
        "",
        paste("Original database error:", e$message),
        sep = "\n"
      ),
      call. = FALSE
    )

    NULL
  }
)


#### 2. Define input paths ####

# Folder where processed point-count data are stored.
point_count_dir <- file.path(path_finder("dataFolder", user = USER), "pointCounts")

# Cached processed breeding bird count file.
breeding_bird_cache <- file.path(
  point_count_dir,
  "breeding_bird_counts_clean.csv"
)

# File or folder containing expert species information.
# Adjust this path if expert data are stored elsewhere.
expert_data_path <- path_finder("expertData", user = USER)


#### 3. Load breeding bird count data ####

if (!file.exists(breeding_bird_cache)) {

  message("Cached breeding bird data not found.")
  message("Retrieving breeding bird routes and counts from database...")

  breedingBirdData <- retrieve_breedingBirds_routesAndCounts(
    con             = con,
    dataFolder      = point_count_dir,
    overwrite_cache = FALSE,
    save_outputs    = TRUE,
    crs             = modelCrs
  )

} else {

  message("Loading cached breeding bird data from: ", breeding_bird_cache)

  breedingBirdData <- readr::read_csv(
    breeding_bird_cache,
    show_col_types = FALSE
  )
}


#### 4. Load expert species list ####

expertSpeciesData <- retrieve_expertSpeciesList(
  path_to_file = expert_data_path,
  ecosystem    = ecosystemType
)

expertSpecies <- expertSpeciesData$expertSpeciesList$Species

message(
  "Number of expert-selected species for ",
  ecosystemType,
  ": ",
  length(expertSpecies)
)


#### 5. Filter breeding bird data ####

ecosystemBreedingBird <- breedingBirdData %>%
  dplyr::filter(Species %in% expertSpecies) %>%
  dplyr::filter(
    Year >= modelYear[1],
    Year <= modelYear[2]
  )

message(
  "Number of records for ",
  ecosystemType,
  " between ",
  modelYear[1],
  " and ",
  modelYear[2],
  ": ",
  nrow(ecosystemBreedingBird)
)

message(
  "Number of species retained in filtered data: ",
  dplyr::n_distinct(ecosystemBreedingBird$Species)
)

