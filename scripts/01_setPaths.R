# Define a function called 'project'
# This function acts as a centralized directory manager.
# Instead of hard-coding file paths throughout multiple scripts,
# all important project folders are stored in one location.
path_finder <- function(name = NULL, user = "all") {
  # Create a named vector of relative folder paths.
  # The names on the left are aliases that can be used
  # elsewhere in the analysis code.
  if (user %in% "all") {
    paths <- c(
      # Folder containing survey route information and
      # environmental covariates used as predictor variables.
      dataFolder = "../../03_Survey_routes_and_covariates",

      # Folder containing expert-derived bird habitat information,
      # species-habitat associations, or key species lists.
      expertData = "../../02_Key_bird_species_by_habitat",

      # Folder where model outputs, predictions,
      # diagnostics, and summaries will be stored.
      results = "../modelResults",

      # Folder where functions are stores
      functions = "scripts/01_functions",

      ninaGeospatialFiles = "R:/GeoSpatialData"
    )
  } else if (user %in% "server") {
    paths <- c(
      # Folder containing survey route information and
      # environmental covariates used as predictor variables.
      dataFolder = "../03_Survey_routes_and_covariates",

      # Folder containing expert-derived bird habitat information,
      # species-habitat associations, or key species lists.
      expertData = "../02_Key_bird_species_by_habitat",

      # Folder where model outputs, predictions,
      # diagnostics, and summaries will be stored.
      results = "../breedingBirdsModelOutput",

      # Folder where functions are stores
      functions = "scripts/01_functions",

      ninaGeospatialFiles = "/data/R/GeoSpatialData"
    )
  } else if (user %in% "kwaku") {
    paths <- c(
      # Folder containing survey route information and
      # environmental covariates used as predictor variables.
      dataFolder = "~/OneDrive - NINA/112635 - Videreutvikling av indikatorer for hekkefugl - Documents/03_Survey_routes_and_covariates",

      # Folder containing expert-derived bird habitat information,
      # species-habitat associations, or key species lists.
      expertData = "~/OneDrive - NINA/112635 - Videreutvikling av indikatorer for hekkefugl - Documents/02_Key_bird_species_by_habitat",

      # Folder where model outputs, predictions,
      # diagnostics, and summaries will be stored.
      #results = "../../ninaProjects/breedingBirdsModelOutput",
      results = "~/OneDrive - NINA/112635 - Videreutvikling av indikatorer for hekkefugl - Documents/ninaProjects/breedingBirdsModelOutput",

      # Folder where functions are stores
      functions = "scripts/01_functions",

      ninaGeospatialFiles = "R:/GeoSpatialData"
    )
  }
  # If no folder name is supplied,
  # return the complete list of project paths.
  if (is.null(name)) {
    return(paths)
  }

  # Otherwise return only the requested path.
  # Example:
  # project("dataFolder")
  # project("results")
  paths[[name]]
}
