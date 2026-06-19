retrieve_expertSpeciesList <- function(
    path_to_file = path_finder("expertData"),
    ecosystem = c("Forest", "Mountain"),
    file_name = "Artsutvalg_030626.xlsx",
    save_outputs = TRUE
) {

  ###########################################################################
  # Function: retrieve_expertSpeciesList
  #
  # Purpose:
  #   Read the expert-defined bird species list used in FIAR analyses and
  #   return species associated with a selected ecosystem.
  #
  # FAIR compliance:
  #   Findable:
  #     - Uses explicit filenames.
  #     - Creates clear output names.
  #
  #   Accessible:
  #     - Reads from a user-defined folder.
  #     - Optionally saves cleaned outputs to disk.
  #
  #   Interoperable:
  #     - Standardizes column names to English.
  #     - Uses scientific species names.
  #
  #   Reusable:
  #     - Checks required packages and file existence.
  #     - Trims species names.
  #     - Removes duplicate species records.
  #     - Creates field-level metadata.
  ###########################################################################

  ###########################################################################
  # 1. Install and load required packages
  ###########################################################################

  required_pkgs <- c("readxl", "dplyr", "stringr", "readr", "tibble")

  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs, dependencies = TRUE)
  }

  invisible(lapply(required_pkgs, library, character.only = TRUE))

  ###########################################################################
  # 2. Validate function inputs
  ###########################################################################

  ecosystem <- match.arg(ecosystem)

  if (missing(path_to_file) || is.null(path_to_file)) {
    stop("'path_to_file' must be supplied.")
  }

  if (!dir.exists(path_to_file)) {
    stop("Folder does not exist: ", path_to_file)
  }

  if (!is.logical(save_outputs)) {
    stop("'save_outputs' must be TRUE or FALSE.")
  }

  file_path <- file.path(path_to_file, file_name)

  if (!file.exists(file_path)) {
    stop("Expert species file not found: ", file_path)
  }

  ###########################################################################
  # 3. Read and standardize expert species list
  #
  # Source file:
  #   Artsutvalg_030626.xlsx
  #
  # Original columns:
  #   Vitenskapelig navn = scientific species name
  #   Skog               = forest
  #   Fjell              = mountain
  ###########################################################################

  # expertData <- readxl::read_excel(file_path) |>
  #   dplyr::select(
  #     Species = `Vitenskapelig navn`,
  #     Forest = Skog,
  #     Mountain = Fjell
  #   ) |>
  #   dplyr::mutate(
  #     Species = stringr::str_trim(Species),
  #     Forest = dplyr::na_if(Forest, ""),
  #     Mountain = dplyr::na_if(Mountain, "")
  #   ) |>
  #   dplyr::filter(
  #     !is.na(.data[[ecosystem]])
  #   ) |>
  #   dplyr::distinct(
  #     Species,
  #     .keep_all = TRUE
  #   )

  tmp <- readxl::read_excel(file_path)

  species_col <- grep("Vitenskapelig", names(tmp), value = TRUE)
  ecosystem_col <- grep("Flerartsindeks", names(tmp), value = TRUE)

  expertData <- tmp |>
    dplyr::select(
      Species = dplyr::all_of(species_col),
      ecosystemType = dplyr::all_of(ecosystem_col)
    ) |>
    dplyr::filter(ecosystemType != "Nei") |>
    dplyr::mutate(
      # Now recode using the Norwegian name
      ecosystemType = dplyr::recode(
        ecosystemType,
        "FJELL" = "Mountain",
        "SKOG" = "Forest",
        "VÅTMARK" = "Wetland",
        "JORDBRUK" = "Farmland"
      )
    ) |>
    dplyr::filter(
     ecosystemType %in% ecosystem
    ) |>
    dplyr::distinct(
       Species,
        .keep_all = TRUE
      )

  ###########################################################################
  # 4. Quality-control checks
  ###########################################################################

  if (any(is.na(expertData$Species))) {
    warning("Some records have missing scientific species names.")
  }

  if (nrow(expertData) == 0) {
    warning("No species found for ecosystem: ", ecosystem)
  }

  duplicated_species <- expertData$Species[duplicated(expertData$Species)]

  if (length(duplicated_species) > 0) {
    warning(
      "Duplicate species detected after filtering: ",
      paste(unique(duplicated_species), collapse = ", ")
    )
  }

  ###########################################################################
  # 5. Create FAIR metadata
  ###########################################################################

  field_metadata <- tibble::tibble(
    field = names(expertData),
    description = c(
      "Scientific species name.",
      "Expert classification ecosystem."
    ),
    source_field = c(
      "Vitenskapelig navn",
      "Flerartsindeks_Norsk_hekkefuglovervåking"
    )
  )

  dataset_metadata <- tibble::tibble(
    metadata_field = c(
      "dataset_name",
      "source_file",
      "selected_ecosystem",
      "processing_date",
      "processed_by",
      "output_file",
      "script_function"
    ),
    value = c(
      "Expert-defined FIAR bird species list",
      normalizePath(file_path, winslash = "/", mustWork = FALSE),
      ecosystem,
      as.character(Sys.time()),
      Sys.info()[["user"]],
      paste0("expert_species_list_", tolower(ecosystem), ".csv"),
      "retrieve_expertSpeciesList"
    )
  )

  ###########################################################################
  # 6. Save FAIR outputs
  ###########################################################################

  if (save_outputs) {

    readr::write_csv(
      expertData,
      file.path(
        path_to_file,
        "processedData",
        paste0("expert_species_list_", tolower(ecosystem), ".csv")
      )
    )

    readr::write_csv(
      field_metadata,
      file.path(
        path_to_file,
        "processedData",
        paste0("expert_species_list_", tolower(ecosystem), "_field_metadata.csv")
      )
    )

    readr::write_csv(
      dataset_metadata,
      file.path(
        path_to_file,
        "processedData",
        paste0("expert_species_list_", tolower(ecosystem), "_dataset_metadata.csv")
      )
    )
  }

  ###########################################################################
  # 7. Return cleaned data and metadata
  ###########################################################################

  return(
    list(
      expertSpeciesList = expertData,
      field_metadata = field_metadata,
      dataset_metadata = dataset_metadata
    )
  )
}
