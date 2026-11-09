#' Reference quantile for each expert reference category
#'
#' Builds the lookup that maps each expert reference category to the quantile
#' of predicted occupancy used as the species reference value. The defaults
#' are "samme" = 0.9 and "høyere" = 0.99. Supply a named value to change
#' a category, keeping the others at their defaults.
#'
#' @param ... Named numeric values in (0, 1), one per category to override,
#'   e.g. `samme = 0.6`. Names must be existing categories (use backticks for
#'   høyere).
#'
#' @return A data frame with columns `reference_category` and
#'   `reference_quantile`.
#'
#' @examples
#' reference_quantile_lookup()
#' reference_quantile_lookup(samme = 0.6)
#' @export
reference_quantile_lookup <- function(...) {
  quantiles <- c(samme = 0.75, "h\u00F8yere" = 0.95)

  overrides <- list(...)
  if (length(overrides)) {
    if (is.null(names(overrides)) || any(names(overrides) == "")) {
      stop("All arguments must be named by reference category.", call. = FALSE)
    }
    unknown <- setdiff(names(overrides), names(quantiles))
    if (length(unknown)) {
      stop(
        "Unknown reference categories: ",
        paste(unknown, collapse = ", "),
        ". Valid categories: ",
        paste(names(quantiles), collapse = ", "),
        call. = FALSE
      )
    }
    values <- unlist(overrides)
    if (
      !is.numeric(values) || anyNA(values) || any(values <= 0 | values >= 1)
    ) {
      stop(
        "Reference quantiles must be numbers between 0 and 1.",
        call. = FALSE
      )
    }
    quantiles[names(values)] <- values
  }

  data.frame(
    reference_category = names(quantiles),
    reference_quantile = unname(quantiles)
  )
}

#' Get the reference quantile for each expert-assessed species
#'
#' Reads the species list and the expert reference categories for an
#' ecosystem, and converts each category to its fixed reference quantile.
#'
#' @param expert_folder Character. Folder containing the species list
#'   workbook and the `referenceValues/` folder.
#' @param ecosystem_type Character. Ecosystem, `"Forest"` or `"Mountain"`.
#' @param ref_quantiles Data frame with `reference_category` and
#'   `reference_quantile`. Defaults to the fixed `reference_quantile_lookup()`.
#' @param species_file Character. Expert species list workbook
#'   (`Artsutvalg_030626.xlsx`) inside `expert_folder`. It needs a Norwegian
#'   name column (`Norsk ...`), a scientific name column (`Vitenskapelig ...`)
#'   and the ecosystem columns used by [retrieve_expertSpeciesList()]
#'   (`Gammelskog ...` for forest and `Fjell_n\u00F8kkel ...` for mountain,
#'   with the value `"Ja"` for the species of the ecosystem).
#' @param species_sheet Character or integer. Sheet in `species_file`
#'   (`NULL`, the default, is the first sheet).
#' @param expert_file Character. Expert reference workbook inside
#'   `referenceValues/`. Defaults to `<ecosystem>_experts_refs.xlsx`.
#' @param expert_sheet Character. Sheet in `expert_file`.
#'
#' @return A tibble with columns `Species` and `reference_quantile`, one row
#'   per species with an expert reference category.
#'
#' @details Categories are matched after trimming whitespace and lower-casing.
#'   A warning is issued for categories not found in `ref_quantiles` and for
#'   expert species names that are missing from the species list.
#' @export
get_reference_quantiles <- function(
  expert_folder,
  ecosystem_type,
  ref_quantiles = reference_quantile_lookup(),
  species_file = "Artsutvalg_030626.xlsx",
  species_sheet = "Vurderinger",
  expert_file = paste0(tolower(ecosystem_type), "_experts_refs.xlsx"),
  expert_sheet = "Sheet 1"
) {
  ref_dir <- file.path(expert_folder, "referenceValues")

  stopifnot(
    is.data.frame(ref_quantiles),
    all(
      c("reference_category", "reference_quantile") %in% names(ref_quantiles)
    ),
    !anyDuplicated(ref_quantiles$reference_category)
  )

  species_path <- file.path(expert_folder, species_file)
  expert_path <- file.path(ref_dir, expert_file)
  for (p in c(species_path, expert_path)) {
    if (!file.exists(p)) {
      stop("File not found: ", p, call. = FALSE)
    }
  }

  clean <- function(x) tolower(trimws(as.character(x)))

  lookup <- data.frame(
    category_key = clean(ref_quantiles$reference_category),
    reference_quantile = ref_quantiles$reference_quantile
  )

  # ---- Species of the ecosystem from the expert species list -----------------
  species_all <- readxl::read_excel(species_path, sheet = species_sheet)

  find_col <- function(pattern, what) {
    col <- grep(pattern, names(species_all), value = TRUE, ignore.case = TRUE)
    if (length(col) == 0) {
      stop(
        "No ",
        what,
        " column (matching '",
        pattern,
        "') in ",
        species_file,
        ". Columns: ",
        paste(names(species_all), collapse = ", "),
        call. = FALSE
      )
    }
    col[1]
  }

  ecosystem_col <- find_col(
    switch(
      ecosystem_type,
      Forest = "Gammelskog",
      Mountain = "Fjell_n\u00F8kkel",
      stop("`ecosystem_type` must be 'Forest' or 'Mountain'.", call. = FALSE)
    ),
    paste(ecosystem_type, "ecosystem")
  )

  species_selection <- species_all |>
    dplyr::select(
      norwegian_name = dplyr::all_of(find_col(
        "^Populærnavn",
        "Norwegian name"
      )),
      Species = dplyr::all_of(find_col("Vitenskapelig", "scientific name")),
      in_ecosystem = dplyr::all_of(ecosystem_col)
    ) |>
    dplyr::mutate(
      norwegian_name = clean(norwegian_name),
      Species = trimws(Species)
    ) |>
    dplyr::filter(in_ecosystem %in% "Ja") |>
    dplyr::distinct(norwegian_name, Species)

  # ---- Expert reference categories ---------------------------------------
  expert_refs <- readxl::read_excel(expert_path, sheet = expert_sheet) |>
    dplyr::select(
      norwegian_name = `norsk navn`,
      reference_category = Referanse
    ) |>
    dplyr::mutate(
      norwegian_name = clean(norwegian_name),
      category_key = clean(reference_category)
    ) |>
    # Species without an expert category ("NA" text or empty cell)
    dplyr::filter(!is.na(category_key), !category_key %in% c("na", ""))

  unknown_categories <- setdiff(expert_refs$category_key, lookup$category_key)
  if (length(unknown_categories)) {
    warning(
      "Reference categories not in `ref_quantiles` (dropped): ",
      paste(unknown_categories, collapse = ", "),
      call. = FALSE
    )
  }

  missing_species <- setdiff(
    expert_refs$norwegian_name,
    species_selection$norwegian_name
  )
  if (length(missing_species)) {
    warning(
      "Expert species not found in the species list for ",
      ecosystem_type,
      " in ",
      species_file,
      " (dropped): ",
      paste(missing_species, collapse = ", "),
      call. = FALSE
    )
  }

  expert_refs |>
    dplyr::inner_join(species_selection, by = "norwegian_name") |>
    dplyr::inner_join(lookup, by = "category_key") |>
    dplyr::select(Species, reference_quantile)
}
