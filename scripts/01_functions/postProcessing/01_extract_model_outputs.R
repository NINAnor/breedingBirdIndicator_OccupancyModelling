#' Extract model outputs from a directory of RDS result files
#'
#' This function reads a set of model result `.rds` files from a directory,
#' extracts standardized model outputs, and returns them as a structured list
#' of data frames with consistent metadata attached.
#'
#' The function supports both single-species and multi-species model outputs.
#' It automatically detects model type per file via `res$single_species` and
#' selects the appropriate summary tables.
#'
#' @param path_to_results Character. Path to directory containing `.rds` result files.
#'
#' @return A named list containing:
#' \describe{
#'   \item{bayes_P}{Data frame of posterior predictive checks (species-level probabilities).}
#'   \item{beta_summary}{Data frame of beta parameter summaries (species or community level).}
#'   \item{alpha_summary}{Data frame of alpha parameter summaries (species or community level).}
#'   \item{occ_effects}{Data frame of occupancy-related effects.}
#'   \item{det_effects}{Data frame of detection-related effects.}
#'   \item{meta}{Data frame containing file-level metadata (file name, model type, spatial flag).}
#' }
#'
#' @details
#' Each `.rds` file is expected to contain a list with at least:
#' \itemize{
#'   \item `single_species` (logical)
#'   \item `fit_spatial_model` (logical or scalar indicator)
#'   \item `model_products` (list of model output tables)
#' }
#'
#' For single-species models, the function extracts:
#' \itemize{
#'   \item `beta_summary`
#'   \item `alpha_summary`
#' }
#'
#' For multi-species models, it extracts:
#' \itemize{
#'   \item `beta_comm_summary`
#'   \item `alpha_comm_summary`
#' }
#'
#' All outputs are row-bound across files using `dplyr::bind_rows()`.
#'
#' @importFrom dplyr select mutate bind_rows
#' @importFrom purrr map compact map_chr map_lgl
#'
#' @export
extract_model_outputs <- function(paths_to_results) {

  # ------------------------------------------------------------
  # 1. Expand directories into files
  # ------------------------------------------------------------
  result_files <- unlist(
    lapply(paths_to_results, list.files, full.names = TRUE),
    use.names = FALSE
  )

  # ------------------------------------------------------------
  # 2. Process each file
  # ------------------------------------------------------------
  all_results <- purrr::map(result_files, function(f) {

    res <- readRDS(f)

    if (is.null(res$model_products)) return(NULL)

    is_single <- isTRUE(res$single_species)
    run_path <- dirname(f)
    run_file <- basename(f)

    # ------------------------------------------------------------
    # helper: attach metadata to any table
    # ------------------------------------------------------------
    tag <- function(df) {
      if (is.null(df) || nrow(df) == 0) return(NULL)

      df$file <- run_file
      df$path <- run_path
      df$single_species <- is_single
      df$spatial <- res$fit_spatial_model

      df
    }

    # ------------------------------------------------------------
    # 3. Extract tables safely
    # ------------------------------------------------------------

    bayes_P <- tag(
      res$model_products$ppc_species %>%
        dplyr::mutate(Mean = bayesian_p) %>%
        dplyr::select(Species, Mean)
    )

    beta_summary <- if (is_single) {
      tag(res$model_products$beta_summary)
    } else {
      tag(
        res$model_products$beta_comm_summary %>%
          dplyr::mutate(covariate = row.names(.))
      )
    }

    alpha_summary <- if (is_single) {
      tag(res$model_products$alpha_summary)
    } else {
      tag(
        res$model_products$alpha_comm_summary %>%
          dplyr::mutate(covariate = row.names(.))
      )
    }

    occ_effects <- tag(res$model_products$occ_effects)
    det_effects <- tag(res$model_products$det_effects)

    # ------------------------------------------------------------
    # 4. Return structured object per file
    # ------------------------------------------------------------
    list(
      file = run_file,
      path = run_path,
      single_species = is_single,
      spatial = res$fit_spatial_model,

      bayes_P = bayes_P,
      beta_summary = beta_summary,
      alpha_summary = alpha_summary,
      occ_effects = occ_effects,
      det_effects = det_effects
    )
  })

  all_results <- purrr::compact(all_results)

  # ------------------------------------------------------------
  # 5. Safe binder
  # ------------------------------------------------------------
  safe_bind <- function(x, key) {
    dplyr::bind_rows(purrr::map(x, key))
  }

  # ------------------------------------------------------------
  # 6. Final output
  # ------------------------------------------------------------
  output <- list(
    bayes_P = safe_bind(all_results, "bayes_P"),
    beta_summary = safe_bind(all_results, "beta_summary"),
    alpha_summary = safe_bind(all_results, "alpha_summary"),
    occ_effects = safe_bind(all_results, "occ_effects"),
    det_effects = safe_bind(all_results, "det_effects"),

    meta = dplyr::tibble(
      file = purrr::map_chr(all_results, "file"),
      path = purrr::map_chr(all_results, "path"),
      spatial = purrr::map_lgl(all_results, "spatial"),
      single_species = purrr::map_lgl(all_results, "single_species")
    )
  )

  class(output) <- "model_outputs"

  output
}
