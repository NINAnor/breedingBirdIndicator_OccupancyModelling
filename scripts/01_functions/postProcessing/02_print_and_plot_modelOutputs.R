#' Print Method for Model Outputs Objects
#'
#' Provides a compact summary of a \code{model_outputs} object, including
#' the number of processed files, model configurations, and available components.
#'
#' This method is primarily intended for quick inspection of extracted model
#' results and supports both single-model and multi-model workflows.
#'
#' @param x An object of class \code{model_outputs}, typically created by
#'   \code{extract_model_outputs()}.
#'
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object \code{x}, allowing the function
#'   to be used in pipelines without interruption.
#'
#' @details
#' The print method summarises:
#' \itemize{
#'   \item Number of processed model files
#'   \item Number of unique model paths
#'   \item Distribution of single-species vs multi-species models
#'   \item Presence of spatial model configurations
#'   \item Availability and size of each extracted component
#' }
#'
#' It also prints the unique model paths included in the object, which is
#' useful for verifying that all intended model configurations were loaded.
#'
#' @section Components Reported:
#' The following components are summarised if present:
#' \itemize{
#'   \item \code{bayes_P}: Posterior predictive checks
#'   \item \code{beta_summary}: Community or species-level beta parameters
#'   \item \code{alpha_summary}: Community or species-level alpha parameters
#'   \item \code{occ_effects}: Occupancy covariate effects
#'   \item \code{det_effects}: Detection covariate effects
#' }
#'
#' Missing components are explicitly flagged as \code{not available}.
#'
#' @seealso
#' \code{\link{extract_model_outputs}}, \code{\link{plot.model_outputs}}
#'
#' @examples
#' \dontrun{
#'
#' outputs <- extract_model_outputs(path_to_results)
#'
#' print(outputs)
#'
#' # or simply
#' outputs
#'
#' }
#'
#' @export
print.model_outputs <- function(x, ...) {

  cat("\n<Model Outputs Object>\n")
  cat("----------------------\n\n")

  cat("Files processed:", nrow(x$meta), "\n")
  cat("Unique paths:", length(unique(x$meta$path)), "\n")
  cat("Single species:", sum(x$meta$single_species), "\n")
  cat("Multi species:", sum(!x$meta$single_species), "\n")
  cat("Spatial models:", sum(x$meta$spatial, na.rm = TRUE), "\n\n")

  cat("Components:\n")

  components <- c("bayes_P", "beta_summary", "alpha_summary",
                  "occ_effects", "det_effects")

  for (comp in components) {
    if (!is.null(x[[comp]])) {
      cat(" -", comp, ":", nrow(x[[comp]]), "rows\n")
    } else {
      cat(" -", comp, ": not available\n")
    }
  }

  cat("\nPaths:\n")
  print(unique(x$meta$path))

  invisible(x)
}
#' Plot Model Outputs from Extracted Occupancy Model Fits
#'
#' S3 plotting method for objects of class \code{model_outputs}. This function
#' provides standardized visualizations of model-derived quantities across
#' one or more fitted models, including occupancy, detection, and community-level
#' summaries.
#'
#' The function supports multiple plot types:
#' \itemize{
#'   \item \code{"error"}: Posterior means with 95\% credible intervals
#'   \item \code{"violin"}: Distribution of posterior estimates
#'   \item \code{"rank"}: Ranked effect sizes by magnitude
#' }
#'
#' When multiple model runs are present, results are automatically grouped
#' and color-coded by run configuration (e.g., species structure and spatial model).
#'
#' @param object An object of class \code{model_outputs}, typically returned
#'   by \code{extract_model_outputs()}.
#'
#' @param type Character. Type of plot to generate. One of:
#' \code{"error"}, \code{"violin"}, or \code{"rank"}.
#'
#' @param components Character vector specifying which model components to include
#'   in the plot. If \code{NULL}, all available components are used. Possible values include:
#'   \code{"occ_effects"}, \code{"det_effects"}, \code{"bayes_P"},
#'   \code{"beta_summary"}, and \code{"alpha_summary"} (when available).
#'
#' @return A \code{ggplot} object. The returned plot can be further modified
#'   using standard \pkg{ggplot2} syntax.
#'
#' @details
#' This function automatically:
#' \itemize{
#'   \item Harmonises column names across different model outputs
#'   \item Extracts run labels from model file paths using \code{extract_run_label()}
#'   \item Groups results by run configuration (species structure + spatial model)
#'   \item Applies consistent faceting across components and parameters
#' }
#'
#' For multi-model comparisons, color and dodge positions are enabled
#' automatically when more than one run is detected.
#'
#' The function is designed for hierarchical occupancy/detection models
#' where outputs may include both species-level and community-level summaries.
#'
#' @section Plot Types:
#' \describe{
#'   \item{\code{error}}{
#'     Displays posterior mean estimates with 95\% credible intervals.
#'     Useful for comparing effect sizes across species, covariates, or model runs.
#'   }
#'
#'   \item{\code{violin}}{
#'     Shows the distribution of posterior estimates across runs or parameters.
#'     Useful for assessing variability and uncertainty structure.
#'   }
#'
#'   \item{\code{rank}}{
#'     Ranks parameters by mean absolute effect size.
#'     Useful for identifying the most influential covariates.
#'   }
#' }
#'
#' @seealso
#' \code{\link{extract_model_outputs}}, \code{\link{extract_run_label}},
#' \code{\link{build_model_paths}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbarh geom_violin
#' @importFrom ggplot2 geom_boxplot geom_col facet_grid facet_wrap
#' @importFrom dplyr mutate filter group_by summarise arrange
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' outputs <- extract_model_outputs(path_to_results)
#'
#' # Occupancy and detection effects
#' plot(outputs, type = "error", components = c("occ_effects", "det_effects"))
#'
#' # Posterior predictive checks
#' plot(outputs, type = "error", components = "bayes_P")
#'
#' # Community summaries
#' plot(outputs, type = "error",
#'      components = c("beta_summary", "alpha_summary"))
#'
#' # Violin plot
#' plot(outputs, type = "violin")
#'
#' # Rank plot
#' plot(outputs, type = "rank")
#'
#' }
plot.model_outputs <- function(
    object,
    type = c("error", "violin", "rank"),
    components = NULL
) {

  type <- match.arg(type)

  meta <- object$meta
  is_multi <- any(meta$single_species == FALSE)

  # ------------------------------------------------------------
  # 1. COMPONENTS
  # ------------------------------------------------------------
  available_components <- list(
    occ_effects = object$occ_effects,
    det_effects = object$det_effects,
    bayes_P = object$bayes_P
  )

  if (is_multi) {
    available_components$beta_summary <- object$beta_summary
    available_components$alpha_summary <- object$alpha_summary
  }

  if (is.null(components)) {
    components <- names(available_components)
  }

  # ------------------------------------------------------------
  # 2. BUILD DATAFRAME (FIXED PATH LOGIC)
  # ------------------------------------------------------------



  df <- purrr::imap_dfr(
        available_components[components],
        function(d, name) {

          if (is.null(d) || nrow(d) == 0) return(NULL)

          d$component <- name

          # CRITICAL: preserve run identity
          #d$path <- run_path
          d$run_label <- c(sapply(d$path, extract_run_label))

          d
        }
      )


  if (is.null(df) || nrow(df) == 0) {
    stop("No data available for plotting")
  }

  # ------------------------------------------------------------
  # 3. STANDARDISE COLUMNS
  # ------------------------------------------------------------
  standardise <- function(df) {

    if ("Mean" %in% names(df) && !"estimate" %in% names(df)) {
      df$estimate <- df$Mean
    }

    if ("2.5%.2.5%" %in% names(df) && !"lower" %in% names(df)) {
      df$lower <- df$`2.5%.2.5%`
    }

    if ("97.5%.97.5%" %in% names(df) && !"upper" %in% names(df)) {
      df$upper <- df$`97.5%.97.5%`
    }

    required <- c("estimate", "lower", "upper", "covariate", "level", "Species")

    for (col in required) {
      if (!col %in% names(df)) df[[col]] <- NA
    }

    df$estimate <- as.numeric(df$estimate)
    df$lower <- as.numeric(df$lower)
    df$upper <- as.numeric(df$upper)

    df
  }

  df <- standardise(df)

  # ------------------------------------------------------------
  # 4. PATH LABEL (FIXED + SIMPLE)
  # ------------------------------------------------------------
  df$path_label <- df$run_label

  # Only collapse if truly single run
  if (length(unique(df$path)) == 1) {
    df$path_label <- "run"
  }

  df$path_label <- factor(df$path_label)

  # ------------------------------------------------------------
  # 5. PARAMETER CONSTRUCTION
  # ------------------------------------------------------------
  df$parameter <- dplyr::case_when(
    !is.na(df$level) ~ paste0(df$covariate, " = ", df$level),
    !is.na(df$covariate) ~ df$covariate,
    TRUE ~ df$component
  )

  # ------------------------------------------------------------
  # 6. PATH LOGIC
  # ------------------------------------------------------------
  n_paths <- length(levels(df$path_label))
  use_colors <- n_paths > 1
  use_dodge <- n_paths > 1

  pd <- ggplot2::position_dodge(width = 0.6)

  if (use_colors) {
    color_scale <- ggplot2::scale_color_brewer(palette = "Set1")
  } else {
    color_scale <- ggplot2::scale_color_manual(values = "black")
  }

  # ============================================================
  # ERROR PLOT
  # ============================================================
  if (type == "error") {

    is_summary <- any(df$component %in% c("alpha_summary", "beta_summary"))

    if(any(meta$single_species == TRUE) && is_summary){
      df <- df %>%
        filter(single_species == FALSE)
    }

    facet_formula <- if (is_summary) {
      stats::as.formula("~ component")
    } else {
      stats::as.formula("component ~ parameter")
    }

    y_var <- if (is_summary) "covariate" else "Species"

    p <- ggplot2::ggplot(df, ggplot2::aes(
      y = .data[[y_var]],
      x = estimate,
      color = path_label,
      group = interaction(path_label, Species)
    ))

    if (use_dodge) {
      p <- p +
        ggplot2::geom_point(position = pd) +
        ggplot2::geom_errorbarh(
          ggplot2::aes(xmin = lower, xmax = upper),
          position = pd,
          height = 0.2
        )
    } else {
      p <- p +
        ggplot2::geom_point() +
        ggplot2::geom_errorbarh(
          ggplot2::aes(xmin = lower, xmax = upper),
          height = 0.2,
          color = "black"
        )
    }

    return(
      p +
        ggplot2::facet_grid(facet_formula) +
        color_scale +
        ggplot2::labs(
          y = if (is_summary) "Covariate" else "Species",
          x = "Posterior mean (95% CI)",
          color = if (use_colors) "Run" else NULL
        ) +
        ggplot2::guides(color = if (use_colors) "legend" else "none")
    )
  }

  # ============================================================
  # VIOLIN PLOT
  # ============================================================
  if (type == "violin") {

    ggplot2::ggplot(df, ggplot2::aes(
      x = component,
      y = estimate,
      fill = path_label
    )) +
      ggplot2::geom_violin(
        position = ggplot2::position_dodge(0.8),
        trim = FALSE
      ) +
      ggplot2::geom_boxplot(
        width = 0.15,
        position = ggplot2::position_dodge(0.8),
        outlier.size = 0.5
      ) +
      ggplot2::facet_wrap(~ parameter, scales = "free_y") +
      color_scale +
      ggplot2::theme_minimal()
  }

  # ============================================================
  # RANK PLOT
  # ============================================================
  if (type == "rank") {

    df_rank <- df %>%
      dplyr::group_by(parameter, component, path_label) %>%
      dplyr::summarise(
        mean_est = mean(estimate, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(abs(mean_est))) %>%
      dplyr::mutate(
        parameter = factor(parameter, levels = unique(parameter))
      )

    ggplot2::ggplot(df_rank, ggplot2::aes(
      x = parameter,
      y = mean_est,
      fill = path_label
    )) +
      ggplot2::geom_col(
        position = if (use_dodge) pd else "identity"
      ) +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~ component, scales = "free_y") +
      color_scale +
      ggplot2::theme_minimal()
  }
}
