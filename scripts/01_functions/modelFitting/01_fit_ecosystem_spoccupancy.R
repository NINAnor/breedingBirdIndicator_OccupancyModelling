#' Fit Occupancy Models by Ecosystem and Species or Guild
#'
#' Fits Bayesian occupancy models using \pkg{spOccupancy} for either
#' single species or groups of species (guilds), with support for both
#' spatial and non-spatial formulations. The function constructs a
#' consistent detection/non-detection dataset, builds model inputs,
#' fits models, extracts posterior summaries, and returns structured outputs.
#'
#' This function supports FAIR principles:
#' \itemize{
#'   \item \strong{Findable}: Outputs are organized by ecosystem, guild, and species.
#'   \item \strong{Accessible}: Results are returned as structured R objects and
#'         optionally saved as `.rds` files.
#'   \item \strong{Interoperable}: Uses standard modeling structures compatible
#'         with \pkg{spOccupancy}, \pkg{coda}, and common ecological workflows.
#'   \item \strong{Reusable}: Ensures consistent survey effort definition, includes
#'         metadata, diagnostics, and derived ecological metrics.
#' }
#'
#' @param ecosystemBreedingBird Data frame containing bird survey data.
#' @param ecosystemType Character. Name of ecosystem (used for labeling outputs).
#' @param maxSpeciesInGrps Integer. Maximum number of species per group if
#'   automatic grouping is used.
#' @param occ_covariates Character vector of occupancy covariate names.
#' @param det_covariates Character vector of detection covariate names.
#' @param species_groupings Optional data frame with columns \code{Species} and
#'   \code{broad_guild} to define species groups.
#' @param site_id_col,species_col,year_col,count_col Character. Column names in
#'   input data defining site ID, species, year, and observation counts.
#' @param lon_col,lat_col Character. Column names for coordinates (longitude, latitude).
#' @param route_id_col,point_id_col Optional character. Explicit route and point identifiers.
#' @param use_route_random_effect Logical. Include random intercept for routes.
#' @param covariates_to_scale Character vector of covariates to scale.
#'   Defaults to all occupancy and detection covariates.
#' @param remove_routes_with_missing_coords Logical. Remove all routes with missing coordinates.
#'
#' @param single_species Logical. If `TRUE`, fits one model per species.
#'   If `FALSE`, fits multi-species models.
#' @param fit_spatial_model Logical. If `TRUE`, fits spatial NNGP models.
#'
#' @param n.batch,batch.length,n.burn,n.thin,n.chains Numeric. MCMC tuning parameters.
#' @param n.neighbors Integer. Number of neighbors for spatial models.
#' @param n.omp.threads Integer. Number of OpenMP threads.
#' @param cov.model Character. Covariance model for spatial models.
#'
#' @param fit_models Logical. If `FALSE`, skips model fitting and returns prepared data.
#' @param save_results Logical. If `TRUE`, saves model outputs to disk.
#' @param verbose Logical. Print fitting progress.
#' @param output_dir Character. Directory for saved outputs.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{ecosystem}: Ecosystem name
#'   \item \code{single_species}: Model type flag
#'   \item \code{fit_spatial_model}: Spatial model flag
#'   \item \code{guild_results}: Nested results per guild or species
#'   \item \code{combined}: Combined summary tables across all groups
#' }
#'
#' The \code{combined} output includes (depending on model type):
#' \itemize{
#'   \item Occupancy probability summaries (\code{psi})
#'   \item Latent occupancy states (\code{z})
#'   \item Spatial random effects (\code{w})
#'   \item Species richness estimates (multi-species only)
#'   \item Posterior predictive checks (PPC)
#'   \item Spatial parameters (\code{sigma_sq}, \code{phi})
#'   \item Covariate effect summaries (occupancy and detection)
#' }
#'
#' @details
#' The function implements a unified modeling workflow:
#' \enumerate{
#'   \item Validates inputs and prepares data.
#'   \item Constructs a global detection/non-detection array
#'         (species x site x year) with consistent survey effort.
#'   \item Removes sites or routes with missing coordinates if specified.
#'   \item Builds occupancy and detection covariate matrices.
#'   \item Optionally scales numeric covariates.
#'   \item Defines model formulas for occupancy and detection processes.
#'   \item Splits species into groups (guilds) if required.
#'   \item Fits models based on selected configuration:
#'     \itemize{
#'       \item Single-species + spatial: \code{spPGOcc()}
#'       \item Single-species + non-spatial: \code{PGOcc()}
#'       \item Multi-species + spatial: \code{spMsPGOcc()}
#'       \item Multi-species + non-spatial: \code{msPGOcc()}
#'     }
#'   \item Extracts posterior summaries and diagnostics:
#'     \itemize{
#'       \item Occupancy probability (\eqn{\\psi})
#'       \item Latent occupancy state (\eqn{z})
#'       \item Spatial random effects (\eqn{w})
#'       \item Species richness (multi-species models)
#'       \item Posterior predictive checks
#'       \item WAIC and convergence diagnostics (R-hat)
#'     }
#'   \item Combines results across species or groups into unified outputs.
#' }
#'
#' Key ecological interpretation:
#' \itemize{
#'   \item \eqn{\\psi}: Probability of occupancy at a site
#'   \item \eqn{z}: Latent true occupancy state
#'   \item \eqn{w}: Spatial variation not explained by covariates
#'   \item Richness: Number of species per site (posterior expectation)
#' }
#'
#' @examples
#' \dontrun{
#' results <- fit_ecosystem_spoccupancy(
#'   ecosystemBreedingBird = bird_data,
#'   ecosystemType = "Forest",
#'   occ_covariates = c("elevation", "forest_cover"),
#'   det_covariates = c("observer_effort"),
#'   single_species = FALSE,
#'   fit_spatial_model = TRUE
#' )
#'
#' results$combined$psi
#' results$combined$richness
#' }
#'
#' @seealso
#' [spOccupancy::spPGOcc()], [spOccupancy::spMsPGOcc()],
#' [spOccupancy::PGOcc()], [spOccupancy::msPGOcc()],
#' [coda::gelman.diag()]
#'
#' @export

# =============================================================================
# GENERAL OCCUPANCY FUNCTION BY ECOSYSTEM AND SPECIES GROUP
# Supports: single-species / multi-species  x  spatial / non-spatial
#
# spOccupancy function routing:
#   single_species = FALSE, fit_spatial_model = TRUE  -> spMsPGOcc  (original)
#   single_species = FALSE, fit_spatial_model = FALSE -> msPGOcc    (original)
#   single_species = TRUE,  fit_spatial_model = TRUE  -> spPGOcc
#   single_species = TRUE,  fit_spatial_model = FALSE -> PGOcc
# =============================================================================

fit_ecosystem_spoccupancy <- function(
    ecosystemBreedingBird = NULL,
    ecosystemType         = NULL,
    maxSpeciesInGrps      = 10,

    occ_covariates = NULL,
    det_covariates = NULL,
    species_groupings = NULL,

    site_id_col  = "Route_point_id",
    species_col  = "Species",
    year_col     = "Year",
    count_col    = "TotalCount",
    lon_col      = "Longitude",
    lat_col      = "Latitude",

    route_id_col = NULL,
    point_id_col = NULL,

    use_route_random_effect = FALSE,

    covariates_to_scale = NULL,
    remove_routes_with_missing_coords = TRUE,

    # --- Model type switches ---
    single_species    = FALSE,   # TRUE -> PGOcc / spPGOcc per species
    fit_spatial_model = TRUE,    # TRUE -> spatial NNGP variant

    n.batch       = 100,
    batch.length  = 5,
    n.burn        = 100,
    n.thin        = 2,
    n.chains      = 2,
    n.neighbors   = 3,
    n.omp.threads = 2,
    cov.model     = "exponential",

    fit_models   = TRUE,
    save_results = FALSE,
    verbose = FALSE,
    output_dir   = getwd()
) {

  # ---------------------------------------------------------------------------
  # Packages
  # ---------------------------------------------------------------------------

  pkgs <- c("spOccupancy", "dplyr", "tidyr", "coda", "readxl")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  }
  invisible(lapply(pkgs, library, character.only = TRUE))

  # ---------------------------------------------------------------------------
  # Basic checks
  # ---------------------------------------------------------------------------

  required_cols <- c(
    site_id_col, species_col, year_col, count_col,
    lon_col, lat_col, occ_covariates, det_covariates
  )

  missing_cols <- setdiff(required_cols, names(ecosystemBreedingBird))
  if (length(missing_cols) > 0) {
    stop("These columns are missing from ecosystemBreedingBird: ",
         paste(missing_cols, collapse = ", "))
  }

  if (is.null(covariates_to_scale)) {
    covariates_to_scale <- unique(c(occ_covariates, det_covariates))
  }

  # ---------------------------------------------------------------------------
  # Species groupings
  # ---------------------------------------------------------------------------

  if (!is.null(species_groupings)) {
    guilds <- sort(unique(species_groupings$broad_guild))
    guilds <- guilds[!is.na(guilds)]

  } else {
    speciesCnts <- length(unique(ecosystemBreedingBird[[species_col]]))

    if (single_species) {
      # One guild per species when fitting individually
      species_groupings <- data.frame(
        Species     = unique(ecosystemBreedingBird[[species_col]]),
        broad_guild = unique(ecosystemBreedingBird[[species_col]]),
        stringsAsFactors = FALSE
      )

    } else if (speciesCnts > maxSpeciesInGrps) {
      species    <- unique(ecosystemBreedingBird[[species_col]])
      n_groups   <- ceiling(length(species) / maxSpeciesInGrps)

      cat("Splitting", speciesCnts, "species into", n_groups,
          "groups (max", maxSpeciesInGrps, "each)\n")

      species_groups <- split(
        species,
        cut(seq_along(species), breaks = n_groups, labels = FALSE)
      )

      species_groupings <- do.call(rbind, lapply(seq_along(species_groups), function(i) {
        data.frame(Species     = species_groups[[i]],
                   broad_guild = paste("Species group", i),
                   stringsAsFactors = FALSE)
      }))

    } else {
      species_groupings <- data.frame(
        Species     = unique(ecosystemBreedingBird[[species_col]]),
        broad_guild = "Species group 1",
        stringsAsFactors = FALSE
      )
    }

    guilds <- sort(unique(species_groupings$broad_guild))
    guilds <- guilds[!is.na(guilds)]
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  # ---------------------------------------------------------------------------
  # Helper functions
  # ---------------------------------------------------------------------------

  make_rhs <- function(covariates, use_route_random_effect = FALSE) {
    rhs <- paste(covariates, collapse = " + ")
    if (use_route_random_effect) rhs <- paste(rhs, "+ (1 | route_id)")
    rhs
  }

  post_summary <- function(samples) {
    # Handles both matrix (multi-species) and vector/1-col (single-species)
    if (is.null(dim(samples))) samples <- matrix(samples, ncol = 1)
    out <- t(apply(samples, 2, function(x) {
      c(Mean     = mean(x, na.rm = TRUE),
        SD       = sd(x, na.rm = TRUE),
        `2.5%`   = quantile(x, 0.025, na.rm = TRUE),
        `97.5%`  = quantile(x, 0.975, na.rm = TRUE),
        `Pr(>0)` = mean(x > 0, na.rm = TRUE))
    }))
    as.data.frame(out)
  }

  extract_named_effect <- function(samples, effect_name) {
    if (is.null(dim(samples))) samples <- matrix(samples, ncol = 1,
                                                 dimnames = list(NULL, effect_name))
    sample_names <- colnames(samples)
    cols <- sample_names[
      sample_names == effect_name |
        startsWith(sample_names, paste0(effect_name, "-")) |
        startsWith(sample_names, effect_name)
    ]
    if (length(cols) == 0) return(NULL)

    out         <- post_summary(samples[, cols, drop = FALSE])
    raw_effect  <- sub("-[^-]+$", "", cols)
    species     <- sub("^.*-", "", cols)
    level       <- ifelse(raw_effect == effect_name, NA_character_,
                          sub(paste0("^", effect_name), "", raw_effect))
    out$Species   <- species
    out$effect    <- raw_effect
    out$covariate <- effect_name
    out$level     <- level
    rownames(out) <- NULL
    out
  }

  summarize_all_named_effects <- function(samples, covariates, process,
                                          species_name = NULL) {
    # For single-species models, samples are vectors/matrices without species suffix
    if (!is.null(species_name) && !is.null(dim(samples))) {
      colnames(samples) <- paste0(colnames(samples), "-", species_name)
    } else if (!is.null(species_name) && is.null(dim(samples))) {
      samples <- matrix(samples, ncol = length(covariates) + 1)
      colnames(samples) <- paste0(c("(Intercept)", covariates), "-", species_name)
    }

    effect_names <- colnames(samples)

    tabs <- lapply(covariates, function(v) {
      matched_names <- effect_names[
        effect_names == v |
          startsWith(effect_names, paste0(v, "-")) |
          startsWith(effect_names, v)
      ]
      if (length(matched_names) == 0) return(NULL)

      dplyr::bind_rows(lapply(matched_names, function(effect_name) {
        z <- extract_named_effect(samples, effect_name)
        if (is.null(z)) return(NULL)
        sp         <- sub("^.*-", "", effect_name)
        raw_effect <- sub("-[^-]+$", "", effect_name)
        level      <- ifelse(raw_effect == v, NA_character_,
                             sub(paste0("^", v), "", raw_effect))
        z$process   <- process
        z$covariate <- v
        z$level     <- level
        z$Species   <- sp
        z$effect    <- effect_name
        z
      }))
    })
    dplyr::bind_rows(tabs)
  }

  check_rhat <- function(samples, parnames, n.chains) {
    if (is.null(dim(samples))) samples <- matrix(samples, ncol = 1,
                                                 dimnames = list(NULL, parnames[1]))
    chains <- lapply(seq_len(n.chains), function(ch) {
      idx <- seq(ch, nrow(samples), by = n.chains)
      coda::as.mcmc(samples[idx, , drop = FALSE])
    })
    rhat <- coda::gelman.diag(coda::as.mcmc.list(chains))$psrf
    #rownames(rhat) <- parnames
    rhat
  }

  mode_value <- function(x) {
    x2 <- x[!is.na(x)]
    if (length(x2) == 0) return(NA)
    out <- unique(x2)[which.max(tabulate(match(x2, unique(x2))))]
    if (is.factor(x)) factor(out, levels = levels(x)) else out
  }

  all_results <- list()

  # ---------------------------------------------------------------------------
  # Build survey effort, detection/non-detection array, covariates, and coords
  # BEFORE filtering to guilds/species
  # ---------------------------------------------------------------------------

  # This is the key consistency step. The full data set is first converted into a
  # complete species x site x year detection/non-detection array. For every
  # surveyed point-year, species without a detection record are assigned 0.
  # Guild and single-species filtering happens only after this array has been
  # created, so multi-species and single-species models use the same survey
  # effort definition.

  survey_dat <- ecosystemBreedingBird %>%
    dplyr::mutate(
      Route_point_id = as.character(.data[[site_id_col]]),
      Species        = as.character(.data[[species_col]]),
      year_chr       = as.character(.data[[year_col]]),
      detect         = as.integer(.data[[count_col]] > 0),
      route_id       = if (is.null(route_id_col))
        sub("_.*", "", as.character(.data[[site_id_col]]))
      else as.character(.data[[route_id_col]]),
      point_id       = if (is.null(point_id_col))
        sub(".*_", "", as.character(.data[[site_id_col]]))
      else as.character(.data[[point_id_col]]),
      Lon = as.numeric(.data[[lon_col]]),
      Lat = as.numeric(.data[[lat_col]])
    )

  # ---------------------------------------------------------------------------
  # Remove routes/sites with missing coordinates once, before any species filter
  # ---------------------------------------------------------------------------

  routes_with_missing_coords <- character(0)

  if (remove_routes_with_missing_coords) {
    routes_with_missing_coords <- survey_dat %>%
      dplyr::filter(is.na(Lon) | is.nan(Lon) | is.na(Lat) | is.nan(Lat)) %>%
      dplyr::pull(route_id) %>%
      unique()

    survey_dat <- survey_dat %>%
      dplyr::filter(!route_id %in% routes_with_missing_coords)
  } else {
    survey_dat <- survey_dat %>%
      dplyr::filter(!is.na(Lon), !is.nan(Lon), !is.na(Lat), !is.nan(Lat))
  }

  if (nrow(survey_dat) == 0) {
    stop("All survey-effort data were removed after coordinate filtering.")
  }

  # ---------------------------------------------------------------------------
  # Global indexing
  # ---------------------------------------------------------------------------

  all_points  <- sort(unique(survey_dat$Route_point_id))
  all_routes  <- sort(unique(survey_dat$route_id))
  all_years   <- sort(unique(survey_dat$year_chr))
  all_species <- sort(unique(survey_dat$Species))

  J_all <- length(all_points)
  R_all <- length(all_years)
  N_all <- length(all_species)

  message(sprintf(
    "Global survey frame: Points: %d | Years: %d | Species: %d | Routes: %d",
    J_all, R_all, N_all, length(all_routes)
  ))

  # ---------------------------------------------------------------------------
  # Global detection/non-detection array: species x site x year
  # ---------------------------------------------------------------------------

  y_all <- array(
    NA_integer_,
    dim      = c(N_all, J_all, R_all),
    dimnames = list(all_species, all_points, all_years)
  )

  # Place observed detections/non-detections from available species rows.
  # If repeated records exist for a species-site-year, keep 1 if any record is 1.
  for (i in seq_len(nrow(survey_dat))) {
    sp  <- survey_dat$Species[i]
    pt  <- survey_dat$Route_point_id[i]
    yr  <- survey_dat$year_chr[i]
    val <- survey_dat$detect[i]

    cur <- y_all[sp, pt, yr]
    y_all[sp, pt, yr] <- if (is.na(cur)) val else max(cur, val, na.rm = TRUE)
  }

  # Survey effort is defined independently of species identity.
  # Every species receives a 0 at surveyed point-years where it was not detected.
  visited_point_year_all <- survey_dat %>%
    dplyr::distinct(Route_point_id, year_chr)

  visited_mat_all <- matrix(
    FALSE,
    nrow = J_all,
    ncol = R_all,
    dimnames = list(all_points, all_years)
  )

  visited_idx_all <- visited_point_year_all %>%
    dplyr::transmute(
      point_i = match(Route_point_id, all_points),
      year_i  = match(year_chr, all_years)
    ) %>%
    dplyr::filter(!is.na(point_i), !is.na(year_i))

  if (nrow(visited_idx_all) > 0) {
    visited_mat_all[as.matrix(visited_idx_all)] <- TRUE
  }

  for (pt in all_points) {
    for (yr in all_years) {
      if (visited_mat_all[pt, yr]) {
        y_all[is.na(y_all[, pt, yr]), pt, yr] <- 0L
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Global point-level covariate table
  # ---------------------------------------------------------------------------

  all_covariates <- unique(c(occ_covariates, det_covariates))

  point_info <- survey_dat %>%
    dplyr::group_by(Route_point_id, route_id, point_id, Year) %>%
    dplyr::summarise(
      Lon = mean(Lon, na.rm = TRUE),
      Lat = mean(Lat, na.rm = TRUE),
      dplyr::across(dplyr::all_of(all_covariates) & where(is.numeric),
                    ~ mean(.x, na.rm = TRUE)),
      dplyr::across(dplyr::all_of(all_covariates) &
                      where(~ is.factor(.x) || is.character(.x)), mode_value),
      .groups = "drop"
    ) %>%
    dplyr::arrange(match(Route_point_id, all_points), Year)

  point_info_raw <- point_info

  scale_existing <- intersect(covariates_to_scale, names(point_info))
  for (v in scale_existing) {
    if (is.numeric(point_info[[v]])) point_info[[v]] <- as.numeric(scale(point_info[[v]]))
  }
  point_info$route_id <- factor(point_info$route_id)

  # ---------------------------------------------------------------------------
  # Global occupancy covariates
  # ---------------------------------------------------------------------------

  occ.covs <- point_info %>%
    dplyr::group_by(Route_point_id, route_id, point_id) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(occ_covariates), ~ {
        if (is.numeric(.x)) mean(.x, na.rm = TRUE)
        else if (is.factor(.x) || is.character(.x)) mode_value(.x)
        else NA
      }),
      .groups = "drop"
    ) %>%
    dplyr::arrange(match(Route_point_id, all_points))

  if (use_route_random_effect) {
    occ.covs <- occ.covs %>%
      dplyr::select(dplyr::all_of(occ_covariates), route_id) %>%
      dplyr::mutate(route_id = as.numeric(as.factor(route_id)))
  } else {
    occ.covs <- occ.covs %>% dplyr::select(dplyr::all_of(occ_covariates))
  }

  occ.covs <- as.data.frame(occ.covs)
  rownames(occ.covs) <- all_points

  # ---------------------------------------------------------------------------
  # Global detection covariates
  # ---------------------------------------------------------------------------

  det.covs <- list()
  for (cov_name in det_covariates) {
    tmp <- point_info %>%
      dplyr::select(Route_point_id, Year, !!rlang::sym(cov_name)) %>%
      tidyr::pivot_wider(names_from = Year, values_from = !!rlang::sym(cov_name))

    tmp <- tmp[match(all_points, tmp$Route_point_id), ]
    mat <- as.matrix(tmp[, all_years, drop = FALSE])
    rownames(mat) <- all_points
    det.covs[[cov_name]] <- mat
  }

  # ---------------------------------------------------------------------------
  # Global coordinates
  # ---------------------------------------------------------------------------

  point_info_coords <- survey_dat %>%
    dplyr::group_by(Route_point_id, route_id, point_id) %>%
    dplyr::summarise(
      Lon = mean(Lon, na.rm = TRUE),
      Lat = mean(Lat, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(match(Route_point_id, all_points))

  coords <- as.matrix(point_info_coords[, c("Lon", "Lat")])
  rownames(coords) <- all_points

  # ---------------------------------------------------------------------------
  # Global missingness alignment
  # ---------------------------------------------------------------------------

  visit_na_all <- is.na(y_all[1, , ])

  for (nm in names(det.covs)) det.covs[[nm]][visit_na_all] <- NA

  y_observed_all <- !is.na(y_all[1, , ])
  for (nm in names(det.covs)) {
    bad <- is.na(det.covs[[nm]]) & y_observed_all
    if (any(bad)) det.covs[[nm]][bad] <- mean(det.covs[[nm]], na.rm = TRUE)
  }

  stopifnot(all(sapply(det.covs, function(x) identical(dim(x), dim(y_all)[2:3]))))

  site_never_visited_all <- apply(visit_na_all, 1, all)
  for (nm in names(occ.covs)) occ.covs[site_never_visited_all, nm] <- NA

  site_observed_all <- !site_never_visited_all
  for (nm in names(occ.covs)) {
    bad <- is.na(occ.covs[[nm]]) & site_observed_all
    if (any(bad)) {
      if (is.numeric(occ.covs[[nm]])) {
        occ.covs[[nm]][bad] <- mean(occ.covs[[nm]], na.rm = TRUE)
      } else {
        occ.covs[[nm]][bad] <- names(which.max(table(occ.covs[[nm]])))
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Formulas and parameter counts are global because site/year covariates are
  # now global and identical across model types.
  # ---------------------------------------------------------------------------

  occ_rhs     <- make_rhs(occ_covariates, use_route_random_effect)
  det_rhs     <- make_rhs(det_covariates, FALSE)
  occ_formula <- as.formula(paste("~", occ_rhs))
  det_formula <- as.formula(paste("~", det_rhs))

  dists  <- NULL; max_d <- NULL; min_d <- NULL
  phi_lo <- NULL; phi_hi <- NULL

  if (fit_spatial_model) {
    dists  <- stats::dist(coords)
    max_d  <- max(dists)
    min_d  <- stats::quantile(dists, 0.05)
    phi_lo <- 3 / max_d
    phi_hi <- 3 / min_d
  }

  p_occ <- ncol(model.matrix(occ_formula, data = occ.covs))
  if (use_route_random_effect) p_occ <- p_occ - 1

  det_cov_df <- as.data.frame(lapply(det.covs, as.vector)) %>%
    dplyr::mutate(dplyr::across(where(is.character), as.factor))
  p_det <- ncol(model.matrix(det_formula, data = det_cov_df))

  # ---------------------------------------------------------------------------
  # Loop over guilds / species
  # ---------------------------------------------------------------------------

  for (guild in guilds) {

    message("\n====================================================")
    message("Ecosystem: ", ecosystemType,
            " | ", if (single_species) "Species" else "Guild", ": ", guild)
    message("Model type: ",
            if (single_species) "Single-species" else "Multi-species",
            " | ",
            if (fit_spatial_model) "Spatial" else "Non-spatial")
    message("====================================================")

    species_in_group <- species_groupings$Species[
      species_groupings$broad_guild == guild
    ]

    species <- sort(intersect(species_in_group, dimnames(y_all)[[1]]))

    if (length(species) == 0) {
      warning("No species in the global survey frame for guild: ", guild)
      next
    }

    # Now filter only the already-completed detection/non-detection array.
    # No survey effort is recomputed inside the guild/species loop.
    y_full <- y_all[species, , , drop = FALSE]

    dat <- survey_dat %>%
      dplyr::filter(Species %in% species)

    points <- all_points
    routes <- all_routes
    years  <- all_years

    J <- length(points)
    R <- length(years)
    N <- length(species)

    message(sprintf("Points: %d | Years: %d | Species: %d | Routes: %d",
                    J, R, N, length(routes)))

    # =========================================================================
    # Branch: Single-species vs Multi-species
    # =========================================================================

    if (single_species) {

      # -----------------------------------------------------------------------
      # SINGLE-SPECIES: fit one model per species in the guild
      # -----------------------------------------------------------------------

      guild_species_results <- list()

      for (sp_name in species) {

        message("  Fitting single-species model for: ", sp_name)

        # y is a 2-D matrix (site x year) for a single species
        y_sp <- y_full[sp_name, , ]

        data.list.sp <- list(
          y        = y_sp,
          occ.covs = occ.covs,
          det.covs = det.covs
        )
        if (fit_spatial_model) data.list.sp$coords <- coords

        # Priors
        priors.sp <- list(
          beta.normal  = list(mean = rep(0, p_occ), var = rep(2.72, p_occ)),
          alpha.normal = list(mean = rep(0, p_det), var = rep(2.72, p_det))
        )
        if (fit_spatial_model) {
          # spPGOcc expects plain numeric vectors, not named lists
          priors.sp$sigma.sq.ig <- c(2, 1)
          priors.sp$phi.unif    <- c(phi_lo, phi_hi)
        }
        if (use_route_random_effect) {
          priors.sp$sigma.sq.psi.ig <- list(a = 2, b = 0.5)
        }

        # Initial values
        z.init.sp <- as.integer(apply(y_sp, 1, function(x) any(x == 1, na.rm = TRUE)))

        inits.sp <- list(
          beta  = rep(0, p_occ),
          alpha = rep(0, p_det),
          z     = z.init.sp
        )
        if (fit_spatial_model) {
          inits.sp$sigma.sq <- 1
          inits.sp$phi      <- 3 / (max_d / 2)
        }
        if (use_route_random_effect) inits.sp$sigma.sq.psi <- 1

        out.sp <- NULL

        if (fit_models) {
          if (fit_spatial_model) {
            out.sp <- spOccupancy::spPGOcc(
              occ.formula   = occ_formula,
              det.formula   = det_formula,
              data          = data.list.sp,
              inits         = inits.sp,
              priors        = priors.sp,
              cov.model     = cov.model,
              NNGP          = TRUE,
              n.neighbors   = n.neighbors,
              n.batch       = n.batch,
              batch.length  = batch.length,
              n.burn        = n.burn,
              n.thin        = n.thin,
              n.chains      = n.chains,
              n.omp.threads = n.omp.threads,
              verbose       = verbose,
              n.report      = max(1, floor(n.batch / 10))
            )
          } else {
            out.sp <- spOccupancy::PGOcc(
              occ.formula   = occ_formula,
              det.formula   = det_formula,
              data          = data.list.sp,
              inits         = inits.sp,
              priors        = priors.sp,
              #n.batch       = n.batch,
              #batch.length  = batch.length,
              n.samples     = n.batch * batch.length,
              n.burn        = n.burn,
              n.thin        = n.thin,
              n.chains      = n.chains,
              n.omp.threads = n.omp.threads,
              verbose       = verbose#,
              #n.report      = max(1, floor(n.batch / 10))
            )
          }
        }

        # --- Extract single-species results ---
        # PGOcc / spPGOcc output shapes differ from the multi-species variants:
        #   psi.samples  : matrix  (n_iter  x  J)          -- no species dimension
        #   z.samples    : matrix  (n_iter  x  J)          -- no species dimension
        #   w.samples    : matrix  (n_iter  x  J)          -- spatial only, no species dim
        #   beta.samples : matrix  (n_iter  x  p_occ)
        #   alpha.samples: matrix  (n_iter  x  p_det)
        #   ppcOcc       : no $sp.names field
        sp_products <- list()

        if (!is.null(out.sp)) {

          ppc.sp  <- spOccupancy::ppcOcc(out.sp, fit.stat = "freeman-tukey", group = 1)
          bp_sp   <- mean(ppc.sp$fit.y.rep > ppc.sp$fit.y, na.rm = TRUE)
          waic.sp <- spOccupancy::waicOcc(out.sp)

          # ppc summary — single-species ppcOcc has no $sp.names; inject species name
          ppc_df.sp <- data.frame(
            ecosystem  = ecosystemType,
            guild      = guild,
            Species    = sp_name,
            bayesian_p = bp_sp,
            row.names  = NULL
          )

          # Coefficient summaries
          beta_summary            <- post_summary(out.sp$beta.samples)
          beta_summary$parameter  <- colnames(out.sp$beta.samples)

          alpha_summary           <- post_summary(out.sp$alpha.samples)
          alpha_summary$parameter <- colnames(out.sp$alpha.samples)

          # Name columns to match multi-species "-Species" convention for helpers
          beta_named  <- out.sp$beta.samples
          alpha_named <- out.sp$alpha.samples
          colnames(beta_named)  <- paste0(colnames(out.sp$beta.samples),  "-", sp_name)
          colnames(alpha_named) <- paste0(colnames(out.sp$alpha.samples), "-", sp_name)

          occ_effects <- summarize_all_named_effects(beta_named,  occ_covariates, "occupancy")
          det_effects <- summarize_all_named_effects(alpha_named, det_covariates, "detection")

          # Spatial parameters (sigma^2, phi, effective range)
          spatial_params <- NULL
          if (fit_spatial_model && !is.null(out.sp$theta.samples)) {
            # theta.samples for spPGOcc: columns named "sigma.sq" and "phi"
            theta_m    <- colMeans(out.sp$theta.samples)
            sigma_mean <- theta_m[grep("sigma.sq", names(theta_m))]
            phi_mean   <- theta_m[grep("^phi",     names(theta_m))]
            spatial_params <- data.frame(
              ecosystem = ecosystemType,
              guild     = guild,
              Species   = sp_name,
              sigma_sq  = sigma_mean,
              phi       = phi_mean,
              eff_range = 3 / phi_mean,
              row.names = NULL
            )
          }

          # ------------------------------------------------------------------
          # Occupancy probability (psi): matrix (n_iter x J) for single-species
          # ------------------------------------------------------------------
          psi.mean.sp <- colMeans(out.sp$psi.samples)          # length-J vector
          names(psi.mean.sp) <- points

          psi.ci.sp <- apply(                                   # 2 x J matrix
            out.sp$psi.samples, 2,
            function(x) stats::quantile(x, c(0.025, 0.975))
          )
          colnames(psi.ci.sp) <- points

          psi_long.sp <- data.frame(
            Species   = sp_name,
            Point     = points,
            psi_mean  = psi.mean.sp,
            psi_lwr   = psi.ci.sp["2.5%",  ],
            psi_upr   = psi.ci.sp["97.5%", ],
            ecosystem = ecosystemType,
            guild     = guild,
            row.names = NULL
          ) %>%
            dplyr::left_join(point_info_raw, by = c("Point" = "Route_point_id"))

          # ------------------------------------------------------------------
          # Latent occupancy state (z): matrix (n_iter x J) for single-species
          # Gives posterior probability of true occupancy at each site.
          # ------------------------------------------------------------------
          z.mean.sp <- colMeans(out.sp$z.samples)               # length-J vector
          names(z.mean.sp) <- points

          z.ci.sp <- apply(
            out.sp$z.samples, 2,
            function(x) stats::quantile(x, c(0.025, 0.975))
          )
          colnames(z.ci.sp) <- points

          z_long.sp <- data.frame(
            Species   = sp_name,
            Point     = points,
            z_mean    = z.mean.sp,
            z_lwr     = z.ci.sp["2.5%",  ],
            z_upr     = z.ci.sp["97.5%", ],
            ecosystem = ecosystemType,
            guild     = guild,
            row.names = NULL
          ) %>%
            dplyr::left_join(point_info_raw, by = c("Point" = "Route_point_id"))

          # ------------------------------------------------------------------
          # Spatial random effect (w): matrix (n_iter x J) for spPGOcc,
          # NULL for non-spatial PGOcc.
          # ------------------------------------------------------------------
          w.mean.sp <- NULL
          w_long.sp <- NULL

          if (fit_spatial_model && !is.null(out.sp$w.samples)) {
            # w.samples for spPGOcc: matrix (n_iter x J)
            w.mean.sp        <- colMeans(out.sp$w.samples)
            names(w.mean.sp) <- points

            w_long.sp <- data.frame(
              Species   = sp_name,
              Point     = points,
              w_mean    = w.mean.sp,
              ecosystem = ecosystemType,
              guild     = guild,
              row.names = NULL
            ) %>%
              dplyr::left_join(point_info_raw, by = c("Point" = "Route_point_id"))
          }

          # Convergence diagnostics
          rhat_beta <- tryCatch(
            check_rhat(out.sp$beta.samples,
                       c("(Intercept)", occ_covariates), n.chains),
            error = function(e) NULL
          )
          rhat_alpha <- tryCatch(
            check_rhat(out.sp$alpha.samples,
                       c("(Intercept)", det_covariates), n.chains),
            error = function(e) NULL
          )

          sp_products <- list(
            model          = out.sp,
            ppc            = ppc.sp,
            ppc_species    = ppc_df.sp,       # consistent name with multi-species
            ppc_overall    = bp_sp,
            waic           = waic.sp,
            beta_summary   = beta_summary,
            alpha_summary  = alpha_summary,
            occ_effects    = occ_effects,
            det_effects    = det_effects,
            spatial_params = spatial_params,
            psi_mean       = psi.mean.sp,     # named length-J vector
            psi_ci         = psi.ci.sp,       # 2 x J matrix
            psi_long       = psi_long.sp,     # long df with lwr/upr
            z_mean         = z.mean.sp,       # named length-J vector
            z_ci           = z.ci.sp,         # 2 x J matrix
            z_long         = z_long.sp,       # long df with lwr/upr
            w_mean         = w.mean.sp,       # named length-J vector or NULL
            w_long         = w_long.sp,       # long df or NULL
            rhat_beta      = rhat_beta,
            rhat_alpha     = rhat_alpha
          )
        }

        guild_species_results[[sp_name]] <- list(
          species      = sp_name,
          guild        = guild,
          ecosystem    = ecosystemType,
          points       = points,
          routes       = routes,
          years        = years,
          dims         = list(J = J, R = R),
          formulas     = list(occ_formula = occ_formula, det_formula = det_formula),
          fit_spatial_model = fit_spatial_model,
          single_species    = TRUE,
          data_used    = dat %>% dplyr::filter(Species == sp_name),
          data_list    = data.list.sp,
          coords       = coords,
          priors       = priors.sp,
          inits        = inits.sp,
          model_products = sp_products
        )

        if (save_results) {
          saveRDS(
            guild_species_results[[sp_name]],
            file = file.path(output_dir,
                             paste0("spPGOcc_", ecosystemType, "_",
                                    gsub("[^A-Za-z0-9]+", "_", sp_name), ".rds"))
          )
        }
      } # end species loop

      all_results[[guild]] <- guild_species_results

    } else {

      # =======================================================================
      # MULTI-SPECIES: original path
      # =======================================================================

      data.list <- list(
        y        = y_full,
        occ.covs = occ.covs,
        det.covs = det.covs
      )
      if (fit_spatial_model) data.list$coords <- coords

      priors <- list(
        beta.comm.normal  = list(mean = rep(0, p_occ), var = rep(2.72, p_occ)),
        tau.sq.beta.ig    = list(a = 0.1, b = 0.1),
        alpha.comm.normal = list(mean = rep(0, p_det), var = rep(2.72, p_det)),
        tau.sq.alpha.ig   = list(a = 0.1, b = 0.1)
      )
      if (fit_spatial_model) {
        priors$sigma.sq.ig <- list(a = 2, b = 1)
        priors$phi.unif    <- list(a = phi_lo, b = phi_hi)
      }
      if (use_route_random_effect) priors$sigma.sq.psi.ig <- list(a = 2, b = 0.5)

      z.init <- apply(y_full, c(1, 2), function(x) as.integer(any(x == 1, na.rm = TRUE)))
      inits  <- list(
        beta.comm    = rep(0, p_occ),
        alpha.comm   = rep(0, p_det),
        tau.sq.beta  = rep(1, p_occ),
        tau.sq.alpha = rep(1, p_det),
        z            = z.init
      )
      if (fit_spatial_model) {
        inits$sigma.sq <- rep(1, N)
        inits$phi      <- rep(3 / (max_d / 2), N)
      }
      if (use_route_random_effect) inits$sigma.sq.psi <- 1

      out <- NULL

      if (fit_models) {
        message("Fitting model for guild: ", guild)

        if (fit_spatial_model) {
          out <- spOccupancy::spMsPGOcc(
            occ.formula   = occ_formula,
            det.formula   = det_formula,
            data          = data.list,
            inits         = inits,
            priors        = priors,
            cov.model     = cov.model,
            NNGP          = TRUE,
            n.neighbors   = n.neighbors,
            n.batch       = n.batch,
            batch.length  = batch.length,
            n.burn        = n.burn,
            n.thin        = n.thin,
            n.chains      = n.chains,
            n.omp.threads = n.omp.threads,
            verbose       = verbose,
            n.report      = max(1, floor(n.batch / 10))
          )
        } else {
          out <- spOccupancy::msPGOcc(
            occ.formula   = occ_formula,
            det.formula   = det_formula,
            data          = data.list,
            inits         = inits,
            priors        = priors,
            #n.batch       = n.batch,
            #batch.length  = batch.length,
            n.samples = batch.length * n.batch,
            n.burn        = n.burn,
            n.thin        = n.thin,
            n.chains      = n.chains,
            n.omp.threads = n.omp.threads,
            verbose       = verbose#,
            #n.report      = max(1, floor(n.batch / 10))
          )
        }
      }

      # -----------------------------------------------------------------------
      # Extract multi-species model products (original logic, unchanged)
      # -----------------------------------------------------------------------

      model_products <- list()

      if (!is.null(out)) {

        ppc <- spOccupancy::ppcOcc(out, fit.stat = "freeman-tukey", group = 1)
        bp_species <- colMeans(ppc$fit.y.rep > ppc$fit.y, na.rm = TRUE)
        bp_overall <- mean(ppc$fit.y.rep > ppc$fit.y, na.rm = TRUE)
        waic       <- spOccupancy::waicOcc(out)

        beta_comm_summary           <- post_summary(out$beta.comm.samples)
        beta_comm_summary$parameter <- c("(Intercept)", occ_covariates)

        alpha_comm_summary           <- post_summary(out$alpha.comm.samples)
        alpha_comm_summary$parameter <- c("(Intercept)", det_covariates)

        occ_effects <- summarize_all_named_effects(
          out$beta.samples, occ_covariates, "occupancy")
        det_effects <- summarize_all_named_effects(
          out$alpha.samples, det_covariates, "detection")

        spatial_params <- NULL
        if (fit_spatial_model && !is.null(out$theta.samples)) {
          sigma_cols <- grep("^sigma.sq-", colnames(out$theta.samples), value = TRUE)
          phi_cols   <- grep("^phi-",      colnames(out$theta.samples), value = TRUE)
          if (length(sigma_cols) > 0 && length(phi_cols) > 0) {
            sigma_mean <- colMeans(out$theta.samples[, sigma_cols, drop = FALSE])
            phi_mean   <- colMeans(out$theta.samples[, phi_cols,   drop = FALSE])
            spatial_params <- data.frame(
              ecosystem = ecosystemType,
              guild     = guild,
              Species   = sub("^sigma.sq-", "", sigma_cols),
              sigma_sq  = sigma_mean,
              phi       = phi_mean,
              eff_range = 3 / phi_mean,
              row.names = NULL
            )
          }
        }

        # ------------------------------------------------------------------
        # Occupancy probability (psi): 4-D array (n_iter x N x J x 1) -> use
        # apply over dims 2 & 3 (species x site).  psi.samples for msPGOcc /
        # spMsPGOcc has shape (n_iter x N x J).
        # ------------------------------------------------------------------
        psi.mean <- apply(out$psi.samples, c(2, 3), mean)       # N x J
        dimnames(psi.mean) <- list(species, points)

        psi.lwr <- apply(out$psi.samples, c(2, 3),              # N x J
                         function(x) stats::quantile(x, 0.025))
        psi.upr <- apply(out$psi.samples, c(2, 3),
                         function(x) stats::quantile(x, 0.975))
        dimnames(psi.lwr) <- dimnames(psi.upr) <- list(species, points)

        # Keep the full 3-D CI array for downstream use (2 x N x J)
        psi.ci <- apply(out$psi.samples, c(2, 3),
                        function(x) stats::quantile(x, c(0.025, 0.975)))

        # ------------------------------------------------------------------
        # Spatial random effect (w): 3-D array (n_iter x N x J) for spMsPGOcc
        # ------------------------------------------------------------------
        w.mean <- NULL
        if (fit_spatial_model && !is.null(out$w.samples)) {
          w.mean <- apply(out$w.samples, c(2, 3), mean)          # N x J
          dimnames(w.mean) <- list(species, points)
        }

        # ------------------------------------------------------------------
        # Latent occupancy state (z): 3-D array (n_iter x N x J)
        # z_mean gives the posterior probability of true occupancy per species x site.
        # Summing over species gives species richness.
        # ------------------------------------------------------------------
        z.mean <- apply(out$z.samples, c(2, 3), mean)            # N x J
        dimnames(z.mean) <- list(species, points)

        z.lwr <- apply(out$z.samples, c(2, 3),
                       function(x) stats::quantile(x, 0.025))
        z.upr <- apply(out$z.samples, c(2, 3),
                       function(x) stats::quantile(x, 0.975))
        dimnames(z.lwr) <- dimnames(z.upr) <- list(species, points)

        # Richness: sum occupied species per site across posterior draws
        richness.samples <- apply(out$z.samples, c(1, 3), sum)   # n_iter x J
        richness.mean    <- colMeans(richness.samples)
        richness.ci      <- apply(richness.samples, 2,
                                  function(x) stats::quantile(x, c(0.025, 0.975)))
        names(richness.mean) <- points

        richness_df <- data.frame(
          ecosystem      = ecosystemType,
          guild          = guild,
          Point          = points,
          Route          = point_info_coords$route_id,
          Lon            = point_info_coords$Lon,
          Lat            = point_info_coords$Lat,
          richness_mean  = richness.mean,
          richness_lwr   = richness.ci[1, ],
          richness_upr   = richness.ci[2, ],
          row.names      = NULL
        ) %>%
          dplyr::left_join(point_info_raw, by = c("Point" = "Route_point_id"))

        # psi long: include lower and upper CIs for consistency with single-species
        psi_long <- as.data.frame(as.table(psi.mean)) %>%
          dplyr::rename(Species = Var1, Point = Var2, psi_mean = Freq) %>%
          dplyr::left_join(
            as.data.frame(as.table(psi.lwr)) %>%
              dplyr::rename(Species = Var1, Point = Var2, psi_lwr = Freq),
            by = c("Species", "Point")
          ) %>%
          dplyr::left_join(
            as.data.frame(as.table(psi.upr)) %>%
              dplyr::rename(Species = Var1, Point = Var2, psi_upr = Freq),
            by = c("Species", "Point")
          ) %>%
          dplyr::mutate(ecosystem = ecosystemType, guild = guild) %>%
          dplyr::left_join(point_info_raw, by = c("Point" = "Route_point_id"))

        # z long: posterior mean occupancy state per species x site
        z_long <- as.data.frame(as.table(z.mean)) %>%
          dplyr::rename(Species = Var1, Point = Var2, z_mean = Freq) %>%
          dplyr::left_join(
            as.data.frame(as.table(z.lwr)) %>%
              dplyr::rename(Species = Var1, Point = Var2, z_lwr = Freq),
            by = c("Species", "Point")
          ) %>%
          dplyr::left_join(
            as.data.frame(as.table(z.upr)) %>%
              dplyr::rename(Species = Var1, Point = Var2, z_upr = Freq),
            by = c("Species", "Point")
          ) %>%
          dplyr::mutate(ecosystem = ecosystemType, guild = guild) %>%
          dplyr::left_join(point_info_raw, by = c("Point" = "Route_point_id"))

        w_long <- NULL
        if (!is.null(w.mean)) {
          w_long <- as.data.frame(as.table(w.mean)) %>%
            dplyr::rename(Species = Var1, Point = Var2, w_mean = Freq) %>%
            dplyr::mutate(ecosystem = ecosystemType, guild = guild) %>%
            dplyr::left_join(point_info_raw, by = c("Point" = "Route_point_id"))
        }

        ppc_df <- data.frame(
          ecosystem  = ecosystemType,
          guild      = guild,
          Species    = ppc$sp.names,
          bayesian_p = bp_species,
          row.names  = NULL
        )

        rhat_beta_comm <- tryCatch(
          check_rhat(out$beta.comm.samples, c("(Intercept)", occ_covariates), n.chains),
          error = function(e) NULL)
        rhat_alpha_comm <- tryCatch(
          check_rhat(out$alpha.comm.samples, c("(Intercept)", det_covariates), n.chains),
          error = function(e) NULL)

        model_products <- list(
          model               = out,
          ppc                 = ppc,
          ppc_species         = ppc_df,
          ppc_overall         = bp_overall,
          waic                = waic,
          beta_comm_summary   = beta_comm_summary,
          alpha_comm_summary  = alpha_comm_summary,
          occ_effects         = occ_effects,
          det_effects         = det_effects,
          spatial_params      = spatial_params,
          psi_mean            = psi.mean,       # N x J matrix
          psi_ci              = psi.ci,          # 2 x N x J array
          psi_long            = psi_long,        # long df with psi_lwr/psi_upr
          z_mean              = z.mean,          # N x J matrix
          z_long              = z_long,          # long df with z_lwr/z_upr
          w_mean              = w.mean,          # N x J matrix or NULL
          w_long              = w_long,
          richness_samples    = richness.samples,
          richness_mean       = richness.mean,
          richness_ci         = richness.ci,
          richness_df         = richness_df,
          rhat_beta_comm      = rhat_beta_comm,
          rhat_alpha_comm     = rhat_alpha_comm
        )
      }

      all_results[[guild]] <- list(
        ecosystem                    = ecosystemType,
        guild                        = guild,
        species                      = species,
        points                       = points,
        routes                       = routes,
        years                        = years,
        dims                         = list(N = N, J = J, R = R),
        removed_routes_missing_coords = routes_with_missing_coords,
        formulas    = list(occ_formula = occ_formula, det_formula = det_formula),
        fit_spatial_model = fit_spatial_model,
        single_species    = FALSE,
        data_used         = dat,
        data_list         = data.list,
        point_info_raw    = point_info_raw,
        point_info_model  = point_info,
        coords            = coords,
        priors            = priors,
        inits             = inits,
        model_products    = model_products
      )

      if (save_results) {
        saveRDS(
          all_results[[guild]],
          file = file.path(output_dir,
                           paste0("spMsPGOcc_", ecosystemType, "_",
                                  gsub("[^A-Za-z0-9]+", "_", guild), ".rds"))
        )
      }
    } # end multi-species branch
  } # end guild loop

  # ---------------------------------------------------------------------------
  # Combine across guilds
  # ---------------------------------------------------------------------------

  if (single_species) {
    # Flatten the nested list (guild -> species -> products) into combined tables

    all_sp_results <- unlist(all_results, recursive = FALSE)

    get_ss_prod <- function(x, nm) {
      if (length(x$model_products) == 0) return(NULL)
      x$model_products[[nm]]
    }

    combined <- list(
      psi                    = dplyr::bind_rows(lapply(all_sp_results, get_ss_prod, "psi_long")),
      z                      = dplyr::bind_rows(lapply(all_sp_results, get_ss_prod, "z_long")),
      spatial_random_effects = dplyr::bind_rows(lapply(all_sp_results, get_ss_prod, "w_long")),
      ppc                    = dplyr::bind_rows(lapply(all_sp_results, get_ss_prod, "ppc_species")),
      spatial_params         = dplyr::bind_rows(lapply(all_sp_results, get_ss_prod, "spatial_params")),
      occ_effects            = dplyr::bind_rows(lapply(all_sp_results, get_ss_prod, "occ_effects")),
      det_effects            = dplyr::bind_rows(lapply(all_sp_results, get_ss_prod, "det_effects"))
    )

  } else {
    get_prod <- function(x, nm) {
      if (length(x$model_products) == 0) return(NULL)
      x$model_products[[nm]]
    }

    combined <- list(
      richness               = dplyr::bind_rows(lapply(all_results, get_prod, "richness_df")),
      psi                    = dplyr::bind_rows(lapply(all_results, get_prod, "psi_long")),
      z                      = dplyr::bind_rows(lapply(all_results, get_prod, "z_long")),
      spatial_random_effects = dplyr::bind_rows(lapply(all_results, get_prod, "w_long")),
      ppc                    = dplyr::bind_rows(lapply(all_results, get_prod, "ppc_species")),
      spatial_params         = dplyr::bind_rows(lapply(all_results, get_prod, "spatial_params")),
      occ_effects            = dplyr::bind_rows(lapply(all_results, get_prod, "occ_effects")),
      det_effects            = dplyr::bind_rows(lapply(all_results, get_prod, "det_effects"))
    )
  }

  return(list(
    ecosystem         = ecosystemType,
    single_species    = single_species,
    fit_spatial_model = fit_spatial_model,
    occ_covariates    = occ_covariates,
    det_covariates    = det_covariates,
    guild_results     = all_results,
    combined          = combined
  ))
}
