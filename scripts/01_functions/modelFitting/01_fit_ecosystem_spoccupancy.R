# =============================================================================
# GENERAL SPATIAL MULTISPECIES OCCUPANCY FUNCTION BY ECOSYSTEM AND SPECIES GROUP
# =============================================================================

fit_ecosystem_spoccupancy <- function(
    ecosystemBreedingBird = NULL,
    ecosystemType = NULL,
    maxSpeciesInGrps = 10,

    occ_covariates = NULL,
    det_covariates = NULL,

    # This will be taken
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
    output_dir   = getwd()
) {

  # ---------------------------------------------------------------------------
  # Packages
  # ---------------------------------------------------------------------------

  pkgs <- c("spOccupancy", "dplyr", "tidyr", "coda", "readxl")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p)
    }
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
    stop(
      "These columns are missing from ecosystemBreedingBird: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (is.null(covariates_to_scale)) {
    covariates_to_scale <- unique(c(occ_covariates, det_covariates))
  }


  if(!is.null(species_groupings)){
  guilds <- sort(unique(species_groupings$broad_guild))
  guilds <- guilds[!is.na(guilds)]

  } else {
    speciesCnts <- length(unique(ecosystemBreedingBird$Species))

    # spOccupancy takes enormous amount of time when fitting models with large number of species
    # so if the species counts are greater than 10, just split the species into groups with at most 10 species
    if(speciesCnts > maxSpeciesInGrps){
      species <- unique(ecosystemBreedingBird$Species)
      n_groups <- ceiling(length(species) / maxSpeciesInGrps)

      cat("Splitting the species into", n_groups,"groups with at most", maxSpeciesInGrps, "species in each group")

      species_groups <- split(
      species,
      cut(
        seq_along(species),
        breaks = n_groups,
        labels = FALSE
      )
    )

      species_groupings <- do.call(
        rbind,
        lapply(seq_along(species_groups), function(i) {
          data.frame(
            Species = species_groups[[i]],
            broad_guild = paste("Species group", i),
            stringsAsFactors = FALSE
          )
        })
      )

    } else {
      species_groupings <- data.frame(Species = unique(ecosystemBreedingBird$Species),
                                      broad_guild = paste("Species group", 1))
    }

    # Select the groupings within each group
    guilds <- sort(unique(species_groupings$broad_guild))
    guilds <- guilds[!is.na(guilds)]
  }

  # ---------------------------------------------------------------------------
  # Helper functions
  # ---------------------------------------------------------------------------

  make_rhs <- function(covariates, use_route_random_effect = FALSE) {
    rhs <- paste(covariates, collapse = " + ")

    if (use_route_random_effect) {
      rhs <- paste(rhs, "+ (1 | route_id)")
    }

    rhs
  }

  post_summary <- function(samples) {
    out <- t(apply(samples, 2, function(x) {
      c(
        Mean     = mean(x, na.rm = TRUE),
        SD       = sd(x, na.rm = TRUE),
        `2.5%`   = quantile(x, 0.025, na.rm = TRUE),
        `97.5%`  = quantile(x, 0.975, na.rm = TRUE),
        `Pr(>0)` = mean(x > 0, na.rm = TRUE)
      )
    }))

    as.data.frame(out)
  }

  extract_named_effect <- function(samples, effect_name) {

    sample_names <- colnames(samples)

    cols <- sample_names[
      sample_names == effect_name |
        startsWith(sample_names, paste0(effect_name, "-")) |
        startsWith(sample_names, effect_name)
    ]

    if (length(cols) == 0) {
      return(NULL)
    }

    out <- post_summary(samples[, cols, drop = FALSE])

    raw_effect <- sub("-[^-]+$", "", cols)
    species <- sub("^.*-", "", cols)

    level <- ifelse(
      raw_effect == effect_name,
      NA_character_,
      sub(paste0("^", effect_name), "", raw_effect)
    )

    out$Species <- species
    out$effect <- raw_effect
    out$covariate <- effect_name
    out$level <- level

    rownames(out) <- NULL

    out
  }

  # Summarise the covariates
  summarize_all_named_effects <- function(samples, covariates, process) {

    effect_names <- colnames(samples)

    tabs <- lapply(covariates, function(v) {

      matched_names <- effect_names[
        effect_names == v |
          startsWith(effect_names, paste0(v, "-")) |
          startsWith(effect_names, v)
      ]

      if (length(matched_names) == 0) {
        return(NULL)
      }

      dplyr::bind_rows(lapply(matched_names, function(effect_name) {

        z <- extract_named_effect(samples, effect_name)
        if (is.null(z)) return(NULL)

        species <- sub("^.*-", "", effect_name)
        raw_effect <- sub("-[^-]+$", "", effect_name)

        level <- ifelse(
          raw_effect == v,
          NA_character_,
          sub(paste0("^", v), "", raw_effect)
        )

        z$process <- process
        z$covariate <- v
        z$level <- level
        z$Species <- species
        z$effect <- effect_name

        z
      }))
    })

    dplyr::bind_rows(tabs)
  }

  check_rhat <- function(samples, parnames, n.chains) {
    chains <- lapply(seq_len(n.chains), function(ch) {
      idx <- seq(ch, nrow(samples), by = n.chains)
      coda::as.mcmc(samples[idx, , drop = FALSE])
    })

    rhat <- coda::gelman.diag(coda::as.mcmc.list(chains))$psrf
    rownames(rhat) <- parnames
    rhat
  }

  mode_value <- function(x) {
    x2 <- x[!is.na(x)]

    if (length(x2) == 0) {
      return(NA)
    }

    out <- unique(x2)[which.max(tabulate(match(x2, unique(x2))))]

    if (is.factor(x)) {
      factor(out, levels = levels(x))
    } else {
      out
    }
  }

  all_results <- list()

  # ---------------------------------------------------------------------------
  # Loop over guilds
  # ---------------------------------------------------------------------------

  for (guild in guilds) {

    message("\n====================================================")
    message("Ecosystem: ", ecosystemType, " | Guild: ", guild)
    message("====================================================")

    species_in_group <- species_groupings$Species[
      species_groupings$broad_guild == guild
    ]

    dat <- ecosystemBreedingBird %>%
      dplyr::filter(.data[[species_col]] %in% species_in_group) %>%
      dplyr::mutate(
        Route_point_id = as.character(.data[[site_id_col]]),
        Species        = as.character(.data[[species_col]]),
        year_chr       = as.character(.data[[year_col]]),
        detect         = as.integer(.data[[count_col]] > 0),
        route_id       = if (is.null(route_id_col)) {
          sub("_.*", "", as.character(.data[[site_id_col]]))
        } else {
          as.character(.data[[route_id_col]])
        },
        point_id       = if (is.null(point_id_col)) {
          sub(".*_", "", as.character(.data[[site_id_col]]))
        } else {
          as.character(.data[[point_id_col]])
        },
        Lon            = as.numeric(.data[[lon_col]]),
        Lat            = as.numeric(.data[[lat_col]])
      )

    if (nrow(dat) == 0) {
      warning("No data for guild: ", guild)
      next
    }

    # -------------------------------------------------------------------------
    # Remove routes with missing coordinates
    # -------------------------------------------------------------------------

    routes_with_missing_coords <- character(0)

    if (remove_routes_with_missing_coords) {
      routes_with_missing_coords <- dat %>%
        dplyr::filter(is.na(Lon) | is.nan(Lon) | is.na(Lat) | is.nan(Lat)) %>%
        dplyr::pull(route_id) %>%
        unique()

      dat <- dat %>%
        dplyr::filter(!route_id %in% routes_with_missing_coords)
    } else {
      dat <- dat %>%
        dplyr::filter(!is.na(Lon), !is.nan(Lon), !is.na(Lat), !is.nan(Lat))
    }

    if (nrow(dat) == 0) {
      warning("All data removed after coordinate filtering for guild: ", guild)
      next
    }

    # -------------------------------------------------------------------------
    # Indexing
    # -------------------------------------------------------------------------

    points  <- sort(unique(dat$Route_point_id))
    routes  <- sort(unique(dat$route_id))
    years   <- sort(unique(dat$year_chr))
    species <- sort(unique(dat$Species))

    J <- length(points)
    R <- length(years)
    N <- length(species)

    message(sprintf(
      "Points: %d | Years: %d | Species: %d | Routes: %d",
      J, R, N, length(routes)
    ))

    # -------------------------------------------------------------------------
    # Detection array: species x point x year
    # -------------------------------------------------------------------------

    y <- array(
      NA_integer_,
      dim = c(N, J, R),
      dimnames = list(species, points, years)
    )

    for (i in seq_len(nrow(dat))) {
      sp  <- dat$Species[i]
      pt  <- dat$Route_point_id[i]
      yr  <- dat$year_chr[i]
      val <- dat$detect[i]

      cur <- y[sp, pt, yr]

      y[sp, pt, yr] <- if (is.na(cur)) {
        val
      } else {
        max(cur, val, na.rm = TRUE)
      }
    }

    visited_point_year <- dat %>%
      dplyr::distinct(Route_point_id, year_chr)

    visited_mat <- matrix(
      FALSE,
      nrow = J,
      ncol = R,
      dimnames = list(points, years)
    )

    visited_mat[
      as.matrix(
        visited_point_year %>%
          dplyr::transmute(
            point_i = match(Route_point_id, points),
            year_i  = match(year_chr, years)
          )
      )
    ] <- TRUE

    for (pt in points) {
      for (yr in years) {
        if (visited_mat[pt, yr]) {
          y[is.na(y[, pt, yr]), pt, yr] <- 0L
        }
      }
    }

    # -------------------------------------------------------------------------
    # Point-level covariate table
    # -------------------------------------------------------------------------

    all_covariates <- unique(c(occ_covariates, det_covariates))

    point_info <- dat %>%
      dplyr::group_by(Route_point_id, route_id, point_id, Year) %>%
      dplyr::summarise(
        Lon = mean(Lon, na.rm = TRUE),
        Lat = mean(Lat, na.rm = TRUE),
        dplyr::across(
          dplyr::all_of(all_covariates) & where(is.numeric),
          ~ mean(.x, na.rm = TRUE)
        ),
        dplyr::across(
          dplyr::all_of(all_covariates) & where(~ is.factor(.x) || is.character(.x)),
          mode_value
        ),
        .groups = "drop"
      ) %>%
      dplyr::arrange(match(Route_point_id, points))

    #stopifnot(all(point_info$Route_point_id == points))

    point_info_raw <- point_info

    # Scale requested numeric covariates
    scale_existing <- intersect(covariates_to_scale, names(point_info))

    for (v in scale_existing) {
      if (is.numeric(point_info[[v]])) {
        point_info[[v]] <- as.numeric(scale(point_info[[v]]))
      }
    }

    point_info$route_id <- factor(point_info$route_id)

    # -------------------------------------------------------------------------
    # Occupancy covariates
    # -------------------------------------------------------------------------

    # Occupancy variables are averaged over the years
    occ.covs <- point_info %>%
      dplyr::group_by(Route_point_id, route_id, point_id) %>%
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(occ_covariates),
          ~ {
            if (is.numeric(.x)) {
              mean(.x, na.rm = TRUE)
            } else if (is.factor(.x) || is.character(.x)) {
              mode_value(.x)
            } else {
              NA
            }
          }
        ),
        .groups = "drop"
      ) %>%
      dplyr::arrange(match(Route_point_id, points))





    if (use_route_random_effect) {
      occ.covs <- occ.covs %>%
        dplyr::select(dplyr::all_of(occ_covariates), route_id) %>%
        dplyr::mutate(route_id = as.numeric(as.factor(route_id)))
    } else {
      occ.covs <- occ.covs %>%
        dplyr::select(dplyr::all_of(occ_covariates))
    }

    occ.covs <- as.data.frame(occ.covs)
    rownames(occ.covs) <- points

    # -------------------------------------------------------------------------
    # Detection covariates
    # -------------------------------------------------------------------------

    det.covs <- list()

    for (cov_name in det_covariates) {

      tmp <- point_info %>%
        dplyr::select(
          Route_point_id,
          Year,
          !!rlang::sym(cov_name)
        ) %>%
        tidyr::pivot_wider(
          names_from = Year,
          values_from = !!rlang::sym(cov_name)
        )

      tmp <- tmp[
        match(points, tmp$Route_point_id),
      ]

      mat <- as.matrix(
        tmp[, years, drop = FALSE]
      )

      rownames(mat) <- points

      det.covs[[cov_name]] <- mat
    }
  #  }

    # -------------------------------------------------------------------------
    # Coordinates
    # -------------------------------------------------------------------------
    point_info <- dat %>%
      dplyr::group_by(Route_point_id, route_id, point_id) %>%
      dplyr::summarise(
        Lon = mean(Lon, na.rm = TRUE),
        Lat = mean(Lat, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(match(Route_point_id, points))

    coords <- as.matrix(point_info[, c("Lon", "Lat")])
    rownames(coords) <- points

    # -------------------------------------------------------------------------
    # Missingness alignment
    # -------------------------------------------------------------------------

    visit_na <- is.na(y[1, , ])

    for (nm in names(det.covs)) {
      det.covs[[nm]][visit_na] <- NA
    }

    # for(nm in names(occ.covs)){
    #   occ.covs[visit_na, nm] <- NA
    # }

    y_observed <- !is.na(y[1, , ])

    for (nm in names(det.covs)) {
      bad <- is.na(det.covs[[nm]]) & y_observed

      if (any(bad)) {
        det.covs[[nm]][bad] <- mean(det.covs[[nm]], na.rm = TRUE)
      }
    }

    stopifnot(all(
      sapply(det.covs, function(x) identical(dim(x), dim(y)[2:3]))
    ))

    # site-level indicator
    site_never_visited <- apply(visit_na, 1, all)

    for (nm in names(occ.covs)) {
      occ.covs[site_never_visited, nm] <- NA
    }

    # Impute missing occupancy covariates for sites that remain in the analysis
    site_observed <- !site_never_visited

    for (nm in names(occ.covs)) {

      bad <- is.na(occ.covs[[nm]]) & site_observed

      if (any(bad)) {

        if (is.numeric(occ.covs[[nm]])) {

          occ.covs[[nm]][bad] <-
            mean(occ.covs[[nm]], na.rm = TRUE)

        } else {

          occ.covs[[nm]][bad] <-
            names(which.max(table(occ.covs[[nm]])))

        }
      }
    }

    # -------------------------------------------------------------------------
    # Formulas
    # -------------------------------------------------------------------------

    occ_rhs <- make_rhs(
      covariates = occ_covariates,
      use_route_random_effect = use_route_random_effect
    )

    det_rhs <- make_rhs(
      covariates = det_covariates,
      use_route_random_effect = FALSE
    )

    occ_formula <- as.formula(paste("~", occ_rhs))
    det_formula <- as.formula(paste("~", det_rhs))

    # -------------------------------------------------------------------------
    # Spatial priors
    # -------------------------------------------------------------------------

    dists <- stats::dist(coords)
    max_d <- max(dists)
    min_d <- stats::quantile(dists, 0.05)

    phi_lo <- 3 / max_d
    phi_hi <- 3 / min_d

    # -------------------------------------------------------------------------
    # Data list
    # -------------------------------------------------------------------------

    data.list <- list(
      y        = y,
      occ.covs = occ.covs,
      det.covs = det.covs,
      coords   = coords
    )

    # -------------------------------------------------------------------------
    # Priors and initial values
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # Number of coefficients, accounting for factor variables
    # -------------------------------------------------------------------------

    # Get number of parameters to be estimated in the occupancy model
    p_occ <- ncol(model.matrix(occ_formula, data = occ.covs))
    if(use_route_random_effect){
      p_occ <- p_occ - 1 # Remove the random effect paramether that is added
    }


    # Get number of parameters to be estimated in the detection model
    det_cov_df <- as.data.frame(
      lapply(det.covs, as.vector)
    )

    det_cov_df <- det_cov_df |>
      dplyr::mutate(
        dplyr::across(
          where(is.character),
          as.factor
        )
      )

    p_det <- ncol(model.matrix(det_formula, data = det_cov_df))

    priors <- list(
      beta.comm.normal  = list(mean = rep(0, p_occ), var = rep(2.72, p_occ)),
      tau.sq.beta.ig    = list(a = 0.1, b = 0.1),
      alpha.comm.normal = list(mean = rep(0, p_det), var = rep(2.72, p_det)),
      tau.sq.alpha.ig   = list(a = 0.1, b = 0.1),
      sigma.sq.ig       = list(a = 2, b = 1),
      phi.unif          = list(a = phi_lo, b = phi_hi)
    )

    if (use_route_random_effect) {
      priors$sigma.sq.psi.ig <- list(a = 2, b = 0.5)
    }

    z.init <- apply(y, c(1, 2), function(x) {
      as.integer(any(x == 1, na.rm = TRUE))
    })

    inits <- list(
      beta.comm    = rep(0, p_occ),
      alpha.comm   = rep(0, p_det),
      tau.sq.beta  = rep(1, p_occ),
      tau.sq.alpha = rep(1, p_det),
      sigma.sq     = rep(1, N),
      phi          = rep(3 / (max_d / 2), N),
      z            = z.init
    )

    if (use_route_random_effect) {
      inits$sigma.sq.psi <- 1
    }
    # -------------------------------------------------------------------------
    # Fit model
    # -------------------------------------------------------------------------

    out <- NULL

    if (fit_models) {

      message("Fitting model for guild: ", guild)

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
        verbose       = TRUE,
        n.report      = max(1, floor(n.batch / 10))
      )
    }

    # -------------------------------------------------------------------------
    # Extract model products
    # -------------------------------------------------------------------------

    model_products <- list()

    if (!is.null(out)) {

      ppc <- spOccupancy::ppcOcc(
        out,
        fit.stat = "freeman-tukey",
        group = 1
      )

      bp_species <- colMeans(ppc$fit.y.rep > ppc$fit.y, na.rm = TRUE)
      bp_overall <- mean(ppc$fit.y.rep > ppc$fit.y, na.rm = TRUE)

      waic <- spOccupancy::waicOcc(out)

      beta_comm_summary <- post_summary(out$beta.comm.samples)
      beta_comm_summary$parameter <- c("(Intercept)", occ_covariates)

      alpha_comm_summary <- post_summary(out$alpha.comm.samples)
      alpha_comm_summary$parameter <- c("(Intercept)", det_covariates)

      occ_effects <- summarize_all_named_effects(
        samples = out$beta.samples,
        covariates = occ_covariates,
        process = "occupancy"
      )

      det_effects <- summarize_all_named_effects(
        samples = out$alpha.samples,
        covariates = det_covariates,
        process = "detection"
      )

      sigma_cols <- grep("^sigma.sq-", colnames(out$theta.samples), value = TRUE)
      phi_cols   <- grep("^phi-", colnames(out$theta.samples), value = TRUE)

      spatial_params <- NULL

      if (length(sigma_cols) > 0 && length(phi_cols) > 0) {
        sigma_mean <- colMeans(out$theta.samples[, sigma_cols, drop = FALSE])
        phi_mean   <- colMeans(out$theta.samples[, phi_cols, drop = FALSE])

        spatial_params <- data.frame(
          ecosystem = ecosystemType,
          guild = guild,
          Species = sub("^sigma.sq-", "", sigma_cols),
          sigma_sq = sigma_mean,
          phi = phi_mean,
          eff_range = 3 / phi_mean,
          row.names = NULL
        )
      }

      species <- sort(unique(dat$Species))
      psi.mean <- apply(out$psi.samples, c(2, 3), mean)
      dimnames(psi.mean) <- list(species, points)

      psi.ci <- apply(
        out$psi.samples,
        c(2, 3),
        function(x) stats::quantile(x, c(0.025, 0.975))
      )

      w.mean <- apply(out$w.samples, c(2, 3), mean)
      dimnames(w.mean) <- list(species, points)

      richness.samples <- apply(out$z.samples, c(1, 3), sum)
      richness.mean <- colMeans(richness.samples)

      richness.ci <- apply(
        richness.samples,
        2,
        function(x) stats::quantile(x, c(0.025, 0.975))
      )

      names(richness.mean) <- points

      richness_df <- data.frame(
        ecosystem = ecosystemType,
        guild = guild,
        Point = points,
        Route = point_info$route_id,
        Lon = point_info$Lon,
        Lat = point_info$Lat,
        richness_mean = richness.mean,
        richness_lwr = richness.ci[1, ],
        richness_upr = richness.ci[2, ],
        row.names = NULL
      ) %>%
        dplyr::left_join(
          point_info_raw,
          by = c("Point" = "Route_point_id")
        )

      psi_long <- as.data.frame(as.table(psi.mean)) %>%
        dplyr::rename(
          Species = Var1,
          Point = Var2,
          psi_mean = Freq
        ) %>%
        dplyr::mutate(
          ecosystem = ecosystemType,
          guild = guild
        ) %>%
        dplyr::left_join(
          point_info_raw,
          by = c("Point" = "Route_point_id")
        )

      w_long <- as.data.frame(as.table(w.mean)) %>%
        dplyr::rename(
          Species = Var1,
          Point = Var2,
          w_mean = Freq
        ) %>%
        dplyr::mutate(
          ecosystem = ecosystemType,
          guild = guild
        ) %>%
        dplyr::left_join(
          point_info_raw,
          by = c("Point" = "Route_point_id")
        )

      ppc_df <- data.frame(
        ecosystem = ecosystemType,
        guild = guild,
        Species = ppc$sp.names,
        bayesian_p = bp_species,
        row.names = NULL
      )

      # Optional convergence summaries
      rhat_beta_comm <- tryCatch(
        check_rhat(
          out$beta.comm.samples,
          c("(Intercept)", occ_covariates),
          n.chains
        ),
        error = function(e) NULL
      )

      rhat_alpha_comm <- tryCatch(
        check_rhat(
          out$alpha.comm.samples,
          c("(Intercept)", det_covariates),
          n.chains
        ),
        error = function(e) NULL
      )

      model_products <- list(
        model = out,
        ppc = ppc,
        ppc_species = ppc_df,
        ppc_overall = bp_overall,
        waic = waic,
        beta_comm_summary = beta_comm_summary,
        alpha_comm_summary = alpha_comm_summary,
        occ_effects = occ_effects,
        det_effects = det_effects,
        spatial_params = spatial_params,
        psi_mean = psi.mean,
        psi_ci = psi.ci,
        psi_long = psi_long,
        w_mean = w.mean,
        w_long = w_long,
        richness_samples = richness.samples,
        richness_mean = richness.mean,
        richness_ci = richness.ci,
        richness_df = richness_df,
        rhat_beta_comm = rhat_beta_comm,
        rhat_alpha_comm = rhat_alpha_comm
      )
    }

    all_results[[guild]] <- list(
      ecosystem = ecosystemType,
      guild = guild,
      species = species,
      points = points,
      routes = routes,
      years = years,
      dims = list(N = N, J = J, R = R),
      removed_routes_missing_coords = routes_with_missing_coords,
      formulas = list(
        occ_formula = occ_formula,
        det_formula = det_formula
      ),
      data_used = dat,
      data_list = data.list,
      point_info_raw = point_info_raw,
      point_info_model = point_info,
      coords = coords,
      priors = priors,
      inits = inits,
      model_products = model_products
    )

    if (save_results) {
      saveRDS(
        all_results[[guild]],
        file = file.path(
          output_dir,
          paste0(
            "spMsPGOcc_",
            ecosystemType,
            "_",
            gsub("[^A-Za-z0-9]+", "_", guild),
            ".rds"
          )
        )
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Combine across guilds
  # ---------------------------------------------------------------------------

  get_prod <- function(x, nm) {
    if (length(x$model_products) == 0) return(NULL)
    x$model_products[[nm]]
  }

  combined <- list(
    richness = dplyr::bind_rows(lapply(all_results, get_prod, "richness_df")),
    psi = dplyr::bind_rows(lapply(all_results, get_prod, "psi_long")),
    spatial_random_effects = dplyr::bind_rows(lapply(all_results, get_prod, "w_long")),
    ppc = dplyr::bind_rows(lapply(all_results, get_prod, "ppc_species")),
    spatial_params = dplyr::bind_rows(lapply(all_results, get_prod, "spatial_params")),
    occ_effects = dplyr::bind_rows(lapply(all_results, get_prod, "occ_effects")),
    det_effects = dplyr::bind_rows(lapply(all_results, get_prod, "det_effects"))
  )

  return(list(
    ecosystem = ecosystemType,
    occ_covariates = occ_covariates,
    det_covariates = det_covariates,
    guild_results = all_results,
    combined = combined
  ))
}

