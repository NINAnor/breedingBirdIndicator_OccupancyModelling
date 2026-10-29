# =============================================================================
# Regional and national species predictions: medians, interval widths,
# model choice and naive occupancy
# =============================================================================

library(sf)
library(dplyr)
library(readr)

# Set to TRUE to ignore cached results and recompute everything
recompute <- FALSE

regions <- c("Central", "East", "North", "South", "West")


# =============================================================================
# Helper functions
# =============================================================================

short_path <- function(path) {
  result <- system(
    sprintf(
      'cmd /c for %%I in ("%s") do @echo %%~sI',
      path
    ),
    intern = TRUE
  )

  result <- trimws(result)

  if (length(result) != 1L || !nzchar(result)) {
    stop(
      "Could not convert path to Windows short path:\n",
      path,
      call. = FALSE
    )
  }

  result
}


# Find prediction RDS files for one model
get_prediction_files <- function(model_path, species = NA_character_) {
  search_path <- if (length(species) == 1L && is.na(species)) {
    model_path
  } else {
    file.path(model_path, species)
  }

  if (!dir.exists(search_path)) {
    warning(
      "Prediction directory does not exist:\n",
      search_path,
      call. = FALSE
    )
    return(character(0))
  }

  list.files(
    search_path,
    pattern = "\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )
}


# -----------------------------------------------------------------------------
# Nature region of each prediction point.
# All species share the same prediction grid, so the (slow) spatial join is done
# once and reused for every file whose coordinates are identical.
# -----------------------------------------------------------------------------

region_cache <- new.env(parent = emptyenv())

assign_regions <- function(points, bio_regions) {
  if (sf::st_crs(bio_regions) != sf::st_crs(points)) {
    bio_regions <- sf::st_transform(bio_regions, sf::st_crs(points))
  }
  bio_regions <- sf::st_make_valid(bio_regions)

  hits <- sf::st_intersects(points, bio_regions)
  # A point on a shared border matches two regions: keep the first match
  first <- vapply(
    hits,
    function(z) if (length(z)) z[[1L]] else NA_integer_,
    integer(1)
  )
  as.character(bio_regions$nature_region)[first]
}


# Read one prediction file and return a plain data frame (no geometry)
read_prediction_file <- function(rds_file, bio_regions) {
  dat <- readRDS(short_path(rds_file))

  if (!inherits(dat, "sf")) {
    stop("Expected an sf object in: ", rds_file, call. = FALSE)
  }

  needed <- c("species", "median", "lower", "upper")
  if (!all(needed %in% names(dat))) {
    stop(
      "Missing column(s) in ",
      rds_file,
      ": ",
      paste(setdiff(needed, names(dat)), collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(dat) == 0L) {
    stop("Prediction file contains no rows: ", rds_file, call. = FALSE)
  }

  xy <- sf::st_coordinates(dat)

  same_grid <- !is.null(region_cache$region) &&
    nrow(xy) == length(region_cache$region) &&
    sf::st_crs(dat) == region_cache$crs &&
    all(xy[, 1] == region_cache$x & xy[, 2] == region_cache$y)

  if (same_grid) {
    region <- region_cache$region
  } else {
    if (!is.null(region_cache$region)) {
      message("  Grid differs from the previous file; recomputing regions.")
    }
    region <- assign_regions(sf::st_geometry(dat), bio_regions)
    region_cache$x <- xy[, 1]
    region_cache$y <- xy[, 2]
    region_cache$crs <- sf::st_crs(dat)
    region_cache$region <- region
  }

  n_na <- sum(is.na(region))
  if (n_na > 0L) {
    warning(
      n_na,
      " points fall outside all bio-regions in: ",
      rds_file,
      call. = FALSE
    )
  }

  data.frame(
    species = as.character(dat$species),
    nature_region = region,
    p_median = dat$median,
    p_width = dat$upper - dat$lower,
    stringsAsFactors = FALSE
  )
}


# -----------------------------------------------------------------------------
# Estimators
# -----------------------------------------------------------------------------

sample_idx <- function(n, max_n, seed) {
  if (n <= max_n) {
    return(seq_len(n))
  }
  set.seed(seed)
  sample.int(n, max_n)
}

# Sample up to max_n rows in total, spread evenly over the regions
stratified_idx <- function(region, max_n, seed) {
  groups <- split(seq_along(region), region)
  per <- ceiling(max_n / length(groups))
  set.seed(seed)
  unlist(
    lapply(groups, function(i) if (length(i) > per) sample(i, per) else i),
    use.names = FALSE
  )
}

# Intercept-only beta regression; returns plogis(intercept).
# region = NULL     -> y ~ 1
# region = a vector -> y ~ 1 + (1 | region)  (typical region, each region equal)
fit_beta_mean <- function(y, region = NULL) {
  ok <- !is.na(y)
  y <- y[ok]
  if (!is.null(region)) {
    region <- region[ok]
  }
  n <- length(y)
  if (n < 10L) {
    return(NA_real_)
  }

  # Beta needs values strictly inside (0, 1): Smithson & Verkuilen squeeze
  d <- data.frame(y = (y * (n - 1) + 0.5) / n)
  frm <- y ~ 1
  if (!is.null(region)) {
    d$region <- factor(region)
    frm <- y ~ 1 + (1 | region)
  }

  fit <- tryCatch(
    suppressWarnings(
      glmmTMB::glmmTMB(frm, data = d, family = glmmTMB::beta_family())
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || !isTRUE(fit$sdr$pdHess)) {
    return(NA_real_)
  }
  stats::plogis(glmmTMB::fixef(fit)$cond[["(Intercept)"]])
}


# One prediction file -> one row of regional and national estimates
summarise_prediction <- function(
  d,
  method = c("beta", "median"),
  max_n = 20000L,
  seed = 1L
) {
  method <- match.arg(method)

  sp <- unique(d$species)
  if (length(sp) != 1L) {
    stop(
      "Expected one species per prediction file, found: ",
      paste(sp, collapse = ", "),
      call. = FALSE
    )
  }

  # Estimate for the median and the width of the same set of points
  fit_pair <- function(med, wid, region = NULL) {
    if (method == "median") {
      return(c(
        stats::median(med, na.rm = TRUE),
        stats::median(wid, na.rm = TRUE)
      ))
    }
    idx <- if (is.null(region)) {
      sample_idx(length(med), max_n, seed)
    } else {
      stratified_idx(region, max_n, seed)
    }
    c(
      fit_beta_mean(med[idx], region[idx]),
      fit_beta_mean(wid[idx], region[idx])
    )
  }

  out <- list(species = sp, method = method)

  for (r in regions) {
    rows <- which(d$nature_region == r)
    est <- if (length(rows)) {
      fit_pair(d$p_median[rows], d$p_width[rows])
    } else {
      c(NA_real_, NA_real_)
    }
    out[[paste0(r, "_median")]] <- est[1]
    out[[paste0(r, "_width")]] <- est[2]
  }

  if (method == "median") {
    # Median over all points, including those without a region
    est <- fit_pair(d$p_median, d$p_width)
  } else {
    # Random intercept for region
    rows <- which(!is.na(d$nature_region))
    est <- fit_pair(d$p_median[rows], d$p_width[rows], d$nature_region[rows])
  }
  out$National_median <- est[1]
  out$National_width <- est[2]

  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}


# All prediction files of one model: read, summarise, free memory, one at a time
process_model_files <- function(
  rds_files,
  model_type,
  method,
  bio_regions,
  max_n = 20000L,
  seed = 1L
) {
  res <- lapply(seq_along(rds_files), function(k) {
    f <- rds_files[k]
    cat(sprintf("  [%d/%d] %s\n", k, length(rds_files), basename(f)))
    d <- read_prediction_file(f, bio_regions)
    out <- summarise_prediction(d, method, max_n, seed)
    out$model_type <- model_type
    rm(d)
    gc(verbose = FALSE)
    out
  })
  dplyr::bind_rows(res)
}


# Save a data frame as an Excel workbook in <output_root>/median_vals
save_xlsx <- function(x, name) {
  f <- file.path(
    output_root,
    "median_vals",
    paste0(ecosystemType, "_", name, ".xlsx")
  )
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  openxlsx::write.xlsx(x, file = f, overwrite = TRUE, na.string = "")
  cat(
    "\nSaved: ",
    normalizePath(f, winslash = "/", mustWork = FALSE),
    " (",
    nrow(x),
    " rows)\n",
    sep = ""
  )
  invisible(f)
}


# =============================================================================
# Inputs
# =============================================================================

best_models <- readr::read_csv(
  file.path(
    output_root,
    "postProcessing_outputs",
    "multipleOutputs",
    ecosystemType,
    "best_models.csv"
  ),
  show_col_types = FALSE
)

best_models$prediction_path <- gsub(
  "modelFit",
  "prediction_outputs",
  best_models$path,
  fixed = TRUE
)

best_models$prediction_files <- lapply(
  seq_len(nrow(best_models)),
  function(i) {
    get_prediction_files(
      model_path = best_models$prediction_path[i],
      species = best_models$Species[i]
    )
  }
)

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

# Model outputs (Bayesian p-values, etc.)
load(file.path(path_to_postProcessing, "model_outputs.RData"))


# =============================================================================
# Regional and national estimates for every best model
# (cached, so the slow step is skipped on re-runs)
# =============================================================================

dir.create(
  file.path(output_root, "median_vals"),
  recursive = TRUE,
  showWarnings = FALSE
)

results_cache <- file.path(
  output_root,
  "median_vals",
  paste0(ecosystemType, "_median_results.rds")
)

if (file.exists(results_cache) && !recompute) {
  cat("Reading cached estimates:", results_cache, "\n")
  median_results <- readRDS(results_cache)
} else {
  median_results <- dplyr::bind_rows(lapply(
    seq_len(nrow(best_models)),
    function(i) {
      rds_files <- best_models$prediction_files[[i]]
      model_type <- if (isTRUE(best_models$single_species[i])) {
        "single_species"
      } else {
        "multi_species"
      }

      if (length(rds_files) == 0L) {
        warning(
          "No prediction RDS files found for model row ",
          i,
          " (",
          model_type,
          ").",
          call. = FALSE
        )
        return(NULL)
      }

      cat(
        "\n",
        model_type,
        ": row ",
        i,
        " of ",
        nrow(best_models),
        "\n",
        sep = ""
      )
      process_model_files(
        rds_files,
        model_type = model_type,
        method = "beta",
        bio_regions = bio_regions
      )
    }
  ))
  saveRDS(median_results, results_cache)
}

if (nrow(median_results) == 0L) {
  stop("No median results were generated. Check the prediction files.")
}

# One row per species and model type is expected
dup <- median_results |>
  dplyr::count(model_type, species) |>
  dplyr::filter(n > 1L)
if (nrow(dup) > 0L) {
  stop(
    "More than one prediction file per species and model type: ",
    paste(dup$species, collapse = ", "),
    ". Keep one file per species before summarising."
  )
}


# =============================================================================
# Medians and widths as separate tables
# =============================================================================

split_stat <- function(x, stat) {
  suffix <- paste0("_", stat)
  x |>
    dplyr::filter(method %in% "beta") |>
    dplyr::select(model_type, species, method, dplyr::ends_with(suffix)) |>
    dplyr::arrange(model_type, species) |>
    dplyr::rename_with(
      ~ sub(paste0(suffix, "$"), "", .x),
      dplyr::ends_with(suffix)
    )
}

all_species_medians <- split_stat(median_results, "median")
all_species_width <- split_stat(median_results, "width")

save_xlsx(all_species_medians, "medians")
save_xlsx(all_species_width, "width")

print(all_species_medians)
print(all_species_width)


# =============================================================================
# Bayesian p-values and naive occupancy
# =============================================================================

bayes_P_best <- outputs$bayes_P |>
  dplyr::semi_join(best_models, by = c("path", "file"))

nrow(bayes_P_best) # expect 32
dplyr::count(bayes_P_best, Species) # expect 2 rows per species
dplyr::count(bayes_P_best, single_species, spatial)


naive_summary <- function(rds_file, route_pattern = "_.*$") {
  output <- readRDS(short_path(rds_file))
  y <- output$data_list$y
  model <- output$model_products$model

  if (is.null(y)) {
    stop("No `output$data_list$y` found in: ", rds_file, call. = FALSE)
  }

  is_single <- length(dim(y)) == 2L

  if (is_single) {
    # Single species: sites x replicates (years). Keep site and year names.
    sp <- if (!is.null(output$species)) {
      output$species
    } else {
      sub(
        "^(spPGOcc|PGOcc)_[^_]+_",
        "",
        tools::file_path_sans_ext(basename(rds_file))
      )
    }
    y <- array(
      y,
      dim = c(1L, dim(y)),
      dimnames = list(sp, rownames(y), colnames(y))
    )
  }
  # Now y is species x sites x replicates

  # Route id = site name without the point number, e.g. "1001_15" -> "1001"
  sites <- dimnames(y)[[2]]
  route <- if (!is.null(sites)) sub(route_pattern, "", sites) else NULL

  stat_names <- c(
    "n_sites",
    "n_sites_detected",
    "naive_occupancy",
    "naive_detection",
    "n_routes",
    "n_routes_detected",
    "naive_occupancy_route"
  )

  per_species <- vapply(
    seq_len(dim(y)[1]),
    function(i) {
      yi <- matrix(y[i, , ], nrow = dim(y)[2]) # sites x replicates

      surveyed <- rowSums(!is.na(yi)) > 0
      detected <- surveyed & rowSums(yi > 0, na.rm = TRUE) > 0
      seen <- yi[detected, , drop = FALSE]

      n_routes <- n_routes_detected <- NA_real_
      if (!is.null(route)) {
        n_routes <- length(unique(route[surveyed]))
        n_routes_detected <- length(unique(route[detected]))
      }

      c(
        n_sites = sum(surveyed),
        n_sites_detected = sum(detected),
        naive_occupancy = sum(detected) / sum(surveyed),
        naive_detection = sum(seen > 0, na.rm = TRUE) / sum(!is.na(seen)),
        n_routes = n_routes,
        n_routes_detected = n_routes_detected,
        naive_occupancy_route = n_routes_detected / n_routes
      )
    },
    FUN.VALUE = stats::setNames(numeric(length(stat_names)), stat_names)
  )

  species <- if (!is.null(dimnames(y)[[1]])) {
    dimnames(y)[[1]]
  } else {
    model$sp.names
  }

  data.frame(
    file = basename(rds_file),
    species = gsub("_", " ", species),
    single_species = is_single,
    spatial = inherits(model, c("spPGOcc", "spMsPGOcc")),
    t(per_species),
    row.names = NULL
  )
}

naive_occupancy_by_species <- function(rds_files, route_pattern = "_.*$") {
  res <- lapply(rds_files, function(f) {
    tryCatch(
      naive_summary(f, route_pattern),
      error = function(e) {
        warning(
          "Skipped ",
          basename(f),
          ": ",
          conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )
  })
  dplyr::bind_rows(res) |>
    dplyr::arrange(species, file)
}

# Naive occupancy for every selected model file (cached: the model files are
# large and only the data matrix is needed)
naive_cache <- file.path(
  output_root,
  "median_vals",
  paste0(ecosystemType, "_naive.rds")
)

if (file.exists(naive_cache) && !recompute) {
  naive <- readRDS(naive_cache)
} else {
  model_files <- unique(file.path(bayes_P_best$path, bayes_P_best$file))
  naive <- naive_occupancy_by_species(model_files)
  saveRDS(naive, naive_cache)
}


# =============================================================================
# Choose the better model for each species
# =============================================================================

pick_model <- function(d) {
  d <- d[order(d$national_width), ]
  if (sum(d$fit_ok) == 1L) {
    return(d[d$fit_ok, ])
  }
  if (sum(d$fit_ok) == 2L) {
    if (abs(diff(d$national_width)) < 0.05) {
      d <- d[order(d$dist_half), ]
    }
    return(d[1, ])
  }
  d[order(d$dist_half, d$national_width), ][1, ]
}

cand <- bayes_P_best |>
  dplyr::transmute(
    species = Species,
    file,
    model_type = ifelse(single_species, "single_species", "multi_species"),
    bayes_p = Mean
  ) |>
  dplyr::left_join(
    all_species_width |>
      dplyr::select(model_type, species, national_width = National),
    by = c("model_type", "species")
  ) |>
  dplyr::left_join(
    all_species_medians |>
      dplyr::select(model_type, species, national_median = National),
    by = c("model_type", "species")
  ) |>
  dplyr::mutate(
    fit_ok = bayes_p >= 0.1 & bayes_p <= 0.9,
    dist_half = abs(bayes_p - 0.5)
  ) |>
  dplyr::group_by(species) |>
  dplyr::mutate(models_disagree = diff(range(national_median)) > 0.2) |>
  dplyr::ungroup()

choice <- cand |>
  dplyr::group_by(species) |>
  dplyr::group_modify(~ pick_model(.x)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    status = dplyr::case_when(
      !fit_ok ~ "Poor fit",
      national_width > 0.5 ~ "Wide uncertainty",
      TRUE ~ "OK"
    )
  )

final_medians <- choice |>
  dplyr::select(
    species,
    file,
    model_type,
    bayes_p,
    national_width,
    models_disagree,
    status
  ) |>
  dplyr::left_join(all_species_medians, by = c("species", "model_type")) |>
  # Naive occupancy from the same model file that was chosen
  dplyr::left_join(
    naive |>
      dplyr::select(
        species,
        file,
        n_sites_detected,
        naive_occupancy,
        naive_detection,
        naive_occupancy_route
      ),
    by = c("species", "file")
  ) |>
  dplyr::select(
    species,
    Central,
    East,
    North,
    South,
    West,
    National,
    status,
    models_disagree,
    bayes_p,
    model_type,
    national_width,
    n_sites_detected,
    naive_occupancy,
    naive_detection,
    naive_occupancy_route
  )

save_xlsx(final_medians, "final_medians")
print(final_medians)
