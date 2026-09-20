################################################################################
# Prepare environmental covariates for occupancy modelling
#
# This script:
#   1. Loads all environmental covariates required by the model.
#   2. Creates a common raster template to ensure consistent resolution and CRS.
#   3. Builds a covariate raster stack.
#   4. Extracts covariate values at bird survey locations.
#   5. Extracts covariate values at prediction locations.
#
# All covariates are aligned to the same grid before extraction to ensure
# consistency between model fitting and prediction.
################################################################################

################################################################################
# Define paths
################################################################################

# Folder containing all raster covariates used in the model.
# Keeping this path centralized improves portability and reproducibility.
cov_folder <- file.path(
  dataFolder,
  "Covariates"
)


################################################################################
# Check that model covariates have been defined
################################################################################

# Detection and occupancy covariates should have been specified in the master
# script. Stop execution if they are unavailable.

required_objects <- c("det_covs", "occ_covs")

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    envir = .GlobalEnv
  )
]

if (length(missing_objects) > 0) {
  stop(
    paste0(
      "The following objects are missing from the global environment:\n",
      paste(missing_objects, collapse = ", "),
      "\n\nPlease run the master script before executing this section."
    ),
    call. = FALSE
  )
}


################################################################################
# Combine detection and occupancy covariates
################################################################################

# Create a unique list of all environmental variables needed for the model.

if (!exists("det_covs") || !exists("occ_covs")) {
  stop(
    "Objects 'det_covs' and 'occ_covs' must be defined. ",
    "Please run the master script first.",
    call. = FALSE
  )
}

all_covariates <- unique(c(det_covs, occ_covs))


################################################################################
# Build the covariate raster stack
################################################################################

# Load existing raster layers or recreate them if their geometry or variable
# names do not match the template raster.

cov_stack <- put_covariates_together(
  x = ecosystemBreedingBird,
  modelCrs = modelCrs,
  all_covariates = all_covariates,
  template = baseRaster,
  output_file = file.path(
    dataFolder,
    "processedData",
    "covariate_stack",
    paste0("covariate_stack_", ecosystemType, ".tif")
  ),
  cov_folder = cov_folder,
  data_folder = dataFolder,
  norway_boundary = norway_boundary,
  modelYear = modelYear,
  overwrite = FALSE
)


################################################################################
# Extract covariates at bird survey locations
################################################################################

# Extract environmental values at each survey point.
# These data will be used for model fitting.

breeding_bird_data_for_modelling <-
  extract_model_covariates(
    x = ecosystemBreedingBird %>%
      dplyr::select(-DEM),
    modelCrs = modelCrs,
    all_covariates = all_covariates,
    lon_col = "Longitude",
    lat_col = "Latitude",
    input_crs = 4326,
    covariate_stack = cov_stack,
    return = "data.frame"
  ) %>%

  # Split Route_point_id into route and point identifiers
  dplyr::mutate(
    route_id = sub("_.*", "", Route_point_id),
    point_id = sub(".*_", "", Route_point_id)
  )

# --------------------------------------------------------
# Apply categorical transformations only when requested
# --------------------------------------------------------

if ("dominantTree" %in% all_covariates) {
  breeding_bird_data_for_modelling <-
    breeding_bird_data_for_modelling %>%
    dplyr::mutate(
      dominantTree = dplyr::if_else(
        is.na(dominantTree),
        4,
        dominantTree
      ),
      dominantTree = factor(
        dominantTree
      )
    )
}

if ("glad_landcover" %in% all_covariates) {
  breeding_bird_data_for_modelling <-
    breeding_bird_data_for_modelling %>%
    dplyr::mutate(
      glad_landcover = factor(
        glad_landcover
      )
    )
}


if ("forestType" %in% all_covariates) {
  breeding_bird_data_for_modelling <-
    breeding_bird_data_for_modelling %>%
    dplyr::mutate(
      forestType = factor(
        forestType
      )
    )
}

# Check if there are missing values

n_total <- nrow(breeding_bird_data_for_modelling)
n_complete <- sum(complete.cases(
  breeding_bird_data_for_modelling %>%
    dplyr::select(
      dplyr::matches(
        paste0("^(", paste(all_covariates, collapse = "|"), ")")
      )
    )
))

dat <- data.frame(
  complete = n_complete,
  incomplete = n_total - n_complete,
  proportion_complete = n_complete / n_total,
  proportion_incomplete = 1 - n_complete / n_total
)

if (dat$complete != n_total) {
  message(paste("There are ", dat$incomplete, "missing records"))

  tt <- breeding_bird_data_for_modelling %>%
    dplyr::group_by(route_id, Year) %>%
    dplyr::mutate(
      dplyr::across(
        where(is.numeric),
        ~ if_else(
          is.na(.x),
          mean(.x, na.rm = TRUE),
          .x
        )
      )
    ) %>%
    dplyr::ungroup()
}

# ============================================================
# Check correlation of variables
# ============================================================

cor_data <- breeding_bird_data_for_modelling %>%

  # Keep one record per Route_point_id
  distinct(Route_point_id, .keep_all = TRUE) %>%

  # Calculate multi-year means
  mutate(
    mean_temperature_avg = rowMeans(
      dplyr::select(
        .,
        starts_with("mean_temperature_")
      ),
      na.rm = TRUE
    ),

    mean_precipitation_avg = rowMeans(
      dplyr::select(
        .,
        starts_with("mean_precipitation_")
      ),
      na.rm = TRUE
    ),

    snowmelt_avg = rowMeans(
      dplyr::select(
        .,
        starts_with("snowmelt_")
      ),
      na.rm = TRUE
    )
  ) %>%

  # Select all original covariates + the new averages
  dplyr::select(
    Route_point_id,
    dplyr::matches(
      paste0(
        "^(",
        paste(all_covariates, collapse = "|"),
        ")"
      )
    ),
    mean_temperature_avg,
    mean_precipitation_avg,
    snowmelt_avg
  )


# ============================================================
# 2. Calculate correlation matrix
# ============================================================

cor_matrix <- cor(
  cor_data %>%
    dplyr::select(-c(Route_point_id, forestType, dominantTree)),
  use = "pairwise.complete.obs",
  method = "pearson"
)

round(cor_matrix, 2)


# ============================================================
# 3. Convert correlation matrix to long format
# ============================================================

cor_long <- as.data.frame(cor_matrix) %>%
  tibble::rownames_to_column("Variable1") %>%
  tidyr::pivot_longer(
    cols = -Variable1,
    names_to = "Variable2",
    values_to = "Correlation"
  )


# ============================================================
# 4. Plot correlation matrix
# ============================================================

ggplot(
  cor_long,
  aes(
    x = Variable1,
    y = Variable2,
    fill = Correlation
  )
) +
  geom_tile(
    colour = "white"
  ) +
  geom_text(
    aes(
      label = sprintf("%.2f", Correlation)
    ),
    size = 3
  ) +
  scale_fill_gradient2(
    limits = c(-1, 1),
    midpoint = 0
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Pearson's r",
    title = "Correlation among environmental covariates"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid = element_blank()
  )


################################################################################
# Extract covariates at prediction locations
################################################################################

# Environmental values are extracted at prediction points so that occupancy
# probabilities can be estimated across the landscape.

# Environmental values are extracted at prediction points so that
# occupancy probabilities can be estimated across the landscape.

predictionPoints_for_modelling <-
  extract_model_covariates(
    x = predictionPoints %>%
      dplyr::select(-dplyr::any_of(all_covariates)),
    modelCrs = modelCrs,
    all_covariates = all_covariates,
    lon_col = "Longitude",
    lat_col = "Latitude",
    input_crs = 4326,
    covariate_stack = cov_stack,
    return = "data.frame"
  )

# --------------------------------------------------------
# Apply dominantTree transformation only if requested
# --------------------------------------------------------

if ("dominantTree" %in% all_covariates) {
  predictionPoints_for_modelling <-
    predictionPoints_for_modelling %>%
    dplyr::mutate(
      dominantTree = dplyr::if_else(
        is.na(dominantTree),
        4,
        dominantTree
      ),
      dominantTree = factor(
        dominantTree,
        levels = levels(
          breeding_bird_data_for_modelling$dominantTree
        )
      )
    )
}

# --------------------------------------------------------
# Apply GLAD land-cover transformation only if requested
# --------------------------------------------------------

if ("glad_landcover" %in% all_covariates) {
  predictionPoints_for_modelling <-
    predictionPoints_for_modelling %>%
    dplyr::mutate(
      glad_landcover = factor(
        glad_landcover,
        levels = levels(
          breeding_bird_data_for_modelling$glad_landcover
        )
      )
    )
}


# --------------------------------------------------------
# Apply forest type transformation only if requested
# --------------------------------------------------------
if ("forestType" %in% all_covariates) {
  predictionPoints_for_modelling <-
    predictionPoints_for_modelling %>%
    dplyr::mutate(
      glad_landcover = factor(
        forestType,
        levels = levels(
          breeding_bird_data_for_modelling$forestType
        )
      )
    )
}
