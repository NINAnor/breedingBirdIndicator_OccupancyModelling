###############################################################################
# Package Setup Script
#
# Purpose:
#   Ensure all packages required for the FIAR modeling workflow are available
#   and loaded into the current R session.
#
# Typical uses:
#   - Database connections (PostgreSQL/PostGIS)
#   - Spatial data processing
#   - Data wrangling and analysis
#   - Report generation
#   - Interactive mapping
###############################################################################

# Retrieve information on packages currently installed in R.
installedPackages <- installed.packages()

###############################################################################
# List of packages required by the project
###############################################################################

necessaryPackages <- c(
  "RPostgreSQL",  # Legacy PostgreSQL database connectivity
  "RPostgres",    # Modern PostgreSQL database connectivity
  "odbc",         # ODBC database connections
  "DBI",          # Database interface used by PostgreSQL and ODBC packages
  "rpostgis",     # Functions for working with PostGIS spatial databases
  "lubridate",    # Date and time manipulation
  "xtable",       # Export tables to HTML/LaTeX
  "rmarkdown",    # Automated report generation
  "tidyverse",    # Data manipulation, visualization, and analysis tools
  "mapview",       # Interactive spatial mapping
  "sf",
  "arrow",
  "readxl",
  "dplyr",
  "stringr",
  "readr",
  "tibble",
  "spOccupancy",
  "coda",
  "geodata",
  "terra",
  "sfarrow",
  "tidyterra"
)

###############################################################################
# Identify packages that are not currently installed
###############################################################################

uninstalledPackages <- necessaryPackages[
  !(necessaryPackages %in% rownames(installedPackages))
]

###############################################################################
# Install any missing packages
#
# This allows new users to run the workflow without manually installing
# dependencies first.
###############################################################################

if (length(uninstalledPackages) > 0) {

  message(
    "Installing missing packages:\n",
    paste(uninstalledPackages, collapse = ", ")
  )

  install.packages(
    uninstalledPackages,
    dependencies = TRUE
  )
}

###############################################################################
# Function: load_required_packages
#
# Purpose:
#   Verify that all required packages are installed and load them into
#   the current R session.
#
# Arguments:
#   pkgs = character vector of package names
###############################################################################

load_required_packages <- function(pkgs) {

  # Check for any packages that still cannot be found
  missing_pkgs <- pkgs[
    !vapply(
      pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  # Stop execution if installation failed
  if (length(missing_pkgs) > 0) {

    stop(
      "The following packages could not be loaded:\n",
      paste(missing_pkgs, collapse = ", "),
      "\nPlease investigate installation issues."
    )
  }

  # Load all packages into the current session
  invisible(
    lapply(
      pkgs,
      library,
      character.only = TRUE
    )
  )

  message("All required packages loaded successfully.")
}

###############################################################################
# Load all required packages
###############################################################################

load_required_packages(necessaryPackages)
