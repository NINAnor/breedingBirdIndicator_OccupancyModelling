#' Extract model outputs, reading archived model files when needed
#'
#' The files of the models that are not selected as best models are moved to
#' an archive folder (see `02_predictions.R`). This wrapper around
#' [extract_model_outputs()] reads a model folder from the archive when it has
#' no model files in `modelFit`, so archived models can still be used without
#' moving them back.
#'
#' The `path` stored in the outputs is changed back to the `modelFit` folder, so
#' that the rest of the pipeline (best-model selection, prediction folders)
#' sees the usual paths.
#'
#' @param paths Character vector of model folders in `modelFit`.
#' @param modelChecks Logical. Passed to [extract_model_outputs()].
#' @param output_root Folder that contains `modelFit`.
#' @param archive_root Archive folder (same structure as `modelFit`).
#'
#' @return A `model_outputs` object, as from [extract_model_outputs()].
#' @export
extract_model_outputs_with_archive <- function(
  paths,
  modelChecks,
  output_root = get("output_root", envir = globalenv()),
  archive_root = get("archive_root", envir = globalenv())
) {
  norm <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

  fit_root <- norm(file.path(output_root, "modelFit"))
  arch_root <- norm(archive_root)

  has_models <- function(d) {
    dir.exists(d) &&
      length(list.files(d, pattern = "\\.rds$", recursive = TRUE)) > 0
  }

  resolved <- vapply(
    paths,
    function(p) {
      if (has_models(p)) {
        return(p)
      }

      alt <- sub(fit_root, arch_root, norm(p), fixed = TRUE)

      if (has_models(alt)) {
        message("Reading archived model files from: ", alt)
        alt
      } else {
        p
      }
    },
    character(1),
    USE.NAMES = FALSE
  )

  outputs <- extract_model_outputs(resolved, modelChecks)

  # Show the modelFit folder in the outputs for the files read from the archive
  to_fit <- function(x) {
    nx <- norm(x)
    in_archive <- startsWith(nx, arch_root)
    x[in_archive] <- paste0(
      fit_root,
      substring(nx[in_archive], nchar(arch_root) + 1)
    )
    x
  }

  for (nm in names(outputs)) {
    comp <- outputs[[nm]]
    if (is.data.frame(comp) && "path" %in% names(comp)) {
      comp$path <- to_fit(comp$path)
      outputs[[nm]] <- comp
    }
  }

  outputs
}
