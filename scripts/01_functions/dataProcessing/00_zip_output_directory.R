#' Compress an output directory
#'
#' @param output_dir Directory to compress.
#' @param zip_file Optional path to the zip file.
#' @param overwrite Logical. Overwrite existing zip file?
#'
#' @return Path to the zip file.
#'
#' @export
zip_output_directory <- function(
    output_dir,
    zip_file = NULL,
    overwrite = TRUE
) {

  stopifnot(dir.exists(output_dir))

  if (is.null(zip_file)) {
    zip_file <- paste0(
      normalizePath(output_dir),
      ".zip"
    )
  }

  if (file.exists(zip_file)) {

    if (!overwrite) {
      stop(
        zip_file,
        " already exists."
      )
    }

    file.remove(zip_file)
  }

  oldwd <- setwd(output_dir)
  on.exit(setwd(oldwd), add = TRUE)

  zip::zipr(
    zipfile = zip_file,
    files = list.files(
      ".",
      recursive = TRUE,
      all.files = FALSE
    )
  )

  zip_file
}
