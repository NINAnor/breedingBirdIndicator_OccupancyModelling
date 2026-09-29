short_path <- function(path) {
  path <- normalizePath(
    path,
    winslash = "\\",
    mustWork = TRUE
  )

  short <- system2(
    "cmd",
    c("/c", "for %I in (", shQuote(path), ") do @echo %~sI"),
    stdout = TRUE
  )

  short <- trimws(short)

  if (length(short) != 1 || !nzchar(short)) {
    stop("Could not convert path to Windows short path:\n", path)
  }

  short
}

short_path <- function(path) {
  result <- system(
    sprintf(
      'cmd /c for %%I in ("%s") do @echo %%~sI',
      path
    ),
    intern = TRUE
  )

  trimws(result)
}
