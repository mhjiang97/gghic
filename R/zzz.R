.onLoad <- function(libname, pkgname) {
  ensure_dir(dir_cache)

  myops <- list(FALSE)
  names(myops) <- glue::glue("{name_pkg}.clean_cache")
  op <- options()
  toset <- !(names(myops) %in% names(op))
  if (any(toset)) options(myops[toset])

  invisible()
}

.onUnload <- function(libpath) {
  if (getOption(glue::glue("{name_pkg}.clean_cache"))) {
    unlink(dir_cache, recursive = TRUE)
  }
}
