name_pkg <- "gghic"

dir_cache <- rappdirs::user_cache_dir(appname = name_pkg)

ensure_dir <- function(paths) {
  for (path in paths) {
    if (!dir.exists(path)) dir.create(path)
  }
}

stop_if_null <- function(object, message) {
  if (is.null(object)) stop(message)
}
