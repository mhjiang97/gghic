scale_data <- function(
    hicexperiment, score_column = "balanced", scale_method = log10
) {
  gis <- InteractionSet::interactions(hicexperiment)
  gis$score <- method(S4Vectors::mcols(gis)[, score_column])

  x <- tibble::as_tibble(gis)
  scores <- x$score[
    InteractionSet::pairdist(gis) != 0 &
      !is.na(InteractionSet::pairdist(gis) != 0)
  ]
  scores <- scores[!is.na(scores) & !is.infinite(scores)]
  M <- max(scores)
  m <- min(scores)

  x$score <- scales::oob_squish(x$score, c(m, M))

  x
}

ensure_data <- function(data) {
  if (is(data, "HiCExperiment")) {
    x <- scale_data(
      data, score_column = score_column, scale_method = scale_method
    )
  } else if (tibble::is_tibble(data) || is(data, "data.frame")) {
    cols_required <- c(
      "seqnames1", "seqnames2", "start1", "end1", "start2", "end2", "score"
    )
    cols_missing <- setdiff(cols_required, colnames(x))
    if (length(cols_missing) > 0) {
      stop(
        "data must have the following columns: ",
        paste(cols_missing, collapse=', ')
      )
    }
    x <- data
  } else {
    stop("data must be a HiCExperiment object or a tibble/data.frame")
  }
  x
}

#' gghic
#'
#' @description A ggplot2 wrapper for HiC data.
#' @param data A HiCExperiment object or a tibble/data.frame
#' @param score_column The column name of which the score is calculated.
#'   Default is `"balanced"`.
#' @param scale_method The function to scale the score. Default is `log10`.
#' @param ideogram Whether to add ideogram or not. Default is `TRUE`.
#' @param annotation Whether to add annotation or not. Default is `TRUE`.
#' @param ... Other parameters passed to [geom_ideogram()], [geom_annotation()].
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' library(HiCExperiment)
#' }
#' @export gghic
gghic <- function(
  data = NULL, score_column = "balanced", scale_method = log10,
  ideogram = TRUE, annotation = TRUE,
  ...
) {
  p <- data |>
    ensure_data() |>
    tidyr::drop_na(score) |>
    ggplot2::ggplot(
      ggplot2::aes(
        seqnames1 = seqnames1, start1 = start1, end1 = end1,
        seqnames2 = seqnames2, start2 = start2, end2 = end2,
        fill = score
      )
    ) +
      geom_hic() +
      theme_hic()

  if (ideogram) p <- p + geom_ideogram(...)


  if (annotation) p <- p + geom_annotation(...)

  p
}
