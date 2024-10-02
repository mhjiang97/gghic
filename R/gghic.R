#' @importFrom S4Vectors mcols
#' @importFrom InteractionSet pairdist interactions
#' @importFrom tibble as_tibble
#' @importFrom scales oob_squish
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

#' gghic
#'
#' @description A ggplot2 wrapper for HiC data.
#' @param data A HiCExperiment object or a tibble/data.frame
#' @param score_column The column name of which the score is calculated.
#'   Default is `"balanced"`.
#' @param scale_method The function to scale the score. Default is `log10`.
#' @param ideogram Whether to add ideogram or not. Default is `TRUE`.
#' @param ... Other parameters passed to [geom_ideogram()].
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' library(HiCExperiment)
#' }
#' @importFrom ggplot2 ggplot aes coord_fixed scale_x_continuous
#' @importFrom tibble is_tibble
#' @export gghic
gghic <- function(
  data = NULL, score_column = "balanced", scale_method = log10,
  ideogram = TRUE,
  ...
) {
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

  # params <- list(...)
  # names_param_ideogram <- c("genome", "highlight", "fontsize")

  p <- x |>
    ggplot2::ggplot(
      ggplot2::aes(
        seqnames1 = seqnames1, start1 = start1, end1 = end1,
        seqnames2 = seqnames2, start2 = start2, end2 = end2,
        fill = score
      )
    ) +
      geom_hic() +
      theme_hic()

  if (ideogram) {
    # p <- p + do.call(
    #   geom_ideogram, params[names(params) %in% names_param_ideogram]
    # )
    p <- p + geom_ideogram(...)
  }

  p
}
