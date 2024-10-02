#' @importFrom ggplot2 ggproto Stat
#' @importFrom dplyr mutate
StatHic <- ggplot2::ggproto(
  "StatHic",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1",
    "seqnames2", "start2", "end2",
    "fill"
  ),
  compute_group = function(data, scales) {
    data |>
      dplyr::mutate(
        x = (end1 + start2) / 2,
        xmin = (start1 + start2) / 2,
        xmax = (start1 + end2) / 2,
        xend = (end1 + end2) / 2,
        y = (x - end1),
        ymin = (xmin - start1),
        ymax = (xmax - start1),
        yend = (xend - end1)
      )
  }
)

#' @importFrom ggplot2 ggproto Geom draw_key_polygon
#' @importFrom grid polygonGrob gpar
GeomHic <- ggplot2::ggproto(
  "GeomHic",
  ggplot2::Geom,
  required_aes = c("x", "xmin", "xmax", "xend", "y", "ymin", "ymax", "yend", "fill"),
  draw_key = ggplot2::draw_key_polygon,
  draw_panel = function(data, panel_params, coord) {
    coords <- coord$transform(data, panel_params)
    grid::polygonGrob(
      x = c(coords$x, coords$xmin, coords$xmax, coords$xend),
      y = c(coords$y, coords$ymin, coords$ymax, coords$yend),
      id = rep(seq_len(nrow(coords)), 4),
      default.units = "native",
      gp = grid::gpar(fill = coords$fill, col = NA)
    )
  }
)

#' geom_hic
#'
#' @description A ggplot2 geom for Hi-C data.
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()].
#' @details
#' Requires the following aesthetics:
#' * seqnames1
#' * start1
#' * end1
#' * seqnames2
#' * start2
#' * end2
#' * fill
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' }
#' @export geom_hic
geom_hic <- function(
  mapping = NULL, data = NULL,
  stat = StatHic, position = "identity",
  na.rm = FALSE, show.legend = NA,
  inherit.aes = TRUE, ...
) {
  ggplot2::layer(
    geom = GeomHic, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
