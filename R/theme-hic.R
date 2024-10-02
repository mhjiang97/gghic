#' @importFrom ggplot2 scale_fill_gradientn
colors_hic <- function() {
  ggplot2::scale_fill_gradientn(
    colors = c(
      "#FFFEF9", "#FCF9CE", "#FFF2A9", "#FDE188", "#FFCA67", "#FAAA4B",
      "#F78E40", "#F15C34", "#ED3024", "#D42027", "#B01F29", "#7A1128",
      "#1A0A10"
    ),
    na.value = "#FFFFFF"
  )
}

#' @importFrom ggplot2 element_blank element_line theme theme_bw %+replace%
#' @importFrom scales unit_format
theme_hic <- function(hide_y = TRUE) {
  `%+replace%` <- ggplot2::`%+replace%`

  t <- ggplot2::theme_bw() %+replace%
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank()
    )

  if (hide_y) {
    t <- t %+replace%
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.line.x.bottom = ggplot2::element_line(color = "black")
      )
  }

  list(
    t,
    colors_hic(),
    ggplot2::coord_fixed(),
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      labels = scales::unit_format(unit = "M", scale = 1e-6)
    )
  )
}
