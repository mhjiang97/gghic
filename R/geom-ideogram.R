ensure_cytoband <- function(genome) {
  path_cytoband <- glue::glue("{dir_cache}/cytoBand.{genome}.rds")
  if (!file.exists(path_cytoband)) {
    session <- rtracklayer::browserSession()
    rtracklayer::genome(session) <- genome
    query <- rtracklayer::ucscTableQuery(session, table = "cytoBandIdeo")
    bands_all <- rtracklayer::getTable(query)
    saveRDS(bands_all, path_cytoband)
  } else {
    bands_all <- readRDS(path_cytoband)
  }
  stop_if_null(bands_all, "cytoBand could not be created or loaded.")
  bands_all
}

retrive_cytoband <- function(data, genome) {
  bands_all <- ensure_cytoband(genome = genome)

  chroms <- unique(c(data$seqnames1, data$seqnames2)) |>
    as.character()
  bands <- bands_all |>
    dplyr::filter(chrom %in% chroms)

  bnames <- bands$name
  indices <- is.na(bnames)
  if (any(indices)) {
    bnames[indices] <- glue::glue("band_na_{seq_len(sum(indices))}")
  }
  if (any(bnames == "")) {
    bnames[bnames == ""] <- glue::glue("band_null_{which(bnames == '')}")
  }

  cols_all <- biovizBase::getBioColor("CYTOBAND")
  cols <- c(
    cols_all[c("gneg", "stalk", "acen")],
    gpos = unname(cols_all["gpos100"]),
    gvar = unname(cols_all["gpos100"])
  )
  gpcols <- unique(grep("gpos", bands$gieStain, value = TRUE))
  crmp <- colorRampPalette(c(cols["gneg"], cols["gpos"]))(100)
  posCols <- setNames(crmp[as.integer(gsub("gpos", "", gpcols))], gpcols)
  cols <- c(cols, posCols) |>
    as.data.frame() |>
    tibble::rownames_to_column("type") |>
    dplyr::rename(band_fill = `c(cols, posCols)`)

  dat <- bands |>
    dplyr::mutate(chromStart = chromStart + 1, name = bnames) |>
    dplyr::rename(
      type = gieStain, seqname = chrom, start = chromStart, end = chromEnd
    ) |>
    dplyr::left_join(cols, by = "type")

  dat
}

StatIdeogram <- ggplot2::ggproto(
  "StatIdeogram",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2"
  ),
  extra_params = c(
    ggplot2::Stat$extra_params,
    "genome", "highlight", "width_ratio", "length_ratio"
  ),
  dropped_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "fill"
  ),
  compute_group = function(
    data, scales,
    genome = "hg19", highlight = TRUE, width_ratio = 1 / 30, length_ratio = 0.8
  ) {
    bands <- retrive_cytoband(data, genome = genome)
    max_y <- max(((data$start1 + data$end2) / 2) - data$start1, na.rm = TRUE)
    max_x <- max((data$end1 + data$end2) / 2, na.rm = TRUE)
    min_x <- min((data$start1 + data$start2) / 2, na.rm = TRUE)
    res <- data$end1[1] - data$start1[1] + 1
    .scale <- ((max_x - min_x) * length_ratio) / max(bands$end)
    .height <- max_y * width_ratio
    # ======================================================================== #
    #   ^                                                                      #
    #   | +--------------------------------------------------------+           #
    #   | | (x, ymax)                                 (xmax, ymax) |           #
    #   | |                                                        |           #
    #   | | (x, y)                                       (xmax, y) |           #
    #   | +--------------------------------------------------------+           #
    #   |                                                                      #
    #   | (min_x, max_y)                                                       #
    #   |                              [HiC plot]                              #
    # --+--------------------------------------------------------------------> #
    #   | (0, 0)                                                               #
    # ======================================================================== #
    ys <- bands |>
      dplyr::distinct(seqname) |>
      dplyr::mutate(
        y = (.height * (dplyr::row_number() - 1)) +
          (dplyr::row_number() * (.height / 2)) +
          max_y +
          res / 2
      )
    dat_band <- bands |>
      dplyr::left_join(ys, by = "seqname") |>
      dplyr::mutate(
        x = start,
        xmax = end,
        ymax = y + .height,
        type = "band",
        x_scale = x * .scale,
        xmax_scale = xmax * .scale,
        x = x_scale + min_x - dplyr::first(x_scale),
        xmax = xmax_scale + min_x - dplyr::first(x_scale)
      ) |>
      dplyr::select(x, y, xmax, ymax, type, seqname, band_fill)
    dat_text <- dat_band |>
      dplyr::group_by(seqname) |>
      dplyr::filter(xmax == max(xmax)) |>
      dplyr::mutate(
        y = (ymax + y) / 2,
        x = xmax + (max(dat_band$xmax) - min(dat_band$x)) / 50,
        type = "text",
        band_fill = "black"
      )
    dat_boundary <- NULL
    if (highlight) {
      dat_boundary <- dat_band |>
        dplyr::distinct(seqname, .keep_all = TRUE) |>
        dplyr::mutate(
          x = min(data$start1) * .scale + min_x,
          xmax = max(data$end1) * .scale + min_x,
          ymax = ymax + .height / 10, y = y - .height / 10,
          type = "boundary",
          band_fill = "red"
        )
    }
    dat <- dplyr::bind_rows(dat_boundary, dat_band, dat_text)
    dat
  }
)

GeomIdeogram <- ggplot2::ggproto(
  "GeomIdeogram",
  ggplot2::Geom,
  required_aes = c("x", "y", "xmax", "ymax", "type", "seqname", "band_fill"),
  extra_params = c(ggplot2::Geom$extra_params, "fontsize"),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(
    data, panel_params, coord,
    fontsize = 10, colour = "red", fill = "#FFE3E680"
  ) {
    coords <- coord$transform(data, panel_params)
    coords_boundary <- coords |>
      dplyr::filter(type == "boundary")
    grob_boundary <- grid::nullGrob()
    if (nrow(coords_boundary) > 0) {
      grob_boundary <- grid::polygonGrob(
        x = c(
          coords_boundary$x, coords_boundary$x,
          coords_boundary$xmax, coords_boundary$xmax
        ),
        y = c(
          coords_boundary$ymax, coords_boundary$y,
          coords_boundary$y, coords_boundary$ymax
        ),
        id = rep(seq_len(nrow(coords_boundary)), 4),
        gp = grid::gpar(col = colour, fill = fill),
        default.units = "native"
      )
    }
    coords_band <- coords |>
      dplyr::filter(type == "band")
    grob_band <- grid::polygonGrob(
      x = c(coords_band$x, coords_band$x, coords_band$xmax, coords_band$xmax),
      y = c(coords_band$ymax, coords_band$y, coords_band$y, coords_band$ymax),
      id = rep(seq_len(nrow(coords_band)), 4),
      gp = grid::gpar(col = "black", fill = coords_band$band_fill),
      default.units = "native"
    )
    coords_text <- coords |>
      dplyr::filter(type == "text")
    grob_text <- grid::textGrob(
      label = coords_text$seqname,
      x = coords_text$x, y = coords_text$y,
      just = c("left", "center"),
      gp = grid::gpar(col = "black", fontsize = fontsize),
      default.units = "native"
    )
    grid::gList(grob_band, grob_boundary, grob_text)
  }
)

#' geom_ideogram
#'
#' @description A ggplot2 geom for chromosome ideogram.
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()].
#' @param ... Parameters including:
#'  * `genome`: The genome name. Default is `"hg19"`.
#'  * `highlight`: Whether to highlight the boundary of the chromosome. Default is `TRUE`.
#'  * `fontsize`: The font size of the chromosome name. Default is `10`.
#'  * `width_ratio`: The ratio of the width of each chromosome ideogram relative to the height of the Hi-C plot. Default is `1/30`.
#'  * `length_ratio`: The ratio of the length of each chromosome ideogram relative to the width of the Hi-C plot. Default is `0.8`.
#' @details
#' Requires the following aesthetics:
#' * seqnames1
#' * start1
#' * end1
#' * seqnames2
#' * start2
#' * end2
#' @return A ggplot object.
#' @examples
#' \dontrun{
#' }
#' @export geom_ideogram
geom_ideogram <- function(
  mapping = NULL, data = NULL, stat = StatIdeogram, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, check.param = FALSE, ...
) {
  ggplot2::layer(
    geom = GeomIdeogram, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE, params = list(na.rm = na.rm, ...)
  )
}

