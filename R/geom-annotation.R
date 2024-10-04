extract_trs <- function(grs, genes, txdump, include_ncrna) {
  chrom <- as.character(grs@seqnames)
  start <- grs@ranges@start
  end <- start + grs@ranges@width - 1

  r <- range(
    genes[
      S4Vectors::subjectHits(
        GenomicRanges::findOverlaps(grs, genes, ignore.strand = TRUE)
      )
    ],
    ignore.strand = TRUE
  )

  .txdump <- txdump
  .txdump$chrominfo <- .txdump$chrominfo[
    .txdump$chrominfo$chrom == chrom, ,
    drop = FALSE
  ]
  .txdump$transcripts <- .txdump$transcripts[
    .txdump$transcripts$tx_chrom == chrom &
      .txdump$transcripts$tx_start < end(r) &
      .txdump$transcripts$tx_end > start(r), ,
    drop = FALSE
  ]
  .txdump$genes <- .txdump$genes[
    .txdump$genes$tx_id %in% .txdump$transcripts$tx_id, ,
    drop = FALSE
  ]
  .txdump$splicings <- .txdump$splicings[
    .txdump$splicings$tx_id %in% .txdump$transcripts$tx_id, ,
    drop = FALSE
  ]
  txdb_subset <- do.call(txdbmaker::makeTxDb, .txdump)

  exons_gviz <- Gviz::GeneRegionTrack(
    txdb_subset,
    chromosome = chrom, start = start(r), end = end(r), strand = "*"
  )
  exons <- exons_gviz@range

  if (!include_ncrna) {
    ncrnas <- unique(exons$transcript[exons$feature == "ncRNA"])
    exons <- exons[!exons$transcript %in% ncrnas]
  }

  trs <- GenomicRanges::split(exons, as.character(exons$transcript))

  trs_exon_intron <- purrr::map(trs, function(tr) {
    introns <- GenomicRanges::gaps(tr)[-1]
    if (length(introns) > 0) {
      introns$gene <- unique(tr$gene)
      introns$feature <- "intron"
    }
    tmp <- GenomicRanges::pintersect(c(tr, introns), grs, ignore.strand = TRUE)
    tmp[width(tmp) > 0]
  })

  GenomicRanges::GRangesList(trs_exon_intron)
}

calculate_lines <- function(genes_span, maxgap) {
  names_gene <- names_gene_copy <- names(genes_span)
  dat_line <- tibble::tibble(gene_id = names_gene, line = 0)

  while (length(names_gene_copy) > 1) {
    n1 <- names_gene_copy[1]
    n1o <- setdiff(names_gene_copy, n1)
    genes_n1 <- genes_span[[n1]]
    genes_n1o <- GenomicRanges::GRangesList(genes_span[n1o]) |>
      unlist()
    genes_n1o_nonol <- genes_n1o[
      !IRanges::overlapsAny(
        genes_n1o, genes_n1, maxgap = maxgap, ignore.strand = TRUE
      )
    ]
    n1o_nonol <- names(genes_n1o_nonol)
    if (length(n1o_nonol) == 0) {
      dat_line$line[dat_line$gene_id == n1] <- max(dat_line$line) + 1
      names_gene_copy <- setdiff(names_gene_copy, n1)
      next
    }
    dat_line$line[dat_line$gene_id %in% c(n1, n1o_nonol)] <- max(dat_line$line) + 1
    names_gene_copy <- setdiff(names_gene_copy, n1o_nonol)
  }
  if (sum(dat_line$line == 0)) {
    dat_line$line[dat_line$line == 0] <- max(dat_line$line) + 1
  }

  dat_line
}

ensure_txdb <- function(txdb, gtf_path) {
  if (is.null(txdb) && !is.null(gtf_path)) {
    path_txdb <- glue::glue("{dir_cache}/{basename(gtf_path)}.sqlite")
    if (!file.exists(path_txdb)) {
      txdb <- txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf")
      AnnotationDbi::saveDb(txdb, path_txdb)
    } else {
      txdb <- suppressMessages(AnnotationDbi::loadDb(path_txdb))
    }
  }
  stop_if_null(txdb, "TxDb could not be created or loaded.")
  txdb
}

ensure_txdump <- function(txdb, gtf_path) {
  path_txdump <- glue::glue("{dir_cache}/{basename(gtf_path)}.txdump.rds")
  if (!file.exists(path_txdump)) {
    txdump <- AnnotationDbi::as.list(txdb)
    saveRDS(txdump, path_txdump)
  } else {
    txdump <- readRDS(path_txdump)
  }
  txdump
}

ensure_tx2gene <- function(tx2gene, gtf_path) {
  if (is.null(tx2gene) && !is.null(gtf_path)) {
    path_tx2gene <- glue::glue("{dir_cache}/{basename(gtf_path)}.tx2gene.rds")
    tsv <- glue::glue("{dir_cache}/{basename(gtf_path)}.tx2gene.tsv")
    if (!file.exists(path_tx2gene)) {
      module_py <- reticulate::import_from_path(
        module = "tx2gene", path = system.file("python", package = name_pkg)
      )
      module_py$generate_tx2gene(gtf_path, tsv)
      vroom::vroom(
        tsv, col_names = c(
          "chrom", "gene_id", "gene_symbol", "tx_id", "tx_name",
          "gene_type", "tx_type"
        ),
        col_types = vroom::cols()
      ) |>
        saveRDS(path_tx2gene)
      file.remove(tsv)
    } else {
      tx2gene <- readRDS(path_tx2gene)
    }
  }
  stop_if_null(tx2gene, "tx2gene data could not be created or loaded.")
  tx2gene
}

retrive_genes <- function(
  data, txdb, tx2gene, gtf_path, maxgap, include_ncrna
) {
  txdb <- ensure_txdb(txdb, gtf_path)
  tx2gene <- ensure_tx2gene(tx2gene, gtf_path)
  txdump <- ensure_txdump(txdb, gtf_path)

  grs <- purrr::map(
    unique(data$seqnames1), function(.x) {
      GenomicRanges::GRanges(
        seqnames = .x,
        ranges = IRanges::IRanges(
          start = min(data$start1[data$seqnames1 == .x]),
          end = max(data$end1[data$seqnames1 == .x])
        )
      )
    }
  )

  genes <- purrr::map(
    grs, function(.x) {
      GenomicFeatures::genes(
        txdb, columns = c("exon_id"), filter = list(tx_chrom = .x@seqnames)
      )
    }
  )

  trs <- purrr::map2(
    grs, genes, extract_trs,
    txdump = txdump, include_ncrna = include_ncrna
  )

  dat <- purrr::pmap(
    list(trs), function(.x) {
      genes <- GenomicRanges::split(unlist(.x), unlist(.x)$gene)

      genes_reduced <- purrr::map(
        genes, function(g) {
          features <- GenomicRanges::split(g, g$feature)
          features_reduced <- purrr::map(
            features, function(f) {
              tmp <- GenomicRanges::reduce(f)
              tmp$feature <- unique(f$feature)
              tmp$gene_id <- unique(g$gene)
              tmp
            }
          )
        }
      )

      genes_span <- purrr::map(
        genes_reduced, function(.x) {
          tmp <- unlist(GenomicRanges::GRangesList(.x))
          GenomicRanges::GRanges(
            seqnames = tmp@seqnames[1],
            ranges = IRanges::IRanges(
              start = min(start(tmp)),
              end = max(end(tmp))
            ),
            strand = tmp@strand[1]
          )
        }
      )

      dat_line <- calculate_lines(genes_span, maxgap)

      dat_gene <- purrr::map_df(
        genes_reduced, function(.x) {
          tmp <- GenomicRanges::GRangesList(.x) |>
            unlist()
          tibble::as_tibble(setNames(tmp, NULL))
        }
      ) |>
        dplyr::left_join(dat_line, by = "gene_id") |>
        dplyr::left_join(
          dplyr::select(
            dplyr::distinct(tx2gene, gene_id, .keep_all = TRUE),
            gene_id, gene_symbol
          ), by = "gene_id"
        )

      dat_gene
    }
  )

  dplyr::bind_rows(dat)
}

StatAnnotation <- ggplot2::ggproto(
  "StatAnnotation",
  ggplot2::Stat,
  required_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2"
  ),
  extra_params = c(
    ggplot2::Stat$extra_params, "txdb", "tx2gene", "gtf_path", "width_ratio",
    "maxgap", "include_ncrna", "style"
  ),
  dropped_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "fill"
  ),
  compute_group = function(
    data, scales,
    txdb = NULL, tx2gene = NULL, gtf_path = NULL, width_ratio = 1 / 50,
    maxgap = -1, include_ncrna = TRUE, style = c("basic", "arrow")
  ) {
    genes <- retrive_genes(data, txdb, tx2gene, gtf_path, maxgap, include_ncrna)
    max_y <- max(((data$start1 + data$end2) / 2) - data$start1, na.rm = TRUE)
    max_x <- max((data$end1 + data$end2) / 2, na.rm = TRUE)
    min_x <- min((data$start1 + data$start2) / 2, na.rm = TRUE)
    res <- data$end1[1] - data$start1[1] + 1
    .height <- max_y * width_ratio
    .height_cds <- .height * 0.8
    ys <- genes |>
      dplyr::distinct(line) |>
      dplyr::mutate(
        y = ((-1 * .height) * (dplyr::row_number() - 1)) -
          (dplyr::row_number() * (.height / 3)) -
          res / 2
      )
    dat_text <- genes |>
      dplyr::left_join(ys, by = "line") |>
      dplyr::group_by(gene_symbol) |>
      dplyr::summarise(
        gene_symbol = dplyr::first(gene_symbol),
        gene_id = dplyr::first(gene_id),
        line = dplyr::first(line),
        seqnames = dplyr::first(seqnames),
        strand = dplyr::first(strand),
        start = min(start),
        end = max(end),
        width = end - start,
        x = mean(c(start, end)),
        y = dplyr::first(y),
        ymin = y - .height_cds - (.height_cds / 5),
        feature = "text",
        xmax = end
      )
    # ======================================================================== #
    #   ^                                                                      #
    #   | (min_x, max_y)                                                       #
    #   |                              [HiC plot]                              #
    # --+--------------------------------------------------------------------> #
    #   | (0, 0)                                                               #
    #   |                      +------------------------------------+          #
    #   | +-----------+        | (x, y)                   (xmax, y) |          #
    #   | |           |--<--<--|                                    |          #
    #   | +-----------+        | (x, ymin)             (xmax, ymin) |          #
    #   |                      +------------------------------------+          #
    #   |                          [gene name]                                 #
    #   |                                                                      #
    #   |    , -----------------------------------------------------+          #
    #   |  ,   (x, y)                                     (xmax, y) |          #
    #   | + (xend, yend)                                            |          #
    #   |  `   (x, ymin)                               (xmax, ymin) |          #
    #   |    ` -----------------------------------------------------+          #
    #   |                          [gene name]                                 #
    #   | +--------------------------------------------------- `               #
    #   | | (x, y)                                      (xmax, y) `            #
    #   | |                                            (xend, yend) +          #
    #   | | (x, ymin)                                (xmax, ymin) ,            #
    #   | +---------------------------------------------------- ,              #
    #   |                          [gene name]                                 #
    # ======================================================================== #
    if (style[1] == "basic") {
      dat_cds <- genes |>
        dplyr::filter(feature == "CDS") |>
        dplyr::left_join(ys, by = "line") |>
        dplyr::mutate(
          x = start,
          xmax = end,
          ymin = y - .height_cds
        )
      dat_intron <- genes |>
        dplyr::filter(feature == "intron") |>
        dplyr::left_join(ys, by = "line") |>
        dplyr::mutate(
          x = start,
          xmax = end,
          ymin = y - .height_cds,
          y = 0.5 * (y + ymin)
        )
      dat_others <- genes |>
        dplyr::filter(feature %in% c("utr3", "utr5", "ncRNA")) |>
        dplyr::left_join(ys, by = "line") |>
        dplyr::mutate(
          x = start,
          xmax = end,
          ymin = y - .height_cds,
          y = y - (.height_cds / 8),
          ymin = ymin + (.height_cds / 8)
        )
      dat <- dplyr::bind_rows(dat_cds, dat_intron, dat_others, dat_text) |>
        dplyr::mutate(
          xmin = min_x,
          xend = max_x
        )
    }
    if (style[1] == "arrow") {
      dat <- genes |>
        dplyr::left_join(ys, by = "line") |>
        dplyr::group_by(gene_symbol) |>
        dplyr::summarise(
          gene_symbol = dplyr::first(gene_symbol),
          gene_id = dplyr::first(gene_id),
          line = dplyr::first(line),
          seqnames = dplyr::first(seqnames),
          strand = dplyr::first(strand),
          y = dplyr::first(y),
          xend_p = max(end),
          xend_n = min(start),
          xmax_n = max(end),
          x_p = min(start),
          x_n = (xmax_n - xend_n) * 0.1 + xend_n,
          xmax_p = xend_p - ((xend_p - x_p) * 0.1),
          xend = dplyr::case_when(
            strand == "+" ~ xend_p,
            strand == "-" ~ xend_n
          ),
          xmax = dplyr::case_when(
            strand == "+" ~ xmax_p,
            strand == "-" ~ xmax_n
          ),
          x = dplyr::case_when(
            strand == "+" ~ x_p,
            strand == "-" ~ x_n
          ),
          ymin = y - .height_cds,
          yend = (y + ymin) / 2,
          feature = "arrow"
        ) |>
        dplyr::bind_rows(dat_text)
    }

    dat
  }
)

GeomAnnotation <- ggplot2::ggproto(
  "GeomAnnotation",
  ggplot2::Geom,
  required_aes = c(
    "x", "y", "xmax", "ymin", "feature", "gene_symbol", "strand"
  ),
  extra_params = c(ggplot2::Geom$extra_params, "fontsize", "style"),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(
    data, panel_params, coord,
    fontsize = 10, style = c("basic", "arrow"),
    colour = "blue", fill = "blue"
  ) {
    coords <- coord$transform(data, panel_params)
    coords_text <- coords |>
      dplyr::filter(feature == "text")
    grob_text <- grid::textGrob(
      label = coords_text$gene_symbol,
      x = coords_text$x, y = coords_text$ymin,
      just = c("center", "top"),
      gp = grid::gpar(col = "black", fontsize = fontsize),
      default.units = "native"
    )
    if (style[1] == "basic") {
      coords_exon <- coords |>
        dplyr::filter(feature %in% c("CDS", "utr3", "utr5", "ncRNA"))
      grob_exon <- grid::polygonGrob(
        x = c(coords_exon$x, coords_exon$x, coords_exon$xmax, coords_exon$xmax),
        y = c(coords_exon$y, coords_exon$ymin, coords_exon$ymin, coords_exon$y),
        id = rep(seq_len(nrow(coords_exon)), 4),
        gp = grid::gpar(col = colour, fill = fill),
        default.units = "native"
      )
      coords_intron <- coords |>
        dplyr::filter(feature == "intron")
      ends <- ifelse(coords_intron$strand == "+", "last", "first")
      lengths_arrow <- rep(
        (coords_intron$xend[1] - coords_intron$xmin[1]) / 80,
        nrow(coords_intron)
      )
      lengths_intron <- coords_intron$xmax - coords_intron$x
      lengths_arrow[lengths_intron < lengths_arrow] <- 0
      grob_intron <- grid::segmentsGrob(
        x0 = coords_intron$x, x1 = coords_intron$xmax,
        y0 = coords_intron$y, y1 = coords_intron$y,
        arrow = grid::arrow(
          type = "open",
          length = grid::unit(lengths_arrow, "native"),
          ends = ends
        ),
        gp = grid::gpar(col = colour),
        default.units = "native"
      )
      grids <- grid::gList(grob_exon, grob_intron, grob_text)
    }
    if (style[1] == "arrow") {
      coords_arrow_p <- coords |>
        dplyr::filter(feature == "arrow", strand == "+")
      grob_arrow_p <- grid::polygonGrob(
        x = c(
          coords_arrow_p$x, coords_arrow_p$x, coords_arrow_p$xmax,
          coords_arrow_p$xend, coords_arrow_p$xmax
        ),
        y = c(
          coords_arrow_p$y, coords_arrow_p$ymin, coords_arrow_p$ymin,
          coords_arrow_p$yend, coords_arrow_p$y
        ),
        id = rep(seq_len(nrow(coords_arrow_p)), 5),
        gp = grid::gpar(col = colour, fill = fill),
        default.units = "native"
      )
      coords_arrow_n <- coords |>
        dplyr::filter(feature == "arrow", strand == "-")
      grob_arrow_n <- grid::polygonGrob(
        x = c(
          coords_arrow_n$x, coords_arrow_n$xend, coords_arrow_n$x,
          coords_arrow_n$xmax, coords_arrow_n$xmax
        ),
        y = c(
          coords_arrow_n$y, coords_arrow_n$yend, coords_arrow_n$ymin,
          coords_arrow_n$ymin, coords_arrow_n$y
        ),
        id = rep(seq_len(nrow(coords_arrow_n)), 5),
        gp = grid::gpar(col = colour, fill = fill),
        default.units = "native"
      )
      grids <- grid::gList(grob_arrow_p, grob_arrow_n, grob_text)
    }
    grids
  }
)

#' geom_annotation
#'
#' @description A ggplot2 geom for gene model tracks.
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()].
#' @param ... Parameters including:
#'  * `txdb`: The TxDb object. Default is `NULL`.
#'  * tx2gene: An optional data frame or tibble that maps transcript information to gene information. It should include the following columns:
#'      * chrom: Chromosome number or name.
#'      * gene_id: Entrez gene ID.
#'      * gene_symbol: Common symbol or name of the gene.
#'      * tx_id: Entrez transcript ID.
#'      * tx_name: Name of the transcript.
#'      * gene_type: Type or classification of the gene.
#'      * tx_type: Type or classification of the transcript.
#'  * `gtf_path`: The path to the GTF file, which is used to generate `txdb` and `tx2gene`. Generated files are saved in the cache directory. Default is `NULL`.
#'  * `width_ratio`: The ratio of the width of each gene model track relative to the height of the Hi-C plot. Default is `0.02`.
#'  * `maxgap`: The maximum gap between genes to be drawn in the same line. Default is `-1`.
#'  * `include_ncrna`: Whether to include ncRNA or not. Default is `TRUE`.
#'  * `style`: The style of the gene model track, which can be `"basic"` or `"arrow"`. Default is `"basic"`.
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
#' @export geom_annotation
geom_annotation <- function(
  mapping = NULL, data = NULL, stat = StatAnnotation, position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, check.param = FALSE, ...
) {
  ggplot2::layer(
    geom = GeomAnnotation, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    check.param = FALSE, params = list(na.rm = na.rm, ...)
  )
}
