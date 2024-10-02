extract_trs <- function(grs, genes, txdb) {
  chrom <- as.character(grs@seqnames)
  start <- grs@ranges@start
  end <- start + grs@ranges@width - 1

  r <- range(
    genes[
      S4Vectors::subjectHits(
        GenomicRanges::findOverlaps(grs, genes, ignore.strand = TRUE)
      )
    ]
  )

  txdump <- AnnotationDbi::as.list(txdb)
  txdump$chrominfo <- txdump$chrominfo[
    txdump$chrominfo$chrom == chrom, ,
    drop = FALSE
  ]
  txdump$transcripts <- txdump$transcripts[
    txdump$transcripts$tx_chrom == chrom &
      txdump$transcripts$tx_start < end(r) &
      txdump$transcripts$tx_end > start(r), ,
    drop = FALSE
  ]
  txdump$genes <- txdump$genes[
    txdump$genes$tx_id %in% txdump$transcripts$tx_id, ,
    drop = FALSE
  ]
  txdump$splicings <- txdump$splicings[
    txdump$splicings$tx_id %in% txdump$transcripts$tx_id, ,
    drop = FALSE
  ]
  txdb_subset <- do.call(txdbmaker::makeTxDb, txdump)

  exons_gviz <- Gviz::GeneRegionTrack(
    txdb_subset,
    chromosome = chrom, start = start(r), end = end(r), strand = "*"
  )
  exons <- exons_gviz@range

  trs <- GenomicRanges::split(exons, as.character(exons$transcript))

  trs[IRanges::overlapsAny(trs, grs, type = "any")]
}

calculate_lines <- function(genes_span, maxgap) {
  names_gene <- names_gene_copy <- names(genes_span)
  dat_line <- tibble::tibble(
    gene_id = names_gene,
    line = 0
  )

  while (length(names_gene_copy) > 1) {
    n1 <- names_gene_copy[1]
    n1o <- setdiff(names_gene_copy, n1)
    genes_n1 <- genes_span[[n1]]
    genes_n1o <- GenomicRanges::GRangesList(genes_span[n1o]) |>
      unlist()
    genes_n1o_nonol <- genes_n1o[
      !IRanges::overlapsAny(genes_n1o, genes_n1, maxgap = maxgap, ignore.strand = TRUE)
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

retrive_genes <- function(data, txdb, tx2gene, gtf_file, maxgap) {
  if (is.null(txdb) && is.null(gtf_file)) {
    stop("txdb or gtf_file must be provided.")
  }
  if (is.null(tx2gene) && is.null(gtf_file)) {
    stop("tx2gene or gtf_file must be provided.")
  }

  if (is.null(txdb) && !is.null(gtf_file)) {
    if (!file.exists(glue::glue("{dir_cache}/{basename(gtf_file)}.sqlite"))) {
      txdb <- txdbmaker::makeTxDbFromGFF(gtf, format = "gtf")
      AnnotationDbi::saveDb(
        txdb, glue::glue("{dir_cache}/{basename(gtf_file)}.sqlite")
      )
    } else {
      txdb <- suppressMessages(
        AnnotationDbi::loadDb(glue::glue("{dir_cache}/{basename(gtf_file)}.sqlite"))
      )
    }
  }

  if (is.null(tx2gene) && !is.null(gtf_file)) {
    if (!file.exists(glue::glue("{dir_cache}/{basename(gtf_file)}.tx2gene.tsv"))) {
      module_py <- reticulate::import_from_path(
        module = "tx2gene", path = system.file("python", package = "gghic")
      )
      module_py$generate_tx2gene(
        gtf_file, glue::glue("{dir_cache}/{basename(gtf_file)}.tx2gene.tsv")
      )
    }
    tx2gene <- vroom::vroom(
      glue::glue("{dir_cache}/{basename(gtf_file)}.tx2gene.tsv"),
      col_names = c(
        "chrom", "gene_id", "gene_symbol", "tx_id",
        "tx_name", "gene_type", "tx_type"
      ),
      col_types = vroom::cols()
    )
  }

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

  trs <- purrr::map2(grs, genes, extract_trs, txdb = txdb)

  dat <- purrr::pmap(
    list(trs), function(.x) {
      genes <- GenomicRanges::split(unlist(.x), unlist(.x)$gene)

      genes_reduced <- purrr::map(
        genes, function(..x) {
          features <- GenomicRanges::split(..x, ..x$feature)
          features_reduced <- purrr::map(
            features, function(...x) {
              tmp <- GenomicRanges::reduce(...x)
              tmp$feature <- unique(...x$feature)
              tmp$gene_id <- unique(..x$gene)
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
          introns <- GenomicRanges::gaps(tmp)[-1]
          if (length(introns) > 0) {
            introns$gene_id <- unique(tmp$gene_id)
            introns$feature <- "intron"
          }
          as.data.frame(c(setNames(tmp, NULL), introns))
        }
      ) |>
        dplyr::left_join(dat_line, by = "gene_id") |>
        dplyr::left_join(
          dplyr::select(
            dplyr::distinct(tx2gene, gene_id, .keep_all = TRUE), gene_id, gene_symbol
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
    "seqnames1", "start1", "end1",
    "seqnames2", "start2", "end2"
  ),
  extra_params = c(
    ggplot2::Stat$extra_params, "txdb", "tx2gene", "gtf_file", "ratio_height",
    "fontsize", "maxgap"
  ),
  dropped_aes = c(
    "seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "fill"
  ),
  compute_group = function(
    data, scales, txdb = NULL, tx2gene = NULL, gtf_file = NULL,
    ratio_height = 1 / 50, fontsize = 10, maxgap = -1
  ) {
    genes <- retrive_genes(
      data,
      txdb = txdb, tx2gene = tx2gene, gtf_file = gtf_file, maxgap = maxgap
    )
    max_y <- max(
      ((data$start1 + data$end2) / 2) - data$start1, na.rm = TRUE
    )
    max_x <- max((data$end1 + data$end2) / 2, na.rm = TRUE)
    min_x <- min((data$start1 + data$start2) / 2, na.rm = TRUE)
    .height <- max_y * ratio_height
    .height_cds <- .height * 0.8
    # ======================================================================== #
    #   ^                                                                      #
    #   | (min_x, max_y)                                                       #
    #   |                              [HiC plot]                              #
    # --+--------------------------------------------------------------------> #
    #   | (0, 0)                                                               #
    #   |                   +------------------------------------+             #
    #   | +-----------+     | (x, y)                   (xmax, y) |             #
    #   | |           |-----|                                    |             #
    #   | +-----------+     | (x, ymin)             (xmax, ymin) |             #
    #   |                   +------------------------------------+             #
    #   |                       [gene name]                                    #
    # ======================================================================== #
    ys <- genes |>
      dplyr::distinct(line) |>
      dplyr::mutate(
        y = ((-1 * .height) * (dplyr::row_number() - 1)) -
          (dplyr::row_number() * (.height / 3)) -
          .height / 2
      )
    dat_cds <- genes |>
      dplyr::filter(feature == "CDS") |>
      dplyr::left_join(ys, by = "line") |>
      dplyr::mutate(
        x = start,
        xmax = end,
        ymin = y - .height_cds
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
        ymin = y - .height_cds - (.height_cds / 20),
        feature = "text",
        xmax = end
      )
    dat <- dplyr::bind_rows(dat_cds, dat_others, dat_text) |>
      dplyr::mutate(fontsize = fontsize)
    # dat
    dat[dat$gene_symbol == dat$gene_symbol[1:2] & dat$feature != "ncRNA", ]
  }
)

GeomAnnotation <- ggplot2::ggproto(
  "GeomAnnotation",
  ggplot2::Geom,
  required_aes = c("x", "y", "xmax", "ymin", "feature", "gene_symbol", "fontsize"),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(data, panel_params, coord) {
    coords <- coord$transform(data, panel_params)
    coords_exon <- coords |>
      dplyr::filter(feature %in% c("CDS", "utr3", "utr5", "ncRNA"))
    grob_exon <- grid::polygonGrob(
      x = c(coords_exon$x, coords_exon$x, coords_exon$xmax, coords_exon$xmax),
      y = c(coords_exon$y, coords_exon$ymin, coords_exon$ymin, coords_exon$y),
      id = rep(seq_len(nrow(coords_exon)), 4),
      gp = grid::gpar(col = "blue" ,fill = "blue"),
      default.units = "native"
    )
    coords_text <- coords |>
      dplyr::filter(feature == "text")
    grob_text <- grid::textGrob(
      label = coords_text$gene_symbol,
      x = coords_text$x, y = coords_text$ymin,
      just = c("center", "top"),
      gp = grid::gpar(col = "black", fontsize = coords_text$fontsize),
      default.units = "native"
    )
    grid::gList(grob_exon, grob_text)
  }
)

geom_annotation <- function(
    mapping = NULL, data = NULL,
    stat = StatAnnotation, position = "identity",
    na.rm = FALSE, show.legend = NA,
    inherit.aes = TRUE, ...
) {
  ggplot2::layer(
    geom = GeomAnnotation, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
