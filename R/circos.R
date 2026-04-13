#' Visualize Host-Virus Integration Events
#'
#' Draw a host-virus circos plot on the current graphics device from an
#' integration table and a chromosome-size table. No intermediate VCF is
#' created and no output file is opened by this function.
#'
#' @param input_file Path to a tab-delimited integration table with columns
#'   \code{chr}, \code{host_loc}, \code{viral_loc}, \code{reads},
#'   \code{sample}, \code{viral_strand}, and \code{method}.
#' @param chrom_file Path to a tab-delimited chromosome-size table with columns
#'   \code{chr}, \code{start}, and \code{end}. It must include the host
#'   chromosomes and the virus sequence.
#' @param layout_list A list of track definitions. Supported \code{type} values
#'   are \code{"ideogram"}, \code{"scatter"}, and \code{"links"}.
#' @param color_file Optional path to a two-column color table.
#' @param visual_ratio Visual proportion assigned to the virus sector.
#' @param clear Logical. Whether to clear the existing circlize plot before
#'   drawing.
#'
#' @return Invisibly returns a list with \code{cfg}, \code{gi}, and
#'   \code{data}.
#' @export
visualize_viral_integration <- function(input_file,
                                        chrom_file,
                                        layout_list = create_host_virus_layout(...),
                                        color_file = NULL,
                                        visual_ratio = 0.1,
                                        clear = TRUE) {
  validate_layout_list(layout_list)

  cfg <- create_config(
    chrom_file = chrom_file,
    visual_ratio = visual_ratio,
    color_file = color_file
  )
  gi <- create_gi_from_table(input_file = input_file, cfg = cfg)
  plot_df <- as.data.frame(gi)

  if (clear) {
    circlize::circos.clear()
  }

  n_sectors <- nrow(cfg$data)
  gaps <- rep(1, n_sectors)
  virus_idx <- which(cfg$data$chr == cfg$virus_name)
  if (length(virus_idx) > 0) {
    gaps[virus_idx] <- 15
  }

  circlize::circos.par(
    start.degree = 90,
    gap.degree = gaps,
    cell.padding = c(0, 0, 0, 0),
    points.overflow.warning = FALSE
  )

  circos_df <- cfg$data[, c("chr", "start", "end"), drop = FALSE]
  circlize::circos.genomicInitialize(
    circos_df,
    plotType = NULL,
    sector.width = cfg$widths
  )

  for (i in seq_along(layout_list)) {
    task <- layout_list[[i]]
    height <- if (is.null(task$height)) 0.05 else task$height

    if (task$type == "ideogram") {
      draw_ideogram(height = height, cfg = cfg)
    } else if (task$type == "scatter") {
      if (!is.null(task$sample_label)) {
        sub_df <- plot_df[plot_df$Label %in% task$sample_label, , drop = FALSE]
      } else {
        sub_df <- plot_df
      }
      draw_scatter(data = sub_df, height = height, cfg = cfg)
    } else if (task$type == "links") {
      draw_link(link_data = plot_df, cfg = cfg, radius = task$radius)
    } else {
      warning("Skipping unsupported track type: ", task$type)
    }
  }

  invisible(list(cfg = cfg, gi = gi, data = plot_df))
}

#' Create a Default Host-Virus Layout
#'
#' @param include_ideogram Logical. Whether to include the ideogram track.
#' @param include_links Logical. Whether to include the link layer.
#' @param sample_label Optional sample labels for the scatter track.
#' @param ideogram_height Numeric. Ideogram track height.
#' @param scatter_height Numeric. Scatter track height.
#' @return A layout list for \code{visualize_viral_integration()}.
#' @export
create_host_virus_layout <- function(include_ideogram = TRUE,
                                     include_links = TRUE,
                                     sample_label = NULL,
                                     ideogram_height = 0.08,
                                     scatter_height = 0.15) {
  layout <- list()
  if (include_ideogram) {
    layout[[length(layout) + 1]] <- list(type = "ideogram", height = ideogram_height)
  }
  layout[[length(layout) + 1]] <- list(
    type = "scatter",
    sample_label = sample_label,
    height = scatter_height
  )
  if (include_links) {
    layout[[length(layout) + 1]] <- list(type = "links")
  }
  layout
}

#' Read Chromosome Sizes
#'
#' @param path Path to a tab-delimited chromosome-size table.
#' @return A data frame with columns \code{chr}, \code{start}, and \code{end}.
#' @keywords internal
read_chrom_sizes <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  df <- utils::read.table(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c("chr", "start", "end")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("chrom_sizes.txt must contain columns: chr, start, end")
  }

  df$chr <- as.character(df$chr)
  df$start <- as.numeric(df$start)
  df$end <- as.numeric(df$end)
  df <- df[!is.na(df$chr) & !is.na(df$start) & !is.na(df$end), , drop = FALSE]

  if (nrow(df) == 0) {
    stop("No valid chromosome sizes were found.")
  }

  df
}

#' Detect Virus Sequence
#'
#' @param chrom_sizes_df Data frame returned by \code{read_chrom_sizes()}.
#' @return The inferred virus sequence name.
#' @keywords internal
get_virus_name <- function(chrom_sizes_df) {
  is_mito <- grepl("^(chr)?(m|mt)$", chrom_sizes_df$chr, ignore.case = TRUE)
  candidates <- chrom_sizes_df[!is_mito, , drop = FALSE]
  candidates <- candidates[order(candidates$end - candidates$start), , drop = FALSE]

  if (nrow(candidates) == 0) {
    stop("No valid chromosome or virus sequence was found.")
  }

  candidates$chr[1]
}

#' Create Plotting Configuration
#'
#' @param chrom_file Path to a chromosome-size table.
#' @param visual_ratio Visual proportion assigned to the virus sector.
#' @param color_file Optional path to a two-column color table.
#' @return A plotting configuration list.
#' @keywords internal
create_config <- function(chrom_file,
                          visual_ratio = 0.1,
                          color_file = NULL) {
  chrom_df <- read_chrom_sizes(chrom_file)
  virus_name <- get_virus_name(chrom_df)

  if (!is.numeric(visual_ratio) || length(visual_ratio) != 1 ||
      is.na(visual_ratio) || visual_ratio <= 0 || visual_ratio >= 1) {
    stop("visual_ratio must be a single number between 0 and 1.")
  }

  widths <- numeric(nrow(chrom_df))
  names(widths) <- chrom_df$chr
  widths[virus_name] <- visual_ratio

  host_mask <- chrom_df$chr != virus_name
  host_len <- chrom_df$end[host_mask] - chrom_df$start[host_mask]
  widths[host_mask] <- host_len / sum(host_len) * (1 - visual_ratio)

  if (!is.null(color_file) && file.exists(color_file)) {
    color_df <- utils::read.table(
      color_file,
      header = FALSE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
    grid_col <- stats::setNames(color_df[[2]], color_df[[1]])
  } else {
    grid_col <- stats::setNames(grDevices::rainbow(nrow(chrom_df)), chrom_df$chr)
  }

  if (!virus_name %in% names(grid_col)) {
    grid_col[virus_name] <- "grey90"
  }

  list(
    data = chrom_df,
    widths = widths,
    virus_name = virus_name,
    grid_col = grid_col
  )
}

#' Convert an Integration Table to GInteractions
#'
#' @param input_file Path to an integration table.
#' @param cfg Plotting configuration returned by \code{create_config()}.
#' @return A \code{GInteractions} object.
#' @export
create_gi_from_table <- function(input_file, cfg) {
  if (!file.exists(input_file)) {
    stop("File not found: ", input_file)
  }

  raw <- utils::read.table(
    input_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c("chr", "host_loc", "viral_loc", "reads", "sample", "viral_strand", "method")
  missing_cols <- setdiff(required_cols, colnames(raw))
  if (length(missing_cols) > 0) {
    stop("data.txt must contain columns: ", paste(required_cols, collapse = ", "))
  }

  host_chr <- match_chr_style(raw$chr, cfg$data$chr)
  virus_chr <- rep(cfg$virus_name, nrow(raw))
  host_pos <- as.numeric(raw$host_loc)
  virus_pos <- as.numeric(raw$viral_loc)

  keep <- !is.na(host_chr) & !is.na(host_pos) & !is.na(virus_pos)
  if (!all(keep)) {
    warning(sum(!keep), " rows were skipped because of missing chromosome or position values.")
  }

  raw <- raw[keep, , drop = FALSE]
  host_chr <- host_chr[keep]
  virus_chr <- virus_chr[keep]
  host_pos <- host_pos[keep]
  virus_pos <- virus_pos[keep]

  if (nrow(raw) == 0) {
    stop("No valid integration records were found in data.txt.")
  }

  gr_host <- GenomicRanges::GRanges(
    seqnames = host_chr,
    ranges = IRanges::IRanges(host_pos, width = 1)
  )
  gr_virus <- GenomicRanges::GRanges(
    seqnames = virus_chr,
    ranges = IRanges::IRanges(virus_pos, width = 1),
    strand = raw$viral_strand
  )

  suppressWarnings(
    InteractionSet::GInteractions(
      gr_host,
      gr_virus,
      mode = "strict",
      Depth = as.numeric(raw$reads),
      Label = raw$sample,
      Source = raw$method,
      ViralStrand = raw$viral_strand
    )
  )
}

#' Draw Ideogram Track in Genomic Mode
#'
#' @param height Numeric. Track height.
#' @param cfg Plotting configuration list.
#' @return Invisibly returns \code{NULL}.
#' @export
draw_ideogram <- function(height, cfg) {
  ideogram_df <- cfg$data[, c("chr", "start", "end"), drop = FALSE]

  circlize::circos.genomicTrackPlotRegion(
    data = ideogram_df,
    ylim = c(-0.2, 1.5),
    bg.border = NA,
    track.height = height,
    panel.fun = function(region, value, ...) {
      chr <- circlize::CELL_META$sector.index
      xlim <- circlize::CELL_META$cell.xlim
      fill_col <- cfg$grid_col[chr]

      if (length(fill_col) == 0 || is.na(fill_col)) {
        fill_col <- "grey90"
      }

      circlize::circos.genomicRect(
        region = region,
        value = value,
        ybottom = 0,
        ytop = 1,
        col = fill_col,
        border = "white"
      )

      circlize::circos.text(
        mean(xlim), 1.05, chr,
        facing = "bending.inside",
        niceFacing = TRUE,
        adj = c(0.5, 0),
        cex = 0.8
      )

      axis_labels <- make_axis_labels(xlim, use_kb = chr == cfg$virus_name)
      circlize::circos.axis(
        h = 1.02,
        major.at = axis_labels$at,
        labels = axis_labels$labels,
        labels.cex = 0.4,
        direction = "outside",
        major.tick.length = 0.02
      )
    }
  )

  invisible(NULL)
}

#' Draw Scatter Track in Genomic Mode
#'
#' @param data Data frame converted from a \code{GInteractions} object.
#' @param height Numeric. Track height.
#' @param cfg Plotting configuration list.
#' @param label_col Label column name. Kept for API compatibility.
#' @param point_color Default point color.
#' @return Invisibly returns \code{NULL}.
#' @export
draw_scatter <- function(data, height, cfg, label_col = "Label", point_color = "blue") {
  if (is.null(data) || nrow(data) == 0) {
    return(invisible(NULL))
  }

  max_depth <- if ("Depth" %in% colnames(data)) max(log10(data$Depth + 1), na.rm = TRUE) else 1
  if (!is.finite(max_depth) || max_depth <= 0) {
    max_depth <- 1
  }

  depth <- if ("Depth" %in% colnames(data)) data$Depth else rep(10, nrow(data))
  source_col <- if ("Source" %in% colnames(data)) data$Source else rep("integration", nrow(data))
  col_map <- c("Pos" = "#FB9A9999", "PB" = "#56B4E999")
  pt_col <- col_map[source_col]
  pt_col[is.na(pt_col)] <- point_color

  scatter_df <- data.frame(
    chr = as.character(data$seqnames1),
    start = data$start1,
    end = data$end1,
    y = 0.5,
    cex = log10(depth + 1) / max_depth * 1.5 + 0.5,
    pt_col = pt_col,
    stringsAsFactors = FALSE
  )

  circlize::circos.genomicTrackPlotRegion(
    data = scatter_df,
    ylim = c(-0.2, 1.2),
    track.height = height,
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      circlize::circos.lines(circlize::CELL_META$cell.xlim, c(0.5, 0.5), col = "grey90")

      if (nrow(region) > 0) {
        circlize::circos.genomicPoints(
          region,
          value,
          numeric.column = "y",
          pch = 19,
          col = value$pt_col,
          cex = value$cex
        )
      }
    }
  )

  invisible(NULL)
}

#' Draw Link Layer in Genomic Mode
#'
#' @param link_data Data frame converted from a \code{GInteractions} object.
#' @param cfg Plotting configuration list.
#' @param radius Optional link radius.
#' @return Invisibly returns \code{NULL}.
#' @export
draw_link <- function(link_data, cfg, radius = NULL) {
  if (is.null(link_data) || nrow(link_data) == 0) {
    return(invisible(NULL))
  }

  source_col <- if ("Source" %in% colnames(link_data)) link_data$Source else rep("integration", nrow(link_data))
  col_map <- c("Pos" = "#FB9A9999", "PB" = "#56B4E999")

  for (i in seq_len(nrow(link_data))) {
    link_col <- col_map[source_col[i]]
    if (is.na(link_col)) {
      link_col <- "grey"
    }

    region1 <- data.frame(
      chr = as.character(link_data$seqnames1[i]),
      start = link_data$start1[i],
      end = link_data$end1[i]
    )
    region2 <- data.frame(
      chr = as.character(link_data$seqnames2[i]),
      start = link_data$start2[i],
      end = link_data$end2[i]
    )

    if (is.null(radius)) {
      circlize::circos.genomicLink(region1, region2, col = link_col, border = link_col)
    } else {
      circlize::circos.genomicLink(region1, region2, rou = radius, col = link_col, border = link_col)
    }
  }

  invisible(NULL)
}

validate_layout_list <- function(layout_list) {
  if (is.null(layout_list) || length(layout_list) == 0) {
    stop("layout_list must contain at least one track definition.")
  }

  invalid_idx <- which(vapply(layout_list, function(x) is.null(x$type) || !nzchar(x$type), logical(1)))
  if (length(invalid_idx) > 0) {
    stop("Every layout item must include a non-empty type. Invalid index: ",
         paste(invalid_idx, collapse = ", "))
  }

  invisible(NULL)
}

make_axis_labels <- function(xlim, use_kb = FALSE) {
  axis_unit <- if (use_kb) "kb" else "Mb"
  axis_scale <- if (use_kb) 1000 else 1e6
  axis_len <- diff(xlim) / axis_scale
  tick_step <- if (use_kb) {
    if (axis_len <= 10) 1 else if (axis_len <= 50) 10 else 50
  } else {
    50
  }

  label_values <- seq(0, axis_len, by = tick_step)
  major_values <- label_values
  major_labels <- paste0(label_values, " ", axis_unit)

  if (length(major_values) == 0 || tail(major_values, 1) < axis_len) {
    major_values <- c(major_values, axis_len)
    major_labels <- c(major_labels, "")
  }

  list(
    at = xlim[1] + major_values * axis_scale,
    labels = major_labels,
    unit = axis_unit
  )
}

match_chr_style <- function(chr, valid_chr) {
  chr <- as.character(chr)
  valid_chr <- as.character(valid_chr)

  stripped <- sub("^chr", "", chr, ignore.case = TRUE)
  prefixed <- ifelse(grepl("^chr", chr, ignore.case = TRUE), chr, paste0("chr", chr))

  out <- chr
  use_stripped <- !out %in% valid_chr & stripped %in% valid_chr
  out[use_stripped] <- stripped[use_stripped]

  use_prefixed <- !out %in% valid_chr & prefixed %in% valid_chr
  out[use_prefixed] <- prefixed[use_prefixed]

  out
}