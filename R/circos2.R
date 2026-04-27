#' Visualize Host-Virus Integration Events
#'
#' Draw a host-virus circos plot on the current graphics device from an
#' integration table and host chromosome-size information. The host genome can
#' be supplied through a standard host name such as \code{"human"} or
#' \code{"mouse"}, which is resolved online to a UCSC assembly and fetched as
#' chromosome sizes, or through a custom host chromosome-size table. The virus
#' sequence name and length are supplied explicitly. No intermediate VCF is
#' created and no output file is opened by this function.
#'
#' @param input_file Path to a tab-delimited integration table with columns
#'   \code{chr}, \code{host_loc}, \code{viral_loc}, \code{reads},
#'   \code{sample}, \code{viral_strand}, and \code{method}.
#' @param host Character scalar specifying a standard host genome name or a
#'   UCSC assembly name. Common aliases such as \code{"human"} and
#'   \code{"mouse"} are supported.
#' @param chrom_file Optional path to a host chromosome-size table with columns
#'   \code{chr}, \code{start}, and \code{end}. Use this when \code{host} is not
#'   supplied or when a custom host genome is needed. This file should contain
#'   host chromosomes only; the virus sequence is added separately from
#'   \code{virus_name} and \code{virus_length}.
#' @param virus_name Name of the virus sequence.
#' @param virus_length Length of the virus sequence in base pairs.
#' @param layout_list A list of track definitions, for example
#'   \code{list(list(type = "ideogram", height = 0.08), list(type = "scatter",
#'   sample_label = "T", height = 0.15), list(type = "scatter",
#'   sample_label = "N", height = 0.15), list(type = "links"))}.
#'   Supported \code{type} values are \code{"ideogram"}, \code{"scatter"},
#'   \code{"histogram"}, and \code{"links"}.
#' @param visual_ratio Visual proportion assigned to the virus sector.
#' @param clear Logical. Whether to clear the existing circlize plot before
#'   drawing.
#'
#' @return Invisibly returns a list with \code{cfg}, \code{gi}, and
#'   \code{data}.
#' @export
visualize_viral_integration <- function(input_file,
                                        host = NULL,
                                        chrom_file = NULL,
                                        virus_name,
                                        virus_length,
                                        layout_list,
                                        visual_ratio = 0.1,
                                        clear = TRUE) {
  validate_layout_list(layout_list)

  cfg <- create_config(
    host = host,
    chrom_file = chrom_file,
    virus_name = virus_name,
    virus_length = virus_length,
    visual_ratio = visual_ratio
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
    gaps[virus_idx] <- 10
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
      draw_ideogram(
        height = height,
        cfg = cfg,
        grid_col = task$grid_col,
        border_col = task$border_col
      )
    } else if (task$type == "scatter") {
      sub_df <- filter_plot_data(plot_df, sample_label = task$sample_label)
      draw_scatter(
        data = sub_df,
        height = height,
        cfg = cfg,
        method_col = task$method_col,
        point_color = task$point_color,
        baseline_col = task$baseline_col
      )
    } else if (task$type == "histogram") {
      sub_df <- filter_plot_data(plot_df, sample_label = task$sample_label)
      draw_histogram(
        data = sub_df,
        height = height,
        cfg = cfg,
        bins = task$bins,
        col = task$col
      )
    } else if (task$type == "links") {
      draw_link(
        link_data = plot_df,
        cfg = cfg,
        radius = task$radius,
        method_col = task$method_col,
        default_col = task$default_col
      )
    } else {
      warning("Skipping unsupported track type: ", task$type)
    }
  }

  legend_spec <- get_method_legend_spec(layout_list, plot_df)
  if (!is.null(legend_spec)) {
    draw_method_legend(legend_spec)
  }

  invisible(list(cfg = cfg, gi = gi, data = plot_df))
}

ensure_host_chrom_utils_loaded <- function() {
  needed <- c(
    "read_chrom_sizes",
    "resolve_ucsc_assembly",
    "fetch_ucsc_chrom_sizes",
    "resolve_host_chrom_sizes"
  )

  if (all(vapply(needed, exists, logical(1), mode = "function", inherits = TRUE))) {
    return(invisible(TRUE))
  }

  helper_path <- file.path(getwd(), "host_chrom_utils.R")
  if (!file.exists(helper_path)) {
    stop("host_chrom_utils.R was not found in the current working directory.", call. = FALSE)
  }

  sys.source(helper_path, envir = .GlobalEnv)

  if (!all(vapply(needed, exists, logical(1), mode = "function", inherits = TRUE))) {
    stop("Failed to load host chromosome helper functions from host_chrom_utils.R.", call. = FALSE)
  }

  invisible(TRUE)
}

validate_virus_info <- function(virus_name, virus_length, host_chr) {
  if (missing(virus_name) || is.null(virus_name) || length(virus_name) != 1 ||
      is.na(virus_name) || !nzchar(virus_name)) {
    stop("virus_name must be a single non-empty sequence name.")
  }

  virus_name <- trimws(gsub("(?i)^chr", "", as.character(virus_name)))

  if (!is.numeric(virus_length) || length(virus_length) != 1 ||
      is.na(virus_length) || virus_length <= 0) {
    stop("virus_length must be a single positive number.", call. = FALSE)
  }

  virus_length <- as.numeric(virus_length)

  if (virus_name %in% host_chr) {
    stop("virus_name conflicts with a host chromosome name: ", virus_name, call. = FALSE)
  }

  list(
    virus_name = virus_name,
    virus_length = virus_length
  )
}

#' Create Plotting Configuration
#'
#' @param host Character scalar specifying a built-in host genome.
#' @param chrom_file Optional path to a host chromosome-size table.
#' @param virus_name Name of the virus sequence.
#' @param virus_length Length of the virus sequence in base pairs.
#' @param visual_ratio Visual proportion assigned to the virus sector.
#' @return A plotting configuration list.
#' @keywords internal
create_config <- function(host = NULL,
                          chrom_file = NULL,
                          virus_name,
                          virus_length,
                          visual_ratio = 0.1) {
  ensure_host_chrom_utils_loaded()
  host_df <- resolve_host_chrom_sizes(host = host, chrom_file = chrom_file)
  virus_info <- validate_virus_info(virus_name, virus_length, host_df$chr)
  virus_name <- virus_info$virus_name
  virus_length <- virus_info$virus_length

  if (!is.numeric(visual_ratio) || length(visual_ratio) != 1 ||
      is.na(visual_ratio) || visual_ratio <= 0 || visual_ratio >= 1) {
    stop("visual_ratio must be a single number between 0 and 1.")
  }

  chrom_df <- rbind(
    host_df,
    data.frame(chr = virus_name, start = 0, end = virus_length, stringsAsFactors = FALSE)
  )

  widths <- numeric(nrow(chrom_df))
  names(widths) <- chrom_df$chr
  widths[virus_name] <- visual_ratio

  host_mask <- chrom_df$chr != virus_name
  host_len <- chrom_df$end[host_mask] - chrom_df$start[host_mask]
  widths[host_mask] <- host_len / sum(host_len) * (1 - visual_ratio)

  list(
    data = chrom_df,
    widths = widths,
    virus_name = virus_name
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

  raw$chr <- trimws(gsub("(?i)^chr", "", as.character(raw$chr), perl = TRUE))
  host_chr <- raw$chr
  virus_chr <- rep(cfg$virus_name, nrow(raw))
  host_pos <- as.numeric(raw$host_loc)
  virus_pos <- as.numeric(raw$viral_loc)

  keep <- !is.na(host_chr) & host_chr %in% cfg$data$chr & !is.na(host_pos) & !is.na(virus_pos)
  if (!all(keep)) {
    warning(sum(!keep), " rows were skipped because chromosome names or position values were invalid.")
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
#' @param grid_col Optional named vector of sector colors.
#' @param border_col Border color for ideogram rectangles.
#' @return Invisibly returns \code{NULL}.
#' @export
draw_ideogram <- function(height, cfg, grid_col = NULL, border_col = "white") {
  ideogram_df <- cfg$data[, c("chr", "start", "end"), drop = FALSE]
  grid_col <- resolve_ideogram_colors(cfg, grid_col = grid_col)

  circlize::circos.genomicTrackPlotRegion(
    data = ideogram_df,
    ylim = c(-0.2, 1.5),
    bg.border = NA,
    track.height = height,
    panel.fun = function(region, value, ...) {
      chr <- circlize::CELL_META$sector.index
      xlim <- circlize::CELL_META$cell.xlim
      fill_col <- grid_col[chr]

      circlize::circos.genomicRect(
        region = region,
        value = value,
        ybottom = 0,
        ytop = 1,
        col = fill_col,
        border = border_col
      )

      circlize::circos.text(
        mean(xlim), 1.2, chr,
        facing = "bending.inside",
        niceFacing = TRUE,
        adj = c(0.5, -0.4),
        cex = 0.95
      )

      axis_labels <- make_axis_labels(xlim, use_kb = chr == cfg$virus_name)
      circlize::circos.axis(
        h = 1.1,
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
#' @param method_col Optional named vector of colors for different methods.
#' @param point_color Default point color.
#' @param baseline_col Baseline color for the scatter track.
#' @return Invisibly returns \code{NULL}.
#' @export
draw_scatter <- function(data, height, cfg, method_col = NULL,
                         point_color = "blue", baseline_col = "grey90") {
  if (is.null(data) || nrow(data) == 0) {
    return(invisible(NULL))
  }

  max_depth <- if ("Depth" %in% colnames(data)) max(log10(data$Depth + 1), na.rm = TRUE) else 1
  if (!is.finite(max_depth) || max_depth <= 0) {
    max_depth <- 1
  }

  depth <- if ("Depth" %in% colnames(data)) data$Depth else rep(10, nrow(data))
  source_col <- if ("Source" %in% colnames(data)) as.character(data$Source) else rep("integration", nrow(data))
  col_map <- resolve_method_colors(source_col, method_col = method_col)
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
    ylim = c(0, 1),
    track.height = height,
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      circlize::circos.lines(circlize::CELL_META$cell.xlim, c(0.5, 0.5), col = baseline_col)

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
#' @param method_col Optional named vector of colors for different methods.
#' @param default_col Fallback link color.
#' @return Invisibly returns \code{NULL}.
#' @export
draw_link <- function(link_data, cfg, radius = NULL, method_col = NULL,
                      default_col = "grey") {
  if (is.null(link_data) || nrow(link_data) == 0) {
    return(invisible(NULL))
  }

  source_col <- if ("Source" %in% colnames(link_data)) as.character(link_data$Source) else rep("integration", nrow(link_data))
  col_map <- resolve_method_colors(source_col, method_col = method_col)

  for (i in seq_len(nrow(link_data))) {
    link_col <- col_map[source_col[i]]
    if (is.na(link_col)) {
      link_col <- default_col
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

#' Draw Histogram Track in Genomic Mode
#'
#' @param data Data frame converted from a \code{GInteractions} object.
#' @param height Numeric. Track height.
#' @param cfg Plotting configuration list.
#' @param bins Number of genomic bins per sector.
#' @param col Histogram fill color.
#' @return Invisibly returns \code{NULL}.
draw_histogram <- function(data, height, cfg, bins = NULL, col = "grey70") {
  if (is.null(data) || nrow(data) == 0) {
    return(invisible(NULL))
  }

  bins <- if (is.null(bins)) 50 else bins
  if (!is.numeric(bins) || length(bins) != 1 || is.na(bins) || bins < 1) {
    stop("bins must be a positive integer.")
  }
  bins <- as.integer(bins)

  hist_df <- make_histogram_data(data = data, cfg = cfg, bins = bins)
  if (nrow(hist_df) == 0) {
    return(invisible(NULL))
  }

  max_count <- max(hist_df$count, na.rm = TRUE)
  if (!is.finite(max_count) || max_count <= 0) {
    max_count <- 1
  }

  circlize::circos.genomicTrackPlotRegion(
    data = hist_df,
    ylim = c(0, max_count),
    track.height = height,
    bg.border = NA,
    panel.fun = function(region, value, ...) {
      if (nrow(region) > 0) {
        circlize::circos.genomicRect(
          region,
          value,
          ybottom = 0,
          ytop = value$count,
          col = col,
          border = NA
        )
      }
    }
  )

  invisible(NULL)
}

make_histogram_data <- function(data, cfg, bins = 50) {
  out <- vector("list", nrow(cfg$data))

  for (i in seq_len(nrow(cfg$data))) {
    chr <- cfg$data$chr[i]
    chr_start <- cfg$data$start[i]
    chr_end <- cfg$data$end[i]
    chr_pos <- data$start1[as.character(data$seqnames1) == chr]
    chr_pos <- chr_pos[!is.na(chr_pos)]

    breaks <- seq(chr_start, chr_end, length.out = bins + 1)
    if (length(unique(breaks)) < 2) {
      next
    }

    counts <- tabulate(
      findInterval(chr_pos, breaks, rightmost.closed = TRUE, all.inside = TRUE),
      nbins = bins
    )

    out[[i]] <- data.frame(
      chr = chr,
      start = breaks[-length(breaks)],
      end = breaks[-1],
      count = counts,
      stringsAsFactors = FALSE
    )
  }

  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) {
    return(data.frame(chr = character(), start = numeric(), end = numeric(), count = numeric()))
  }

  do.call(rbind, out)
}

draw_method_legend <- function(method_col) {
  if (is.null(method_col) || length(method_col) == 0) {
    return(invisible(NULL))
  }

  graphics::legend(
    "topright",
    legend = names(method_col),
    col = unname(method_col),
    pch = 19,
    title = "Methods",
    bty = "n",
    cex = 0.8,
    inset = 0.02
  )

  invisible(NULL)
}

get_method_legend_spec <- function(layout_list, plot_df) {
  for (task in layout_list) {
    if (task$type %in% c("scatter", "links")) {
      methods <- if ("Source" %in% colnames(plot_df)) as.character(plot_df$Source) else "integration"
      return(resolve_method_colors(methods, method_col = task$method_col))
    }
  }

  NULL
}

filter_plot_data <- function(plot_df, sample_label = NULL) {
  if (is.null(sample_label)) {
    return(plot_df)
  }

  plot_df[plot_df$Label %in% sample_label, , drop = FALSE]
}

resolve_ideogram_colors <- function(cfg, grid_col = NULL) {
  chr <- as.character(cfg$data$chr)
  if (is.null(grid_col)) {
    grid_col <- stats::setNames(grDevices::rainbow(length(chr)), chr)
  } else {
    grid_col <- normalize_named_colors(grid_col, chr)
  }

  if (!cfg$virus_name %in% names(grid_col) || is.na(grid_col[cfg$virus_name])) {
    grid_col[cfg$virus_name] <- "grey70"
  }

  grid_col
}

resolve_method_colors <- function(methods, method_col = NULL) {
  methods <- sort(unique(as.character(methods)))
  methods <- methods[!is.na(methods) & nzchar(methods)]
  if (length(methods) == 0) {
    methods <- "integration"
  }

  if (is.null(method_col)) {
    base_col <- c(
      Pos = "#D55E00",
      PB = "#0072B2",
      integration = "#009E73"
    )
    method_col <- base_col[methods]
    missing <- is.na(method_col)
    if (any(missing)) {
      method_col[missing] <- grDevices::hcl.colors(sum(missing), "Dark 3")
    }
    return(stats::setNames(unname(method_col), methods))
  }

  normalize_named_colors(method_col, methods)
}

normalize_named_colors <- function(colors, keys) {
  if (is.null(colors)) {
    return(stats::setNames(rep(NA_character_, length(keys)), keys))
  }

  if (!is.null(names(colors))) {
    matched <- colors[keys]
    missing <- is.na(matched)
    if (any(missing)) {
      matched[missing] <- rep(unname(colors)[1], sum(missing))
    }
    return(stats::setNames(unname(matched), keys))
  }

  if (length(colors) == 1) {
    return(stats::setNames(rep(colors, length(keys)), keys))
  }

  stats::setNames(rep(colors, length.out = length(keys)), keys)
}

validate_layout_list <- function(layout_list) {
  if (is.null(layout_list) || length(layout_list) == 0) {
    stop("layout_list must contain at least one track definition.")
  }

  invalid_obj_idx <- which(!vapply(layout_list, is.list, logical(1)))
  if (length(invalid_obj_idx) > 0) {
    stop("Every layout item must be a list. Invalid index: ",
         paste(invalid_obj_idx, collapse = ", "))
  }

  invalid_idx <- which(vapply(layout_list, function(x) is.null(x$type) || !nzchar(x$type), logical(1)))
  if (length(invalid_idx) > 0) {
    stop("Every layout item must include a non-empty type. Invalid index: ",
         paste(invalid_idx, collapse = ", "))
  }

  valid_types <- c("ideogram", "scatter", "histogram", "links")
  invalid_type_idx <- which(vapply(layout_list, function(x) !x$type %in% valid_types, logical(1)))
  if (length(invalid_type_idx) > 0) {
    stop("Unsupported layout type at index: ",
         paste(invalid_type_idx, collapse = ", "),
         ". Supported types: ", paste(valid_types, collapse = ", "),
         call. = FALSE)
  }

  scatter_missing_label <- which(vapply(
    layout_list,
    function(x) identical(x$type, "scatter") && (is.null(x$sample_label) || !nzchar(x$sample_label)),
    logical(1)
  ))
  if (length(scatter_missing_label) > 0) {
    stop("Each scatter layout item must include a non-empty sample_label. Invalid index: ",
         paste(scatter_missing_label, collapse = ", "),
         call. = FALSE)
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