#' Visualize Viral Integration Events
#' Generates a Circos plot to visualize viral integration sites on host chromosomes.
#' @param input_file Path to the viral integration results file (Tab-delimited).
#' @param host_file Path to the chromosome sizes file (Must contain host chromosomes AND the virus info).
#' @param output_file Path to the output PDF file. Default is "Result.pdf".
#' @param color_file Optional path to a CSV file specifying chromosome colors (No header, Col1=Name, Col2=Color).
#' @param visual_ratio The visual proportion of the circle occupied by the virus sector. Default is 0.1 (10%).
#'
#' @return None. A PDF file is generated at the specified path.
#' @export
#' @importFrom grDevices pdf dev.off
#' @examples
#' \dontrun{
#' # 获取包内的示例数据路径
#' data_path <- system.file("extdata", "data.txt", package = "ViHostPlot")
#' host_path <- system.file("extdata", "chrom_sizes.txt", package = "ViHostPlot")
#'
#' # 定义布局
#' layout_demo <- list(
#'   list(type = "ideogram", height = 0.05),
#'   list(type = "scatter", sample_label = "T", height = 0.15),
#'   list(type = "links")
#' )
#'
#' # 运行绘图
#' visualize_viral_integration(
#'   input_file = data_path,
#'   host_file = host_path,
#'   layout_list = layout_demo,
#'   output_file = "Example_Result.pdf"
#' )
#' }
visualize_viral_integration <- function(input_file, host_file, layout_list, output_file="Result.pdf",
                                        color_file=NULL, visual_ratio=0.1) {

  temp_vcf <- tempfile(pattern = "integration_", fileext = ".vcf")
  on.exit(if(file.exists(temp_vcf)) file.remove(temp_vcf))

  # Extract Virus Info
  chrom_sizes <- read_chrom_sizes(host_file)
  virus_info  <- get_virus_name(chrom_sizes)

  # Convert to VCF
  raw_to_vcf_simple(input_file, temp_vcf, host_file, virus_info)

  # Create Config & Object
  cfg <- create_config(host_file, virus_info, visual_ratio, color_file)
  gi_all <- create_gi(temp_vcf)

  all_data_df <- as.data.frame(gi_all)

  # plot
  if (!is.null(output_file)) {
    grDevices::pdf(output_file, width=10, height=10)
  }

  circlize::circos.clear()

  n_sectors <- nrow(cfg$data)
  gaps <- rep(1, n_sectors)
  v_idx <- which(cfg$data$chr == cfg$virus_name)
  if(length(v_idx)>0) {
    #if(v_idx > 1) gaps[v_idx-1] <- 15 else gaps[n_sectors] <- 15
    gaps[v_idx] <- 15
  }

  circlize::circos.par(start.degree = 90, gap.degree = gaps,
                       cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE)
  circlize::circos.initialize(factors=cfg$data$chr, xlim=cbind(0, cfg$data$end), sector.width=cfg$widths)

  for (i in seq_along(layout_list)) {
    task <- layout_list[[i]]
    type <- task$type

    # default height
    h <- if(is.null(task$height)) 0.05 else task$height

    if (type == "ideogram") {
      draw_ideogram(height = h, cfg = cfg)

    } else if (type == "scatter") {
      if(is.null(task$sample_label)) {
        warning(paste("Skipping scatter track", i, ": missing 'sample_label'"))
        next
      }
      sub_data <- all_data_df[all_data_df$Label %in% task$sample_label, ]
      draw_scatter(data = sub_data, height = h, cfg = cfg)

    } else if (type == "links") {
      draw_link(link_data = all_data_df, cfg = cfg, radius = task$radius)
    }
  }
  if (!is.null(output_file)) {
    grDevices::dev.off()

  }
}


#' Read Chromosome Sizes File
#' @param path Path to the chromosome sizes text file (two columns: name, length).
#' @return A data frame with columns 'chr' and 'end'.
#' @keywords internal
read_chrom_sizes <- function(path) {
  if (!file.exists(path)) stop("Error: File not found: ", path)
  df <- utils::read.table(path, header = FALSE, stringsAsFactors = FALSE)
  df <- df[, 1:2]
  colnames(df) <- c("chr", "end")
  return(df)
}


#' Detect Virus Name from Chromosome Sizes
#' logic: Exclusion of mitochondria, then finding the shortest sequence.
#' @param chrom_sizes_df Data frame with 'chr' and 'end'.
#' @return A list containing 'name' and 'length' of the virus.
#' @keywords internal
get_virus_name <- function(chrom_sizes_df) {
  # Exclude mitochondria
  is_mito <- grepl("^(chr)?(m|mt)$", chrom_sizes_df$chr, ignore.case = TRUE)
  candidates <- chrom_sizes_df[!is_mito, ]
  if(nrow(candidates) == 0) stop("Error: No valid chromosomes found after excluding mitochondria.")

  # Find the shortest sequence
  candidates <- candidates[order(candidates$end), ]
  target_virus <- candidates$chr[1]
  target_len   <- candidates$end[1]
  return(list(name = target_virus, length = target_len))
}


#' Convert Raw Table to Temporary VCF
#' @param input_file Path to input tabular data.
#' @param output_vcf Path to output temporary VCF.
#' @param host_file Path to host chromosome sizes file.
#' @param virus_info List containing virus name and length.
#' @importFrom dplyr filter transmute bind_rows mutate select
#' @importFrom utils read.table write.table
#' @importFrom stringr str_extract
#' @keywords internal
raw_to_vcf_simple <- function(input_file, output_vcf, host_file, virus_info) {

  chrom_sizes <- read_chrom_sizes(host_file)
  raw <- utils::read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  # Normalize Chromosome Names (Remove 'chr')
  if("Chr" %in% colnames(raw)) {
    raw$Chr <- gsub("^chr", "", raw$Chr, ignore.case = TRUE)
  }

  virus_name <- virus_info$name

  # Clean Data
  d1_filtered <- dplyr::filter(raw, !is.na(.data$Pos), !is.na(.data$Pos_HBV))
  d1 <- dplyr::transmute(d1_filtered, C=.data$Chr, Pos=.data$Pos, D=.data$Pos_depth, V_Str=.data$Pos_HBV, Lbl=.data$sample_label, Src="Pos")
  d2_filtered <- dplyr::filter(raw, !is.na(.data$PB), !is.na(.data$PB_HBV))
  d2 <- dplyr::transmute(d2_filtered, C=.data$Chr, Pos=.data$PB, D=.data$PB_depth, V_Str=.data$PB_HBV, Lbl=.data$sample_label, Src="PB")
  tidy <- dplyr::bind_rows(d1, d2)

  if(nrow(tidy) == 0) stop("Error: No valid data found in input file.")

  # Parse Virus Info
  tidy$V_Pos <- as.numeric(stringr::str_extract(tidy$V_Str, "^\\d+"))
  tidy$V_Strand <- ifelse(grepl("-", tidy$V_Str), "-", "+")

  if (all(is.na(tidy$V_Pos))) stop("Error: Failed to parse virus position.")

  # Build Header
  header_lines <- c(
    "##fileformat=VCFv4.2",
    "##source=ViHostPlot_Simple",
    sprintf("##contig=<ID=%s,length=%d>", chrom_sizes$chr, chrom_sizes$end),
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">',
    '##INFO=<ID=LBL,Number=1,Type=String,Description="Label">',
    '##INFO=<ID=SRC,Number=1,Type=String,Description="Source">',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
  )

  # Build Body
  alt_col <- ifelse(tidy$V_Strand == "+",
                    paste0("N[", virus_name, ":", tidy$V_Pos, "["),
                    paste0("N]", virus_name, ":", tidy$V_Pos, "]"))

  dp_safe <- ifelse(is.na(tidy$D), 0, tidy$D)
  body <- data.frame(
    CHROM  = tidy$C,
    POS    = as.integer(tidy$Pos),
    ID     = paste0("evt_", seq_len(nrow(tidy))),
    REF    = "N",
    ALT    = alt_col,
    QUAL   = ".",
    FILTER = "PASS",
    INFO   = paste0("SVTYPE=BND;LBL=", tidy$Lbl, ";SRC=", tidy$Src, ";DP=", dp_safe),
    FORMAT = "GT",
    SAMPLE = "0/1",
    stringsAsFactors = FALSE
  )

  writeLines(header_lines, output_vcf)
  utils::write.table(body, output_vcf, append=TRUE, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}


#' Create GInteractions Object from VCF
#' @param vcf_file Path to the VCF file.
#' @importFrom VariantAnnotation readVcf info alt
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomicRanges GRanges mcols mcols<- strand<-
#' @importFrom IRanges IRanges
#' @importFrom InteractionSet GInteractions
#' @importFrom stringr str_match
#' @return A GInteractions object.
#' @keywords internal
create_gi <- function(vcf_file) {
  vcf <- VariantAnnotation::readVcf(vcf_file)
  gr_host <- SummarizedExperiment::rowRanges(vcf)
  vcf_info <- VariantAnnotation::info(vcf)

  alt_obj <- VariantAnnotation::alt(vcf)
  alt_strings <- as.character(unlist(alt_obj))

  matches <- stringr::str_match(alt_strings, "[\\[\\]]([^:]+):(\\d+)[\\[\\]]")

  gr_virus <- GenomicRanges::GRanges(seqnames = matches[,2], ranges = IRanges::IRanges(as.numeric(matches[,3]), width=1))
  GenomicRanges::strand(gr_virus) <- ifelse(grepl("\\[", alt_strings), "+", "-")

  gi <- suppressWarnings(InteractionSet::GInteractions(gr_host, gr_virus, mode="strict"))
  GenomicRanges::mcols(gi)$Depth <- vcf_info$DP
  GenomicRanges::mcols(gi)$Label <- vcf_info$LBL
  GenomicRanges::mcols(gi)$Source <- vcf_info$SRC
  return(gi)
}


#' Create Plotting Configuration
#' @param host_file Path to host file.
#' @param virus_info List of virus name and length.
#' @param visual_ratio Visual ratio for the virus sector.
#' @param color_file Optional path to color CSV.
#' @importFrom stats setNames
#' @importFrom grDevices rainbow
#' @importFrom utils read.csv
#' @return A list containing plotting configuration.
#' @keywords internal
create_config <- function(host_file, virus_info, visual_ratio, color_file) {

  full_data <- read_chrom_sizes(host_file)
  virus_name <- virus_info$name

  full_data$start <- 0

  # Calculate widths
  widths <- numeric(nrow(full_data)); names(widths) <- full_data$chr
  widths[virus_name] <- visual_ratio

  host_mask <- full_data$chr != virus_name
  widths[host_mask] <- (full_data$end[host_mask] / sum(full_data$end[host_mask])) * (1 - visual_ratio)

  # Colors
  if (!is.null(color_file) && file.exists(color_file)) {
    c_df <- utils::read.csv(color_file, header=FALSE, stringsAsFactors=FALSE)
    grid_col <- stats::setNames(c_df[,2], c_df[,1])
    if(!virus_name %in% names(grid_col)) grid_col[virus_name] <- "grey90"
  } else {
    grid_col <- stats::setNames(grDevices::rainbow(nrow(full_data)), full_data$chr)
  }

  list(data=full_data, widths=widths, virus_name=virus_name, grid_col=grid_col)
}

#' Plan Track Layout
#' @param labels Vector of sample labels.
#' @param space Total available space.
#' @return Named vector of track heights.
#' @keywords internal
plan_layout <- function(labels, space) {
  stats::setNames(rep(space/length(labels), length(labels)), labels)
}


#' Draw Ideogram Track
#' @param height Numeric. The height of the track
#' @param cfg List. The configuration list containing `grid_col` (colors).
#' @importFrom circlize circos.track circos.rect circos.text CELL_META
#' @export
draw_ideogram <- function(height, cfg) {

  circlize::circos.track(ylim=c(0, 1), bg.border=NA, track.height=height,
                         panel.fun=function(x, y) {
                           chr <- circlize::CELL_META$sector.index
                           xlim <- circlize::CELL_META$cell.xlim

                           circlize::circos.rect(xlim[1], 0, xlim[2], 1, col=cfg$grid_col[chr], border="white")

                           circlize::circos.text(mean(xlim), 1, chr,
                                                 facing="bending.inside", niceFacing=TRUE,
                                                 adj = c(0.5, -1.5), cex=0.8)

                           circlize::circos.axis(h = "top", labels.cex = 0.4, direction = "outside")
                         })
}

#' Draw Link Layer (Arcs)
#' @param link_data Data frame. Must contain columns: seqnames1, start1, seqnames2, start2, Source(optional).
#' @param cfg List. Configuration list containing data about chromosomes.
#' @param radius Numeric (Optional). The radius to draw links on. If NULL, defaults to the bottom of the current track.
#' @importFrom circlize circos.link get.cell.meta.data
#' @export
draw_link <- function(link_data, cfg, radius = NULL) {
  if(is.null(link_data) || nrow(link_data) == 0) return()

  if (is.null(radius)) {
    radius <- circlize::get.cell.meta.data("cell.bottom.radius", sector.index = cfg$data$chr[1])
  }

  # 颜色映射
  col_map <- c("Pos"="#FB9A9999", "PB"="#56B4E999")

  for(i in seq_len(nrow(link_data))) {
    src <- if("Source" %in% colnames(link_data)) link_data$Source[i] else "Pos"
    lnk_col <- if(src %in% names(col_map)) col_map[src] else "grey"

    circlize::circos.link(
      link_data$seqnames1[i], link_data$start1[i],
      link_data$seqnames2[i], link_data$start2[i],
      col = lnk_col,
      h.ratio = 0.5, # 弧度高度
      rou = radius   # 关键：控制连线画在多大的圈上
    )
  }
}

#' Draw Scatter Plot Track
#' @param data Data frame. The full dataset (must be converted from GI object).
#'             Must contain columns: seqnames1, start1, Depth.
#' @param height Numeric. The height of the track.
#' @param cfg List. Configuration list.
#' @param label_col String. The column name in `data` to use for identifying the sample.
#' @param point_color String. Color of the points. Default is "blue".
#' @importFrom circlize circos.track circos.lines circos.points circos.text CELL_META
#' @export
draw_scatter <- function(data, height, cfg, label_col="Label", point_color="blue") {

  if(is.null(data) || nrow(data) == 0) return()

  # 计算全局的最大深度，用于统一缩放点的大小
  max_depth_all <- if("Depth" %in% colnames(data)) max(log10(data$Depth + 1), na.rm=TRUE) else 1

  col_map <- c("Pos"="#FB9A9999", "PB"="#56B4E999")

  circlize::circos.track(ylim=c(0,1), track.height=height, bg.border=NA, panel.fun=function(x,y) {
    # 画一条灰色的中心线或底线
    circlize::circos.lines(circlize::CELL_META$cell.xlim, c(0,0), col="grey90")

    # 筛选属于当前扇区(sector)的点
    curr_chr <- circlize::CELL_META$sector.index
    sub_df <- data[data$seqnames1 %in% curr_chr, , drop=FALSE]

    if(nrow(sub_df) > 0) {
      val <- if("Depth" %in% colnames(sub_df)) sub_df$Depth else rep(10, nrow(sub_df))
      cex_val <- log10(val + 1) / max_depth_all * 1.5 + 0.5

      pt_cols <- rep(point_color, nrow(sub_df))
      if("Source" %in% colnames(sub_df)) {
        mapped_cols <- col_map[sub_df$Source]
        valid_idx <- !is.na(mapped_cols)
        pt_cols[valid_idx] <- mapped_cols[valid_idx]
      }

      circlize::circos.points(sub_df$start1, rep(0.5, nrow(sub_df)),
                              pch=19, col=pt_cols, cex=cex_val)
    }
  })
}


