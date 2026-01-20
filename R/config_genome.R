#' Create Genome Configuration for Circos
#' @param host_file Path to chromosome sizes file.
#' @param virus_name Name of the virus (e.g., "HBV").
#' @param virus_length Length of the viral genome in bp.
#' @param virus_visual_ratio The ratio of the circle occupied by the virus.
#' @param color_file Optional path to a CSV with chromosome colors.
#' @importFrom grDevices rainbow
#' @export
create_genome_config <- function(host_file, virus_name, virus_length,
                                 virus_visual_ratio = 0.1, color_file = NULL) {
  # 注意：在包内应使用系统路径或用户提供路径，取消 setwd()
  host_data <- read.table(host_file, header = FALSE)[,1:2]
  colnames(host_data) <- c("chr", "end")

  virus_data <- data.frame(chr = virus_name, end = virus_length)
  full_data <- rbind(host_data, virus_data)
  full_data$start <- 0

  host_total_bp <- sum(host_data$end)
  host_widths <- (host_data$end / host_total_bp) * (1 - virus_visual_ratio)
  widths_vec <- c(host_widths, virus_visual_ratio)
  names(widths_vec) <- full_data$chr

  # 默认颜色逻辑
  if (is.null(color_file)) {
    grid_col <- setNames(grDevices::rainbow(nrow(full_data)), full_data$chr)
  } else {
    c_df <- read.csv(color_file, header = FALSE)
    grid_col <- setNames(c_df[,2], c_df[,1])
  }

  structure(list(data = full_data, widths = widths_vec,
                 virus_name = virus_name, grid_col = grid_col),
            class = "GenomeConfig")
}
