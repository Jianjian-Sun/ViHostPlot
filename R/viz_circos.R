#' Global Style Configuration for ViHostPlot
#' @param ... Style parameters to override defaults.
#' @export
create_global_style <- function(...) {
  default_style <- list(
    start_degree = 90, gap_degree_host = 1, gap_degree_virus = 15,
    color_pos = "#FB9A99", color_pb = "#56B4E9",
    point_size_min = 0.5, point_size_max = 2.0, alpha_links = "CC"
    # ... 其余参数保持一致 ...
  )
  utils::modifyList(default_style, list(...))
}

#' Draw Viral Integration Tracks
#' @param gi Normalized GenomicInteractions object.
#' @param config GenomeConfig object.
#' @param style Global style list.
#' @import circlize
#' @export
plot_virus_integration <- function(gi, config, style = create_global_style()) {
  layout_info <- calculate_layout_metrics(gi, radius_start = 0.95, radius_end = 0.40)

  # 初始化逻辑
  initialize_circos(config, style) # 这里假设您已将 initialize_circos 放入包中

  # 循环绘制轨道
  for (lbl in layout_info$labels) {
    draw_single_track(gi = gi, label = lbl, config = config,
                      height = layout_info$track_height, style = style)
  }

  # 绘制最内层连线
  draw_links_separately(gi = gi, labels_to_plot = layout_info$labels, style = style)

  message("ViHostPlot Circos drawing completed.")
}
