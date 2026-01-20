#' Parse Viral Position String
#' @description Extracts numeric position and strand from strings like "826 (-)".
#' @param pos_strs A character vector of viral position strings.
#' @return A list with numeric positions and strand characters.
#' @export
parse_viral_string <- function(pos_strs) {
  vals <- as.numeric(stringr::str_extract(pos_strs, "\\d+"))
  strands <- ifelse(grepl("-", pos_strs), "-", "+")
  return(list(pos = vals, strand = strands))
}


#' Calculate Circular Layout Metrics
#' @description Automatically computes track heights based on sample count.
#' @param gi A GenomicInteractions object containing a 'Label' column.
#' @param radius_start Outer radius starting point, default 0.95.
#' @param radius_end Inner radius ending point, default 0.30.
#' @return A list of layout parameters including track_height and labels.
#' @export
calculate_layout_metrics <- function(gi, radius_start = 0.95, radius_end = 0.30) {
  if (!"Label" %in% names(S4Vectors::mcols(gi))) stop("GI object missing Label column")

  target_labels <- unique(S4Vectors::mcols(gi)$Label)
  num_tracks <- length(target_labels)
  if (num_tracks == 0) stop("No sample labels detected")

  available_height <- radius_start - radius_end
  single_height <- available_height / num_tracks

  return(list(
    track_height = single_height,
    labels = target_labels,
    num_tracks = num_tracks,
    radius_start = radius_start,
    radius_end = radius_end
  ))
}
