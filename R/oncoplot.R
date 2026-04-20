#' Generate a oncoplot with clinical tracks
#'
#' @param maf MAF object.
#' @param genes the gene names or the number, default is 20.
#' @param clinical_vars A character vector or list defining clinical variables to display as bottom tracks.
#' @param n_breaks Integer. Number of breaks for continuous clinical variables. Default is 5.
#' @param colors A named list of customized colors for clinical tracks.
#' @param clinical_palettes A named list or character vector of RColorBrewer palettes for clinical tracks.
#'
#' @returns A combined plot object (class depends on \code{aplot} output).
#' @importFrom aplot insert_top insert_right insert_bottom
#' @export
#' @examples
#' \donttest{
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml.clin <- system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')
#' laml <- maftools::read.maf(maf = laml.maf, clinicalData = laml.clin)
#' var_names <- c("FAB_classification", "days_to_last_followup")
#' oncoplot(maf = laml, genes = 20, clinical_vars = var_names)
#' }
oncoplot <- function(maf, genes = 20, clinical_vars = NULL, n_breaks = 5,
                     colors = NULL, clinical_palettes = NULL) {

  p_main <- oncoplot_main(maf, genes)
  p_top <- aplotExtra:::oncoplot_sample(maf, genes)
  p_right <- aplotExtra:::oncoplot_gene(maf, genes, ylab = 'percentage')
  p_spacer <- ggplot() + ggfun::theme_transparent()

  # Base assembly using aplot
  pp <- p_main |>
    aplot::insert_top(p_spacer, height = 0.02) |>
    aplot::insert_top(p_top, height = 0.2) |>
    aplot::insert_right(p_right, width = 0.2)

  tracks <- oncoplot_clinical_track(maf, genes = genes, clinical_vars = clinical_vars,
                                    n_breaks = n_breaks, colors = colors,
                                    clinical_palettes = clinical_palettes)

  if (length(tracks) > 0) {
    for (var in names(tracks)) {
      pp <- aplot::insert_bottom(pp, tracks[[var]], height = 0.05)
    }
  }

  return(pp)
}


#' @importFrom ggplot2 geom_tile
#' @importFrom rlang .data
oncoplot_main <- function(maf, genes = 20) {
  d <- aplotExtra:::oncoplot_tidy_onco_matrix(maf, genes)

  p <- ggplot(d, aes(x = .data$Sample, y = .data$Gene, fill = .data$Type)) +
    geom_tile(colour="white", linewidth=.01) + 
    oncoplot_setting(continuous = FALSE, fill_name = "Mutation Type") +
    theme(
        legend.position = "right", 
        axis.text.y.left = element_text(face='italic')
    )
}

#' @importFrom maftools getClinicalData
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom rlang .data
#' @importFrom dplyr select mutate transmute all_of
oncoplot_clinical_track <- function(maf, genes = 20, clinical_vars, n_breaks,
                                    colors = NULL, clinical_palettes = NULL) {

  if (is.null(clinical_vars)) return(list())

  clinical_data <- as.data.frame(maftools::getClinicalData(maf))
  var_names <- intersect(
    if(is.list(clinical_vars)) names(clinical_vars) else clinical_vars,
    colnames(clinical_data)
  )

  onco_matrix <- aplotExtra:::oncoplot_tidy_onco_matrix(maf, genes)
  sample_order <- unique(as.character(onco_matrix$Sample))

  df <- clinical_data |>
    dplyr::select(.data$Tumor_Sample_Barcode, dplyr::all_of(var_names)) |>
    dplyr::mutate(Tumor_Sample_Barcode = factor(.data$Tumor_Sample_Barcode,
                                                levels = sample_order))

  plot_list <- list()

  for (var in var_names) {
    tmp <- dplyr::transmute(
      df,
      Tumor_Sample_Barcode = .data$Tumor_Sample_Barcode,
      value = .data[[var]],
      y_label = var
    )
    is_continuous <- is.numeric(tmp$value)

    tmp$value <- if (is_continuous) {
      binning_numeric(tmp$value, n_breaks = n_breaks)
    } else {
      as.factor(tmp$value)
    }

    # Used .data$ prefix here to prevent R CMD check "no visible binding" notes
    p <- ggplot(tmp, aes(x = .data$Tumor_Sample_Barcode, y = .data$y_label, fill = .data$value)) +
      geom_tile(color = "white", linewidth = 0.01) +
      theme_minimal() +
      theme(
        legend.position = "right",
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = 'italic', color = "black", size = 10),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(t = 5, r = 2, b = 0, l = 5)
      ) +
      labs(fill = var)

    if (!is.null(colors) && var %in% names(colors)) {
      p <- p + scale_fill_manual(values = colors[[var]], na.value = "grey90")
    } else {
      track_colors <- clinical_track_colors(
        tmp$value,
        continuous = is_continuous,
        var = var,
        clinical_palettes = clinical_palettes
      )
      if (!is.null(track_colors)) {
        p <- p + scale_fill_manual(values = track_colors, na.value = "grey90")
      }
    }
    plot_list[[var]] <- p
  }

  return(plot_list)
}

clinical_track_colors <- function(x, continuous = FALSE, var = NULL, clinical_palettes = NULL) {
  lvls <- levels(as.factor(x))
  n <- length(lvls)

  if (n == 0) {
    return(NULL)
  }

  palette <- clinical_track_palette(
    var = var,
    continuous = continuous,
    clinical_palettes = clinical_palettes
  )

  if (length(palette) > 1) {
    base_colors <- palette
    max_colors <- length(base_colors)
  } else {
    max_colors <- clinical_brewer_maxcolors(palette)
    base_n <- min(max(n, 3), max_colors)
    base_colors <- brewer.pal(base_n, palette)
  }

  if (n > max_colors) {
    values <- colorRampPalette(base_colors)(n)
  } else {
    values <- base_colors[seq_len(n)]
  }

  stats::setNames(values, lvls)
}

clinical_track_palette <- function(var = NULL, continuous = FALSE, clinical_palettes = NULL) {
  default_palette <- if (continuous) "YlGnBu" else "Set3"

  if (is.null(clinical_palettes)) {
    return(default_palette)
  }

  if (!is.null(var) && var %in% names(clinical_palettes)) {
    return(clinical_palettes[[var]])
  }

  type_names <- if (continuous) {
    c("continuous", "numeric", "default_continuous")
  } else {
    c("discrete", "categorical", "factor", "default_discrete")
  }
  matched_type <- intersect(type_names, names(clinical_palettes))

  if (length(matched_type) > 0) {
    return(clinical_palettes[[matched_type[1]]])
  }

  if (is.character(clinical_palettes) && is.null(names(clinical_palettes)) &&
      length(clinical_palettes) == 1) {
    return(clinical_palettes)
  }

  default_palette
}

clinical_brewer_maxcolors <- function(palette) {
  max_colors <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]

  if (length(max_colors) == 0 || is.na(max_colors)) {
    stop("Unknown RColorBrewer palette: ", palette, call. = FALSE)
  }

  max_colors
}

oncoplot_setting <- function(noxaxis = TRUE, continuous = TRUE, scale = 'y',
                             fill_name = NULL) {  
    list(
        theme_minimal(),
        if (noxaxis) ggfun::theme_noxaxis(),
        theme(legend.position = "none", panel.grid.major = element_blank()),
        aplotExtra:::oncoplot_scale(continuous = continuous, scale = scale),
        oncoplot_fill(name = fill_name),
        xlab(NULL),
        ylab(NULL)
    )
}

oncoplot_fill <- function(breaks = NULL, values = NULL, name = NULL, na.value = "#bdbdbd") {
    vc_col <- aplotExtra:::get_vcColors(websafe = FALSE)

    if (is.null(values)) {
        values <- vc_col
    } 

    if (is.null(breaks)) {
        vc_lev <- names(vc_col)
        breaks <- rev(vc_lev)
    }

    scale_fill_manual(
        name = name,
        breaks = breaks,
        values = values,
        na.value = na.value
    )
}


#' @param x A numeric vector
#' @importFrom stats quantile
#' @noRd
binning_numeric <- function(x, n_breaks = 5) {
  if (!is.numeric(x)) return(as.factor(x))

  probs <- seq(0, 1, length.out = n_breaks + 1)
  brks <- stats::quantile(x, probs = probs, na.rm = TRUE)
  brks <- unique(brks)

  if (length(brks) < 2) {
    return(as.factor(x))
  }

  binned <- cut(x, breaks = brks, include.lowest = TRUE)

  lvls <- levels(binned)
  lvls <- gsub("\\.0", "", lvls)
  lvls <- gsub("\\(|\\[|\\]", "", lvls)
  lvls <- gsub(",", "–", lvls)
  levels(binned) <- lvls

  return(binned)
}
