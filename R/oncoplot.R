#' Generate a oncoplot with clinical tracks
#'
#' @param maf MAF object.
#' @param genes the gene names or the number, default is 20.
#' @param clinical_vars A character vector or list defining clinical variables to display as bottom tracks.
#' @param n_breaks Integer. Number of breaks for continuous clinical variables. Default is 5.
#' @param colors A named list of customized colors for clinical tracks.
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
oncoplot <- function(maf, genes = 20, clinical_vars = NULL, n_breaks = 5, colors = NULL) {

  p_main <- oncoplot_main(maf, genes)
  p_top <- aplotExtra:::oncoplot_sample(maf, genes)
  p_right <- aplotExtra:::oncoplot_gene(maf, genes, ylab = 'percentage')
  p_spacer <- ggplot2::ggplot() + ggfun::theme_transparent()

  # Base assembly using aplot
  pp <- p_main |>
    aplot::insert_top(p_spacer, height = 0.02) |>
    aplot::insert_top(p_top, height = 0.2) |>
    aplot::insert_right(p_right, width = 0.2)

  tracks <- oncoplot_clinical_track(maf, genes = genes, clinical_vars = clinical_vars,
                                    n_breaks = n_breaks, colors = colors)

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

  ggplot2::ggplot(d, ggplot2::aes(x = .data$Sample, y = .data$Gene, fill = .data$Type)) +
    ggplot2::geom_tile(colour = "white", linewidth = .01) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      axis.text.y.left = ggplot2::element_text(face = 'italic'),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 5, r = 2, b = 0, l = 5)
    ) +
    ggplot2::labs(fill = "Mutation Type")
}

#' @importFrom maftools getClinicalData
#' @importFrom rlang .data
oncoplot_clinical_track <- function(maf, genes = 20, clinical_vars, n_breaks, colors = NULL) {

  if (is.null(clinical_vars)) return(list())

  clinical_data <- as.data.frame(maftools::getClinicalData(maf))
  var_names <- intersect(
    if(is.list(clinical_vars)) names(clinical_vars) else clinical_vars,
    colnames(clinical_data)
  )

  onco_matrix <- aplotExtra:::oncoplot_tidy_onco_matrix(maf, genes)
  sample_order <- unique(as.character(onco_matrix$Sample))

  df <- clinical_data[, c("Tumor_Sample_Barcode", var_names), drop = FALSE]
  df$Tumor_Sample_Barcode <- factor(df$Tumor_Sample_Barcode, levels = sample_order)

  plot_list <- list()

  for (var in var_names) {
    tmp <- df[, c("Tumor_Sample_Barcode", var)]
    colnames(tmp)[2] <- "value"
    tmp$y_label <- var

    tmp$value <- if (is.numeric(tmp$value)) {
      binning_numeric(tmp$value, n_breaks = n_breaks)
    } else {
      as.factor(tmp$value)
    }

    # Used .data$ prefix here to prevent R CMD check "no visible binding" notes
    p <- ggplot2::ggplot(tmp, ggplot2::aes(x = .data$Tumor_Sample_Barcode, y = .data$y_label, fill = .data$value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.01) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "right",
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(face = 'italic', color = "black", size = 10),
        axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 5, r = 2, b = 0, l = 5)
      ) +
      ggplot2::labs(fill = var)

    if (!is.null(colors) && var %in% names(colors)) {
      p <- p + ggplot2::scale_fill_manual(values = colors[[var]], na.value = "grey90")
    }
    plot_list[[var]] <- p
  }

  return(plot_list)
}


#' @param x A numeric vector
#' @importFrom stats quantile
#' @noRd
binning_numeric <- function(x, n_breaks = 5) {
  if (!is.numeric(x)) return(as.factor(x))

  probs <- seq(0, 1, length.out = n_breaks + 1)
  brks <- stats::quantile(x, probs = probs, na.rm = TRUE)
  brks <- unique(brks)

  binned <- cut(x, breaks = brks, include.lowest = TRUE)

  lvls <- levels(binned)
  lvls <- gsub("\\.0", "", lvls)
  lvls <- gsub("\\(|\\]|\\]", "", lvls)
  lvls <- gsub(",", "–", lvls)
  levels(binned) <- lvls

  return(binned)
}
