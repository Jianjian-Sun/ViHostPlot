#' Generate a oncoplot with clinical tracks
#'
#' @param maf MAF object.
#' @param genes the gene names or the number, default is 20.
#' @param clinical_vars A character vector or list defining clinical variables to display as bottom tracks.
#' @param n_breaks Integer. Number of breaks for continuous clinical variables. Default is 5.
#' @param colors A named list of customized colors for clinical tracks.
#'   If omitted for a clinical variable, a scientific fallback palette is used:
#'   Okabe-Ito colors for discrete variables and Viridis-style colors for
#'   binned continuous variables.
#' @param sort_by Character. Clinical variable used to sort samples left-to-right.
#'   When \code{NULL}, sample order follows the oncoplot matrix default.
#' @param clinical_order Named list of character vectors giving custom level order
#'   per clinical variable. Used when \code{sort_by} matches a list name; if omitted
#'   for the sorting variable, samples are ordered by increasing values (numeric)
#'   or alphabetically (categorical).
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
#' oncoplot(maf = laml, genes = 20, clinical_vars = var_names,
#'           sort_by = "FAB_classification")
#' oncoplot(maf = laml, genes = 20, clinical_vars = var_names,
#'           sort_by = "FAB_classification",
#'           clinical_order = list(FAB_classification = c("M0", "M1", "M2")))
#' }
utils::globalVariables(".data")

oncoplot <- function(maf, genes = 20, clinical_vars = NULL, n_breaks = 5,
                     colors = NULL, sort_by = NULL,
                     clinical_order = NULL) {

  sample_order <- oncoplot_resolve_sample_order(
    maf, genes, sort_by = sort_by, clinical_order = clinical_order
  )

  p_main <- oncoplot_main(maf, genes, sample_order = sample_order)
  p_top <- oncoplot_apply_sample_order(
    aplotExtra:::oncoplot_sample(maf, genes), sample_order
  )
  p_right <- aplotExtra:::oncoplot_gene(maf, genes, ylab = 'percentage')
  p_spacer <- ggplot2::ggplot() + ggfun::theme_transparent()

  # Base assembly using aplot
  pp <- p_main |>
    aplot::insert_top(p_spacer, height = 0.02) |>
    aplot::insert_top(p_top, height = 0.2) |>
    aplot::insert_right(p_right, width = 0.2)

  tracks <- oncoplot_clinical_track(maf, genes = genes, clinical_vars = clinical_vars,
                                    n_breaks = n_breaks, colors = colors,
                                    sample_order = sample_order)

  if (length(tracks) > 0) {
    for (var in names(tracks)) {
      pp <- aplot::insert_bottom(pp, tracks[[var]], height = 0.05)
    }
  }

  return(pp)
}


#' @importFrom ggplot2 geom_tile
#' @importFrom rlang .data
oncoplot_main <- function(maf, genes = 20, sample_order = NULL) {
  d <- aplotExtra:::oncoplot_tidy_onco_matrix(maf, genes)
  if (!is.null(sample_order)) {
    d$Sample <- factor(as.character(d$Sample), levels = sample_order)
  }

  .data <- rlang::.data
  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = .data$Sample, y = .data$Gene, fill = .data$Type)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = .01) +
    oncoplot_setting(continuous = FALSE, fill_name = "Mutation Type") +
    ggplot2::theme(
        legend.position = "right",
        axis.text.y.left = ggplot2::element_text(face = 'italic')
    )

  p
}

#' @importFrom maftools getClinicalData
#' @importFrom rlang .data
#' @importFrom dplyr select mutate transmute all_of
oncoplot_clinical_track <- function(maf, genes = 20, clinical_vars, n_breaks,
                                    colors = NULL, sample_order = NULL) {

  if (is.null(clinical_vars)) return(list())

  clinical_data <- as.data.frame(maftools::getClinicalData(maf))
  var_names <- intersect(
    if(is.list(clinical_vars)) names(clinical_vars) else clinical_vars,
    colnames(clinical_data)
  )

  if (is.null(sample_order)) {
    onco_matrix <- aplotExtra:::oncoplot_tidy_onco_matrix(maf, genes)
    sample_order <- unique(as.character(onco_matrix$Sample))
  }

  df <- clinical_data[, c("Tumor_Sample_Barcode", var_names), drop = FALSE]
  df$Tumor_Sample_Barcode <- factor(
    df$Tumor_Sample_Barcode,
    levels = sample_order
  )

  plot_list <- list()

  for (var in var_names) {
    tmp <- data.frame(
      Tumor_Sample_Barcode = df$Tumor_Sample_Barcode,
      value = df[[var]],
      y_label = var,
      stringsAsFactors = FALSE
    )
    is_continuous <- is.numeric(tmp$value)

    tmp$value <- if (is_continuous) {
      binning_numeric(tmp$value, n_breaks = n_breaks)
    } else {
      as.factor(tmp$value)
    }

    .data <- rlang::.data
    p <- ggplot2::ggplot(
      tmp,
      ggplot2::aes(
        x = .data$Tumor_Sample_Barcode,
        y = .data$y_label,
        fill = .data$value
      )
    ) +
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
    } else {
      lvls <- levels(as.factor(tmp$value))
      default_colors <- if (is_continuous) {
        grDevices::hcl.colors(length(lvls), palette = "Viridis")
      } else {
        okabe_ito <- c(
          "#E69F00", "#56B4E9", "#009E73", "#F0E442",
          "#0072B2", "#D55E00", "#CC79A7", "#999999"
        )
        if (length(lvls) <= length(okabe_ito)) {
          okabe_ito[seq_along(lvls)]
        } else {
          grDevices::hcl.colors(length(lvls), palette = "Dark 3")
        }
      }
      p <- p + ggplot2::scale_fill_manual(
        values = stats::setNames(default_colors, lvls),
        na.value = "grey90"
      )
    }
    plot_list[[var]] <- p
  }

  return(plot_list)
}

#' Resolve sample order for oncoplot panels
#' @noRd
oncoplot_resolve_sample_order <- function(maf, genes, sort_by = NULL,
                                          clinical_order = NULL) {
  onco_matrix <- aplotExtra:::oncoplot_tidy_onco_matrix(maf, genes)
  default_order <- unique(as.character(onco_matrix$Sample))

  if (is.null(sort_by)) {
    return(default_order)
  }

  clinical_data <- as.data.frame(maftools::getClinicalData(maf))
  if (!sort_by %in% colnames(clinical_data)) {
    stop(
      "sort_by variable '", sort_by,
      "' was not found in clinical data.",
      call. = FALSE
    )
  }

  custom_order <- NULL
  if (!is.null(clinical_order) && sort_by %in% names(clinical_order)) {
    custom_order <- clinical_order[[sort_by]]
  }

  sorted_ids <- oncoplot_sort_samples(
    clinical_data[, c("Tumor_Sample_Barcode", sort_by), drop = FALSE],
    var = sort_by,
    custom_order = custom_order
  )

  sorted_ids <- sorted_ids[sorted_ids %in% default_order]
  c(sorted_ids, setdiff(default_order, sorted_ids))
}

#' Sort sample barcodes by a clinical variable
#' @noRd
oncoplot_sort_samples <- function(clinical_data, var, custom_order = NULL) {
  ids <- as.character(clinical_data$Tumor_Sample_Barcode)
  vals <- clinical_data[[var]]

  if (!is.null(custom_order)) {
    vals_ord <- factor(vals, levels = custom_order)
    ord <- order(vals_ord, na.last = TRUE)
  } else if (is.numeric(vals)) {
    ord <- order(vals, na.last = TRUE)
  } else {
    lvls <- sort(unique(as.character(vals[!is.na(vals)])))
    vals_ord <- factor(as.character(vals), levels = lvls)
    ord <- order(vals_ord, na.last = TRUE)
  }

  ids[ord]
}

#' Apply sample order to a ggplot object with a Sample aesthetic
#' @noRd
oncoplot_apply_sample_order <- function(p, sample_order) {
  if (is.null(sample_order) || is.null(p) || is.null(p$data) ||
      !"Sample" %in% names(p$data)) {
    return(p)
  }

  p$data$Sample <- factor(as.character(p$data$Sample), levels = sample_order)
  p
}

oncoplot_setting <- function(noxaxis = TRUE, continuous = TRUE, scale = 'y',
                             fill_name = NULL) {  
    list(
        ggplot2::theme_minimal(),
        if (noxaxis) ggfun::theme_noxaxis(),
        ggplot2::theme(
          legend.position = "none",
          panel.grid.major = ggplot2::element_blank()
        ),
        aplotExtra:::oncoplot_scale(continuous = continuous, scale = scale),
        oncoplot_fill(name = fill_name),
        ggplot2::xlab(NULL),
        ggplot2::ylab(NULL)
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

    ggplot2::scale_fill_manual(
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
