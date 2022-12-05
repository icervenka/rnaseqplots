library(magrittr, include.only = "%>%")

#' Create an annotated volcano plot
#'
#' @param diffexp_data data frame of gene differential expression, must contain
#' column for gene IDs/symbols, log2 fold-changes and a p-value
#' @param x column name in diffexp_data to be plotted on the x axis, supplied
#' as variable. default: log2FoldChange
#' @param y column name in diffexp_data to be plotted on the y axis, supplied
#' as variable. default: pvalue
#' @param id_colname column name of gene IDs/symbols in differential expression
#' data file, supplied as variable. default: SYMBOL
#' @param filter_sig_on usually column name of p-values in differential
#' expression data file values of which will determine significance,
#' supplied as variable. default: padj
#' @param sig_threshold double, significance threshold for filter_sig_on column
#' @param log2fc_threshold
#' @param label_bottom_n integer, number of labels to display for genes with
#' highest y value. default: 5
#' @param label_top_n integer, number of labels to display for genes with
#' lowest y value. default: 5
#' @param gene_list character vector of gene IDs/symbols to display, overrides
#' label_top_n and label_bottom_n settings. default: NULL
#' @param add_vhlines logical, wheteher to add vertical and horizontal lines
#' indicating significance levels for both axes. default: TRUE
#' @param vhline_color character string, color of vertical and horizontal lines.
#' default: grey40
#' @param vhline_type character string, line type of vertical and horizontal
#' lines. Lines types for ggplot geom_line are accepted. default: dashed
#' @param xlab_label character string, label for x axis. default: log2FC
#' @param ylab_label character string, label for x axis. default: -log10FDR
#' @param color_palette character vector of color palette, length has to be
#' equal to 3 for unchaged, downregulated and upregulated genes. default:
#' c("gray70", "steelblue","darkred")
#'
#' @return ggplot volcano plot
#' @export
#'
#' @examples
volcano_plot <- function(diffexp_data,
                         x = log2FoldChange,
                         y = pvalue,
                         id_colname = SYMBOL,
                         filter_sig_on = padj,
                         sig_threshold = 0.05,
                         log2fc_threshold = 0.585, # rougly 1.5x
                         label_bottom_n = 5,
                         label_top_n = 5,
                         gene_list = NULL,
                         add_vhlines = TRUE,
                         vhline_color = "gray40",
                         vhline_type = "dashed",
                         xlab_label = expression(paste("lo", g[2], "FC")),
                         ylab_label = expression(paste("-lo", g[10], "FDR")),
                         color_palette = c(
                           "gray70",
                           "steelblue",
                           "darkred"
                         )) {

  # validate function arguments
  if (length(color_palette) != 3) {
    stop("Color palette requires a vector length 3.")
  }

  # Get data based on user supplied parameters
  df <- diffexp_data %>%
    dplyr::select({{ id_colname }}, {{ x }}, {{ y }}, {{ filter_sig_on }}) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(significant = dplyr::case_when(
      {{ filter_sig_on }} <= sig_threshold &
        {{ x }} > log2fc_threshold ~ "up",
      {{ filter_sig_on }} <= sig_threshold &
        {{ x }} <= (-1) * log2fc_threshold ~ "down",
      TRUE ~ "unchanged"
    )) %>%
    dplyr::mutate(significant = factor(significant,
      levels = c("unchanged", "down", "up")
    ))

  # Create basic plot
  p <- df %>%
    ggplot2::ggplot(ggplot2::aes(x = {{ x }}, y = -log10({{ y }}))) +
    ggplot2::geom_point(ggplot2::aes(color = significant), size = 1)

  # Add vertical and horizontal lines
  if (add_vhlines) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = c(-log2fc_threshold, log2fc_threshold),
        linetype = vhline_type,
        color = vhline_color
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(sig_threshold),
        linetype = vhline_type,
        color = vhline_color
      )
  }

  # Add gene labels for geom_points depending on what type of list is
  # provided by the user
  # TODO clean up the code, lots of repetition
  if (is.null(gene_list)) {
    p <- p +
      ggrepel::geom_label_repel(
        data = df %>%
          dplyr::filter(significant != "unchanged") %>%
          dplyr::top_n(-label_bottom_n, {{ x }}),
        ggplot2::aes(label = {{ id_colname }}),
        min.segment.length = grid::unit(0, "lines")
      ) +
      ggrepel::geom_label_repel(
        data = df %>%
          dplyr::filter(significant != "unchanged") %>%
          dplyr::top_n(label_top_n, {{ x }}),
        ggplot2::aes(label = {{ id_colname }}),
        min.segment.length = grid::unit(0, "lines")
      )
  } else if (is.vector(gene_list, mode = "character")) {
    p <- p +
      ggrepel::geom_label_repel(
        data = df %>%
          dplyr::filter(tolower({{ id_colname }}) %in% tolower(gene_list)),
        ggplot2::aes(label = {{ id_colname }}),
        min.segment.length = grid::unit(0, "lines")
      )
  }

  # Add final styling
  p <- p +
    ggplot2::theme_bw() +
    ggplot2::xlab(xlab_label) +
    ggplot2::ylab(ylab_label) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_color_manual(values = color_palette)
  return(p)
}

#' Create an annotated volcano plot from cuffdiff data
#'
#' @param diffexp_data data frame of gene differential expression, must contain
#' column for gene IDs/symbols, log2 fold-changes and a p-value
#' @param group1 character string, name of experimental group in the data.
#' @param group2 character string, name of group in the data that will be used
#' as a baseline.
#' @param ... other parameters to plot_volcano function
#'
#' @return ggplot volcano plot
#' @export
#'
#' @examples
volcano_plot_cuffdiff <- function(diffexp_data, group1, group2, ...) {
  df <- diffexp_data %>%
    dplyr::filter(sample_1 == group1 & sample_2 == group2)
  return(volcano_plot(df, ...))
}
