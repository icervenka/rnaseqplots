#' Draws MA-plot based on differential expression data
#'
#' Adds labels for genes of interest and colors based on upregulated and
#' downregulated genes.
#'
#' @param diffexp_data data frame of gene differential expression, must contain
#' column for gene IDs/symbols, log2 fold-changes and a p-value
#' @param x column name in diffexp_data to be plotted on the x axis, supplied
#' as variable. For MA plots, log2(x) is plotted. default: baseMean
#' @param y column name in diffexp_data to be plotted on the y axis, supplied
#' as variable. default: log2FoldChange
#' @param id_colname column name of gene IDs/symbols in differential expression
#' data file, supplied as variable
#' @param filter_sig_on usually column name of p-values in differential
#' expression data file values of which will determine significance,
#' supplied as variable. default: padj
#' @param sig_threshold double, significance threshold for filter_sig_on column
#' @param x_threshold double, significance threshold for x column. default: 0
#' @param y_threshold double, significance threshold for y column. default: 0
#' @param label_bottom_n integer, number of labels to display for genes with
#' highest y value. default: 5
#' @param label_top_n integer, number of labels to display for genes with
#' lowest y value. default: 5
#' @param gene_list character vector of gene IDs/symbols to display, overrides
#' label_top_n and label_bottom_n settings. default: NULL
#' @param xlab_label character string, label for x axis. default: log2Expression
#' @param ylab_label character string, label for x axis. default: log2FC
#' @param color_palette character vector of color palette, length has to be
#' equal to 3 for unchaged, downregulated and upregulated genes. default:
#' c("gray70", "steelblue","darkred")
#'
#' @return MA plot with colored significat genes and optionally gene labels
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
plot_ma <- function(diffexp_data,
                    x = baseMean,
                    y = log2FoldChange,
                    id_colname = SYMBOL,
                    filter_sig_on = padj,
                    sig_threshold = 0.05,
                    x_threshold = 0,
                    y_threshold = 0,
                    label_bottom_n = 5,
                    label_top_n = 5,
                    gene_list = NULL,
                    xlab_label = expression(paste("lo", g[2], "Expression")),
                    ylab_label = expression(paste("lo", g[2], "FC")),
                    color_palette = c(
                      "gray70",
                      "steelblue",
                      "darkred"
                    )) {

  # validate function arguments
  if (length(color_palette) != 3) {
    stop("Color palette requires a vector length 3.")
  }

  # filter significant (downregulate, upregulated genes and assign them
  # to groups)
  data_fil <- diffexp_data %>%
    dplyr::select({{ id_colname }}, {{ x }}, {{ y }}, {{ filter_sig_on }}) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(significant = dplyr::case_when(
      {{ filter_sig_on }} <= sig_threshold &
        {{ y }} > y_threshold &
        {{ x }} > x_threshold ~ "up",
      {{ filter_sig_on }} <= sig_threshold &
        {{ y }} <= (-1) * y_threshold &
        {{ x }} > x_threshold ~ "down",
      TRUE ~ "unchanged"
    )) %>%
    dplyr::mutate(significant = factor(significant,
      levels = c("unchanged", "down", "up")
    ))

  # create basic plot
  p <- data_fil %>%
    ggplot2::ggplot(ggplot2::aes(x = log10({{ x }}), y = {{ y }})) +
    ggplot2::geom_point(ggplot2::aes(color = significant), size = 1)

  # add gene list labels based on gene_list vector or top/bottom genes
  if (is.null(gene_list)) {
    p <- p +
      ggrepel::geom_label_repel(
        data = data_fil %>%
          dplyr::filter(significant != "unchanged") %>%
          dplyr::top_n(-label_bottom_n, {{ y }}),
        ggplot2::aes(label = {{ id_colname }}),
        min.segment.length = grid::unit(0, "lines")
      ) +
      ggrepel::geom_label_repel(
        data = data_fil %>%
          dplyr::filter(significant != "unchanged") %>%
          dplyr::top_n(label_top_n, {{ y }}),
        ggplot2::aes(label = {{ id_colname }}),
        min.segment.length = grid::unit(0, "lines")
      )
  } else if (is.vector(gene_list, mode = "character")) {
    p <- p +
      ggrepel::geom_label_repel(
        data = data_fil %>%
          dplyr::filter(tolower({{ id_colname }}) %in% tolower(gene_list)),
        ggplot2::aes(label = {{ id_colname }}),
        min.segment.length = grid::unit(0, "lines")
      )
  }

  # style the plot
  p <- p +
    ggplot2::theme_bw() +
    ggplot2::xlab(xlab_label) +
    ggplot2::ylab(ylab_label) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_color_manual(values = color_palette)
  return(p)
}
