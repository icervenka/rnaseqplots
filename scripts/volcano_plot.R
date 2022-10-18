library(magrittr, include.only = "%>%")

# label_genes if specified overrides top_n/bottom_n arguments,
volcano_plot <- function(results,
                         label = SYMBOL,
                         x = log2FoldChange,
                         y = pvalue,
                         sig_threshold = 0.05,
                         log2fc_threshold = 0.585, # rougly 1.5x
                         filter_sig_on = padj,
                         label_bottom_n = 5,
                         label_top_n = 5,
                         label_genes = NULL,
                         add_vhlines = TRUE,
                         vhline_color = "gray40",
                         vhline_type = "dashed",
                         xlab_label = expression(paste("lo", g[2], "FC")),
                         ylab_label = expression(paste("lo", g[10], "FDR")),
                         color_palette = c(
                           "gray70",
                           "steelblue",
                           "darkred"
                         )) {

  # validate function arguments
  if (length(color_palette) != 3) {
    stop("Color palette requires a vector length 3.")
  }

  results_fil <- results %>%
    dplyr::select({{ label }}, {{ x }}, {{ y }}, {{ filter_sig_on }}) %>%
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

  p <- results_fil %>%
    ggplot2::ggplot(ggplot2::aes(x = {{ x }}, y = -log10({{ y }}))) +
    ggplot2::geom_point(ggplot2::aes(color = significant), size = 1)

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

  if (is.null(label_genes)) {
    p <- p +
      ggrepel::geom_label_repel(
        data = results_fil %>%
          dplyr::filter(significant != "unchanged") %>%
          dplyr::top_n(-label_bottom_n, {{ x }}),
        ggplot2::aes(label = {{ label }}),
        min.segment.length = grid::unit(0, "lines")
      ) +
      ggrepel::geom_label_repel(
        data = results_fil %>%
          dplyr::filter(significant != "unchanged") %>%
          dplyr::top_n(label_top_n, {{ x }}),
        ggplot2::aes(label = {{ label }}),
        min.segment.length = grid::unit(0, "lines")
      )
  } else if (is.vector(label_genes, mode = "character")) {
    p <- p +
      ggrepel::geom_label_repel(
        data = results_fil %>%
          dplyr::filter({{ label }} %in% label_genes),
        ggplot2::aes(label = {{ label }}),
        min.segment.length = grid::unit(0, "lines")
      )
  }

  p <- p +
    ggplot2::theme_bw() +
    ggplot2::xlab(xlab_label) +
    ggplot2::ylab(ylab_label) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_color_manual(values = color_palette)
  p
}

volcano_plot_cuffdiff <- function(results, sample1, sample2, ...) {
  results <- results %>%
    dplyr::filter(sample_1 == sample1 & sample_2 == sample2)
  volcano_plot(results, ...)
}
