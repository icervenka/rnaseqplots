ma_plot <- function(diffexp_data,
                    label = SYMBOL,
                    x = baseMean,
                    y = log2FoldChange,
                    sig_threshold = 0.05,
                    log2fc_threshold = 0,
                    filter_sig_on = padj,
                    label_bottom_n = 5,
                    label_top_n = 5,
                    label_genes = NULL,
                    add_vhlines = TRUE,
                    xlab_label = "mean expression",
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

  data_fil <- diffexp_data %>%
    dplyr::select({{ label }}, {{ x }}, {{ y }}, {{ filter_sig_on }}) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(significant = dplyr::case_when(
      {{ filter_sig_on }} <= sig_threshold &
        {{ y }} > log2fc_threshold ~ "up",
      {{ filter_sig_on }} <= sig_threshold &
        {{ y }} <= (-1) * log2fc_threshold ~ "down",
      TRUE ~ "unchanged"
    )) %>%
    dplyr::mutate(significant = factor(significant,
      levels = c("unchanged", "down", "up")
    ))

  p <- data_fil %>%
    ggplot2::ggplot(ggplot2::aes(x = log10({{ x }}), y = {{ y }})) +
    ggplot2::geom_point(ggplot2::aes(color = significant), size = 1)

  if (is.null(label_genes)) {
    p <- p +
      ggrepel::geom_label_repel(
        data = data_fil %>%
          dplyr::filter(significant != "unchanged") %>%
          dplyr::top_n(-label_bottom_n, {{ y }}),
        ggplot2::aes(label = {{ label }}),
        min.segment.length = grid::unit(0, "lines")
      ) +
      ggrepel::geom_label_repel(
        data = data_fil %>%
          dplyr::filter(significant != "unchanged") %>%
          dplyr::top_n(label_top_n, {{ y }}),
        ggplot2::aes(label = {{ label }}),
        min.segment.length = grid::unit(0, "lines")
      )
  } else if (is.vector(label_genes, mode = "character")) {
    p <- p +
      ggrepel::geom_label_repel(
        data = data_fil %>%
          dplyr::filter(tolower({{ label }}) %in% tolower(label_genes)),
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
  return(p)
}
