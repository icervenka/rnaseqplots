plot_corr <- function(gene_expr, param_values) {
  p <- base::cbind.data.frame(
    gene_expr = gene_expr,
    param = param_values,
    color_by = param_values
  ) %>%
    ggplot2::ggplot(ggplot2::aes(x = gene_expr, y = param, color = color_by)) +
    ggplot2::geom_point()
  return(p)
}

plot_lfc_scatter <- function(lfc_data) {
  p <- lfc_data %>%
    dplyr::filter(p_chi < 0.1) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = log2_fold_change.x,
      y = log2_fold_change.y,
      text = gene
    )) +
    ggplot2::geom_point(ggplot2::aes(size = -log10(p_chi), color = err_sq)) +
    colorspace::scale_colour_continuous_sequential("viridis", rev = F) +
    ggplot2::theme_bw() +
    ggplot2::labs(size = "-log10(p-value)", colour = "squared\nerror")
  return(p)
}
