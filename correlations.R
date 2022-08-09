plot_corr = function(gene_expr, param_values) {
 p = cbind.data.frame(expr = gene_expr,
                       param = param_values,
                       color_by = param_values) %>%
    ggplot2::ggplot(aes(x = expr, y = param, color = color_by)) +
    geom_point()
 return(p)
}

plot_lfc_scatter = function(lfc_data) {
  p = lfc_data %>%
    dplyr::filter(p_chi < 0.1) %>%
    ggplot2::ggplot(aes(x = log2_fold_change.x,
                        y = log2_fold_change.y,
                        text = gene)) +
    geom_point(aes(size = -log10(p_chi), color = err_sq)) +
    colorspace::scale_colour_continuous_sequential("viridis", rev = F) +
    theme_bw() +
    labs(size="-log10(p-value)", colour="squared\nerror")
  return(p)
}
