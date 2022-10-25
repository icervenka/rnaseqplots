library(magrittr, include.only = "%>%")

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

plot_lfc_scatter <- function(diffexp_data_1,
                             diffexp_data_2,
                             data_colnames_1 = c(id = "ENSEMBL", log2fc = "log2FoldChange", pvalue = "padj")
                             data_colnames_2 = c(id = "ENSEMBL", log2fc = "log2FoldChange", pvalue = "padj")
                             pvalue_threshold,
                             colo_quadrants = TRUE) {

  same_lfc_colname = data_colnames_1[["lof2fc"]] == data_colnames_2[["lof2fc"]]
  same_pvalue_colname = data_colnames_1[["pvalue"]] == data_colnames_2[["pvalue"]]

  common_genes <- intersect(
    diffexp_data_1 %>%
      dplyr::filter(matches(data_colnames_1[["id"]]) < pvalue_threshold) %>%
      dplyr::pull(data_colnames_1[["id"]]),
    diffexp_data_2 %>%
      dplyr::filter(matches(data_colnames_2[["id"]]) < pvalue_threshold) %>%
      dplyr::pull(data_colnames_2[["id"]])
  )

  lfc_data <- diffexp_data_1 %>%
    dplyr::select(matches(data_colnames_1)) %>%
    dplyr::filter(matches(data_colnames_1[["id"]]) %in% common_genes) %>%
    dplyr::left_join(diffexp_data_2 %>%
     dplyr::select(matches(data_colnames_2)), 
    by = c(data_colnames_1[["id"]] = data_colnames_2[["id"]]), 
    suffix = c("_x", "_y"))

  if (same_lfc_colname) {
    data_colnames_1[["log2fc"]] = c(data_colnames_1[["log2fc"]], "_x")
    data_colnames_2[["log2fc"]] = c(data_colnames_2[["log2fc"]], "_y")
  }

  if (same_pvalue_colname) {
    data_colnames_1[["pvalue"]] = c(data_colnames_1[["pvalue"]], "_x")
    data_colnames_2[["pvalue"]] = c(data_colnames_2[["pvalue"]], "_y")
  }

dplyr::mutate(chi_pcomb = -2*(log(q_value.x)+log(q_value.y))) %>% # Fisher method of combining p-values
    dplyr::mutate(p_chi = pchisq(chi_pcomb, df=4, lower.tail=FALSE)) %>% # calculate significance of combined p-values
    dplyr::mutate(err_sq = (log2_fold_change.y - log2_fold_change.x)^2) # calculate squared error of genes for ranking
  

  lfc_data = lfc_data %>%
    dplyr::mutate(chi_pcomb = -2*(log(matches(data_colnames_1[["pvalue"]])) +
    log(matches(data_colnames_2[["pvalue"]])))) %>% # Fisher method of combining p-values
    dplyr::mutate(p_chi = pchisq(chi_pcomb, df=4, lower.tail=FALSE, log.p = TRUE)) %>% # calculate significance of combined p-values
    dplyr::mutate(err_sq = (matches(data_colnames_1[["log2fc"]]) - 
    matches(data_colnames_2[["log2fc"]]))^2) # calculate squared error of genes for ranking
 

  p <- lfc_data %>%
    dplyr::filter(p_chi < pvalue_threshold) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = !!as.symbol(data_colnames_1[["log2fc"]]),
      y = !!as.symbol(data_colnames_2[["log2fc"]]),
      text = gene
    )) +
    ggplot2::geom_point(ggplot2::aes(size = -log10(p_chi), color = err_sq)) +
    colorspace::scale_colour_continuous_sequential("viridis", rev = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::labs(size = "-log10(p-value)", colour = "squared\nerror")
  return(p)
}
