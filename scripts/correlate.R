library(magrittr, include.only = "%>%")

cor_test_df <- function(df, x, y) {
  cor.test(
    df[[deparse(substitute(x))]],
    df[[deparse(substitute(y))]],
    method = "pearson"
  ) %>%
    tidy()
}

plot_param_corr <- function(expression_data,
                            metadata,
                            gene_list,
                            params,
                            id_colname = SYMBOL,
                            sample_colname = sample,
                            show_regression = TRUE) {
  sample_colname_str <- deparse(substitute(sample_colname))

  filtered_expression <- expression_data %>%
    dplyr::filter({{ id_colname }} %in% gene_list) %>%
    dplyr::select({{ id_colname }}, matches(metadata[[sample_colname_str]])) %>%
    tidyr::pivot_longer(-{{ id_colname }},
      names_to = sample_colname_str,
      values_to = "gene_expression"
    )

  filtered_metadata <- metadata %>%
    dplyr::select({{ sample_colname }}, matches(params)) %>%
    tidyr::pivot_longer(-{{ sample_colname }},
      names_to = "param",
      values_to = "param_value"
    )

  param_df <- filtered_expression %>%
    dplyr::left_join(filtered_metadata, by = sample_colname_str)

  # calculates pearson correlation, results are currently unused
  # kept in case of future features
  cor_df <- param_df %>%
    group_by({{ id_colname }}, param) %>%
    do(cor_test_df((.), gene_expression, param_value))

  p <- param_df %>%
    ggplot2::ggplot(ggplot2::aes(x = param_value, y = gene_expression)) +
    ggplot2::geom_point() +
    ggplot2::stat_smooth(
      method = "lm",
      se = FALSE,
      size = 0.5,
      color = "steelblue4"
    ) +
    ggplot2::facet_grid(SYMBOL ~ param, scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::xlab("parameters") +
    ggplot2::ylab("normalized gene expression ")
  return(p)
}

# TODO color quadrants
plot_lfc_scatter <- function(diffexp_data_1,
                             diffexp_data_2,
                             data_colnames_1 = c(
                               id = "ENSEMBL",
                               log2fc = "log2FoldChange",
                               pvalue = "padj"
                             ),
                             data_colnames_2 = c(
                               id = "ENSEMBL",
                               log2fc = "log2FoldChange",
                               pvalue = "padj"
                             ),
                             pvalue_threshold = 0.05,
                             color_quadrants = TRUE) {
  same_lfc_colname <- data_colnames_1["log2fc"] == data_colnames_2["log2fc"]
  same_pvalue_colname <- data_colnames_1["pvalue"] == data_colnames_2["pvalue"]

  common_genes <- intersect(
    diffexp_data_1 %>%
      dplyr::filter(!!as.symbol(data_colnames_1["pvalue"]) < pvalue_threshold) %>%
      dplyr::pull(data_colnames_1["id"]),
    diffexp_data_2 %>%
      dplyr::filter(!!as.symbol(data_colnames_2["pvalue"]) < pvalue_threshold) %>%
      dplyr::pull(data_colnames_2["id"])
  )
  
  by_x = data_colnames_1["id"] %>% unname()
  by_y = data_colnames_2["id"] %>% unname()
  
  print(by_x)
  print(by_y)
  
  lfc_data <- diffexp_data_1 %>%
    dplyr::select(dplyr::all_of(data_colnames_1 %>% unname)) %>%
    dplyr::filter(!!as.symbol(data_colnames_1["id"]) %in% common_genes) %>%
    dplyr::left_join(diffexp_data_2 %>%
                       dplyr::select(dplyr::all_of(data_colnames_2 %>% unname)),
                     by = "ENSEMBL",
                     suffix = c("_x", "_y")
                     )
  
  if (same_lfc_colname) {
    data_colnames_1["log2fc"] <- paste0(data_colnames_1["log2fc"], "_x")
    data_colnames_2["log2fc"] <- paste0(data_colnames_2["log2fc"], "_y")
  }

  if (same_pvalue_colname) {
    data_colnames_1["pvalue"] <- paste0(data_colnames_1["pvalue"], "_x")
    data_colnames_2["pvalue"] <- paste0(data_colnames_2["pvalue"], "_y")
  }

  lfc_data <- lfc_data %>%
    # Fisher method of combining p-values
    dplyr::mutate(chi_pcomb = -2 * (log(!!as.symbol(data_colnames_1["pvalue"])) +
      log(!!as.symbol(data_colnames_2["pvalue"])))) %>%
    # calculate significance of combined p-values
    # TODO check the degrees of freedom
    dplyr::mutate(p_chi = pchisq(chi_pcomb,
      df = 2,
      lower.tail = FALSE
    )) %>%
    # calculate squared error of genes for ranking
    dplyr::mutate(err_sq = (!!as.symbol(data_colnames_1["log2fc"]) -
                              !!as.symbol(data_colnames_2["log2fc"]))^2)

  print(lfc_data)
  p <- lfc_data %>%
    dplyr::filter(p_chi < pvalue_threshold) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = !!as.symbol(data_colnames_1["log2fc"]),
      y = !!as.symbol(data_colnames_2["log2fc"]),
      text = !!as.symbol(data_colnames_1["id"])
    )) +
    ggplot2::geom_point(ggplot2::aes(size = -log10(p_chi), color = err_sq)) +
    colorspace::scale_colour_continuous_sequential("viridis", rev = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::labs(size = "-log10(p-value)", colour = "squared\nerror")
  return(p)
}
