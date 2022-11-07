library(magrittr, include.only = "%>%")

plot_param_corr <- function(expression_data,
                            metadata,
                            gene_list,
                            params,
                            id_colname = SYMBOL,
                            group_colname = group,
                            sample_colname = sample,
                            palette = NULL,
                            show_regression = TRUE) {
  sample_colname_str <- deparse(substitute(sample_colname))

  filtered_expression <- expression_data %>%
    dplyr::filter(tolower({{ id_colname }}) %in% tolower(gene_list)) %>%
    dplyr::select({{ id_colname }}, matches(metadata[[sample_colname_str]])) %>%
    tidyr::pivot_longer(-{{ id_colname }},
      names_to = sample_colname_str,
      values_to = "gene_expression"
    )

  filtered_metadata <- metadata %>%
    dplyr::select({{ sample_colname }}, {{ group_colname }}, matches(params)) %>%
    tidyr::pivot_longer(-c({{ sample_colname }}, {{ group_colname }}),
      names_to = "param",
      values_to = "param_value"
    )

  print(filtered_metadata)
  param_df <- filtered_expression %>%
    dplyr::left_join(filtered_metadata, by = sample_colname_str)

  # calculates pearson correlation, results are currently unused
  # kept in case of future features
  cor_df <- param_df %>%
    group_by({{ id_colname }}, param) %>%
    do(cor_test_df((.), gene_expression, param_value))

  p <- param_df %>%
    ggplot2::ggplot(ggplot2::aes(x = param_value, y = gene_expression)) +
    ggplot2::geom_point(aes(color = {{ group_colname }})) +
    { if (!is.null(palette)) scale_color_manual(values = palette) } + 
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
                             color_quadrants = TRUE,
                             alpha = 1) {
  same_lfc_colname <- data_colnames_1["log2fc"] == data_colnames_2["log2fc"]
  same_pvalue_colname <- data_colnames_1["pvalue"] == data_colnames_2["pvalue"]

  common_genes <- intersect(
    diffexp_data_1 %>%
      dplyr::filter(
        !!as.symbol(data_colnames_1["pvalue"]) < pvalue_threshold
      ) %>%
      dplyr::pull(data_colnames_1["id"]),
    diffexp_data_2 %>%
      dplyr::filter(
        !!as.symbol(data_colnames_2["pvalue"]) < pvalue_threshold
      ) %>%
      dplyr::pull(data_colnames_2["id"])
  )

  lfc_data <- diffexp_data_1 %>%
    dplyr::select(dplyr::all_of(unname(data_colnames_1))) %>%
    dplyr::filter(!!as.symbol(data_colnames_1["id"]) %in% common_genes) %>%
    dplyr::left_join(diffexp_data_2 %>%
      dplyr::select(dplyr::all_of(unname(data_colnames_2))),
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
    dplyr::mutate(
      chi_pcomb = -2 * (log(!!as.symbol(data_colnames_1["pvalue"])) +
        log(!!as.symbol(data_colnames_2["pvalue"])))
    ) %>%
    # calculate significance of combined p-values
    # TODO check the degrees of freedom
    dplyr::mutate(p_chi = pchisq(chi_pcomb,
      df = 1,
      lower.tail = FALSE
    )) %>%
    # calculate squared error of genes for ranking
    dplyr::mutate(err_sq = (!!as.symbol(data_colnames_1["log2fc"]) -
      !!as.symbol(data_colnames_2["log2fc"]))^2)

  lfc_data <- lfc_data %>%
    dplyr::mutate(
      sign_1 = sign(!!as.symbol(data_colnames_1["log2fc"])),
      sign_2 = sign(!!as.symbol(data_colnames_2["log2fc"]))
    ) %>%
    dplyr::mutate(quadrant = case_when(
      sign_1 > 0 & sign_2 > 0 ~ 1,
      sign_1 > 0 & sign_2 < 0 ~ 2,
      sign_1 < 0 & sign_2 < 0 ~ 3,
      sign_1 < 0 & sign_2 > 0 ~ 4
    )) %>%
    dplyr::mutate(quadrant = factor(quadrant))

  p <- lfc_data %>%
    dplyr::filter(p_chi < pvalue_threshold) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = !!as.symbol(data_colnames_1["log2fc"]),
      y = !!as.symbol(data_colnames_2["log2fc"]),
      text = !!as.symbol(data_colnames_1["id"])
    ))

  if (color_quadrants) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(size = -log10(p_chi), color = quadrant),
                          alpha = alpha) + 
      ggplot2::labs(colour = "quadrant")
  } else {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(size = -log10(p_chi), color = err_sq),
                          alpha = alpha) +
      colorspace::scale_colour_continuous_sequential("viridis", rev = FALSE) + 
      ggplot2::labs(colour = "squared\nerror")
  }

  p <- p +
    ggplot2::theme_bw() +
    ggplot2::labs(size = "-log10(p-value)")
  return(p)
}
