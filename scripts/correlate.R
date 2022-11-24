library(magrittr, include.only = "%>%")

#' Create a correlation graph of gene expression and sample parameter
#'
#' @param expression_data data frame of normalized gene expression with samples
#' as columns and an id column for gene names/symbols
#' @param metadata data frame with column of samples names, column for group
#' names and optionally other columns describing the sample features that will
#' be correlated with gene expression
#' @param gene_list character vector of gene IDs/symbols, which will be used
#' to  plot the correlations
#' @param params character vector stating which parameter columns to choose
#' from metadata to plot the  gene expression
#' @param id_colname column name of gene IDs/symbols in expression data file,
#' supplied as variable
#' @param sample_colname column name of samples in metadata file, supplied as
#' variable. default: sample
#' @param group_colname column name of experimental groups in metadata file,
#' supplied as variable. default: group
#' @param palette character vector of color palette, length has to be equal to
#' the amount of groups in metadata file
#' @param show_regression logical, whether to show regression line on the graph
#'
#' @return ggplot facet graph of correlation between specified sample parameters
#' and gene expression of selected genes
#' @export
#'
#' @importFrom dplyr filter select left_join group_by
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual stat_smooth
#' facet_grid theme_bw xlab ylab
#'
#' @examples
plot_param_corr <- function(expression_data,
                            metadata,
                            gene_list,
                            params,
                            id_colname = SYMBOL,
                            sample_colname = sample,
                            group_colname = group,
                            palette = NULL,
                            show_regression = TRUE) {
  sample_colname_str <- deparse(substitute(sample_colname))

  # filter genes supplied in an gene list and samples supplied in metadata and
  # transform to a long data frame
  filtered_expression <- expression_data %>%
    dplyr::filter(tolower({{ id_colname }}) %in% tolower(gene_list)) %>%
    dplyr::select({{ id_colname }}, matches(metadata[[sample_colname_str]])) %>%
    tidyr::pivot_longer(-{{ id_colname }},
      names_to = sample_colname_str,
      values_to = "gene_expression"
    )

  # transform metadata into long form
  filtered_metadata <- metadata %>%
    dplyr::select(
      {{ sample_colname }},
      {{ group_colname }},
      matches(params)
    ) %>%
    tidyr::pivot_longer(-c({{ sample_colname }}, {{ group_colname }}),
      names_to = "param",
      values_to = "param_value"
    )

  # join the gene expression and metadata together
  param_df <- filtered_expression %>%
    dplyr::left_join(filtered_metadata, by = sample_colname_str)

  # calculates pearson correlation, results are currently unused
  # kept in case of future features
  cor_df <- param_df %>%
    dplyr::group_by({{ id_colname }}, param) %>%
    do(cor_test_df((.), gene_expression, param_value))

  # create ggplot with custom palette and facet based on gene ~ param
  p <- param_df %>%
    ggplot2::ggplot(ggplot2::aes(x = param_value, y = gene_expression)) +
    ggplot2::geom_point(aes(color = {{ group_colname }})) +
    { # nolint
      if (!is.null(palette)) ggplot2::scale_color_manual(values = palette)
    } +
    ggplot2::stat_smooth(
      method = "lm",
      se = FALSE,
      size = 0.5,
      color = "steelblue4"
    ) +
    # TODO replace symbol with proper column name
    ggplot2::facet_grid(SYMBOL ~ param, scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::xlab("parameters") +
    ggplot2::ylab("normalized gene expression ")
  return(p)
}

#' Create scatter plot of log2 fold-changes from two data sets
#'
#' @param diffexp_data_1 data frame of differential expression of first data set
#' @param diffexp_data_2 data frame of differential expression of second data set
#' @param data_colnames_1 named character vector of length 3 with names 'id',
#' 'log2fc' and 'pvalue' corresponding to the names of these columns in first
#' data set
#' @param data_colnames_2 named character vector of length 3 with names 'id',
#' 'log2fc' and 'pvalue' corresponding to the names of these columns in second
#' data set
#' @param pvalue_threshold double, maximal p-value threshold of the gene to be
#' displayed on graph based on combined p-values from two data sets
#' @param color_quadrants logical value denoting whether points should be
#' colored based on graph quadrant, squared error is used otherwise
#' @param alpha double, transparency of points in the graph
#'
#' @return ggplot2 scatter plot of log2 fold changes in x and y axes and point
#' size based on combined p-values
#' @export
#'
#' @importFrom dplyr filter pull select left_join mutate case_when
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw
#' @importFrom colorspace scale_colour_continuous_sequential
#' @importFrom stats pchisq
#'
#' @examples
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
  # check if the column names for pval and log2fc are the same
  # this will be needed later as the names will get suffixes because of join
  same_lfc_colname <- data_colnames_1["log2fc"] == data_colnames_2["log2fc"]
  same_pvalue_colname <- data_colnames_1["pvalue"] == data_colnames_2["pvalue"]

  # extract common genes from the datasets
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

  # create a combined data frame with all the values necessary for plotting
  lfc_data <- diffexp_data_1 %>%
    dplyr::select(dplyr::all_of(unname(data_colnames_1))) %>%
    dplyr::filter(!!as.symbol(data_colnames_1["id"]) %in% common_genes) %>%
    dplyr::left_join(diffexp_data_2 %>%
      dplyr::select(dplyr::all_of(unname(data_colnames_2))),
    by = "ENSEMBL",
    suffix = c("_x", "_y")
    )

  # update the column names character vectors with join suffixes in case they
  # were the same
  if (same_lfc_colname) {
    data_colnames_1["log2fc"] <- paste0(data_colnames_1["log2fc"], "_x")
    data_colnames_2["log2fc"] <- paste0(data_colnames_2["log2fc"], "_y")
  }

  if (same_pvalue_colname) {
    data_colnames_1["pvalue"] <- paste0(data_colnames_1["pvalue"], "_x")
    data_colnames_2["pvalue"] <- paste0(data_colnames_2["pvalue"], "_y")
  }

  # update final data frame with combined p-values and squared errors
  lfc_data <- lfc_data %>%
    # Fisher method of combining p-values
    dplyr::mutate(
      chi_pcomb = -2 * (log(!!as.symbol(data_colnames_1["pvalue"])) +
        log(!!as.symbol(data_colnames_2["pvalue"])))
    ) %>%
    # calculate significance of combined p-values
    # TODO check the degrees of freedom
    dplyr::mutate(p_chi = stats::pchisq(chi_pcomb,
      df = 1,
      lower.tail = FALSE
    )) %>%
    # calculate squared error of genes for ranking
    dplyr::mutate(err_sq = (!!as.symbol(data_colnames_1["log2fc"]) -
      !!as.symbol(data_colnames_2["log2fc"]))^2)

  # in case quadrant coloring is specified, create the coloring variable
  lfc_data <- lfc_data %>%
    dplyr::mutate(
      sign_1 = sign(!!as.symbol(data_colnames_1["log2fc"])),
      sign_2 = sign(!!as.symbol(data_colnames_2["log2fc"]))
    ) %>%
    dplyr::mutate(quadrant = dplyr::case_when(
      sign_1 > 0 & sign_2 > 0 ~ 1,
      sign_1 > 0 & sign_2 < 0 ~ 2,
      sign_1 < 0 & sign_2 < 0 ~ 3,
      sign_1 < 0 & sign_2 > 0 ~ 4
    )) %>%
    dplyr::mutate(quadrant = factor(quadrant))

  # filter data based on p-value threshold and create basic graph
  p <- lfc_data %>%
    dplyr::filter(p_chi < pvalue_threshold) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = !!as.symbol(data_colnames_1["log2fc"]),
      y = !!as.symbol(data_colnames_2["log2fc"]),
      text = !!as.symbol(data_colnames_1["id"])
    ))

  # add point colors (quadrant or squared error)
  if (color_quadrants) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(size = -log10(p_chi), color = quadrant),
        alpha = alpha
      ) +
      ggplot2::labs(colour = "quadrant")
  } else {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(size = -log10(p_chi), color = err_sq),
        alpha = alpha
      ) +
      colorspace::scale_colour_continuous_sequential("viridis", rev = FALSE) +
      ggplot2::labs(colour = "squared\nerror")
  }

  # graph styling
  p <- p +
    ggplot2::theme_bw() +
    ggplot2::labs(size = "-log10(p-value)")
  return(p)
}
