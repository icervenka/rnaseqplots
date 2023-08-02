#' Wrapper for cor.test for data frames
#'
#'
#' @param df data frame with numeric data to correlate, usually grouped by
#' another column
#' @param x numeric data frame colum, supplied as variable
#' @param y numeric data frame colum, supplied as variable
#'
#' @return tidy data frame of correlations
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
cor_test_df <- function(df, x, y) {
  df <- stats::cor.test(
    df[[deparse(substitute(x))]],
    df[[deparse(substitute(y))]],
    method = "pearson"
  ) %>%
    generics::tidy()
  return(df)
}

#' Filter DESeq2 result object
#'
#' Converts to data frame and adds a column of values based on rownames
#'
#' @param res DESEeq2 result object
#' @param colname column to filter, supplied as variable. default: padj
#' @param threshold double, filtering threshold. Will be applied on absolute
#' values in colname. default: 0.05
#' @param rownames_to character string, name of the column that will be created
#' from rownames of DESeq2 result object. default ensembl_gene_id
#' @param comparison function, which comparison function to use to compare the
#' values to threshold. default `<`
#'
#' @return data frame
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
filter_deseq_results <- function(res,
                                 colname = padj,
                                 threshold = 0.05,
                                 comparison = `<`,
                                 rownames_to = "ensembl_gene_id") {
  # fix global variable binding notes
  padj <- NULL
  
  if (!(rownames_to %in% names(res))) {
    res <- res %>%
      data.frame() %>%
      tibble::rownames_to_column(rownames_to)
  }

  if (!is.null(colname)) {
    res <- res %>%
      dplyr::filter(comparison(abs(colname), threshold))
  }
  return(res)
}

#' Filter column values based on regex
#'
#' @param df data frame to filter
#' @param colname character string column name of values to filter
#' @param str regex to filter by
#'
#' @return data frame
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
filter_string <- function(df, colname, str) {
  df <- df %>%
    dplyr::filter(grepl(str, colname))
  return(df)
}

#' Get significantly regulated genes in character vector
#'
#' @param diffexp_data data frame of gene differential expression, must contain
#' column for gene IDs/symbols and a p-value
#' @param id_colname column name of gene IDs/symbols in differential expression
#' data file, supplied as variable
#' @param filter_sig_on usually column name of p-values in differential
#' expression data file values of which will determine significance,
#' supplied as variable. default: padj
#' @param sig_threshold double, significance threshold for filter_sig_on column
#'
#' @return character vector of significant gene IDs/symbols
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
get_sig_genes <- function(diffexp_data,
                          id_colname = SYMBOL,
                          filter_sig_on = padj,
                          sig_threshold = 0.05) {
  # fix global variable binding notes
  SYMBOL <- padj <- NULL
                            
  gene_vec <- diffexp_data %>%
    dplyr::filter({{ filter_sig_on }} < sig_threshold) %>%
    dplyr::pull({{ id_colname }})
  return(gene_vec)
}

#' Get gene ID mapping from biomart
#'
#' An alternative to clusterProfiler::bitr function
#'
#' @param mart biomart mart to use
#' @param id_attribute character string, source attribute name
#' @param mapped_attribute character vector, names of target attributes
#' @param ... other parameters to biomaRt getBM function
#'
#' @return data frame of id_attributed IDs mapped to mapped_attribute IDs
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
get_biomart_gene_mapping <- function(mart,
                                     id_attribute,
                                     mapped_attribute,
                                     ...) {
  # fix global variable binding notes
  id <- NULL
  
  if (!all(c(
    id_attribute %in% list(...)[["attributes"]],
    mapped_attribute %in% list(...)[["attributes"]]
  ))) {
    stop("Id attribute or mapped attribute not present in attribute list.")
  }

  mapping <- biomaRt::getBM(mart = mart, ...)
  all_values <- data.frame(
    id = list(...)[["values"]],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(mapping, by = c("id" = id_attribute)) %>%
    dplyr::mutate(!!as.symbol(mapped_attribute) := dplyr::case_when(
      is.na(!!as.symbol(mapped_attribute)) ~ id,
      !!as.symbol(mapped_attribute) == "" ~ id,
      TRUE ~ !!as.symbol(mapped_attribute)
    )) %>%
    dplyr::group_by(id) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(!!as.symbol(id_attribute) := id, .before = id) %>%
    dplyr::select(-id)

  return(all_values)
}
