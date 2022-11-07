library(magrittr, include.only = "%>%")

cor_test_df <- function(df, x, y) {
  cor.test(
    df[[deparse(substitute(x))]],
    df[[deparse(substitute(y))]],
    method = "pearson"
  ) %>%
    tidy()
}

merge_cuffdiff_fc <- function(x, y, annot, by = "gene") {
  lfcs <- merge(x %>%
    dplyr::select(gene, gene_id, log2_fold_change, q_value),
  y %>%
    dplyr::select(gene, gene_id, log2_fold_change, q_value),
  by = by
  )
  lfcs <- lfcs %>%
    # Fisher method of combining p-values
    dplyr::mutate(chi_pcomb = -2 * (log(q_value.x) + log(q_value.y))) %>%
    # calculate significance of combined p-values
    dplyr::mutate(p_chi = pchisq(chi_pcomb, df = 4, lower.tail = FALSE)) %>%
    # calculate squared error of genes for ranking
    dplyr::mutate(err_sq = (log2_fold_change.y - log2_fold_change.x)^2)

  lfcs <- lfcs %>%
    dplyr::mutate(err_sq = (log2_fold_change.y - log2_fold_change.x)^2)
  return(lfcs)
}

filter_deseq_results <- function(res, padj_threshold = 0.05) {
  filtered_res <- res %>%
    data.frame() %>%
    tibble::rownames_to_column("ensembl_gene_id")

  if (!is.null(padj_threshold)) {
    filtered_res <- filtered_res %>%
      filter(padj < padj_threshold)
  }
  return(filtered_res)
}

filter_file <- function(df, str) {
  df %>%
    dplyr::filter(grepl(str, file))
}

filter_name <- function(df, str) {
  df %>%
    dplyr::filter(grepl(str, name))
}

get_sig_genes <- function(diffexp_data,
                          colname = SYMBOL,
                          sig_threshold = 0.05,
                          filter_sig_on = padj) {
  gene_vec <- diffexp_data %>%
    dplyr::filter({{ filter_sig_on }} < sig_threshold) %>%
    pull({{ colname }})
  return(gene_vec)
}

get_biomart_gene_mapping <- function(mart,
                                     id_attribute,
                                     mapped_attribute,
                                     ...) {
  if (!all(c(
    id_attribute %in% list(...)[["attributes"]],
    mapped_attribute %in% list(...)[["attributes"]]
  ))) {
    stop("Id attribute or mapped attribute not present in attribute list.")
  }

  mapping <- biomart::getBM(mart = mart, ...)
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