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

  lfcs <- lfcs %>% dplyr::mutate(err_sq = (log2_fold_change.y - log2_fold_change.x)^2)
  return(lfcs)
}
