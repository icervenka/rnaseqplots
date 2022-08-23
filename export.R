library(egg)
ggsave_fixed = function(file, plot = last_plot(), units = "mm", margin = 1,
                        plot_width = 4, plot_height = 4,
                        width = round(dev.size()[1], digits = 1),
                        height = round(dev.size()[1], digits = 1)) {
  pf = egg::set_panel_size(p = plot,
                           file = NULL,
                           margin = unit(margin, units),
                           width = unit(plot_width, units),
                           height = unit(plot_height, units))
  ggplot2::ggsave(file, plot = pf, units = units, width = width, height = height)
}

create_gsea_normalized = function(data,
                                  metadata,
                                  data_name_col,
                                  data_description_col = NULL,
                                  metadata_sample_col = sample,
                                  out = NULL) {
  
  samples = metadata[[deparse(substitute(metadata_sample_col))]]
  data_cols = names(data)[names(data) %in% samples]
  
  out_df = data %>%
    dplyr::select(NAME = {{ data_name_col }}, any_of(data_cols))
  
  if(is.null(data_description_col)) {
    out_df = out_df %>%
      dplyr::mutate(DESCRIPTION = "") %>%
      dplyr::relocate(DESCRIPTION, .after = NAME)
  } else {
    out_df = out_df %>%
      dplyr::mutate(DESCRIPTION = {{ data_description_col }}) %>%
      dplyr::relocate(DESCRIPTION, .after = NAME)
  }
  
  if(!is.null(out)) {
    write.table(out_df,
                out,
                sep = "\t",
                row.names = F,
                quote = F)
  } else {
    return(out_df)
  }
  
}

create_gsea_cls = function(data,
                           metadata,
                           out = NULL,
                           sample_col = sample,
                           group_col = group) {
  metadata_filtered = metadata %>%
    dplyr::select({{ sample_col }}, {{ group_col }})
  
  samples = metadata_filtered[[deparse(substitute(sample_col))]]
  groups = metadata_filtered[[deparse(substitute(group_col))]] %>%
    stringr::str_replace_all("[:space:]+", "")
  data_cols = names(data)[names(data) %in% samples]
  sample_indices = match(samples, data_cols)
  
  text_lines = c(paste(samples %>% unique %>% length,
                       groups %>% unique %>% length,
                       "1"),
                 paste("#", paste(groups %>% unique, collapse = " ")),
                 paste(groups[order(sample_indices)], collapse = " ")
  )
  
  if(!is.null(out)) {
    file_conn = file(out)
    writeLines(text_lines, file_conn)
    close(file_conn)
  } else {
    return(text_lines)
  }
}

create_gsea_rank = function(data,
                            out = NULL,
                            ranking_equation = ~ -log10(pvalue)*sign(log2FoldChange),
                            gene_id_column = ensembl_gene_id,
                            .inf = "replace") {
  out_df = data %>%
    dplyr::mutate(rank = !!lazyeval::f_rhs(ranking_equation)) %>%
    dplyr::select({{ gene_id_column }}, rank) %>%
    dplyr::arrange(desc(rank)) %>%
    tidyr::drop_na() 
  
  if(.inf == "drop") {
    out_df = out_df %>%
      dplyr::filter(is.finite(rank))
  } else if(.inf == "replace") {
    min_val = out_df %>%
      dplyr::filter(is.finite(rank)) %>%
      pull(rank) %>%
      min()
    max_val = out_df %>%
      dplyr::filter(is.finite(rank)) %>%
      pull(rank) %>%
      max()
    
    out_df$rank[which(out_df$rank == Inf)] = max_val + 0.01 * max_val
    out_df$rank[which(out_df$rank == -Inf)] = min_val + 0.01 * min_val
  }
  
  if(!is.null(out)) {
    write.table(out_df,
                out,
                sep = "\t",
                row.names = F,
                quote = F)
  } else {
    return(out_df)
  }
}