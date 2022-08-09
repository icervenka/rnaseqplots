collate_pathways = function(pathway_files_basepath, pattern, padj_threshold = 0.05) {
  if(stringr::str_sub(path_to_pathway_files, -1) != "/") {
    pathway_files_basepath = paste0(pathway_files_basepath, "/")
  }

  pathway_files = list.files(pathway_files_basepath,
                             pattern = pattern,
                             full.names = F)

  cp_unified_colnames = c("ID",	"Description",	"GeneRatio/NES",	"pvalue",
                          "p.adjust",	"SYMBOL",	"ENTREZID",	"log2FoldChange")

  diffexp_pathways = purrr::map_dfr(pathway_files, function(x) {
    readr::read_delim(paste0(pathway_files_basepath, x),
                      delim = "\t", escape_double = FALSE,
                      trim_ws = TRUE,
                      show_col_types = F)  %>%
      setNames(cp_unified_colnames) %>%
      # TODO needs to be changed
      dplyr::mutate(source = gsub("_contrast.txt", "", x),
                    ID = as.character(ID)) %>%
      # might not reflect actual columns
      dplyr::select(-SYMBOL, -ENTREZID, -log2FoldChange) %>%
      unique()
  }) %>%
    filter(p.adjust < padj_threshold)
  return(diffexp_pathways)
}


plot_pathways_meta = function(df, top_pathways = 30) {
  pathways_summary = df %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(n_pathways = n(),
                     mean_enrichment = mean(`GeneRatio/NES`))

  p1 = pathways_summary %>%
    ggplot2::ggplot() +
    geom_histogram(aes(x = n_pathways),
                   bins = 100,
                   fill = "steelblue") +
    theme_bw()

  # top pathway contributors
  p2 = pathways_summary %>%
    dplyr::top_n(30, n_pathways) %>%
    dplyr::arrange(n_pathways) %>%
    dplyr::mutate(source = factor(source, levels = source)) %>%
    ggplot2::ggplot(aes(x = n_pathways, y = source)) +
    geom_bar(stat="identity", fill = "steelblue") +
    theme_bw()

  # bottom pathway contributors
  p3 = pathways_summary %>%
    dplyr::top_n(30, -n_pathways) %>%
    dplyr::arrange(n_pathways) %>%
    dplyr::mutate(source = factor(source, levels = source)) %>%
    ggplot2::ggplot(aes(x = n_pathways, y = source)) +
    geom_bar(stat="identity", fill = "steelblue") +
    theme_bw()

  return(plot_grid(p1, p2, p3, nrow = 1))
}

plot_pathway_bargraph = function(df, pathway_source, top_n = 20, truncate_desc = 80) {
  df %>%
    dplyr::filter(grepl(pathway_source, source)) %>%
    dplyr::arrange(-abs(`GeneRatio/NES`)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::arrange(`GeneRatio/NES`) %>%
    dplyr::mutate(Description = stringr::str_trunc(Description, truncate_desc)) %>%
    dplyr::mutate(Description = factor(Description, levels = (.) %>% pull(Description))) %>%
    ggplot2::ggplot(aes(x = `GeneRatio/NES`,
               y = Description,
               color = -log10(p.adjust),
               fill = -log10(p.adjust))) +
    geom_bar(stat="identity") +
    ylab("") +
    theme_bw()
}
