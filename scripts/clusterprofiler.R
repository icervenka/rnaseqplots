library(magrittr, include.only = "%>%")

plot_pathways_meta <- function(df, top_pathways = 30) {
  pathways_summary <- df %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(
      n_pathways = dplyr::n(),
      mean_enrichment = mean(`GeneRatio/NES`)
    )

  p1 <- pathways_summary %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x = n_pathways),
      bins = 100,
      fill = "steelblue"
    ) +
    ggplot2::theme_bw()

  # top pathway contributors
  p2 <- pathways_summary %>%
    dplyr::top_n(30, n_pathways) %>%
    dplyr::arrange(n_pathways) %>%
    dplyr::mutate(source = factor(source, levels = source)) %>%
    ggplot2::ggplot(ggplot2::aes(x = n_pathways, y = source)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::theme_bw()

  # bottom pathway contributors
  p3 <- pathways_summary %>%
    dplyr::top_n(30, -n_pathways) %>%
    dplyr::arrange(n_pathways) %>%
    dplyr::mutate(source = factor(source, levels = source)) %>%
    ggplot2::ggplot(ggplot2::aes(x = n_pathways, y = source)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::theme_bw()

  return(cowplot::plot_grid(p1, p2, p3, nrow = 1))
}

plot_pathway_bargraph <- function(df,
                                  pathway_source,
                                  top_n = 20,
                                  truncate_desc = 80) {
  df %>%
    dplyr::filter(grepl(pathway_source, source)) %>%
    dplyr::arrange(-abs(`GeneRatio/NES`)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::arrange(`GeneRatio/NES`) %>%
    dplyr::mutate(
      Description =
        stringr::str_trunc(Description, truncate_desc)
    ) %>%
    dplyr::mutate(Description = factor(Description,
      levels = (.) %>% dplyr::pull(Description)
    )) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = `GeneRatio/NES`,
      y = Description,
      color = -log10(p.adjust),
      fill = -log10(p.adjust)
    )) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::ylab("") +
    ggplot2::theme_bw()
}
