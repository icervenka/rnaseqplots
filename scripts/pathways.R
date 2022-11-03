library(magrittr, include.only = "%>%")

plot_pathways_rank <- function(pathway_data,
                               rank_include,
                               x_axis,
                               y_axis,
                               bar_fill,
                               plot_x_text = TRUE,
                               plot_legend = TRUE,
                               ylab_text = NULL,
                               legend_text = NULL,
                               label_offset = 0.4) {
  if (max(rank_include) > nrow(pathway_data)) {
    stop("Some ranks not found in filtered dataset.")
  }

  pathway_data <- pathway_data %>%
    # mutate(pathway_rank = row_number()) %>%
    arrange(desc({{ x_axis }})) %>%
    mutate(rank_of = paste0(pathway_rank, "/", nrow(.))) %>%
    filter(pathway_rank %in% rank_include) %>%
    mutate({{ y_axis }} := stringr::str_replace_all({{ y_axis }}, "_", " ")) %>%
    mutate({{ y_axis }} := tools::toTitleCase(tolower({{ y_axis }}))) %>%
    mutate({{ y_axis }} := factor({{ y_axis }}, levels = rev({{ y_axis }})))

  p <- pathway_data %>%
    ggplot(aes(x = {{ y_axis }}, y = {{ x_axis }})) +
    geom_bar(stat = "identity", aes(fill = -log10({{ bar_fill }})), size = 1) +
    geom_label(aes(label = rank_of),
      nudge_y = (-1) *
        sign(pathway_data[[deparse(substitute(x_axis))]]) *
        label_offset,
      size = 2.5
    ) +
    coord_flip() +
    theme_minimal() +
    scale_fill_viridis() +
    xlab("")

  if (!is.null(ylab_text)) {
    p <- p + ylab(ylab_text)
  }
  if (!is.null(legend_text)) {
    p <- p + ylab(legend_text)
  }
  if (plot_legend == FALSE) {
    p <- p + theme(legend.position = "none")
  }
  if (plot_x_text == FALSE) {
    cat(as.character(pathway_data[[deparse(substitute(y_axis))]]), sep = "\n")
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}

plot_pathways_volcano <- function(pathway_data,
                                  x_axis,
                                  y_axis,
                                  label = name,
                                  x_axis_threshold = 0,
                                  y_axis_threshold = 0,
                                  label_pathways = NULL,
                                  xlab_text = "Normalized Enrichment Score",
                                  ylab_text = "-log10(FDR)",
                                  add_vhlines = TRUE,
                                  vhline_color = "gray40",
                                  vhline_type = "dashed",
                                  alpha = 0.2,
                                  color_palette = c(
                                    "gray70",
                                    "steelblue4",
                                    "darkred"
                                  )) {
  pathway_data <- pathway_data %>%
    dplyr::mutate(significant_x = dplyr::case_when(
      abs({{ x_axis }}) > x_axis_threshold ~ 1,
      TRUE ~ 0
    )) %>%
    dplyr::mutate(significant_y = dplyr::case_when(
      abs({{ y_axis }}) > y_axis_threshold ~ 1,
      TRUE ~ 0
    )) %>%
    dplyr::mutate(significant = significant_x + significant_y) %>%
    dplyr::mutate(significant = as.factor(significant))

  color_palette <- color_palette[1:max(as.numeric(pathway_data$significant))]

  p <- pathway_data %>%
    ggplot2::ggplot(ggplot2::aes(x = {{ x_axis }}, y = -log10({{ y_axis }}))) +
    ggplot2::geom_point(ggplot2::aes(color = significant),
      size = 3,
      alpha = alpha
    ) +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab(xlab_text) +
    ggplot2::ylab(ylab_text)

  if (add_vhlines) {
    if (x_axis_threshold > 0) {
      p <- p + ggplot2::geom_vline(
        xintercept = {{ x_axis_threshold }},
        color = vhline_color,
        linetype = vhline_type
      )
    }
    if (y_axis_threshold > 0) {
      p <- p + ggplot2::geom_hline(
        yintercept = -log10({{ y_axis_threshold }}),
        color = vhline_color,
        linetype = vhline_type
      )
    }
  }

  if (is.vector(label_pathways, mode = "character")) {
    p <- p +
      ggrepel::geom_label_repel(
        data = pathway_data %>%
          dplyr::filter({{ label }} %in% label_pathways),
        ggplot2::aes(label = {{ label }}),
        min.segment.length = grid::unit(0, "lines")
      )
  } else if (is.vector(label_pathways, mode = "numeric")) {
    p <- p +
      ggrepel::geom_label_repel(
        data = pathway_data %>%
          dplyr::filter(pathway_rank %in% label_pathways),
        ggplot2::aes(label = {{ label }}),
        min.segment.length = grid::unit(0, "lines")
      )
  }
  return(p)
}

plot_cp_pathways_meta <- function(pathway_data, top_pathways = 30) {
  pathways_summary <- pathway_data %>%
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
    ggplot2::theme_bw() +
    xlab("# of significant pathways") +
    ggtitle("Significant pathways histogram")

  # top pathway contributors
  p2 <- pathways_summary %>%
    dplyr::top_n(top_pathways, n_pathways) %>%
    dplyr::arrange(n_pathways) %>%
    dplyr::mutate(source = factor(source, levels = source)) %>%
    ggplot2::ggplot(ggplot2::aes(x = n_pathways, y = source)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::theme_bw() +
    xlab("# of significant pathways") +
    ggtitle("Analyses with highest\n# of significant pathways")

  # bottom pathway contributors
  p3 <- pathways_summary %>%
    dplyr::top_n(top_pathways, -n_pathways) %>%
    dplyr::arrange(n_pathways) %>%
    dplyr::mutate(source = factor(source, levels = source)) %>%
    ggplot2::ggplot(ggplot2::aes(x = n_pathways, y = source)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::theme_bw() +
    xlab("# of significant pathways") +
    ggtitle("Analyses with lowest\n# of significant pathways")

  return(cowplot::plot_grid(p1, p2, p3, nrow = 1))
}

plot_cp_pathways_bargraph <- function(pathway_data,
                                      pathway_source_pattern,
                                      top_pathways = 20,
                                      truncate_description = 80) {
  p <- pathway_data %>%
    dplyr::filter(grepl(pathway_source_pattern, source)) %>%
    dplyr::arrange(-abs(`GeneRatio/NES`)) %>%
    dplyr::slice_head(n = top_pathways) %>%
    dplyr::arrange(`GeneRatio/NES`) %>%
    dplyr::mutate(
      Description =
        stringr::str_trunc(Description, truncate_description)
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
  return(p)
}
