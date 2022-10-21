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
      nudge_y = (-1) * sign(pathway_data[[deparse(substitute(x_axis))]]) * label_offset,
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
  if (plot_legend == F) {
    p <- p + theme(legend.position = "none")
  }

  if (plot_x_text == F) {
    cat(as.character(pathway_data[[deparse(substitute(y_axis))]]), sep = "\n")
    p <- p +
      # theme(axis.text.y = element_blank(), legend.position = "none")
      theme(axis.text.y = element_blank())
  }
  p
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

  color_palette = color_palette[1:max(as.numeric(pathway_data$significant))]
  
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
