#' Simple barplot of differentally regulated pathways
#'
#' Uses data from clusterprofiler_reports_snakemake analysis.
#'
#' @param pathway_data data frame with data for pathway analysis. Should include
#' a columns for enrichment, p-value, pathway name and sorted based on their
#' rank.
#' @param pathway_source_pattern character string, filter for a specific pathway
#' source such as Gene Ontology, WikiPathways. See
#' clusterprofiler_reports_snakemake output file to see the possible values.
#' @param top_pathways integer, top n pathways to display. default: 20
#' @param truncate_name integer, maximum number of characters after which the
#' pathway names will be truncated. default: 80
#'
#' @return ggplot barplot
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
plot_cp_pathways_bargraph <- function(pathway_data,
                                      pathway_source_pattern,
                                      top_pathways = 20,
                                      truncate_name = 80) {
  p <- pathway_data %>%
    dplyr::filter(grepl(pathway_source_pattern, source)) %>%
    dplyr::arrange(-abs(`GeneRatio/NES`)) %>%
    dplyr::slice_head(n = top_pathways) %>%
    dplyr::arrange(`GeneRatio/NES`) %>%
    dplyr::mutate(
      Description =
        stringr::str_trunc(Description, truncate_name)
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

# #' Title
# #'
# #' @param json_file
# #' @param top_pathways
# #'
# #' @return
# #' @export
# #'
# #' @examples
# plot_basic_amigo = function(json_file, top_pathways = 15) {
#   df = read_pathways_amigo(fromJSON(file = json_file))
#   df %>%
#     filter(level == 0) %>%
#     arrange(fold_enrichment) %>%
#     slice_tail(n = top_n) %>%
#     mutate(label = tools::toTitleCase(label)) %>%
#     mutate(label = factor(label, levels = label)) %>%
#     ggplot(aes(x = label, y = fold_enrichment)) +
#     geom_bar(aes(fill = -log10(pvalue)), stat = "identity") +
#     #geom_point(aes(color = -log10(pvalue)), size = 3) +
#     coord_flip() +
#     theme_bw() +
#     scale_fill_viridis()
# }

#' Create barplot for pathways with indicated rank
#'
#' Ranking is displayed as 'x/y' where x is the pathway rank and y is the total
#' number of idenfied pathways. Ranking labels are displayed inside the bars.
#'
#' @param pathway_data data frame with data for pathway analysis. Should include
#' a columns for enrichment, p-value, pathway name and sorted based on their
#' rank.
#' @param rank_include charater or integer vector. Which pathways to include in
#' the graph. If character vector, pathways are selected by names, if integer
#' vector, pathways are selected by their rank.
#' @param x What to plot on x axis, supplied as variable present as
#' a column in pathway_data variable. Note, if the axes  are flipped for
#' legibility, this will become the vertical axis.
#' @param y What to plot on y axis, supplied as variable present as
#' a column in pathway_data variable. Note, if the axes  are flipped for
#' legibility, this will become the horizontal axis.
#' @param bar_fill What to plot as bar fill, supplied as variable.
#' @param truncate_name integer, maximum number of characters after which the
#' pathway names will be truncated. default: 80
#' @param flip_axis logical, whether axes should be flipped for legibilty.
#' default: TRUE
#' @param plot_x_text logical, whether to plot x axis labels (usually text).
#' Provided in case user wants to ensure equal sizing of graph area or make sure
#' that the axis labels are an editable text in graphing processing software of
#' choice.
#' @param plot_legend logical, whether to plot legend. See plot_x_text.
#' @param ylab_text character string, text label for y axis.
#' @param label_offset double, offset of rank labels from the end of bar for
#' each pathway. Optimal offset will depend on the final size of the graph.
#' Will probably require some tweaking.
#'
#' @return ggplot barplot with pathway ranking
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
plot_pathways_rank <- function(pathway_data,
                               rank_include,
                               x,
                               y,
                               bar_fill,
                               truncate_name = 80,
                               flip_axis = TRUE,
                               plot_x_text = TRUE,
                               plot_legend = TRUE,
                               ylab_text = NULL,
                               label_offset = 0.4) {
  if (max(rank_include) > nrow(pathway_data)) {
    stop("Some ranks not found in filtered dataset.")
  }

  # process pathways data
  pathway_data <- pathway_data %>%
    # mutate(pathway_rank = row_number()) %>%
    dplyr::arrange(dplyr::desc({{ y }})) %>%
    dplyr::mutate(rank_of = paste0(pathway_rank, "/", nrow(.))) %>%
    dplyr::filter(pathway_rank %in% rank_include) %>%
    dplyr::mutate({{ x }} := stringr::str_replace_all({{ x }}, "_", " ")) %>%
    dplyr::mutate({{ x }} := tools::toTitleCase(tolower({{ x }}))) %>%
    dplyr::mutate(
      {{ x }} :=
        stringr::str_trunc({{ x }}, truncate_name)
    ) %>%
    dplyr::mutate({{ x }} := factor({{ x }}, levels = rev({{ x }})))

  # Create graph
  p <- pathway_data %>%
    ggplot2::ggplot(aes(x = {{ x }}, y = {{ y }})) +
    ggplot2::geom_bar(
      stat = "identity",
      ggplot2::aes(fill = -log10({{ bar_fill }})),
      size = 1
    ) +
    ggplot2::geom_label(ggplot2::aes(label = rank_of),
      nudge_y = (-1) *
        sign(pathway_data[[deparse(substitute(y))]]) *
        label_offset,
      size = 2.5
    ) +
    ggplot2::theme_minimal() +
    { # nolint
      if (!is.null(ylab_text)) ggplot2::ylab(ylab_text)
    } + # nolint
    # { if (!is.null(legend_text)) ylab(legend_text) } + # nolint
    { # nolint
      if (plot_legend == FALSE) ggplot2::theme(legend.position = "none")
    } + # nolint
    viridis::scale_fill_viridis() +
    ggplot2::xlab("")

  # If the x axis text is not plotted, it is written to the console so the user
  # can copy-paste it to the program of choice
  if (plot_x_text == FALSE) {
    cat(as.character(pathway_data[[deparse(substitute(x))]]), sep = "\n")
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }

  if (flip_axis) {
    p <- p + ggplot2::coord_flip()
  }

  return(p)
}

#' Create volcano plot for identified pathways
#'
#' Works with pathway analysis that supplies a measure of overrepresentation
#' with directionality, such as GSEA.
#'
#' @param pathway_data data frame with data for pathway analysis. Should include
#' a columns for enrichment, p-value and a pathway name or ID.
#' @param x  What to plot on x axis, supplied as variable present as
#' a column in pathway_data variable.
#' @param y What to plot on y axis, supplied as variable present as
#' a column in pathway_data variable.
#' @param label Which column in the pathway_data contains names, supplied as
#' variable.
#' @param x_threshold double, significance threshold for x axis.
#' @param y_threshold double, significance threshold for y axis.
#' @param label_pathways charater or integer vector. Which pathways to include
#' in the graph. If character vector, pathways are selected by names, if integer
#' vector, pathways are selected by their rank.
#' @param xlab_text character string, text label for x axis. default:
#' "Normalized Enrichment Score"
#' @param ylab_text character string, text label for y axis. default:
#' "-log10(FDR)"
#' @param add_vhlines logical, whether to include horizontal a vertical lines
#' that show the significance threholds. default: TRUE
#' @param vhline_color character string, color of the vertical and horizontal
#' lines. default: grey40
#' @param vhline_type character string, line type of vertical and horizonal.
#' default: dashed
#' @param alpha double, transparency of the geom_point.
#' @param color_palette character vector of length 3, colors for the unchaged,
#' down- and upregulated pathways.
#'
#' @return ggplot pathway volcano plot
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
plot_pathways_volcano <- function(pathway_data,
                                  x,
                                  y,
                                  label = name,
                                  x_threshold = 0,
                                  y_threshold = 0,
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
  # prepare the pathway data for plotting, setting thresholds
  pathway_data <- pathway_data %>%
    dplyr::mutate(significant_x = dplyr::case_when(
      abs({{ x }}) > x_threshold ~ 1,
      TRUE ~ 0
    )) %>%
    dplyr::mutate(significant_y = dplyr::case_when(
      abs({{ y }}) > y_threshold ~ 1,
      TRUE ~ 0
    )) %>%
    dplyr::mutate(significant = significant_x + significant_y) %>%
    dplyr::mutate(significant = as.factor(significant))

  # subset the color palette if there are only up or downregulated pathways
  color_palette <- color_palette[1:max(as.numeric(pathway_data$significant))]

  # create basic plot
  p <- pathway_data %>%
    ggplot2::ggplot(ggplot2::aes(x = {{ x }}, y = -log10({{ y }}))) +
    ggplot2::geom_point(ggplot2::aes(color = significant),
      size = 3,
      alpha = alpha
    ) +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab(xlab_text) +
    ggplot2::ylab(ylab_text)

  # add vertical and horizontal lines
  if (add_vhlines) {
    if (x_threshold > 0) {
      p <- p + ggplot2::geom_vline(
        xintercept = {{ x_threshold }},
        color = vhline_color,
        linetype = vhline_type
      )
    }
    if (y_threshold > 0) {
      p <- p + ggplot2::geom_hline(
        yintercept = -log10({{ y_threshold }}),
        color = vhline_color,
        linetype = vhline_type
      )
    }
  }

  # add pathways labels depending on whether character or integer vector is
  # specified by the user
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

# TODO add normalized contributions based
#' Meta analysis of pathway data
#'
#' This function is used to plot the overview of differentially regulated
#' pathway sets from the  clusterprofiler_pipline workflow. Creates a grid of
#' 3 plots: i) histogram of number of differentially regulated pathways from
#' pathways databases, ii) databases with highest amount of differentially
#' regulated pathways, iii) databases with lowest amount of differentially
#' regulated pathways.
#'
#' @param pathway_data data frame with data for pathway analysis. Should include
#' a columns for enrichment, p-value, pathway name and sorted based on their
#' rank. Usually it would be the whole collated output from
#' clusterprofiler_pipline.
#' @param top_pathways integer, number of top and bottom patways to plot.
#' default: 30
#'
#' @return grid of 3 ggplot plots
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
plot_cp_pathways_meta <- function(pathway_data, top_pathways = 30) {
  pathways_summary <- pathway_data %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(
      n_pathways = dplyr::n(),
      mean_enrichment = mean(`GeneRatio/NES`)
    )

  # plot histogram of significant pathways identified by different pathway
  # databases
  p1 <- pathways_summary %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x = n_pathways),
      bins = 100,
      fill = "steelblue"
    ) +
    ggplot2::theme_bw() +
    ggplot2::xlab("# of significant pathways") +
    ggplot2::ggtitle("Significant pathways histogram") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

  # top pathway contributors
  p2 <- pathways_summary %>%
    dplyr::top_n(top_pathways, n_pathways) %>%
    dplyr::arrange(n_pathways) %>%
    dplyr::mutate(source = factor(source, levels = source)) %>%
    ggplot2::ggplot(ggplot2::aes(x = n_pathways, y = source)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::theme_bw() +
    ggplot2::xlab("# of significant pathways") +
    ggplot2::ggtitle("Analyses with highest\n# of significant pathways") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

  # bottom pathway contributors
  p3 <- pathways_summary %>%
    dplyr::top_n(top_pathways, -n_pathways) %>%
    dplyr::arrange(n_pathways) %>%
    dplyr::mutate(source = factor(source, levels = source)) %>%
    ggplot2::ggplot(ggplot2::aes(x = n_pathways, y = source)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::theme_bw() +
    ggplot2::xlab("# of significant pathways") +
    ggplot2::ggtitle("Analyses with lowest\n# of significant pathways") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

  return(cowplot::plot_grid(p1, p2, p3, nrow = 1))
}
