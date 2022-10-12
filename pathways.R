library(magrittr, include.only = "%>%")

plot_pathways_rank = function(pathway_data,
                              rank_include,
                              x_axis,
                              y_axis,
                              bar_fill,
                              plot_x_text = T,
                              plot_legend = T,
                              ylab_text = NULL,
                              legend_text = NULL,
                              label_offset = 0.4) {

  if(max(rank_include) > nrow(pathway_data)) stop("Some ranks not found in filtered dataset.")

  pathway_data = pathway_data %>%
    #mutate(pathway_rank = row_number()) %>%
    arrange(desc({{ x_axis }})) %>%
    mutate(rank_of = paste0(pathway_rank, "/", nrow(.))) %>%
    filter(pathway_rank %in% rank_include) %>%
    mutate({{ y_axis }} := stringr::str_replace_all({{ y_axis }}, "_", " ")) %>%
    mutate({{ y_axis }} := tools::toTitleCase(tolower({{ y_axis }}))) %>%
    mutate({{ y_axis }} := factor({{ y_axis }}, levels = rev({{ y_axis }})))

  p = pathway_data %>%
    ggplot(aes(x = {{ y_axis }}, y = {{ x_axis }})) +
    geom_bar(stat = "identity", aes(fill = -log10({{ bar_fill }})), size = 1) +
    geom_label(aes(label = rank_of),
               nudge_y = (-1)*sign(pathway_data[[deparse(substitute(x_axis))]])*label_offset,
               size = 2.5) +
    coord_flip() +
    theme_minimal() +
    scale_fill_viridis() +
    xlab("")

  if(!is.null(ylab_text)) { p = p + ylab(ylab_text) }
  if(!is.null(legend_text)) { p = p + ylab(legend_text) }
  if(plot_legend == F) { p = p + theme(legend.position = "none") }

  if(plot_x_text == F) {
    cat(as.character(pathway_data[[deparse(substitute(y_axis))]]), sep = "\n")
    p = p +
      # theme(axis.text.y = element_blank(), legend.position = "none")
      theme(axis.text.y = element_blank())
  }
  p
}