library(magrittr, include.only = "%>%")

plot_dire <- function(df, dot_size = 2.5, color = "steelblue4") {
  p <- df %>%
    ggplot2::ggplot(ggplot2::aes(x = Occurrence, y = Importance)) +
    ggplot2::geom_point(color = color, size = dot_size) +
    ggplot2::theme_bw()
  return(p)
}

plot_dire_labeled <- function(df,
                              occurrence_threshold = 0.05,
                              importance_threshold = 0.05,
                              label_genes = NULL,
                              combine_thresholds = `&`,
                              ...) {
  p <- df %>%
    plot_dire(...)

  if (is.null(label_genes)) {
    label_data <- df %>%
      dplyr::filter(combine_thresholds(
        Occurrence > occurrence_threshold,
        Importance > importance_threshold)
      )
  } else if (is.vector(label_genes, mode = "character")) {
    label_data <- df %>%
      dplyr::filter(`Transcription Factor` %in% label_genes)
  }

  p <- p +
    ggrepel::geom_label_repel(
      data = label_data,
      ggplot2::aes(label = `Transcription Factor`),
      min.segment.length = grid::unit(0, "lines")
    )

  return(p)
}
