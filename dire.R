plot_dire <- function(df, color = "steelblue") {
  p <- df %>%
    ggplot2::ggplot(ggplot2::aes(x = Occurrence, y = Importance)) +
    ggplot2::geom_point(color = color) +
    ggplot2::theme_bw()
  return(p)
}

plot_dire_labeled <- function(df,
                              occurrence_threshold = 0.05,
                              importance_threshold = 0.05,
                              dot_size = 3.5,
                              color = "steelblue") {
  p <- df %>%
    plot_dire(color) +
    ggrepel::geom_label_repel(
      data = . %>%
        dplyr::filter(Occurrence > occurrence_threshold | Importance > importance_threshold),
      ggplot2::aes(x = Occurrence, y = Importance, label = `Transcription Factor`),
      size = dot_size
    )
  return(p)
}
