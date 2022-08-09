plot_dire = function(df, color = "steelblue") {
  p = df %>%
    ggplot2::ggplot(aes(x = Occurrence, y = Importance)) +
    geom_point(color = color) +
    theme_bw()
  return(p)
}

plot_dire_labeled = function(df,
                             occurrence_threshold = 0.05,
                             importance_threshold = 0.05,
                             dot_size = 3.5,
                             color = "steelblue") {

  p = df %>%
    plot_dire(color) +
    geom_label_repel(data = . %>%
                       filter(Occurrence > occurrence_threshold | Importance > importance_threshold),
                     aes(x = Occurrence, y = Importance, label = `Transcription Factor`),
                     size = dot_size)
  return(p)
}
