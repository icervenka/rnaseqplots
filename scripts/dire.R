library(magrittr, include.only = "%>%")

#' Creates plot from dire.dcode.com upstream transcription factor data
#'
#' @param df data frame containing column of transcription factor names,
#' Occurence and Importance columns
#' @param dot_size size of scatter plot symbols
#' @param color color of scatter plot symbols
#'
#' @return ggplot scatter plot
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point theme_bw
#'
#' @examples
plot_dire <- function(df, dot_size = 2.5, color = "steelblue4") {
  p <- df %>%
    ggplot2::ggplot(ggplot2::aes(x = Occurrence, y = Importance)) +
    ggplot2::geom_point(color = color, size = dot_size) +
    ggplot2::theme_bw()
  return(p)
}

#' Creates labeled plot from dire.dcode.com upstream transcription factor data
#'
#' @param df data frame containing column of transcription factor names,
#' Occurence and Importance columns
#' @param occurrence_threshold double, minimum occurrence value to plot gene
#' labels
#' @param importance_threshold double, minimum importance value to plot gene
#' labels
#' @param gene_list character vector, label specific set of genes, overrides
#' occurrence and importance thresholds
#' @param combine_thresholds binary logical operator, how to combine the
#' occurrence and importance thresholds, default '&'
#' @param ... other parameters to plot_dire function
#'
#' @return ggplot scatter plot with labeled transcription factors of interest
#' @export
#'
#' @importFrom dplyr filter
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggplot2 aes
#' @importFrom grid unit
#'
#' @examples
plot_dire_labeled <- function(df,
                              occurrence_threshold = 0.05,
                              importance_threshold = 0.05,
                              gene_list = NULL,
                              combine_thresholds = `&`,
                              ...) {
  p <- df %>%
    plot_dire(...)

  if (is.null(gene_list)) {
    label_data <- df %>%
      dplyr::filter(combine_thresholds(
        Occurrence > occurrence_threshold,
        Importance > importance_threshold
      ))
  } else if (is.vector(gene_list, mode = "character")) {
    label_data <- df %>%
      dplyr::filter(`Transcription Factor` %in% gene_list)
  }

  p <- p +
    ggrepel::geom_label_repel(
      data = label_data,
      ggplot2::aes(label = `Transcription Factor`),
      min.segment.length = grid::unit(0, "lines")
    )

  return(p)
}
