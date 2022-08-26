library(SummarizedExperiment, include_only = NULL)
library("magrittr")

plot_pca_deseq <- function(dds,
                           group = "group",
                           plot_center = TRUE,
                           linetype = "solid",
                           palette = NULL) {
  vsd <- DESeq2::vst(dds, blind = FALSE)

  pca_data <- DESeq2::plotPCA(vsd, intgroup = c(group), returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  segments <- pca_data %>%
    dplyr::group_by(!!as.symbol(group)) %>%
    dplyr::summarise(xend = mean(PC1), yend = mean(PC2))
  pca_data <- merge(pca_data, segments, by = group)

  no_colors <- pca_data[, group] %>%
    unique() %>%
    length()

  p <- pca_data %>%
    ggplot2::ggplot(ggplot2::aes(
      PC1,
      PC2,
      color = !!as.symbol(group)
    )) +
    ggplot2::geom_point(size = 3)

  if (plot_center == TRUE) {
    p <- p +
      ggplot2::geom_segment(ggplot2::aes(
        x = PC1,
        y = PC2,
        xend = xend,
        yend = yend
      ),
      size = 0.5,
      linetype = linetype
      ) +
      ggplot2::geom_point(
        data = segments,
        ggplot2::aes(x = xend, y = yend),
        size = 1
      )
  }

  p <- p +
    ggplot2::xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ggplot2::ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")

  if (is.null(palette)) {
    p <- p +
      ggplot2::scale_color_manual(values = viridis::viridis(no_colors + 1))
  } else {
    p <- p +
      ggplot2::scale_color_manual(values = palette)
  }
  return(p)
}

plot_pca_common <- function(expression_data,
                            metadata,
                            color_by = "group") {
  pca <- prcomp(
    t(
      expression_data %>%
        dplyr::select(dplyr::matches(metadata$sample))
    ),
    scale = TRUE
  )
  pca_data <- as.data.frame(pca$x) %>%
    tibble::rownames_to_column("id")
  pca_data <- merge(pca_data, metadata, by.x = "id", by.y = "sample")
  segments <- pca_data %>%
    dplyr::group_by(!!as.symbol(color_by)) %>%
    dplyr::summarise(xend = mean(PC1), yend = mean(PC2))
  pca_data <- merge(pca_data, segments, by = color_by)

  percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
  percentage <- paste0(
    colnames(pca),
    "(",
    paste(as.character(percentage),
      "%",
      ")",
      sep = ""
    )
  )

  no_colors <- pca_data[, color_by] %>%
    unique() %>%
    length()

  p <- pca_data %>%
    ggplot2::ggplot(ggplot2::aes(
      PC1,
      PC2,
      color = !!as.symbol(color_by),
      shape = !!as.symbol(color_by)
    )) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_segment(ggplot2::aes(
      x = PC1,
      y = PC2,
      xend = xend,
      yend = yend
    ),
    size = 0.5,
    linetype = "dashed"
    ) +
    ggplot2::geom_point(
      data = segments,
      ggplot2::aes(x = xend, y = yend),
      size = 2
    ) +
    ggplot2::xlab(paste0("PC1: ", percentage[1], " variance")) +
    ggplot2::ylab(paste0("PC2: ", percentage[2], " variance")) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(
      values = viridis::viridis(no_colors + 1)
      [-length(viridis::viridis(no_colors + 1))]
    ) +
    ggplot2::theme(legend.position = "top")
  return(p)
}
