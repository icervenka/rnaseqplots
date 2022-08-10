
library(SummarizedExperiment, include_only = NULL)
plot_pca_deseq <- function(dds, group = "group") {
  vsd <- DESeq2::vst(dds, blind = F)

  pcaData <- DESeq2::plotPCA(vsd, intgroup = c(group), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  segments <- pcaData %>%
    dplyr::group_by(!!as.symbol(group)) %>%
    dplyr::summarise(xend = mean(PC1), yend = mean(PC2))
  pcaData <- merge(pcaData, segments, by = group)

  no_colors <- pcaData[, group] %>%
    unique() %>%
    length()

  p <- pcaData %>%
    ggplot2::ggplot(ggplot2::aes(PC1,
      PC2,
      color = !!as.symbol(group),
      shape = !!as.symbol(group)
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
    ggplot2::geom_point(data = segments, ggplot2::aes(x = xend, y = yend), size = 2) +
    ggplot2::xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ggplot2::ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = viridis::viridis(no_colors + 1)) +
    ggplot2::theme(legend.position = "top")

  return(p)
}

plot_pca_common <- function(expression_data, metadata, color_by = "group") {
  pca <- prcomp(t(expression_data %>% dplyr::select(dplyr::matches(metadata$sample))), scale = T)
  pcaData <- as.data.frame(pca$x) %>% tibble::rownames_to_column("id")
  pcaData <- merge(pcaData, metadata, by.x = "id", by.y = "sample")
  segments <- pcaData %>%
    dplyr::group_by(!!as.symbol(color_by)) %>%
    dplyr::summarise(xend = mean(PC1), yend = mean(PC2))
  pcaData <- merge(pcaData, segments, by = color_by)

  percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
  percentage <- paste0(colnames(pca), "(", paste(as.character(percentage), "%", ")", sep = ""))

  no_colors <- pcaData[, color_by] %>%
    unique() %>%
    length()

  p <- pcaData %>%
    ggplot2::ggplot(ggplot2::aes(PC1,
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
    ggplot2::geom_point(data = segments, ggplot2::aes(x = xend, y = yend), size = 2) +
    ggplot2::xlab(paste0("PC1: ", percentage[1], " variance")) +
    ggplot2::ylab(paste0("PC2: ", percentage[2], " variance")) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = viridis::viridis(no_colors + 1)[-length(viridis::viridis(no_colors + 1))]) +
    ggplot2::theme(legend.position = "top")

  return(p)
}
