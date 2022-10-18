library(magrittr, include.only = "%>%")

plot_pca_deseq <- function(dds,
                           norm = "vst",
                           group = "group",
                           plot_center = TRUE,
                           linetype = "solid",
                           palette = NULL) {
  if (norm == "vst") {
    dds_norm <- DESeq2::vst(dds, blind = FALSE)
  } else if (norm == "rlog") {
    dds_norm <- DESeq2::rlog(dds, blind = FALSE)
  }


  pca_data <- DESeq2::plotPCA(dds_norm, intgroup = c(group), returnData = TRUE) %>%
    tibble::rownames_to_column(var = "id")
  percent_var_labs <- round(100 * attr(pca_data, "percentVar"))
  segments <- pca_data %>%
    dplyr::group_by(!!as.symbol(group)) %>%
    dplyr::summarise(xend = mean(PC1), yend = mean(PC2))
  pca_data <- pca_data %>%
    dplyr::select(dplyr::matches(paste0("^", group, "$")), id, PC1, PC2) %>%
    dplyr::left_join(segments, by = group)

  p <- plot_pca(
    pca_data,
    percent_var_labs,
    group,
    plot_center,
    linetype,
    palette
  )
  return(p)
}

# TODO results from prcomp give very different numbers
plot_pca_common <- function(expression_data,
                            metadata,
                            group = "group",
                            plot_center = TRUE,
                            linetype = "solid",
                            palette = NULL) {
  pca <- prcomp(
    t(
      expression_data %>%
        dplyr::select(dplyr::matches(metadata$sample))
    ),
    scale = TRUE
  )
  pca_data <- as.data.frame(pca$x) %>%
    tibble::rownames_to_column("id") %>%
    dplyr::left_join(metadata, by = c("id" = "sample"))
  segments <- pca_data %>%
    dplyr::group_by(!!as.symbol(group)) %>%
    dplyr::summarise(xend = mean(PC1), yend = mean(PC2))
  pca_data <- pca_data %>%
    dplyr::left_join(segments, by = group)

  percent_var_labs <- round(pca$sdev / sum(pca$sdev) * 100, 2)
  percent_var_labs <- paste0(
    colnames(pca),
    "(",
    paste(as.character(percent_var_labs),
      "%",
      ")",
      sep = ""
    )
  )

  p <- plot_pca(
    pca_data,
    percent_var_labs,
    group,
    plot_center,
    linetype,
    palette
  )
}

plot_pca <- function(pca_data,
                     percent_var_labs,
                     group,
                     plot_center,
                     linetype,
                     palette) {
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
        ggplot2::aes(x = xend, y = yend),
        size = 1
      )
  }

  p <- p +
    ggplot2::xlab(paste0("PC1: ", percent_var_labs[1], "% variance")) +
    ggplot2::ylab(paste0("PC2: ", percent_var_labs[2], "% variance")) +
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
