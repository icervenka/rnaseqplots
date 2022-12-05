library(magrittr, include.only = "%>%")

#' Draw PCA plot from DESeq2 dds data
#'
#' Adds optional group centers and palette.
#'
#' @param dds DESeq2 dds file
#' @param norm character string, one of "vst" or "rlog". Type of count
#' normalization to be performed. default: vst.
#' @param ... other parameters to the plot_pca function, see ?plot_pca for more
#' information
#'
#' @return ggplot PCA plot
#' @export
#'
#' @examples
plot_pca_deseq <- function(dds,
                           norm = "vst",
                           ...) {

  # Normalized the DESeq data
  if (norm == "vst") {
    dds_norm <- DESeq2::vst(dds, blind = FALSE)
  } else if (norm == "rlog") {
    dds_norm <- DESeq2::rlog(dds, blind = FALSE)
  }

  # Calculate PCA, return only values without plot
  pca_data <- DESeq2::plotPCA(
    dds_norm,
    intgroup = c(group),
    returnData = TRUE
  ) %>%
    tibble::rownames_to_column(var = "id")

  # Create numeric values for variance explained to pass to axis labels
  percent_var_labs <- round(100 * attr(pca_data, "percentVar"))

  # Calculate coordinates for segments connecting sample points to group
  # centers
  segments <- pca_data %>%
    dplyr::group_by(!!as.symbol(group)) %>%
    dplyr::summarise(xend = mean(PC1), yend = mean(PC2))

  # Extract PCA data in format for generic PCA plotting function
  pca_data <- pca_data %>%
    dplyr::select(dplyr::matches(paste0("^", group, "$")), id, PC1, PC2) %>%
    dplyr::left_join(segments, by = group)

  # Plot using generic PCA plotting function
  p <- plot_pca(
    pca_data,
    percent_var_labs,
    !!!list(...)
  )
  return(p)
}

#' Draw PCA plot from normalized gene counts
#'
#' @param expression_data data frame of normalized gene expression with samples
#' as columns and an id column for gene names/symbols
#' @param metadata data frame with column of samples names, column for group
#' names and optionally other columns describing the sample features
#' @param ... other parameters to the plot_pca function, see ?plot_pca for more
#' information
#'
#' @return ggplot PCA plot
#' @export
#'
#' @examples
plot_pca_common <- function(expression_data,
                            metadata,
                            ...) {
  # calculate PCA components for samples specified in metadata
  pca <- prcomp(
    t(
      expression_data %>%
        dplyr::select(dplyr::matches(metadata$sample))
    ),
    scale = TRUE
  )

  # Create final PCA data frame to plot
  pca_data <- as.data.frame(pca$x) %>%
    tibble::rownames_to_column("id") %>%
    dplyr::left_join(metadata, by = c("id" = "sample"))

  # Calculate coordinates for segments connecting sample points to group
  # centers
  segments <- pca_data %>%
    dplyr::group_by(!!as.symbol(group)) %>%
    dplyr::summarise(xend = mean(PC1), yend = mean(PC2))

  # Join with the PCA data plot
  pca_data <- pca_data %>%
    dplyr::left_join(segments, by = group)

  # Create numeric values for variance explained to pass to axis labels
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

  # Plot using generic PCA plotting function
  p <- plot_pca(
    pca_data,
    percent_var_labs,
    !!!list(...)
  )
  return(p)
}

#' Generic PCA plotting function
#'
#' @param pca_data data frame of PCA data, with required PC1 and PC2 columns
#' @param percent_var_labs double vector of length 2, numeric values of
#' explained variance for the PCA plot axes
#' @param group
#' @param plot_center logical, whether to plot the center of group in the PCA
#' plot
#' @param linetype character string, a valid linetype for ggplot2. If the
#' plot_center is set to TRUE, this will determine the line type of drawn
#' connections between sample points to group center. default: solid
#' be
#' @param palette character vector of the same length as the number of groups in
#' the data. Custom color palette for groups. If NULL, viridis color palette
#' will be used. default: NULL
#'
#' @return ggplot PCA plot
#' @export
#'
#' @examples
plot_pca <- function(pca_data,
                     percent_var_labs,
                     group = "group",
                     plot_center = TRUE,
                     linetype = "solid",
                     palette = NULL) {

  # get number of colors needed to display
  no_colors <- pca_data[, group] %>%
    unique() %>%
    length()

  # Create basic PCA plot
  p <- pca_data %>%
    ggplot2::ggplot(ggplot2::aes(
      PC1,
      PC2,
      color = !!as.symbol(group)
    )) +
    ggplot2::geom_point(size = 3)

  # Add group center and segments if specified by the user
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

  # Add axes labels and style the graph
  p <- p +
    ggplot2::xlab(paste0("PC1: ", percent_var_labs[1], "% variance")) +
    ggplot2::ylab(paste0("PC2: ", percent_var_labs[2], "% variance")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")

  # Apply custom palette
  if (is.null(palette)) {
    p <- p +
      ggplot2::scale_color_manual(values = viridis::viridis(no_colors + 1))
  } else {
    p <- p +
      ggplot2::scale_color_manual(values = palette)
  }
  return(p)
}
