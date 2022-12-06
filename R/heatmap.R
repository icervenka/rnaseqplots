library(magrittr, include.only = "%>%")

#' Draws a heatmap of gene expression
#'
#' Uses pheatmap package.
#'
#' @param expression_data data frame of normalized gene expression with samples
#' as columns and an id column for gene names/symbols
#' @param metadata data frame with column of samples names, column for group
#' names and optionally other columns describing the sample features
#' @param gene_list can be one of:
#' - character vector of gene IDs/symbols to display
#' - integer, top ranking genes based on gene_ranking_fun will be displayed
#' - NULL, all genes will be displayed
#' @param include_groups character vector, which groups from metadata file to
#' include. If NULL, all groups will be included. default: NULL
#' @param id_colname column name of gene IDs/symbols in expression data file,
#' supplied as variable. default: SYMBOL
#' @param process_gene_dup_fun function, how to treat the expression data for
#'  potential gene ID/symbol duplicates. default: sum
#' @param sample_colname column name of samples in metadata file, supplied as
#' variable. default: sample
#' @param group_colname column name of experimental groups in metadata file,
#' supplied as variable. default: group
#' @param gene_ranking_fun function, in case of displaying top genes, how should
#' the genes in the matrix be ranked. default: matrixStats::rowVars
#' @param cell_dims interger vector of length 2 specifying width and height of
#' individual heatmap cells. default: c(10, 1)
#' @param palette character, which colorbrewer palette to use. default: "RdBu"
#' @param ... other parameters to pheatmap function
#'
#' @return pheatmap plot
#' @export
#'
#' @examples
plot_heatmap <- function(expression_data,
                         metadata,
                         gene_list = NULL,
                         include_groups = NULL,
                         id_colname = SYMBOL,
                         process_gene_dup_fun = sum,
                         sample_colname = sample,
                         group_colname = group,
                         gene_ranking_fun = matrixStats::rowVars,
                         cell_dims = c(10, 1),
                         palette = "RdBu",
                         ...) {
  # Convert column name variables to string, needed for certain functions
  sample_colname_str <- deparse(substitute(sample_colname))
  group_colname_str <- deparse(substitute(metadta_group_colname))
  id_colname_str <- deparse(substitute(id_colname))
  # Filter which groups to include based on metadata
  if (!is.null(include_groups)) {
    metadata <- metadata %>%
      dplyr::filter({{ group_colname }} %in% include_groups)
  }
  # Filter which samples to include based on metadata
  samples <- metadata[[sample_colname_str]]

  # Filter which genes to include
  if (!is.null(gene_list)) {
    # If character list, filter based on gene ID
    if (is.character(gene_list)) {
      expression_data <- expression_data %>%
        dplyr::filter(tolower({{ id_colname }}) %in% tolower(gene_list))
      # If integer, select top genes
    } else if (is.numeric(gene_list)) {
      # Create matrix
      expr_matrix <- expression_data %>%
        dplyr::select(dplyr::matches(samples)) %>%
        as.matrix()
      # Order matrix based on supplied function and select top genes
      expression_data <- expression_data[order(gene_ranking_fun(expr_matrix),
        decreasing = TRUE
      ), ][1:gene_list, ]
    } else {
      stop("Unrecognized type of argument for gene list.")
    }
  }

  # Stop if no genes are present after filtering
  if (nrow(expression_data) == 0) {
    stop("Filtered expression dataset is empty.")
  }

  # Process duplicated gene IDs/symbols. Required because row names of matrix
  # create later need to be unique
  expression_data <- expression_data %>%
    dplyr::select({{ id_colname }}, dplyr::matches(samples)) %>%
    dplyr::group_by({{ id_colname }}) %>%
    dplyr::summarise(
      dplyr::across(dplyr::matches(samples),
        process_gene_dup_fun,
        na.rm = TRUE
      )
    ) %>%
    tibble::column_to_rownames(id_colname_str)

  # Create heatmap
  p <- expression_data %>%
    pheatmap::pheatmap(
      scale = "row",
      color = colorRampPalette(S4Vectors::rev(RColorBrewer::brewer.pal(
        n = 7, name = palette
      )))(100),
      cellwidth = cell_dims[1],
      cellheight = cell_dims[2],
      border_color = "white",
      treeheight_row = 15,
      treeheight_col = 20,
      annotation_col = metadata %>%
        tibble::column_to_rownames(sample_colname_str),
      ...
    )
  return(p)
}

# named vector, column annotation defaults to 'group'

#' Creates a heatmap combining expression and differential expression data
#'
#' Heatmap will have a plot for p-values and log2 fold-changes of selected
#' comparison on the right side. Uses ComplexHeatmap.
#'
#' @param expression_data data frame of normalized gene expression with samples
#' as columns and an id column for gene names/symbols
#' @param diffexp_data data frame of gene differential expression, must contain
#' column for gene IDs/symbols, log2 fold-changes and a p-value
#' @param metadata data frame with column of samples names, column for group
#' names and optionally other columns describing the sample features
#' @param gene_list character vector of gene IDs/symbols to display
#' @param include_groups character vector, which groups from metadata file to
#' include. If NULL, all groups will be included. default: NULL
#' @param id_colname column name of gene IDs/symbols in expression data file,
#' supplied as variable. default: SYMBOL
#' @param fc_colname column name of log fold-changes in the data frame with
#' differential expression data, supplied as variable. default: log2FolcChange
#' @param .fc_colors character vector of length 2, colors that will be mapped to
#' highest and lowest log fold-changes respectively. default:
#' c("darkred", "steelblue")
#' @param pval_colname column name of p-values or adjusted p-values in the
#' data frame with differential expression data, supplied as variable.
#' default: padj
#' @param pval_threshold double, significance threshold for p-values, used for 
#' mapping the lower value from .pval_colors to this value. default: 0.05
#' @param .pval_colors character vector of length 2, colors that will be mapped
#' to highest and lowest p-values respectively. default: c("white", "steelblue")
#' @param sample_colname column name of samples in metadata file, supplied as
#' variable. default: sample
#' @param group_colname column name of experimental groups in metadata file,
#' supplied as variable. default: group
#' @param .col_annot_colors names character vector of length corresponding to
#' the number of displayed groups. Names of vector correspond group names and
#' values to their annotation colors. If NULL, standard ggplot color scheme
#' will be used. default: NULL
#' @param ... extra params to ComplexHeatmap
#'
#' @return ComplexHeatmap heatmap
#' @export
#'
#' @examples
plot_heatmap_fc <- function(expression_data,
                            diffexp_data,
                            metadata,
                            gene_list,
                            include_groups = NULL,
                            id_colname = SYMBOL,
                            fc_colname = log2FoldChange,
                            .fc_colors = c("darkred", "steelblue"),
                            pval_colname = padj,
                            pval_threshold = 0.05,
                            .pval_colors = c("white", "steelblue"),
                            sample_colname = sample,
                            group_colname = group,
                            .col_annot_colors = NULL,
                            ...) {
  # Convert column name variables to string, needed for certain functions
  group_colname_str <- deparse(substitute(group_colname))
  ctr_id <- deparse(rlang::enexpr(id_colname))
  # Filter which groups to include based on metadata
  if (!is.null(include_groups)) {
    metadata <- metadata %>%
      dplyr::filter({{ group_colname }} %in% include_groups)
  }
  # Filter which samples to include based on metadata
  samples <- metadata %>%
    dplyr::pull({{ sample_colname }})

  # Filter which genes to include from expression data file
  expression_data_fil <- expression_data %>%
    dplyr::arrange({{ id_colname }}) %>%
    dplyr::filter(tolower({{ id_colname }}) %in% tolower(gene_list)) %>%
    dplyr::select({{ id_colname }}, dplyr::matches(samples)) %>%
    tibble::column_to_rownames(var = ctr_id)
  # Filter which genes to include from differential expression data file
  diffexp_data_fil <- diffexp_data %>%
    dplyr::arrange({{ id_colname }}) %>%
    dplyr::filter(tolower({{ id_colname }}) %in% tolower(gene_list))
  # Complexheatmap requires the row annotation to be a vector of values
  log2fc_vals <- diffexp_data_fil %>%
    dplyr::pull({{ fc_colname }})
  # Assign corresponding colors to values
  log2fc_colors <- ifelse(log2fc_vals < 0, .fc_colors[2], .fc_colors[1])

  # Complexheatmap requires the row annotation to be a vector of values
  pvals <- diffexp_data_fil %>%
    # min double is added to pval, in cas of pval = 0 it becomes Inf
    dplyr::mutate(log_pval = -log10({{ pval_colname }} +
      .Machine$double.xmin)) %>%
    dplyr::pull(log_pval)
  # Assign corresponding colors to values
  pvals_colors <- circlize::colorRamp2(
    breaks = c(-log10(pval_threshold), max(pvals)),
    colors = .pval_colors
  )
  # Create ComplexHeatmap legend for pvalues
  pvalues_legend <- ComplexHeatmap::Legend(
    col_fun = pvals_colors,
    title = "-log10(p-value)"
  )

  # Combine log fold-changes and p-values into row annotation with colors
  har <- ComplexHeatmap::rowAnnotation(
    pvalue = ComplexHeatmap::anno_simple(pvals,
      col = pvals_colors,
      gp = grid::gpar(col = "black", lwd = 1)
    ),
    log2fc = ComplexHeatmap::anno_barplot(log2fc_vals,
      baseline = 0,
      bar_width = 0.9,
      gp = grid::gpar(fill = log2fc_colors, col = "white")
    ),
    simple_anno_size = grid::unit(0.5, "cm"), width = grid::unit(2, "cm"),
    gap = grid::unit(2, "mm")
  )

  # Create color annotations for groups in the columns
  # If NULL, standard ggplot color scheme will be used, or user supplied colors
  # otherwise
  if (is.null(.col_annot_colors)) {
    selected_groups <- c(group_colname_str)
    col_colors <- metadata[[selected_groups]] %>%
      unique() %>%
      length() %>%
      gg_color_hue() %>%
      setNames(unique(metadata[[selected_groups]]))
    col_list <- list(group_colname_str = col_colors)
  } else {
    selected_groups <- names(.col_annot_colors)
    col_list <- .col_annot_colors
  }

  # Create column annotation with colors
  hat <- ComplexHeatmap::columnAnnotation(
    df = metadata %>%
      dplyr::select(dplyr::contains(selected_groups)) %>%
      as.list(),
    col = col_list,
    border = TRUE
  )

  # Create final heatmap
  ht <- ComplexHeatmap::Heatmap(
    expression_data_fil %>%
      as.matrix() %>%
      t() %>%
      scale() %>%
      t(),
    name = "expression",
    right_annotation = har,
    top_annotation = hat,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    border_gp = grid::gpar(col = "black", lty = 1),
    rect_gp = grid::gpar(col = "black", lwd = 1),
    width = ncol(expression_data_fil) * grid::unit(5, "mm"),
    height = nrow(expression_data_fil) * grid::unit(5, "mm"),
    ...
  )

  ComplexHeatmap::draw(
    ht,
    annotation_legend_list = list(pvalues_legend),
    merge_legends = TRUE
  )
}
