library(magrittr, include.only = "%>%")

plot_heatmap <- function(expression_data,
                         metadata,
                         gene_list = NULL,
                         include_groups = NULL,
                         geneid_colname = SYMBOL,
                         process_gene_dup_fun = sum,
                         metadata_sample_colname = sample,
                         metadata_group_colname = group,
                         gene_ranking_fun = matrixStats::rowVars,
                         cell_dims = c(10, 1),
                         palette = "RdBu",
                         ...) {
  sample_colname_str <- deparse(substitute(metadata_sample_colname))
  group_colname_str <- deparse(substitute(metadta_group_colname))
  geneid_colname_str <- deparse(substitute(geneid_colname))

  if (!is.null(include_groups)) {
    metadata <- metadata %>%
      dplyr::filter({{ metadata_group_colname }} %in% include_groups)
  }

  samples <- metadata[[sample_colname_str]]

  if (!is.null(gene_list)) {
    if (is.character(gene_list)) {
      expression_data <- expression_data %>%
        dplyr::filter(tolower({{ geneid_colname }}) %in% tolower(gene_list))
    } else if (is.numeric(gene_list)) {
      expr_matrix <- expression_data %>%
        dplyr::select(dplyr::matches(samples)) %>%
        as.matrix()
      expression_data <- expression_data[order(gene_ranking_fun(expr_matrix),
        decreasing = TRUE
      ), ][1:gene_list, ]
    } else {
      stop("Unrecognized type of argument for gene list.")
    }
  }

  if (nrow(expression_data) == 0) {
    stop("Filtered expression dataset is empty.")
  }

  expression_data <- expression_data %>%
    dplyr::select({{ geneid_colname }}, dplyr::matches(samples)) %>%
    dplyr::group_by({{ geneid_colname }}) %>%
    dplyr::summarise(
      across(dplyr::matches(samples),
      process_gene_dup_fun,
      na.rm = TRUE
    )) %>%
    tibble::column_to_rownames(geneid_colname_str)

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
# (...) extra params to ComplexHeatmap
plot_heatmap_fc <- function(expression_data,
                            diffexp_data,
                            metadata,
                            gene_list,
                            include_groups = NULL,
                            id_colname = SYMBOL,
                            fc_colname = log2FoldChange,
                            .fc_colors = c("darkred", "steelblue"),
                            padj_colname = padj,
                            .pval_colors = c("white", "steelblue"),
                            metadata_sample_colname = sample,
                            metadata_group_colname = group,
                            .col_annot_colors = NULL,
                            ...) {
  group_colname_str <- deparse(substitute(metadata_group_colname))
  ctr_id <- deparse(rlang::enexpr(id_colname))
  if (!is.null(include_groups)) {
    metadata <- metadata %>%
      dplyr::filter({{ metadata_group_colname }} %in% include_groups)
  }
  samples <- metadata %>%
    dplyr::pull({{ metadata_sample_colname }})


  expression_data_fil <- expression_data %>%
    dplyr::arrange({{ id_colname }}) %>%
    dplyr::filter(tolower({{ id_colname }}) %in% tolower(gene_list)) %>%
    dplyr::select({{ id_colname }}, dplyr::matches(samples)) %>%
    tibble::column_to_rownames(var = ctr_id)

  diffexp_data_fil <- diffexp_data %>%
    dplyr::arrange({{ id_colname }}) %>%
    dplyr::filter(tolower({{ id_colname }}) %in% tolower(gene_list))

  log2fc_vals <- diffexp_data_fil %>%
    dplyr::pull({{ fc_colname }})
  log2fc_colors <- ifelse(log2fc_vals < 0, .fc_colors[2], .fc_colors[1])

  # min double is added to pval, in cas of pval = 0 it becomes Inf
  pvals <- diffexp_data_fil %>%
    dplyr::mutate(log_pval = -log10({{ padj_colname }} +
      .Machine$double.xmin)) %>%
    dplyr::pull(log_pval)
  pvals_colors <- circlize::colorRamp2(
    breaks = c(-log10(0.05), max(pvals)),
    colors = .pval_colors
  )
  pvalues_legend <- ComplexHeatmap::Legend(
    col_fun = pvals_colors,
    title = "-log10(p-value)"
  )

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

  if (is.null(.col_annot_colors)) {
    selected_groups <- c(group_colname_str)
    col_colors <- gg_color_hue(length(unique(metadata[[selected_groups]]))) %>%
      setNames(unique(metadata[[selected_groups]]))
    col_list <- list(group_colname_str = col_colors)
  } else {
    selected_groups <- names(.col_annot_colors)
    col_list <- .col_annot_colors
  }

  hat <- ComplexHeatmap::columnAnnotation(
    df = metadata %>%
      dplyr::select(dplyr::contains(selected_groups)) %>%
      as.list(),
    col = col_list,
    border = TRUE
  )

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
