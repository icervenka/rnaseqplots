plot_heatmap <- function(expression_df, metadata,
                         gene_list = NULL,
                         geneid_colname = SYMBOL,
                         metadata_sample_colname = sample,
                         gene_ranking_fun = matrixStats::rowVars,
                         cell_dims = c(10, 1),
                         palette = "RdBu",
                         ...) {
  samples <- ds(metadata[[metadata_sample_colname]])

  if (!is.null(gene_list)) {
    if (is.character(gene_list)) {
      expression_df <- expression_df %>%
        dplyr::filter({{ geneid_colname }} %in% gene_list)
    } else if (is.numeric(gene_list)) {
      expr_matrix <- expression_df %>%
        dplyr::select(dplyr::matches(samples)) %>%
        as.matrix()
      expression_df <- expression_df[order(gene_ranking_fun(expr_matrix),
        decreasing = T
      ), ][1:gene_list, ]
    } else {
      stop("Unrecognized type of argument for gene list.")
    }
  }

  expression_df <- expression_df %>%
    dplyr::select(dplyr::matches(samples))

  p <- expression_df %>%
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
        tibble::column_to_rownames(sc),
      ...
    )
  return(p)
}

# ISSUE doesn't plot heatmap on mac rstudio
plot_heatmap_fc <- function(expression_data, diffexp_data, metadata, gene_list,
                            id_colname = SYMBOL, fc_colname = log2FoldChange,
                            .fc_colors = c("darkred", "steelblue"),
                            padj_colname = padj, .pval_colors = c("white", "steelblue"),
                            metadata_sample_colname = sample, .col_annot_colors = NULL, # named vector, column annotation defaults to 'group'
                            ) {
  samples = metadata[[ds(metadata_sample_colname)]]
  
  expression_data_fil <- expression_data %>%
    dplyr::arrange({{ id_colname }}) %>%
    dplyr::filter({{ id_colname }} %in% gene_list) %>%
    dplyr::select({{ id_colname }}, dplyr::matches(samples)) %>%
    tibble::column_to_rownames(ds(subid_colname))

  diffexp_data_fil <- diffexp_data %>%
    dplyr::arrange({{ id_colname }}) %>%
    dplyr::filter({{ id_colname }} %in% gene_list)

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
  
  if(is.null(.col_annot_colors)) {
    selected_gropus <- c("group")
    col_colors <- gg_color_hue(length(unique(metadata[[selected_groups]]))) %>%
      setNames(unique(metadata[[selected_groups]]))
    col_list = list("group" = col_colors)
  } else {
    selected_gropus = names(.col_annot_colors)
    col_list = .col_annot_colors
  }
  
  hat <- ComplexHeatmap::columnAnnotation(
    # TODO check if this works properly
    df = metadata %>% select(contains(selected_groups)),
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
    width = ncol(expression_data_fil) * unit(5, "mm"),
    height = nrow(expression_data_fil) * unit(5, "mm")
  )

  ComplexHeatmap::draw(
    ht,
    annotation_legend_list = list(pvalues_legend),
    merge_legends = TRUE
  )
}
