# run master loading script that loads all the necessary, packages, functions
# and data
source("scripts/load_all.R", local = TRUE)

############################# test functions ###################################
## pca plots -------------------------------------------------------------------
# A
plot_pca_deseq(dds)
ggsave_param(
  path_to_output_directory,
  plot_params[["pca"]],
  filename_suffix = "_deseq"
)

# B
plot_pca_common(expression_data, metadata)
ggsave_param_wrapper("pca")

plot_pca_common(
  expression_data,
  metadata,
  group = "group",
  plot_center = TRUE,
  linetype = "dashed",
  palette = c("steelblue", "darkred")
)

## ma plot ---------------------------------------------------------------------
ma_plot(
  diffexp_data$ex_1,
  label_top_n = 8,
  label_bottom_n = 8
)
ggsave_param_wrapper("ma")

## volcano plots ---------------------------------------------------------------
# A
volcano_plot(
  diffexp_data,
  label_bottom_n = 5,
  label_top_n = 10
)
ggsave_param(
  path_to_output_directory,
  list(
    filename = "volcano",
    device = "png",
    dpi = 600,
    unit = "cm",
    width = 16,
    height = 8
  ),
  filename_suffix = "_top"
)

# B
volcano_plot(
  diffexp_data$ex_1,
  label = SYMBOL,
  x = log2FoldChange,
  y = pvalue,
  sig_threshold = 0.1,
  log2fc_threshold = 1,
  filter_sig_on = padj,
  label_genes = gene_lists[["volcano_example"]],
  add_vhlines = TRUE,
  vhline_color = "darkolivegreen4",
  vhline_type = "solid",
  xlab_label = "log2FoldChange",
  ylab_label = "-log10(p-value)",
  color_palette = c("bisque2", "darkorchid4", "cyan4")
)
ggsave_param_wrapper("volcano_plot")

# TODO add cuffdiff volcano plot

## heatmaps  -------------------------------------------------------------------
# A
plot_heatmap(
  expression_data,
  metadata,
  gene_list = 1000,
  show_rownames = FALSE
)

# B
plot_heatmap(
  expression_data,
  metadata,
  gene_list = gene_lists[["metabolism"]],
  geneid_colname = SYMBOL,
  metadata_sample_colname = sample,
  gene_ranking_fun = rowMeans,
  cell_dims = c(9, 9),
  palette = "PuOr"
)
ggsave_param_wrapper("heatmap")

# C
plot_heatmap_fc(
  expression_data,
  diffexp_data$ex_1,
  metadata,
  gene_lists[["tfs"]]
)

# D
plot_heatmap_fc(
  expression_data,
  diffexp_data$ex_1,
  metadata,
  gene_lists[["mito_1"]],
  id_colname = SYMBOL,
  fc_colname = log2FoldChange,
  .fc_colors = c("darkred", "steelblue"),
  padj_colname = padj,
  .pval_colors = c("white", "steelblue"),
  metadata_sample_colname = sample,
  .col_annot_colors = list(group = c(
    "wt" = "gray60",
    "ko" = "steelblue"
  ))
)
ggsave_param_wrapper("heatmap_fc")

## compare two datasets --------------------------------------------------------
# A
plot_param_corr(
  expression_data,
  metadata,
  gene_lists[["mito_1"]],
  c("weight")
)

# B
plot_lfc_scatter(diffexp_data$ex_1, diffexp_data$ex_2)

# C
plot_venn(
  list(
    get_sig_genes(diffexp_data$ex_1),
    get_sig_genes(diffexp_data$ex_1)
  ),
  plot_params[["venn"]],
  font_size = 1
)


## TF dire (dcode) plots -------------------------------------------------------
# A
plot_dire(dire)

# B
plot_dire_labeled(
  dire,
  occurrence_threshold = 0.05,
  importance_threshold = 0.05
)
ggsave_param_wrapper("dire")

## gsea plots  -----------------------------------------------------------------
# TODO write pathway name prettifier
# A
plot_pathways_rank(
  gsea_data %>%
    dplyr::filter(grepl("hallmark", name, ignore.case = TRUE)),
  pathway_lists[["example_1"]],
  x_axis = nes,
  y_axis = name,
  bar_fill = fdr_qval
)
ggsave_param_wrapper("pathways_rank")

# B
plot_pathways_volcano(
  gsea_data %>%
    dplyr::filter(grepl("hallmark", name, ignore.case = TRUE)),
  x_axis = nes,
  y_axis = fdr_qval,
  label_pathways = pathway_lists[["example_volcano"]],
  alpha = 0.8,
  color_palette = c(
    "steelblue4",
    "darkred"
  )
)
ggsave_param_wrapper("pathways_volcano")

## other pathway plots (amigo, clusterprofiler, etc.)  -------------------------
# A
plot_cp_pathways_meta(cp_pathways, top_pathways = 30)

# B
plot_cp_pathways_bargraph(
  cp_pathways,
  pathway_source_pattern = "KEGG_GSEA",
  top_pathways = 10,
  truncate_description = 80
)
ggsave_param(
  path_to_output_directory,
  plot_params[["pathways_bargraph"]],
)
