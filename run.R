# run master loading script that loads all the necessary, packages, functions
# and data
# TODO think about moving config and path scripts to main run script
source("scripts/load_all.R", local = TRUE)

############################# test functions ###################################
## pca plots -------------------------------------------------------------------
# A)
plot_pca_deseq(dds)
ggsave_param(
  output_directory,
  get_plot_params("pca", "config.json"),
  filename_suffix = "_deseq"
)

# B)
# not exported
plot_pca_common(expression_data, metadata)

# C)
# not exported
plot_pca_common(
  expression_data,
  metadata,
  group = "group",
  plot_center = TRUE,
  linetype = "dashed",
  palette = viridis::viridis(5)
)

## ma plot ---------------------------------------------------------------------
# A)
ma_plot(
  diffexp_data$ex_1,
  label_top_n = 5,
  label_bottom_n = 5
)
ggsave_param_wrapper("ma")

## volcano plots ---------------------------------------------------------------
# A)
volcano_plot(
  diffexp_data$ex_1,
  label_bottom_n = 5,
  label_top_n = 6
)
ggsave_param(
  output_directory,
  list(
    filename = "volcano",
    device = "png",
    dpi = 600,
    unit = "cm",
    width = 10,
    height = 10
  ),
  filename_suffix = "_top"
)

# B)
# not exported
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

## heatmaps  -------------------------------------------------------------------
# pheatmap and ComplexHeatmap package can't use last_plot() function,
# created plot has to be piped to the ggsave_param(_wrapper) functions
# A)
# not exported
plot_heatmap(
  expression_data,
  metadata,
  gene_list = 1000,
  show_rownames = FALSE
)

# B)
plot_heatmap(
  expression_data,
  metadata,
  gene_list = gene_lists[["oxphos"]],
  geneid_colname = SYMBOL,
  metadata_sample_colname = sample,
  gene_ranking_fun = rowMeans,
  cell_dims = c(9, 9),
  palette = "PuOr"
) %>%
ggsave_param_wrapper("heatmap", plot = .)

# C)
# not exported
plot_heatmap_fc(
  expression_data,
  diffexp_data$ex_1 %>% dplyr::mutate(padj = tidyr::replace_na(padj, 1)),
  metadata,
  gene_lists$oxphos
)

# D
plot_heatmap_fc(
  expression_data,
  diffexp_data$ex_2,
  metadata,
  gene_lists[["heatmap_fc_example_2"]],
  include_groups = c("ctrl", "ne", "cgrp"),
  id_colname = SYMBOL,
  fc_colname = log2FoldChange,
  .fc_colors = c("darkred", "steelblue"),
  padj_colname = padj,
  .pval_colors = c("white", "steelblue"),
  metadata_sample_colname = sample,
  .col_annot_colors = list(group = c(
    "ctrl" = "gray60",
    "ne" = "steelblue",
    "cgrp" = "darkseagreen",
    "npy" = "tomato4",
    "sp" = "lightgoldenrod2"
  ))
) %>%
ggsave_param_wrapper("heatmap_fc", plot = .)

## compare two datasets --------------------------------------------------------
# A)
plot_param_corr(
  expression_data,
  metadata,
  gene_lists[["mito_1"]],
  c("weight"),
  palette = viridis::viridis(5)
)
ggsave_param_wrapper("param_correlation")

# B)
plot_lfc_scatter(
  diffexp_data$ex_1,
  diffexp_data$ex_2,
  color_quadrants = TRUE,
  alpha = 1
)
ggsave_param_wrapper("lfc_scatter")

# C)
## Venn diagram function takes parameters for plot export as one of the
## function arguments
plot_venn(
  list(
    get_sig_genes(diffexp_data$ex_1),
    get_sig_genes(diffexp_data$ex_2)
  ),
  get_plot_params("venn", "config.json"),
  font_size = 1
)

## TF dire (dcode) plots -------------------------------------------------------
# A)
# not exported
plot_dire(dire %>% filter_file("dire_1"))

# B)
plot_dire_labeled(
  dire %>% filter_file("dire_2"),
  occurrence_threshold = 0.10,
  importance_threshold = 0.15,
  combine_thresholds = `|`)
ggsave_param_wrapper("dire")

## gsea plots  -----------------------------------------------------------------
# A)
plot_pathways_rank(
  gsea_data %>% filter_name("HALLMARK"),
  pathway_lists[["example_1"]],
  x_axis = nes,
  y_axis = name,
  bar_fill = fdr_qval,
  label_offset = 0.8
)
ggsave_param_wrapper("pathways_rank")

# B)
plot_pathways_volcano(
  gsea_data %>% filter_name("KEGG"),
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
# A)
plot_cp_pathways_meta(cp_pathways, top_pathways = 30)
ggsave_param_wrapper("pathways_meta")

# B)
plot_cp_pathways_bargraph(
  cp_pathways,
  pathway_source_pattern = "KEGG_GSEA",
  top_pathways = 10,
  truncate_description = 80
)
ggsave_param(
  output_directory,
  get_plot_params("pathways_bargraph", plot_params),
)
