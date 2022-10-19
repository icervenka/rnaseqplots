# load required packages
source("scripts/load_packages.R", local = TRUE)

# load package functions
source("scripts/load_scripts.R", local = TRUE)

# parses paths.json file for locations of user supplied data
source("scripts/load_paths.R", local = TRUE)

# parses config.json file for parameters
source("scripts/load_config.R", local = TRUE)

# load data and metadata
source("scripts/load_data.R", local = TRUE)

### test functions -------------------------------------------------------------
## pca plots
# A
plot_pca_deseq(dds)
ggsave_param(
  path_to_output_directory,
  get_export_params("pca", config_path),
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

## volcano plots
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
  diffexp_data,
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

## heatmaps
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
  diffexp_data,
  metadata,
  gene_lists[["tfs"]]
)

# D
plot_heatmap_fc(
  expression_data,
  diffexp_data,
  metadata,
  gene_lists[["mito_1"]],
  id_colname = SYMBOL,
  fc_colname = log2FoldChange,
  .fc_colors = c("darkred", "steelblue"),
  padj_colname = padj,
  .pval_colors = c("white", "steelblue"),
  metadata_sample_colname = sample,
  .col_annot_colors = NULL
)
ggsave_param_wrapper("heatmap_fc")

## TF dire (dcode) plots
plot_dire(dire)

plot_dire_labeled(dire, 0.05, 0.05)

## pathway plots
pathways_df <- collate_pathways(path_to_pathway_files, pattern = "_contrast")
plot_pathways_meta(pathways_df, top_pathways = 30)
plot_pathway_bargraph(pathways_df, pathway_source, top_n = 20, truncate_desc = 80)

# gsea plots

# compare two datasets
plot_corr()
plot_lfc_scatter()
plot_venn2()
