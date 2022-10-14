# load required packages
source("scripts/load_packages.R", local = TRUE)

# load package functions
source("scripts/load_scripts.R", local = TRUE)

# parses paths.json file for locations of user supplied data and parameters
source("scripts/load_paths.R", local = TRUE)

# load data and metadata
source("scripts/load_data.R", local = TRUE)

# pathway_source <- "ConsensusPathDB_HumanCyc_GSEA_all"

### test functions -------------------------------------------------------------
# pca plots
plot_pca_deseq(dds)
ggsave_param(
  path_to_output_directory,
  get_export_params("pca", path_to_plot_export_params),
  filename_suffix = "_deseq"
)

plot_pca_common(sample_expression, metadata)
ggsave_param_wrapper("pca")

# volcano plots
# TODO add call to volcano plot with parameters
volcano_plot(diffexp_data)

# heatmaps
plot_heatmap_all(expression_data, metadata)
plot_heatmap_topn(expression_data, metadata, n = 1000)
plot_heatmap_diffexp(expression_data, diffexp_data, metadata, palette = "RdYlBu")
plot_heatmap_fc(expression_data, diffexp_data, metadata, gene_list)

# TF dire (dcode) plots
plot_dire(dire)
plot_dire_labeled(dire, 0.05, 0.05)

# pathway plots
pathways_df <- collate_pathways(path_to_pathway_files, pattern = "_contrast")
plot_pathways_meta(pathways_df, top_pathways = 30)
plot_pathway_bargraph(pathways_df, pathway_source, top_n = 20, truncate_desc = 80)

# compare two datasets
plot_corr()
plot_lfc_scatter()
plot_venn2()
