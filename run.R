# load required packages
source("load_packages.R", local = TRUE)

# load package functions
source("load_scripts.R", local = TRUE)

## supplied by user, contains paths to data and metadata files, will contain
## some of the following variables
# path_to_metadata
# path_to_sample_expression
# path_to_diffexp_data
# path_to_deseq_dds
# path_to_cuffdiff
# path_to_dire
# path_to_pathway_files
# path_to_plot_export_params (json format)
# path_to_output_directory
source("load_paths.R", local = TRUE)

# load data and metadata
source("load_data.R")

# load gene list, pathway lists etc.
gene_list <- c("Ppargc1a", "Mb", "Myog", "Mstn", "ND5", "Cyc1", "Sdha", "Atp5a1")
gene_labels <- c("Mb", "Mstn", "Cyc1", "Sln", "Myh3")

pathway_source <- "ConsensusPathDB_HumanCyc_GSEA_all"

### test functions -------------------------------------------------------------
# pca plots
plot_pca_deseq(dds)
ggsave_param(
  path_to_output_directory,
  get_export_params("pca", path_to_plot_export_params),
  filename_suffix = "_deseq"
)
plot_pca_common(sample_expression, metadata)
ggsave_param(
  path_to_output_directory,
  get_export_params("pca", path_to_plot_export_params)
)

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
