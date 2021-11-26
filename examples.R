# supplied by user
source("paths.R")
source("rnaseq_visualisation_misc.R")

### load data and set global variables -----------------------------------------
metadata = read.table(path_to_metadata, header = T) %>%
  select(sample, group) %>%
  unique() %>%
  arrange(rev(group))
rownames(metadata) = NULL

sample_expression = read_delim(
  path_to_sample_expression,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE)

diffexp_data = read_delim(
  path_to_diffexp_data,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE)

dds = readRDS(path_to_deseq_dds)
cufdiff_diff = read_cuffdiff_diff(path_to_cuffdiff)
dire = read_dire(path_to_dire, "dire_all")

gene_list = c("Ppargc1a", "Mb", "Myog", "Mstn", "ND5", "Cyc1", "Sdha", "Atp5a1"))
gene_labels = c("Mb", "Mstn", "Cyc1", "Sln", "Myh3")

pathway_source = "MeSH_g2p_C@gene2pubmed_GSEA_all"

plot_dimensions_mm = list(
  "pca" = c(100, 100),
  "volcano" = c(100, 100),
  "heatmap" = c(100, 100),
  "heatmap_sample" = c(100, 100),
  "heatmap_fc" = c(100, 100),
  "dire" = c(100, 100),
  "pathways_meta" = c(100, 100),
  "pathways" = c(100, 100),
  "corr" = c(100, 100),
  "fc_scatter" = c(100, 100),
  "venn" = c(100, 100)
)
plot_dpi = 600

### test functions -------------------------------------------------------------
# pca plots
plot_pca_deseq(dds)
plot_pca_common(sample_expression, metadata)

# volcano plots
plot_volcano(diffexp_data)
plot_volcano_labeled(diffexp_data, gene_labels, symbol_colname = "SYMBOL")
#plot_volcano_cuffdiff()

# heatmaps
plot_heatmap(expression_data, metadata)
plot_heatmap_topn(expression_data, metadata, n = 1000)
plot_heatmap_diffexp(expression_data, metadata, pallete = "RdYlBu")
plot_heatmap_fc(expression_data, diffexp_data, metadata, gene_list)

# TF dire (dcode) plots
plot_dire()
plot_dire_labeled()

# pathway plots
pathways_df = collate_pathways(path_to_pathway_files)
plot_pathways_meta(pathways_df, top_pathways = 30)
plot_pathway_bargraph(df, pathway_source, top_n = 20, truncate_desc = 80)

# compare two datasets
plot_corr()
plot_lfc_scatter()
plot_venn2()
