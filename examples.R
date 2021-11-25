# supplied by user
source("paths.R")
source("rnaseq_visualisation_misc.R")

# load data
metadata = read.table(path_to_metadata, header = T) %>%
  select(sample, group) %>%
  unique() %>%
  arrange(rev(group))
rownames(metadata) = NULL

sample_expression = read_delim(path_to_sample_expression,
                               delim = "\t", 
                               escape_double = FALSE,
                               trim_ws = TRUE)

diffexp_data = read_delim(path_to_diffexp_data,
                          delim = "\t", 
                          escape_double = FALSE,
                          trim_ws = TRUE)

dds = readRDS(path_to_deseq_dds)

dire_up = read_dire(path_to_dire)

# plot_heatmap_fc(sample_expression, diffexp_data, metadata,
#                 gene_list = c("Ppargc1a", "Mb", "Myog", "Mstn", "ND5", "Cyc1", "Sdha", "Cox5a", "Atp5a1"))

# test functions
plot_pca_deseq()
plot_pca_common()
plot_corr()
plot_lfc_scatter()
plot_venn2()
plot_dire()
plot_dire_labeled()
plot_heatmap()
plot_sample_heatmap()
plot_diffexp_heatmap()
plot_volcano()
plot_volcano_labeled()
plot_heatmap_fc()
collate_pathways()
plot_pathways_meta()
plot_pathway_bargraph()
