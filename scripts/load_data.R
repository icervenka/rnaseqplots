# load metadata
if (nchar(path_to_metadata) > 0) {
  metadata <- read_data(path_to_metadata) %>%
    dplyr::select(sample, group) %>%
    base::unique() %>%
    dplyr::mutate(group = factor(group, levels = group_levels)) %>%
    dplyr::arrange(group)
  rownames(metadata) <- NULL
} else {
  warning("Path to experiment metadata is missing.")
}

# load expression data
if (nchar(path_to_expression_data) > 0) {
  expression_data <- read_data(path_to_expression_data)
} else {
  warning("Path to experiment expression data is missing.")
}

# load differential expression data
if (nchar(path_to_diffexp_data) > 0) {
  diffexp_data <- read_data(path_to_diffexp_data)
} else {
  warning("Path to experiment expression data containing differential
  gene expression is missing.")
}

# load DESeq2 dds object
if (nchar(path_to_deseq_dds) > 0) {
  dds <- readRDS(path_to_deseq_dds)
}

# load cuffdiff object
if (nchar(path_to_cuffdiff) > 0) {
  cufdiff_diff <- read_cuffdiff_diff(path_to_cuffdiff)
}

# load dire analysis
if (nchar(path_to_dire) > 0) {
  dire <- read_dire_xlsx(path_to_dire, dire_sheet_name)
}

# load gene list json file
if (nchar(path_to_gene_lists) > 0) {
  gene_lists <- rjson::fromJSON(file = path_to_gene_lists)
}

# load pathway ranks json file
if (nchar(path_to_pathway_ranks) > 0) {
  pathway_ranks <- rjson::fromJSON(file = path_to_pathway_ranks)
}

# load gsea pathway data
if (nchar(path_to_gsea_files) > 0) {
  gsea_data <- batch_read_filter_gsea(
    dir = append_dir_slash(path_to_gsea_files),
    pvalue_threshold = 1
  )
}
