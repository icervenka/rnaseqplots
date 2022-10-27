# load metadata
if (nchar(path_to_metadata) > 0) {
  metadata <- read_data(path_to_metadata) %>%
    dplyr::select(-c(fq, lane, read)) %>%
    #dplyr::select(sample, group) %>%
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
  diffexp_data <- purrr::walk(path_to_diffexp_data, function(x) {
    read_data(x)
  }) %>%
  setNames(names(path_to_diffexp_data))

} else {
  warning("Path to experiment expression data containing differential
  gene expression is missing.")
}

# load DESeq2 dds object
if (nchar(path_to_deseq_dds) > 0) {
  dds <- readRDS(path_to_deseq_dds)
}

# load dire analysis
if (nchar(path_to_dire) > 0) {
  if (endsWith(path_to_dire, "xlsx")) {
    dire <- read_dire_xlsx(path_to_dire, dire_sheet_name)
  } else {
    dire <- read_data(path_to_dire)
  }
}

# load cuffdiff object
if (nchar(path_to_cuffdiff) > 0) {
  cufdiff_diff <- read_cuffdiff(append_dir_slash(path_to_cuffdiff))
}

# load clusterprofiler pathway files
if (nchar(path_to_cp_pathway_files) > 0) {
  cp_pathways <- collate_cp_pathways(
    append_dir_slash(path_to_cp_pathway_files),
    pattern = cp_pathways_pattern)
}

# load gsea pathway data
if (nchar(path_to_gsea_files) > 0) {
  gsea_data <- batch_read_filter_gsea(
    dir = append_dir_slash(path_to_gsea_files),
    pvalue_threshold = 1
  )
}

# load gene list json file
if (nchar(path_to_gene_lists) > 0) {
  gene_lists <- rjson::fromJSON(file = path_to_gene_lists)
}

# load pathway ranks json file
if (nchar(path_to_pathway_lists) > 0) {
  pathway_lists <- rjson::fromJSON(file = path_to_pathway_lists)
}