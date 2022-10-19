# TODO fix factor ordering for group
metadata <- read_data(path_to_metadata) %>%
  dplyr::select(sample, group) %>%
  base::unique() %>%
  dplyr::mutate(group = factor(group, levels = group_levels)) %>%
  dplyr::arrange(group)
rownames(metadata) <- NULL

expression_data <- read_data(path_to_expression_data)
diffexp_data <- read_data(path_to_diffexp_data)

dds <- readRDS(path_to_deseq_dds)
# cufdiff_diff <- read_cuffdiff_diff(path_to_cuffdiff)
dire <- read_dire_xlsx(path_to_dire, dire_sheet_name)

gene_lists <- rjson::fromJSON(file = path_to_gene_lists)