# TODO fix factor ordering for group
metadata <- read_data(path_to_metadata) %>%
  dplyr::select(sample, group) %>%
  base::unique()
rownames(metadata) <- NULL

sample_expression <- read_data(path_to_sample_expression)
diffexp_data <- read_data(path_to_diffexp_data)

dds <- readRDS(path_to_deseq_dds)
cufdiff_diff <- read_cuffdiff_diff(path_to_cuffdiff)
dire <- read_dire(path_to_dire, dire_sheet_name)