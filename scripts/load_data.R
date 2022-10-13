metadata <- read.table(path_to_metadata, header = TRUE) %>%
  select(sample, group) %>%
  unique() %>%
  arrange(rev(group))
rownames(metadata) <- NULL

sample_expression <- read_delim(
  path_to_sample_expression,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

diffexp_data <- read_delim(
  path_to_diffexp_data,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

dds <- readRDS(path_to_deseq_dds)
cufdiff_diff <- read_cuffdiff_diff(path_to_cuffdiff)
dire <- read_dire(path_to_dire, "dire_all")