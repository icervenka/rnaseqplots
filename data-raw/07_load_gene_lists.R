gene_lists <- load_gene_lists(
    "../inst/extdata/rnaseq/lists/gene_lists.json",
    list()
)

usethis::use_data(gene_lists, overwrite = TRUE)