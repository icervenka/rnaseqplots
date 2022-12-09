gene_lists <- load_gene_lists(
    "inst/extdata/lists/gene_lists.json",
    list("modify_gene_lists" = "mouseify_gene_symbols")
)

usethis::use_data(gene_lists, overwrite = TRUE)
