pathway_lists <- load_pathway_lists(
    "inst/extdata/lists/pathway_lists.json",
    list()
)

usethis::use_data(pathway_lists, overwrite = TRUE)
