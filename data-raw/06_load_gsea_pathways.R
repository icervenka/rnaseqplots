gsea_pathways <- load_pathways_gsea(
    "inst/extdata/rnaseq/gsea",
    list(gsea_fdr_cutoff = 0.25)
)

usethis::use_data(gsea_pathways, overwrite = TRUE)
