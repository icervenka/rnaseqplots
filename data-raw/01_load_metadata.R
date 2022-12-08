metadata <- load_metadata(
    "../inst/extdata/rnaseq/metadata.tsv",
    list(
        "sample_colname" = "sample",
        "group_colname" = "group"
    )
)

usethis::use_data(metadata, overwrite = TRUE)
