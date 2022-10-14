## supplied by user, contains paths to data and metadata files, will contain
## some of the following variables
# path_to_metadata
# path_to_sample_expression
# path_to_diffexp_data
# path_to_deseq_dds
# path_to_cuffdiff
# path_to_dire
# path_to_pathway_files
# path_to_plot_export_params (json format)
# path_to_output_directory

path_list <- rjson::fromJSON(file = "paths.json")

input_directory <- path_list$input_directory %>%
    dir_exists() %>%
    append_dir_slash()

output_directory <- path_list$output_directory %>%
    dir_exists() %>%
    append_dir_slash()

purrr::walk2(
    c(
        names(path_list$input_data_paths),
        names(path_list$param_paths)
    ),
    c(
        path_list$input_data_paths,
        path_list$param_paths
    ), function(x, y) {
        assign(x, y, envir = .GlobalEnv)
    }
)
