path_list <- rjson::fromJSON(file = "paths.json")

input_directory <- path_list$input_directory %>%
  dir_exists() %>%
  append_dir_slash()

output_directory <- path_list$output_directory %>%
  dir_exists() %>%
  append_dir_slash()

# load paths for input files
purrr::walk2(
  names(path_list$input_data_paths), path_list$input_data_paths,
  function(x, y, input_directory) {
    assign(x, paste0(input_directory, y), envir = .GlobalEnv)
  }, input_directory = input_directory
)

# load paths for params
purrr::walk2(
  names(path_list$param_paths), path_list$param_paths,
  function(x, y) {
    assign(x, y, envir = .GlobalEnv)
  }
)
