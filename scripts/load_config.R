# TODO load all config in a walk function
config_path <- "config.json"

# TODO rename
config <- rjson::fromJSON(file = config_path)[["plot_export_params"]] %>%
  purrr::map_dfr(data.frame)

group_levels <- rjson::fromJSON(file = config_path)[["group_levels"]]

gsea_fdr_cutoff <- rjson::fromJSON(file = config_path)[["gsea_fdr_cutoff "]]