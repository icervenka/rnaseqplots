config_path <- "config.json"
config <- rjson::fromJSON(file = config_path)

plot_params <- config[["plot_export_params"]] %>%
  purrr::map_dfr(data.frame)

config_items <- c(
  "group_levels",
  "gsea_fdr_cutoff",
  "dire_pattern",
  "cp_pathways_pattern",
  "dire_sheet_name"
)
purrr::walk(config_items, function(x, config) {
  if (validate_config(x, config)) {
    assign(x, config[[x]], envir = .GlobalEnv)
  }
}, config = config)
