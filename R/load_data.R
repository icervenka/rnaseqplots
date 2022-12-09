#TODO add amigo and ipa data load

load_json_config <- function(config, into = "list") {
  config_list <- rjson::fromJSON(file = config)

  if (into == "list") {
    out <- config_list
  } else if (into == "env") {
    out <- new.env()
    list2env(config_list, envir = out)
  }
  return(out)
}

# TODO reorganize due to params
load_metadata <- function(path, params) {
  metadata <- read_data(path) %>%
    # dplyr::select(-dplyr::any_of(c("fq", "lane", "read"))) %>%
    dplyr::select(all_of(c(params$sample_colname, params$group_colname))) %>%
    base::unique() %>%
    dplyr::mutate(!!as.symbol(params$group_colname) := factor(group,
      levels = params$group_levels
    )) %>%
    dplyr::arrange(!!as.symbol(params$group_colname))
  rownames(metadata) <- NULL
  return(metadata)
}

load_expression_data <- function(path, params) {
  expression_data <- read_data(path)
  return(expression_data)
}

load_diffexp_data <- function(path, params) {
  diffexp_data <- purrr::map(path, function(x) {
    if (nchar(x) > 0) {
      read_data(x)
    }
  }) %>%
    setNames(names(path))
  return(diffexp_data)
}

load_dire <- function(path, params) {
  dire <- collate_dire_pathways(
    path,
    params$dire_pattern,
    params$dire_sheet_name
  )
  return(dire)
}

load_deseq <- function(path, params) {
  dds <- readRDS(path)
}

load_cuffdiff <- function(path, params) {
  cufdiff_diff <- read_cuffdiff(append_dir_slash(path))
  return(cuffdiff_diff)
}

load_cp_pathways <- function(path, params) {
  cp_pathways <- collate_cp_pathways(
    append_dir_slash(path),
    pattern = params$cp_pathways_pattern
  )
  return(cp_pathways)
}

load_gsea_pathways <- function(path, params) {
  gsea_data <- batch_read_filter_gsea(
    dir = append_dir_slash(path),
    pvalue_threshold = params$pvalue_threshold
  )
  return(gsea_data)
}

load_ipa <- function(path) {

}

load_gene_lists <- function(path, params) {
  if (nchar(params$modify_gene_lists) == 0) {
    modify_gene_lists <- identity
  } else {
    modify_gene_lists <- get(params$modify_gene_lists)
  }

  gene_lists <- rjson::fromJSON(file = path)
  gene_lists_to_read <- gene_lists[purrr::map_int(gene_lists, length) == 1]
  gene_lists_updated <- purrr::map(
    gene_lists_to_read,
    read_gene_list_from_file
  ) %>%
    setNames(names(gene_lists_to_read))
  gene_lists[names(gene_lists_updated)] <- NULL
  gene_lists <- append(gene_lists, gene_lists_updated) %>%
    purrr::map(modify_gene_lists)
  return(gene_lists)
}

load_pathway_lists <- function(path, params) {
  pathway_lists <- rjson::fromJSON(file = path)
}

load_function_map <- list(
  "path_to_metadata" = load_metadata,
  "path_to_expression_data" = load_expression_data,
  "path_to_diffexp_data" = load_diffexp_data,
  "path_to_deseq_dds" = load_deseq,
  "path_to_cuffdiff" = load_cuffdiff,
  "path_to_cp_pathway_files" = load_cp_pathways,
  "path_to_gsea_files" = load_gsea_pathways,
  "path_to_ipa_files" = load_ipa,
 # "path_to_amigo_files" = load_amigo,
  "path_to_dire_files" = load_dire,
  "path_to_gene_lists" = load_gene_lists,
  "path_to_pathway_lists" = load_pathway_lists
)

load_data <- function(input_config,
                      param_config,
                      load_function_map,
                      into = "list") {

  if (is.character(input_config) && endsWith("json")) {
    paths <- load_json_config(input_config)
  } else if (is.environment(input_config)) {
    paths <- as.list(input_config)
  } else if (is.list(input_config)) {
    paths <- input_config
  } else {
    stop("Unrecognized config type.")
  }

  data <- list()
  purrr::walk2(names(input_config), input_config, function(x, y) {
    if (nchar(y) > 0) {
      data[[x]] <- load_function_map$x(y, param_config)
    }
  })

  if (into == "list") {
    out <- data()
  } else {
    list2env(data, envir = out)
  }
  return(out)
}
