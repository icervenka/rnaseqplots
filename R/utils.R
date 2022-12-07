#' Checks is item is present in a configuration dictionary
#'
#' @param entry character string to check
#' @param config_dict named list usually obtained by parsing json configuration
#' file
#'
#' @return TRUE or FALSE depending whether the entry is present as a name in the
#' list
#' @export
#'
#' @examples
validate_config <- function(entry, config_dict) {
  # empty strings are sometimes needed for patterns for file listing
  if (entry %in% names(config_dict)) { # && nchar(config[[item]]) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Adds a slash at the end of path
#'
#' Slash is added only if there is none present.
#'
#' @param dir_string character string representing path to directory
#'
#' @return path too directory with slash as a last character
#'
#' @export
#'
#' @examples
append_dir_slash <- function(dir_string) {
  if (!base::endsWith(dir_string, "/")) {
    dir_string <- paste0(dir_string, "/")
  }
  return(dir_string)
}

#' Check if directory exists
#'
#' @param dir_string character string of path to check
#'
#' @return character string of path if it exists, stops executions otherwise
#'
#' @export
#'
#' @examples
dir_exists <- function(dir_string) {
  if (!file.exists(dir_string)) {
    stop(paste0("Specified directory does not exist: ", dir_string))
  } else {
    return(dir_string)
  }
}

#' Create default ggplot2 color hues
#'
#' @param i integer, number of colors to generate
#'
#' @return character vector with equally spaced color hues
#' @export
#'
#' @examples
gg_color_hue <- function(i) {
  hues <- seq(15, 375, length = i + 1)
  hcl(h = hues, l = 65, c = 100)[1:i]
}

#' Helper switch logical toggle for descending ranking function
#'
#' @param descending logical, whether to rank data frame columns in descending
#' order
#'
#' @return function  for ranking data frame entries
#' @export
#'
#' @examples
rank_how <- function(descending) {
  if (descending == TRUE) {
    return(dplyr::desc)
  } else {
    return(identity)
  }
}

#' Text processing of gene symbols
#'
#' @param gene_list character vector of gene symbols to process
#' @param text_replacements named vector, where names are patterns and values
#' replacement. Will be applied in order it is specified.
#' @param text_modifications list of single argument text processing functions
#' to apply to gene_list. Will be applied from right to left.
#'
#' @return character vector of modified gene names
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
modify_gene_symbols <- function(gene_list,
                                text_replacements,
                                text_modifications) {
  gene_list <- gene_list %>%
    stringr::str_replace_all(regex(text_replacements, ignore_case = TRUE)) %>%
    purrr::map_chr(purrr::compose(!!!text_modifications))
  return(gene_list)
}

#' Replace mouse gene symbols
#'
#' Applies simple text replacements for mitochondrial gene symbols  and converts
#' gene symbols to title case. This is a more specific version of
#' modify_gene_symbols function and a fast way to replace human gene symbols
#' with reasonable accuracy.
#'
#' @param gene_list character vector of gene symbols to process
#'
#' @return character vector of modified gene names
#' @export
#'
#' @examples
mouseify_gene_symbols <- function(gene_list) {
  modify_gene_symbols(gene_list,
    text_replacements = c(
      mito_atp_replacer,
      c("MT-CO" = "MT-COX", "MT-" = "")
    ),
    text_modifications = list(tools::toTitleCase, tolower)
  )
}

load_json_config <- function(config, into = "list") {
  config_list <- rjson::fromJSON(file = config)

  if (into == "list") {
    out <- config_list
  } else if (into == "env") {
    out <- new.env()
    # purrr::walk2(
    #   names(config_list), config_list,
    #   function(x, y) {
    #     assign(x, y, envir = out)
    #   }
    # )
    list2env(config_list, envir = out)
  }
  return(out)
}

# TODO reorganize due to params
load_metadata <- function(path, params) {
  metadata <- read_data(path) %>%
    # dplyr::select(-dplyr::any_of(c("fq", "lane", "read"))) %>%
    dplyr::select({{ sample_colname }}, {{ group_colname }}) %>%
    base::unique() %>%
    dplyr::mutate({{ group_colname }} := factor(group,
      levels = group_levels
    )) %>%
    dplyr::arrange({{ group_colname }})
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
  "path_to_amigo_files" = load_amigo,
  "path_to_dire_files" = load_dire,
  "path_to_gene_lists" = load_gene_lists,
  "path_to_pathway_lists" = load_pathway_lists,
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
