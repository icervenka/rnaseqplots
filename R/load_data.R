#' Helper for load config into list or environment
#'
#' @param x named list, config file parsed from json or other source
#' @param into character string, one of c("list", "env") whether to return
#' loaded config file as list or as an environment
#'
#' @return list or environment of parsed config file
#' @export
#'
#' @examples
load_into <- function(x, into) {
  match.arg(into, c("list", "env"))
  if (into == "env") {
    e <- new.env()
    list2env(x, envir = e)
    return(e)
  } else {
    return(x)
  }
}

#' Parse config file
#'
#' @param x charcter string, list or environment to parse config from.
#' Character string has to be a location to a json file
#' @param into character string, one of c("list", "env") whether to return
#' loaded config file as list or as an environment
#' @param ... other parameters passed to load_config
#'
#' @return list or environment of parsed config file
#' @export
#'
#' @examples
load_config <- function(x, into = "list", ...) {
  UseMethod("load_json_config", x)
}

## S3 method for class 'character'
load_config.character <- function(x, into = "list", ...) {
  stopifnot(endsWith(x, "json"))
  config_list <- rjson::fromJSON(file = x)
  return(load_into(config_list, into))
}

## S3 method for class 'list'
load_config.list <- function(x, into = "list", ...) {
  return(load_into(x, into))
}

## S3 method for class 'environment'
load_config.environment <- function(x, into = "list", ...) {
  return(load_into(as.list(x), into))
}

#' Load metadata information
#'
#' Wrapper function for loading data that exposes common interface
#' consisting of path to file/directory and parameter list parsed from json
#' config file
#'
#' @param path character string, path to metadata file
#' @param params list of parameters parsed from config file, requires
#' sample_colname, group_colname and group_levels to be present
#'
#' @return data frame
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
load_metadata <- function(path, params) {
  df <- read_data(path) %>%
    # dplyr::select(-dplyr::any_of(c("fq", "lane", "read"))) %>%
    dplyr::select(all_of(c(params$sample_colname, params$group_colname))) %>%
    base::unique() %>%
    dplyr::mutate(!!as.symbol(params$group_colname) := factor(group,
      levels = params$group_levels
    )) %>%
    dplyr::arrange(!!as.symbol(params$group_colname))
  rownames(df) <- NULL
  return(df)
}

#' Load normalized expression data
#'
#' Wrapper function for loading data that exposes common interface
#' consisting of path to file/directory and parameter list parsed from json
#' config file
#'
#' @param path character string, path to metadata file
#' @param params list of parameters parsed from config file. Currently unused
#'
#' @return data frame
#' @export
#'
#' @examples
load_data_expression <- function(path, params) {
  df <- read_data(path)
  return(df)
}

#' Load data containing differetial expression
#'
#' Wrapper function for loading data that exposes common interface
#' consisting of path to file/directory and parameter list parsed from json
#' config file
#'
#' @param path named list, names have to be unique, values are character strings
#' corresponding to file paths of data with differential expression
#' @param params list of parameters parsed from config file. Currently unused
#'
#' @return data frame
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
load_data_diffexp <- function(path, params) {
  df <- purrr::map(path, function(x) {
    if (nchar(x) > 0) {
      read_data(x)
    }
  }) %>%
    setNames(names(path))
  return(df)
}

#' Load DESeq2 dds file
#'
#' Wrapper function for loading data that exposes common interface
#' consisting of path to file/directory and parameter list parsed from json
#' config file
#'
#' @param path character string, path to dds file
#' @param params list of parameters parsed from config file. Currently unused
#'
#' @return data frame
#' @export
#'
#' @examples
load_data_deseq <- function(path, params) {
  df <- readRDS(path)
  return(df)
}

#' Load data from cuffdiff analysis
#'
#' Wrapper function for loading data that exposes common interface
#' consisting of path to file/directory and parameter list parsed from json
#' config file
#'
#' @param path character string, path to directory with cuffdiff data
#' @param params list of parameters parsed from config file. Currently unused
#'
#' @return data frame
#' @export
#'
#' @examples
load_data_cuffdiff <- function(path, params) {
  df <- read_data_cuffdiff(append_dir_slash(path))
  return(df)
}

#' Load data from clusterProfiler analysis
#'
#' Wrapper function for loading data that exposes common interface
#' consisting of path to file/directory and parameter list parsed from json
#' config file. Accepts files from
#' https://github.com/icervenka/clusterprofiler_reports_snakemake
#'
#' @param path character string, path to directory with
#' clusterprofiler_reports_snakemake csv output. Will search the path
#' recursively, see collate_pathways_cp for more information
#' @param params list of parameters parsed from config file. Requires
#' cp_pathways_pattern to be present
#'
#' @return data frame
#' @export
#'
#' @examples
load_pathways_cp <- function(path, params) {
  pathway_data <- collate_pathways_cp(
    append_dir_slash(path),
    pattern = params$cp_pathways_pattern
  )
  return(pathway_data)
}

#' Load data from GSEA analysis
#'
#' Wrapper function for loading data that exposes common interface
#' consisting of path to file/directory and parameter list parsed from json
#' config file. Accepts data from GUI GSEA app.
#'
#' @param path character string, path to directory with GSEA analysis.
#' Will search the path recursively, see collate_pathways_gsea for more
#' information
#' @param params list of parameters parsed from config file. Requires
#' gsea_pathways_pattern and pvalue_threshold to be present
#'
#' @return data frame
#' @export
#'
#'
#' @examples
load_pathways_gsea <- function(path, params) {
  pathway_data <- collate_pathways_gsea(
    dir = append_dir_slash(path),
    pattern = params$gsea_pathways_pattern,
    pvalue_threshold = params$gsea_fdr_cutoff
  )
  return(pathway_data)
}

#' Load data from IPA analysis
#'
#' Wrapper function for loading data that exposes common interface
#' consisting of path to file/directory and parameter list parsed from json
#' config file. Accepts data from Qiagen IPA analysis that were processed by
#' github.com/icervenka/ipa_reports_snakemake.
#'
#' @param path character string, path to directory with IPA analysis.
#' Will search the path recursively
#' @param params list of parameters parsed from config file. Requires
#' ipa_pathways_pattern and ipa_rank_by to be present. ipa_rank_by denotes the
#' column name to rank the pathways by
#'
#' @return data frame
#' @export
#'
#' @examples
load_pathways_ipa <- function(path, params) {
  pathway_data <- collate_pathways_ipa(
    append_dir_slash(pathway_data),
    pattern = params$ipa_pathways_pattern
    rank_by = params$ipa_rank_by
  )
  return(ipa_data)
}

#' Load data from AmiGO analysis
#'
#' Wrapper function for loading data that exposes common interface
#' consisting of path to file/directory and parameter list parsed from json
#' config file. Accepts data from AmiGO analysis in the form of json files
#' due to the preservation of hierarchical pathway analysis.
#'
#' @param path character string, path to directory with AmiGO analysis.
#' Will search the path recursively
#' @param params list of parameters parsed from config file. Requires
#' amigo_pathways_pattern and amigo_filter_level_up_to to be present
#'
#' @return data frame
#' @export
#'
#' @examples
load_pathways_amigo <- function(path, params) {
  pathway_data <- collate_pathways_amigo(
    append_dir_slash(path),
    pattern = params$amigo_pathways_pattern,
    filter_level_up_to = params$filter_level_up_to)
}

#' Load dire analysis files
#'
#' Wrapper function for loading data that exposes common interface
#' consisting of path to file/directory and parameter list parsed from json
#' config file
#'
#' @param path character string, path to directory containing files with dire
#' analysis. Will search the path recursively, see collate_pathways_dire for
#' more information
#' @param params list of parameters parsed from config file. Requires
#' dire_pathways_pattern and dire_sheet_name to be present
#'
#' @return data frame
#' @export
#'
#' @examples
load_pathways_dire <- function(path, params) {
  pathway_data <- collate_pathways_dire(
    path,
    pattern = params$dire_pathways_pattern,
    sheet_name = params$dire_sheet_name
  )
  return(pathway_data)
}

#' Load user specified gene lists from json file
#'
#' Wrapper function for loading data that exposes common interface
#' consisting of path to file/directory and parameter list parsed from json
#' config file. Gene list json file has user specified names and values can be
#' either lists of lenght > 1 containing gene symbols or lists of length == 1
#' containing path to text file with gene list to be read. File based gene lists
#' contain single column of gene names without header.
#'
#' @param path character string, path to json gene list file
#' @param params list of parameters parsed from config file. Currently unused
#'
#' @return named list, where names correspond to gene list names supplied by
#' user and values are character vectors of gene symbols/IDs
#' @export
#'
#' @importFrom magrittr %>%
#'
#' @examples
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

#' Load user specified pathway lists from json file
#'
#' Wrapper function for loading data that exposes common interface consisting
#' of path to file/directory and parameter list parsed from json config file.
#' Gene list json file has user specified names and values can be either
#' integer lists corresponding to the rank of the pathway in data or character
#' list with names of pathways.
#'
#' @param path character string, path to json pathway list file
#' @param params list of parameters parsed from config file. Currently unused
#'
#' @return
#' @export named list, where names correspond to gene list names supplied by
#' user and values are integer or character vectors
#'
#' @examples
load_pathway_lists <- function(path, params) {
  pathway_lists <- rjson::fromJSON(file = path)
}

# lookup table mapping path values from input config file to loading functions
load_function_map <- list(
  "path_to_metadata" = load_metadata,
  "path_to_expression_data" = load_data_expression,
  "path_to_diffexp_data" = load_data_diffexp,
  "path_to_deseq_dds" = load_data_deseq,
  "path_to_cuffdiff" = load_data_cuffdiff,
  "path_to_cp_pathway_files" = load_pathways_cp,
  "path_to_gsea_files" = load_pathways_gsea,
  "path_to_ipa_files" = load_pathways_ipa,
  "path_to_amigo_files" = load_pathways_amigo,
  "path_to_dire_files" = load_pathways_dire,
  "path_to_gene_lists" = load_gene_lists,
  "path_to_pathway_lists" = load_pathway_lists
)

#' Loads all user supplied data conforming to parameters specified
#'
#' @param input_config character string, list or environment or path to json
#' file containing locations of data files/directories to load
#' @param param_config character string, list or environment or path to json
#' file containing locations of parameter data file to load
#' @param load_function_map named list mapping input config file paths to data
#' loading function, where names are path_to_* names from input config and
#' values are load_* functions
#' @param into character string, one of c("list", "env") whether to return
#' loaded config file as list or as an environment
#'
#' @return list or environment with loaded data files
#' @export
#'
#' @examples
load_all_data <- function(input_config,
                          param_config,
                          load_function_map,
                          into = "list") {
  match.arg(into, c("list", "env"))

  paths <- load_json_config(input_config)
  params <- load_json_config(param_config)

  data <- list()
  purrr::walk2(names(input_config), input_config, function(x, y) {
    if (nchar(y) > 0) {
      data[[x]] <- load_function_map$x(y, params)
    }
  })

  return(load_into(data, into))
}
