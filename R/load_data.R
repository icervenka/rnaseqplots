#' Parse json config file
#'
#' @param config charcter string, path to json config file to load
#' @param into character string, one of c("list", "env") whether to return
#' loaded config file as list or as an environment
#'
#' @return list or environment of parsed config file
#' @export
#'
#' @examples
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
  expression_data <- read_data(path)
  return(expression_data)
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
  diffexp_data <- purrr::map(path, function(x) {
    if (nchar(x) > 0) {
      read_data(x)
    }
  }) %>%
    setNames(names(path))
  return(diffexp_data)
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
  dds <- readRDS(path)
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
  cufdiff_diff <- read_data_cuffdiff(append_dir_slash(path))
  return(cuffdiff_diff)
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
#' cp_pathways_pattern_to be present
#'
#' @return data frame
#' @export
#'
#' @examples
load_pathways_cp <- function(path, params) {
  cp_pathways <- collate_pathways_cp(
    append_dir_slash(path),
    pattern = params$cp_pathways_pattern
  )
  return(cp_pathways)
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
#' pvalue_threshold to be present
#'
#' @return data frame
#' @export
#'
#'
#' @examples
load_pathways_gsea <- function(path, params) {
  gsea_data <- collate_pathways_gsea(
    dir = append_dir_slash(path),
    pvalue_threshold = params$gsea_fdr_cutoff
  )
  return(gsea_data)
}

#' Title
#'
#' @param path
#' @param params
#'
#' @return
#' @export
#'
#' @examples
load_pathways_ipa <- function(path, params) {
  ipa_data <- collate_pathways_ipa(append_dir_slash(path), params$ipa_rank_by)
  return(ipa_data)
}

#' Title
#'
#' @param path
#' @param params
#'
#' @return
#' @export
#'
#' @examples
load_pathways_amigo <- function(path, params) {

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
#' dire_pattern and dire_sheet_name to be present
#'
#' @return data frame
#' @export
#'
#' @examples
load_pathways_dire <- function(path, params) {
  dire <- collate_pathways_dire(
    path,
    params$dire_pattern,
    params$dire_sheet_name
  )
  return(dire)
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


  # TODO move to helper function
  if (is.character(input_config) && endsWith("json")) {
    paths <- load_json_config(input_config)
  } else if (is.environment(input_config)) {
    paths <- as.list(input_config)
  } else if (is.list(input_config)) {
    paths <- input_config
  } else {
    stop("Unrecognized input config type.")
  }

  if (is.character(param_config) && endsWith("json")) {
    params <- load_json_config(param_config)
  } else if (is.environment(param_config)) {
    params <- as.list(param_config)
  } else if (is.list(param_config)) {
    params <- param_config
  } else {
    stop("Unrecognized input config type.")
  }

  data <- list()
  purrr::walk2(names(input_config), input_config, function(x, y) {
    if (nchar(y) > 0) {
      data[[x]] <- load_function_map$x(y, params)
    }
  })

  if (into == "list") {
    out <- data()
  } else {
    list2env(data, envir = out)
  }
  return(out)
}
