## require the function files
script_folder = "scripts/"
# script_list = c("misc.R",
#                 "clusterprofiler.R",
#                 "correlations.R",
#                 "dire.R",
#                 "export.R",
#                 "heatmap.R",
#                 "ma.R",
#                 "pathways.R",
#                 "pca_plot.R",
#                 "read_data.R",
#                 "transform_data.R",
#                 "venn.R",
#                 "volcano_plot.R")
# 
# purrr::walk(script_list, function(x) { 
#   source(paste0(script_folder, "/", x), local = TRUE)
# })

source(paste0(script_folder, "misc.R"), local = TRUE)

source(paste0(script_folder, "clusterprofiler.R"), local = TRUE)
source(paste0(script_folder, "correlations.R"), local = TRUE)
source(paste0(script_folder, "dire.R"), local = TRUE)
source(paste0(script_folder, "export.R"), local = TRUE)
source(paste0(script_folder, "heatmap.R"), local = TRUE)
source(paste0(script_folder, "ma.R"), local = TRUE)
source(paste0(script_folder, "pathways.R"), local = TRUE)
source(paste0(script_folder, "pca_plot.R"), local = TRUE)
source(paste0(script_folder, "read_data.R"), local = TRUE)
source(paste0(script_folder, "transform_data.R"), local = TRUE)
source(paste0(script_folder, "venn.R"), local = TRUE)
source(paste0(script_folder, "volcano_plot.R"), local = TRUE)