# load required packages
source("scripts/load_packages.R", local = TRUE)

# load package functions
source("scripts/load_scripts.R", local = TRUE)

# parses paths.json file for locations of user supplied data
source("scripts/load_paths.R", local = TRUE)

# load gene name replacers
source("scripts/gene_name_replacers.R", local = TRUE)

# parses config.json file for parameters
source("scripts/load_config.R", local = TRUE)

# load data and metadata
source("scripts/load_data.R", local = TRUE)