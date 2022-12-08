# load data and output parameters ----------------------------------------------
analysis_data <- load_data(
    "config/input_config.json",
    "config/param_config.json",
    load_function_map
)

output_parameters <- rjson::fromJSON(file = "config/output_config.json")

# create and save plots --------------------------------------------------------
plot_pca_deseq(dds)
ggsave_param("deseq")

ma_plot(
  diffexp_data$ex_1,
  label_top_n = 5,
  label_bottom_n = 5
)
ggsave_param_wrapper("ma")