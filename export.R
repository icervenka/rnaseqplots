library(egg)
ggsave_fixed = function(file, plot = last_plot(), units = "mm", margin = 1,
                        plot_width = 4, plot_height = 4,
                        width = round(dev.size()[1], digits = 1),
                        height = round(dev.size()[1], digits = 1)) {
  pf = egg::set_panel_size(p = plot,
                           file = NULL,
                           margin = unit(margin, units),
                           width = unit(plot_width, units),
                           height = unit(plot_height, units))
  ggplot2::ggsave(file, plot = pf, units = units, width = width, height = height)
}
