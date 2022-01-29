# Packages ----------------------------------------------------------
library(tidyverse)
library(patchwork)
library(egg)
library(grid)
library(cowplot)
library(ungeviz)
library(scales)
source(file = "Simulation_functions_parameter_sweep.R")
theme_set(theme_cowplot())

# General plotting customization -------------------------------------------------
Figure_width = 7 # in width for one column figure in PNAS
Figure_height = 9
Figure_label_size = 10
margin_y_axis = 0.2
margin_right = 0.2
margin_left = 0.2
margin_x_axis = 0.1
plot_width = Figure_width - margin_right - margin_left - margin_y_axis
y_axis_label = expression(paste("Conjugation estimate (ml ", 'cfu'^-1, ' hr'^-1, ")"))
x_axis_label = expression("Time (hr)")
text_size_axis_title = 8

# basecase estimates  --------------------------------------------
basecase_frame <- Time_sweep_load_data(treatments_filename = "data/basecase",
                                 zero_value = 1e-16)
basecase_mean <- Time_sweep_mean_calculation(treatments_filename = "data/basecase",
                                             zero_value = 1e-16)

basecase <- Time_box_plot(sim_data_frame = basecase_frame, 
                          mean_dataframe = basecase_mean,
                          Set_conjugation_rate = 1e-14, 
                          y_axis_min = 1e-16, 
                          y_axis_max = 1e-12,
                          mean_line_length = 0.4)

# plasmid_effects estimates  --------------------------------------------

plasmid_effects_frame <- Time_sweep_load_data(treatments_filename = "data/plasmid_effects",
                                       zero_value = 1e-16)
plasmid_effects_mean <- Time_sweep_mean_calculation(treatments_filename = "data/plasmid_effects",
                                             zero_value = 1e-16)

plasmid_effects <- Time_box_plot(sim_data_frame = plasmid_effects_frame, 
                          mean_dataframe = plasmid_effects_mean,
                          Set_conjugation_rate = 1e-14, 
                          y_axis_min = 1e-16, 
                          y_axis_max = 1e-12,
                          mean_line_length = 0.4)
# host_effects estimates  --------------------------------------------

host_effects_frame <- Time_sweep_load_data(treatments_filename = "data/host_effects",
                                              zero_value = 1e-16)
host_effects_mean <- Time_sweep_mean_calculation(treatments_filename = "data/host_effects",
                                                    zero_value = 1e-16)

host_effects <- Time_box_plot(sim_data_frame = host_effects_frame, 
                                 mean_dataframe = host_effects_mean,
                                 Set_conjugation_rate = 1e-14, 
                                 y_axis_min = 1e-16, 
                                 y_axis_max = 1e-12,
                                 mean_line_length = 0.4)
# gammaT estimates  --------------------------------------------
gammaT_frame <- Time_sweep_load_data(treatments_filename = "data/gammaT",
                                           zero_value = 1e-16)
gammaT_mean <- Time_sweep_mean_calculation(treatments_filename = "data/gammaT",
                                                 zero_value = 1e-16)

gammaT <- Time_box_plot(sim_data_frame = gammaT_frame, 
                              mean_dataframe = gammaT_mean,
                              Set_conjugation_rate = 1e-14, 
                              y_axis_min = 1e-16, 
                              y_axis_max = 1e-12,
                              mean_line_length = 0.4)
# tau estimates  --------------------------------------------
tau_frame <- Time_sweep_load_data(treatments_filename = "data/tau",
                                     zero_value = 1e-16)
tau_mean <- Time_sweep_mean_calculation(treatments_filename = "data/tau",
                                           zero_value = 1e-16)

tau <- Time_box_plot(sim_data_frame = tau_frame, 
                        mean_dataframe = tau_mean,
                        Set_conjugation_rate = 1e-14, 
                        y_axis_min = 1e-16, 
                        y_axis_max = 1e-12,
                        mean_line_length = 0.4)



# Fig assemble  ------------------------------------------------------------
Fig <- (basecase/ plasmid_effects/ host_effects/ gammaT/ tau) +
  plot_annotation(tag_levels = list(c('a', 'b', 'c', 'd', 'e'))) & 
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(size = Figure_label_size)) &
  plot_layout(widths = unit(c(plot_width), c('in')))
y.grob <- textGrob(y_axis_label, gp=gpar(fontsize=text_size_axis_title), rot=90)
x.grob <- textGrob(x_axis_label, gp=gpar(fontsize=text_size_axis_title))
Figure <- grid.arrange(y.grob, as_grob(Fig), ncol = 2, widths = unit(c(margin_y_axis, (plot_width + margin_right + margin_left)), c('in', 'in')))
Figure2 <- grid.arrange(as_grob(Figure), x.grob, nrow = 2, heights = unit(c(Figure_height-margin_x_axis, margin_x_axis), c('in', 'in')))
save_plot("SI_Fig5.pdf", plot = Figure2, base_width = Figure_width, base_height = Figure_height)
