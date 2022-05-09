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
Figure_height = 7.5
Figure_label_size = 10
margin_y_axis = 0.2
margin_right = 0.2
margin_left = 0.2
plot_width = Figure_width - margin_right - margin_left - margin_y_axis
y_axis_label = expression(paste("Conjugation estimate (ml ", 'cfu'^-1, ' hr'^-1, ")"))
text_size_axis_title = 8

# psiD sweep data -------------------------------------------------
psiD_frame <- Parameter_sweep_load_data(treatments_filename = "Treatments_master.csv",
                                        simulation_estimates_folder = "data/psiD",
                                        new_zero_value = 1e-9, 
                                        parameter_for_xaxis_values = 'psiD')
psiD_mean_frame <- Parameter_sweep_mean_calculation(treatments_filename = "Treatments_master.csv",
                                                    simulation_estimates_folder = "data/psiD",
                                                    new_zero_value = 1e-9, 
                                                    parameter_for_xaxis_values = 'psiD')
psiD <- Box_plot(sim_data_frame = psiD_frame, 
                 mean_dataframe = psiD_mean_frame,
                 x_axis_label = expression(paste("Donor growth rate (", 'hr'^-1, ")")),
                 Set_conjugation_rate = 1e-6, 
                 y_axis_min = 1e-9, 
                 y_axis_max = 1e-4)

# psiR sweep data -------------------------------------------------
psiR_frame <- Parameter_sweep_load_data(treatments_filename = "Treatments_master.csv",
                                        simulation_estimates_folder = "data/psiR",
                                        new_zero_value = 1e-9, 
                                        parameter_for_xaxis_values = 'psiR')
psiR_mean_frame <- Parameter_sweep_mean_calculation(treatments_filename = "Treatments_master.csv",
                                                    simulation_estimates_folder = "data/psiR",
                                                    new_zero_value = 1e-9, 
                                                    parameter_for_xaxis_values = 'psiR')
psiR <- Box_plot(sim_data_frame = psiR_frame, 
                 mean_dataframe = psiR_mean_frame,
                 x_axis_label = expression(paste("Recipient growth rate (", 'hr'^-1, ")")),
                 Set_conjugation_rate = 1e-6, 
                 y_axis_min = 1e-9, 
                 y_axis_max = 1e-4)

# psiT sweep data -------------------------------------------------
psiT_frame <- Parameter_sweep_load_data(treatments_filename = "Treatments_master.csv",
                                        simulation_estimates_folder = "data/psiT",
                                        new_zero_value = 1e-9, 
                                        parameter_for_xaxis_values = 'psiT')
psiT_mean_frame <- Parameter_sweep_mean_calculation(treatments_filename = "Treatments_master.csv",
                                                    simulation_estimates_folder = "data/psiT",
                                                    new_zero_value = 1e-9, 
                                                    parameter_for_xaxis_values = 'psiT')
psiT <- Box_plot(sim_data_frame = psiT_frame, 
                 mean_dataframe = psiT_mean_frame,
                 x_axis_label = expression(paste("Transconjugant growth rate (", 'hr'^-1, ")")),
                 Set_conjugation_rate = 1e-6, 
                 y_axis_min = 1e-9, 
                 y_axis_max = 1e-2)

# Assemble plot -------------------------------------------------
Fig <- (psiD/ psiR/ psiT) +
  plot_annotation(tag_levels = list(c('a', 'b', 'c'))) & 
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(size = Figure_label_size)) &
  plot_layout(widths = unit(c(plot_width), c('in')))
y.grob <- textGrob(y_axis_label, gp=gpar(fontsize=text_size_axis_title), rot=90)
Figure <- grid.arrange(y.grob, as_grob(Fig), ncol = 2, widths = unit(c(margin_y_axis, (plot_width + margin_right + margin_left)), c('in', 'in')))
save_plot("SI_Fig2.pdf", plot = Figure, base_width = Figure_width, base_height = Figure_height)

