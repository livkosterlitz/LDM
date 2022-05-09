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
Figure_height = 2.5
Figure_label_size = 10
margin_y_axis = 0.2
margin_right = 0.2
margin_left = 0.2
plot_width = Figure_width - margin_right - margin_left - margin_y_axis
y_axis_label = expression(paste("Conjugation estimate (ml ", 'cfu'^-1, ' hr'^-1, ")"))
text_size_axis_title = 8

# gammaTR sweep data -------------------------------------------------
gammaTR_frame <- Parameter_sweep_load_data(treatments_filename = "Treatments_master.csv",
                                           simulation_estimates_folder = "data/gammaTR",
                                           new_zero_value = 1e-8, 
                                           parameter_for_xaxis_values = 'gammaTR')
gammaTR_mean_frame <- Parameter_sweep_mean_calculation(treatments_filename = "Treatments_master.csv",
                                                       simulation_estimates_folder = "data/gammaTR",
                                                       new_zero_value = 1e-8, 
                                                       parameter_for_xaxis_values = 'gammaTR')
gammaTR <- Gamma_box_plot(sim_data_frame = gammaTR_frame, 
                          mean_dataframe = gammaTR_mean_frame,
                          x_axis_label = expression(paste("Transconjugant conjugation rate (ml ", 'cfu'^-1, ' hr'^-1, ")")),
                          Set_conjugation_rate = 1e-6, 
                          y_axis_min = 1e-8, 
                          y_axis_max = 1e-2)

# Assemble plot -------------------------------------------------
Fig <- (gammaTR) +
  plot_annotation(tag_levels = list(c(''))) & 
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(size = Figure_label_size)) &
  plot_layout(widths = unit(c(plot_width), c('in')))
y.grob <- textGrob(y_axis_label, gp=gpar(fontsize=text_size_axis_title), rot=90)
Figure <- grid.arrange(y.grob, as_grob(Fig), ncol = 2, widths = unit(c(margin_y_axis, (plot_width + margin_right + margin_left)), c('in', 'in')))
save_plot("SI_Fig3.pdf", plot = Figure, base_width = Figure_width, base_height = Figure_height)

