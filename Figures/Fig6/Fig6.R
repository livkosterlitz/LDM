# Load packages -----------------------------------------------------------
library(tidyverse)
library(readxl)
library(scales)
library(cowplot)
library(ggpubr)
library(rstatix)  
library(ggpubr)
library(ggpmisc)
theme_set(theme_cowplot())

# General plotting customization -------------------------------------------------------
Figure_width = 3.42 # in width for one column figure in PNAS
plot_width = 3
plot_height = plot_width
Figure_height = 4
Figure_label_size = 10
text_size_axis_title = 8
text_size_axis_tick_labels = 8
axis_line_size = 0.75/(ggplot2::.pt*72.27/96) # pt size for the axis lines
axis_tick_size = 0.75/(ggplot2::.pt*72.27/96)
axis_tick_lengths = 0.03
y_axis_min = 1e-14
y_axis_max = 1e-3
y_axis_number_of_ticks = 8
margins<- c(0,0,0,0)
y_axis_label = expression(paste("Conjugation estimate (ml ", 'cfu'^-1, ' hr'^-1, ")"))
color = c('tan4', c(rgb(205, 102, 0, maxColorValue = 255))) #tan4 and darkorange3
color_fill = c(c(rgb(226, 210, 195, maxColorValue = 255)), c(rgb(255, 229, 204, maxColorValue = 255))) #light tan and light orange
point_jitter_width = 0.15
point_size = 1

# Box plot customization -------------------------------------------------
box_line_size = 0.75/(ggplot2::.pt*72.27/96) 

# Load conjugation data -------------------------------------------------------
Corrected_data <- read.csv(file = "Fig6_experimental_data.csv")
Data_conjugation <- Corrected_data %>%
  select(Method, Time, Set, Corrected, Approximate) %>%
  mutate(method = Method) %>%
  unite("labels", Method:Set) 

# Statistics -------------------------------------------------------
Data_conjugation$labels <- factor(Data_conjugation$labels, levels= c("LDM_4_G", "SIM_24_G", "LDM_1.25_D", "SIM_5_G"), 
       labels = c("cross-spp.\n LDM", "cross-spp.\n SIM", "within-sp.\n LDM", "cross-spp.\n SIM truncated"))

my_comparisons <- list(c("cross-spp.\n LDM", "cross-spp.\n SIM truncated"),
                       c("cross-spp.\n SIM", "cross-spp.\n SIM truncated"),
                       c("cross-spp.\n LDM", "within-sp.\n LDM"),
                       c("cross-spp.\n LDM", "cross-spp.\n SIM"))

# Plot -------------------------------------------------------
p <- ggplot(data = Data_conjugation, 
       aes(x=labels, y=Corrected, color = method, fill = method))+
  geom_boxplot(notch=FALSE, size = box_line_size, outlier.shape = NA)+
  geom_jitter(width = point_jitter_width, size = point_size) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = 'p.signif') +
  ylab(y_axis_label)+
  theme(axis.title.x = element_blank()) +
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text = element_text(size = text_size_axis_tick_labels))+
  theme(legend.position = "none") +
  theme(plot.margin = unit(margins, "in"))+
  scale_y_log10(limits = c(y_axis_min, y_axis_max), breaks = trans_breaks("log10", function(x) 10^x, n = y_axis_number_of_ticks), labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_manual(values = color_fill)+
  scale_color_manual(values = color, guide = "none")+
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank())
p
# Assemble plot -------------------------------------------------------
save_plot("Fig.pdf", plot = p, base_width = Figure_width, base_height = Figure_height)




