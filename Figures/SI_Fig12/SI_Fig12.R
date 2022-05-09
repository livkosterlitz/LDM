# Packages ----------------------------------------------------------
library(tidyverse)
library(patchwork)
library(egg)
library(grid)
library(cowplot)
library(devtools)
#devtools::install_github("wilkelab/ungeviz")
library(ungeviz)
library(scales)
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
text_size_axis_tick_labels = 8
axis_line_size = 0.75/(ggplot2::.pt*72.27/96) # pt size for the axis lines
axis_tick_size = 0.75/(ggplot2::.pt*72.27/96)
axis_tick_lengths = 0.03
margins_top_left <- c(0,0,0,0.0)
margins_bottom_left <- c(0.1,0,0,0.0)
margins_top_right <- c(0,0,0,0.0)
margins_bottom_right <- c(0.35,0,0,0.0)
# Filters and function ------------
new_zero_value = 1e-16
set_rate = 1e-14
y_min = 1e-16
y_max = 1e-12

Time_sweep_load_data_nofilter <- function(treatments_filename,
                                          zero_value){
  csvfiles <- list.files(treatments_filename, full.names = TRUE)
  csvfileslist <- lapply(csvfiles, read.csv, header = TRUE, stringsAsFactors = FALSE)
  treatments <- strsplit(csvfiles,"/")
  treatment_names <- rep(NA, length(treatments)) #This is to get the time points 
  for (i in 1:length(treatments)){
    filename <- treatments[[i]][length(treatments[[i]])]
    treatment_names[i] <- as.numeric(strsplit(filename,'_')[[1]][2])}
  new_zero_value = zero_value
  for (l in 1:length(csvfileslist)) {
    csvfileslist[[l]] <- csvfileslist[[l]][2:101,] #grab the long incubation populations
    csvfileslist[[l]] <- subset(csvfileslist[[l]], select = c('LDM','SIM'))
    csvfileslist[[l]]['time'] <- rep(treatment_names[l], nrow(csvfileslist[[l]]))
    place_holder <- pivot_longer(csvfileslist[[l]], cols=c('LDM', 'SIM'), names_to = "estimate",  values_to = "conjugationestimate")
    place_holder <- do.call(data.frame,lapply(place_holder, function(x) replace(x, is.infinite(x),NA)))
    n <- subset(place_holder, select = c('estimate', 'conjugationestimate')) %>%  # number of estimates at each time point
      group_by(estimate) %>%
      summarise(conjugationestimate = sum(conjugationestimate>=0))
    t <- subset(place_holder, select = c('estimate', 'conjugationestimate')) %>%   # number of NA estimates at each time point
      group_by(estimate) %>%
      summarise(conjugationestimate = sum(is.na(conjugationestimate)))
    z <- subset(place_holder, select = c('estimate', 'conjugationestimate')) %>% # number of zero estimates at each time point
      group_by(estimate) %>%
      summarise(conjugationestimate = sum(conjugationestimate==0))
    for (e in n$estimate) { 
      csvfileslist[[l]]['LDM'][csvfileslist[[l]]['LDM'] == 0] <-new_zero_value
      csvfileslist[[l]]['SIM'][csvfileslist[[l]]['SIM'] == 0] <-new_zero_value
    }
    csvfileslist[[l]] <- pivot_longer(csvfileslist[[l]], cols=c('LDM', 'SIM'), names_to = "estimate", values_to = "conjugationestimate")
  }
  Data_frame <- bind_rows(csvfileslist)
  Data_frame$estimate <- factor(Data_frame$estimate, levels=c("LDM", "SIM"), ordered = TRUE, labels=c("LDM", "SIM"))
  return(Data_frame)
}

Time_sweep_load_EPdata_nofilter <- function(treatments_filename,
                                            zero_value){
  csvfiles <- list.files(treatments_filename, full.names = TRUE)
  csvfileslist <- lapply(csvfiles, read.csv, header = TRUE, stringsAsFactors = FALSE)
  treatments <- strsplit(csvfiles,"/")
  treatment_names <- rep(NA, length(treatments)) #This is to get the time points 
  for (i in 1:length(treatments)){
    filename <- treatments[[i]][length(treatments[[i]])]
    treatment_names[i] <- as.numeric(strsplit(filename,'_')[[1]][3])}
  new_zero_value = zero_value
  for (l in 1:length(csvfileslist)) {
    csvfileslist[[l]] <- csvfileslist[[l]][2:101,] #grab the long incubation populations
    csvfileslist[[l]] <- subset(csvfileslist[[l]], select = c('LDM', 'SIM'))
    csvfileslist[[l]]['time'] <- rep(treatment_names[l], nrow(csvfileslist[[l]]))
    place_holder <- pivot_longer(csvfileslist[[l]], cols=c('LDM', 'SIM'), names_to = "estimate",  values_to = "conjugationestimate")
    place_holder <- do.call(data.frame,lapply(place_holder, function(x) replace(x, is.infinite(x),NA)))
    n <- subset(place_holder, select = c('estimate', 'conjugationestimate')) %>%  # number of estimates at each time point
      group_by(estimate) %>%
      summarise(conjugationestimate = sum(conjugationestimate>=0))
    t <- subset(place_holder, select = c('estimate', 'conjugationestimate')) %>%   # number of NA estimates at each time point
      group_by(estimate) %>%
      summarise(conjugationestimate = sum(is.na(conjugationestimate)))
    z <- subset(place_holder, select = c('estimate', 'conjugationestimate')) %>% # number of zero estimates at each time point
      group_by(estimate) %>%
      summarise(conjugationestimate = sum(conjugationestimate==0))
    for (e in n$estimate) { 3
      csvfileslist[[l]]['LDM'][csvfileslist[[l]]['LDM'] == 0] <-new_zero_value
      csvfileslist[[l]]['SIM'][csvfileslist[[l]]['SIM'] == 0] <-new_zero_value
    }
    csvfileslist[[l]] <- pivot_longer(csvfileslist[[l]], cols=c('LDM', 'SIM'), names_to = "estimate", values_to = "conjugationestimate")
  }
  Data_frame <- bind_rows(csvfileslist)
  Data_frame$estimate <- factor(Data_frame$estimate, levels=c("LDM", "SIM"), ordered = TRUE, labels=c("LDM", "SIM"))
  return(Data_frame)
}


# basecase estimates  --------------------------------------------
basecase_frame <- Time_sweep_load_data_nofilter(treatments_filename = "Time_series/base_case",
                                       zero_value = new_zero_value)

# EP0 estimates  --------------------------------------------

EP0_frame <- Time_sweep_load_EPdata_nofilter(treatments_filename = "Time_series/EP0",
                                    zero_value = new_zero_value)


# EP0.5 estimates  --------------------------------------------

EP0.5_frame <- Time_sweep_load_EPdata_nofilter(treatments_filename = "Time_series/EP0.5",
                                      zero_value = new_zero_value)


# EP0.9 estimates  --------------------------------------------

EP0.9_frame <- Time_sweep_load_EPdata_nofilter(treatments_filename = "Time_series/EP0.9",
                                      zero_value = new_zero_value)


# EP0.99 estimates  --------------------------------------------

EP0.99_frame <- Time_sweep_load_EPdata_nofilter(treatments_filename = "Time_series/EP0.99",
                                       zero_value = new_zero_value)

# tile plot ------------
basecase_frame$EP <- paste("no plating\nor extinction")
EP0_frame$EP <- 0
EP0.5_frame$EP <- 0.5
EP0.9_frame$EP <- 0.9
EP0.99_frame$EP <- 0.99
df <- rbind(basecase_frame, EP0_frame, EP0.5_frame, EP0.9_frame, EP0.99_frame)
df_restructure <- df %>%
  filter(estimate %in% c('LDM', 'SIM')) %>%
  mutate(zero = ifelse(conjugationestimate == new_zero_value | conjugationestimate==0 | is.na(conjugationestimate), 1, 0),
         infinite = ifelse(conjugationestimate == Inf, 1, 0),
         conjugationestimate = ifelse(conjugationestimate == new_zero_value | conjugationestimate == 0 | conjugationestimate == Inf | is.na(conjugationestimate), NA, conjugationestimate),
         finite.nonzero = ifelse(is.na(conjugationestimate), 0, 1)) %>%
  group_by(estimate, EP, time) %>%
  summarise(N = n(), 
            zero = sum(zero, na.rm = T),
            infinite = sum(infinite, na.rm = T),
            finite.nonzero = sum(finite.nonzero, na.rm = T), 
            varestimate = var(conjugationestimate, na.rm = T),
            conjugationestimate = mean(conjugationestimate, na.rm = T)) %>%
  filter(!(is.na(conjugationestimate)),
         finite.nonzero >= 75)

df_restructure$EP <- factor(df_restructure$EP, levels = c(paste("no plating\nor extinction"), 0, 0.5, 0.9, 0.99))

Estimate.tile <- ggplot(data = df_restructure, aes(x = time, y = EP, fill = -1*(-14 - log10(conjugationestimate)))) +
  geom_tile()+
  facet_wrap(~estimate, nrow = 1)+
  ylab('Extinction probability') + 
  xlab(x_axis_label)+
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text = element_text(size = text_size_axis_tick_labels))+
  theme(plot.margin = unit(margins_bottom_left, "in"))+
  labs(fill = paste("Estimate deviation\n(log10)"))+
  theme(legend.title = element_text(size = text_size_axis_title))+
  theme(legend.text  = element_text(size = text_size_axis_title))+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(size = text_size_axis_tick_labels))+
  scale_fill_gradientn(colours = c("blue", "light grey", "red"), values = c(0, 0, 1))
Estimate.tile

Variance.tile <- ggplot(data = df_restructure, aes(x = time, y = EP, fill = log10(varestimate))) +
  geom_tile()+
  facet_wrap(~estimate, nrow = 1)+
  ylab('Extinction probability') + 
  xlab(x_axis_label)+
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text = element_text(size = text_size_axis_tick_labels))+
  theme(plot.margin = unit(margins_bottom_left, "in"))+
  labs(fill = paste("Estimate variation\n(log10)"))+
  theme(legend.title = element_text(size = text_size_axis_title))+
  theme(legend.text  = element_text(size = text_size_axis_title))+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(size = text_size_axis_tick_labels))+
  scale_fill_gradientn(colours = c("black", 'light grey'))
Variance.tile

Infinite.Numbers.tile <- ggplot(data = df_restructure, aes(x = time, y = EP, fill = infinite)) +
  geom_tile()+
  facet_wrap(~estimate, nrow = 1)+
  ylab('Extinction probability') + 
  xlab(x_axis_label)+
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text = element_text(size = text_size_axis_tick_labels))+
  theme(plot.margin = unit(margins_bottom_left, "in"))+
  labs(fill = paste("number of \ninfinite\nestimates",'\n'))+
  theme(legend.title = element_text(size = text_size_axis_title))+
  theme(legend.text  = element_text(size = text_size_axis_title))+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(size = text_size_axis_tick_labels))
Infinite.Numbers.tile

Zero.Numbers.tile <- ggplot(data = df_restructure, aes(x = time, y = EP, fill = zero)) +
  geom_tile()+
  facet_wrap(~estimate, nrow = 1)+
  ylab('Extinction probability') + 
  xlab(x_axis_label)+
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text = element_text(size = text_size_axis_tick_labels))+
  theme(plot.margin = unit(margins_bottom_left, "in"))+
  labs(fill = paste("number of \nzero\nestimates",'\n'))+
  theme(legend.title = element_text(size = text_size_axis_title))+
  theme(legend.text  = element_text(size = text_size_axis_title))+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(size = text_size_axis_tick_labels))
Zero.Numbers.tile

Fig_var <- (Estimate.tile / Variance.tile/ Infinite.Numbers.tile/ Zero.Numbers.tile) + 
  plot_layout(widths = unit(c(4), units = c('in')), heights = unit(c(1.5,1.5,1.5,1.5), c('in', 'in')))
Fig_var <- (Fig_var)  +
  plot_annotation(tag_levels = list(c('a', 'b', 'c', 'd'))) & 
  theme(plot.tag.position = c(0, 1.05), plot.tag = element_text(size = Figure_label_size))
save_plot("Fig_var.pdf", plot = Fig_var, base_width = Figure_width, base_height = 8.8)

