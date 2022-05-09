# Packages ----------------------------------------------------------
library(tidyverse)
library(cowplot)
library(patchwork)
library(egg)
library(grid)
library(ungeviz)
theme_set(theme_cowplot())

# Density data (ODE + stochastic) --------------------------------------------
dat <- read.csv("ODE_data/T55.csv")
T_1 = (dat$time[dat$Transconjugant>1]/3600)[1]
num_of_rows_to_skip = 100 # This skips datapoints so that the ggplot linetypes will work
dat2 <- dat %>%
  filter(row_number() %% num_of_rows_to_skip == 1)
dat1 <- dat2 %>%
  select(-PlasmidFreeDonor, -Total) %>%
  gather(Cell.type, Density, 2:4)
dat1$time <- dat1$time/3600
csvfiles <- list.files("Sim_data/", full.names = TRUE) #stochastic data
csvfileslist <- lapply(csvfiles, read.csv, header = TRUE, stringsAsFactors = FALSE)
treatments <- strsplit(csvfiles,"/")
treatment_names <- rep(NA, length(treatments)) #This is to get the time points 
for (i in 1:length(treatments)){
  filename <- treatments[[i]][length(treatments[[i]])]
  treatment_names[i] <- strsplit(filename,'_')[[1]][2]}
for (l in 1:length(csvfileslist)) {
  csvfileslist[[l]]["time"] <- csvfileslist[[l]]["time"]/3600
  csvfileslist[[l]]['replicate'] <- rep(treatment_names[l], nrow(csvfileslist[[l]]))}  
density_frame <- bind_rows(csvfileslist)
stochastic_frame <- density_frame %>%
  select(-PlasmidFreeDonor) %>%
  gather(Cell.type, Density, 2:4)
# Plot filtering ----------------------------------------------------------
Incubation_time = 8.516944
Set_conjugation_rate = 1e-14 # The simulated conjugation rate

# General plotting customization -------------------------------------------------
Figure_height = 5
Figure_width = 7
Figure_label_size = 10
plot_width = 2.5
text_size_axis_title = 8
text_size_axis_tick_labels = 8
axis_line_size = 0.75/(ggplot2::.pt*72.27/96) # pt size for the axis lines
axis_tick_size = 0.75/(ggplot2::.pt*72.27/96)
axis_tick_lengths = 0.03
x_axis_min = 0
x_axis_max = 9
x_axis_breaks = c(0,3,6,9)
x_axis_label = "Time (hr)"
margins_top_left <- c(0,0,0,0.0)
margins_bottom_left <- c(0.35,0,0,0.0)
margins_top_right <- c(0,0,0,0.0)
margins_bottom_right <- c(0.35,0,0,0.0)

# Density plotting customization -------------------------------------------------
density_cell_colors <- c('firebrick1', 'blue2', 'darkorchid2') # the line colors for donor, recipient, and transconjugant trajectories
density_linetype = 'solid' # the linetype used to plot the density trajectories
density_line_size = 1.5/(ggplot2::.pt*72.27/96) # pt size for the density trajectories
stochastic_line_size = 0.75/(ggplot2::.pt*72.27/96)
density_line_dodge = 0.8 # to jitter the overlapping density trajectories
density_y_axis_min = 1
density_y_axis_max = 1e9
density_y_axis_number_of_ticks = 4
density_height = plot_width*0.75
density_y_axis_label = expression(atop("Cell density", paste("(cfu ", 'ml'^-1, ")")))

# Conjugation rate plotting customization -------------------------------------------------
SIM_line_size = 1.5/(ggplot2::.pt*72.27/96)
SIM_y_axis_number_of_ticks = 4
SIM_y_axis_min = 1e-16
SIM_y_axis_max = 1e-12
SIM_height = plot_width * 0.5
SIM_y_axis_label = expression(atop("Conjugation rate estimate", paste("(ml ", 'cfu'^-1, ' hr'^-1, ")")))
mean_line_size = 1/(ggplot2::.pt*72.27/96) 
mean_line_color = 'black'
box_line_size = 0.75/(ggplot2::.pt*72.27/96) 
set_conjugation_rate_linetype = 'dashed'
basecase_color= c(rgb(169, 169, 169, maxColorValue = 255))
basecase_color_fill = c(rgb(221, 221, 221, maxColorValue = 255))
set_conjugation_rate_color = c(rgb(169, 169, 169, maxColorValue = 255))
set_conjugation_rate_line_size = 0.75/(ggplot2::.pt*72.27/96)
outlier_point_size = 0
outlier_shape = 19
SIM_mean_line_length = 0.34

# Simulation estimates  --------------------------------------------
csvfiles <- list.files("Sim_estimates/", full.names = TRUE)
csvfileslist <- lapply(csvfiles, read.csv, header = TRUE, stringsAsFactors = FALSE)
treatments <- strsplit(csvfiles,"/")
treatment_names <- rep(NA, length(treatments)) #This is to get the time points 
for (i in 1:length(treatments)){
  filename <- treatments[[i]][length(treatments[[i]])]
  treatment_names[i] <- as.numeric(strsplit(filename,'_')[[1]][2])
}
meantable = data.frame()
new_zero_value = 1e-16
for (l in 1:length(csvfileslist)) {
  csvfileslist[[l]] <- csvfileslist[[l]][2:101,] #grab the long incubation populations
  csvfileslist[[l]] <- subset(csvfileslist[[l]], select = c('LDM', 'TDR', 'SIM', 'ASM'))
  csvfileslist[[l]]['time'] <- rep(treatment_names[l], nrow(csvfileslist[[l]]))
  place_holder <- pivot_longer(csvfileslist[[l]], cols=c('LDM', 'TDR', 'SIM', 'ASM'), names_to = "estimate",  values_to = "conjugationestimate")
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
  for (e in n$estimate) { # filter the time points based on the number of estimates (n), the number of NA (t), and the number of zero estimate values (z)
    if (!is.na(n$conjugationestimate[n$estimate == e]) && n$conjugationestimate[n$estimate == e] < 90 || # If there is estimates and there is less than 90
        is.na(n$conjugationestimate[n$estimate == e]) && t$conjugationestimate[t$estimate == e] >= 1  || # If there is estimates and there is any NA estimates
        z$conjugationestimate[z$estimate == e] > 10) # If more than 10 of the estimates are zero 
    {csvfileslist[[l]][e] <- rep(NA, nrow(csvfileslist[[l]])) # If the previous statements are true, fill in all the estimates as NA so that we don't use this time point
    }
  csvfileslist[[l]]['LDM'][csvfileslist[[l]]['LDM'] == 0] <-new_zero_value
  csvfileslist[[l]]['SIM'][csvfileslist[[l]]['SIM'] == 0] <-new_zero_value
  csvfileslist[[l]]['TDR'][csvfileslist[[l]]['TDR'] == 0] <-new_zero_value
  csvfileslist[[l]]['ASM'][csvfileslist[[l]]['ASM'] == 0] <-new_zero_value
  }
  csvfileslist[[l]] <- pivot_longer(csvfileslist[[l]], cols=c('LDM', 'TDR', 'SIM', 'ASM'), names_to = "estimate", values_to = "conjugationestimate")
  if (all(is.na(csvfileslist[[l]]$conjugationestimate)) == FALSE) {
    meanlist <- aggregate(csvfileslist[[l]]$conjugationestimate ~ csvfileslist[[l]]$estimate, FUN = mean)
    colnames(meanlist) <- c("estimate", "conjugationestimate")
    meanlist['time'] <- rep(treatment_names[l], nrow(meanlist))
    meantable <- rbind(meantable, meanlist)
  }
}
tau_TR_frame <- bind_rows(csvfileslist)
tau_TR_frame$estimate <- factor(tau_TR_frame$estimate, levels=c("SIM", "LDM", "TDR",  "ASM"))
meantable$estimate <- factor(meantable$estimate, levels=c("SIM", "LDM", "TDR",  "ASM"))


# a. Density plot for SIM panel ------------------------------------------------------------
top_left <- ggplot(data = dat1 %>% filter(time < Incubation_time, Cell.type == 'Transconjugant'), aes(x = time, y = Density, color = Cell.type)) +
  geom_line(size = density_line_size, linetype = density_linetype) +
  geom_line(data = dat1 %>% filter(time < Incubation_time, Cell.type != 'Transconjugant'), position=position_dodge(width=density_line_dodge), size = density_line_size, linetype = density_linetype)+
  geom_line(data = stochastic_frame %>% filter(replicate == 'E4', time < Incubation_time, Cell.type == 'Transconjugant'), aes(x = time, y = Density, group = interaction(Cell.type, replicate), color = Cell.type), size = stochastic_line_size) +
  scale_color_manual(values = density_cell_colors) +
  scale_x_continuous(limits = c(x_axis_min, x_axis_max), breaks = x_axis_breaks)+
  scale_y_log10(limits = c(density_y_axis_min, density_y_axis_max), breaks = scales::trans_breaks("log10", function(x) 10^x, n = density_y_axis_number_of_ticks), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  ylab(density_y_axis_label) + 
  theme(axis.title.x = element_blank()) +
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text.y = element_text(size = text_size_axis_tick_labels))+
  theme(axis.text.x = element_blank()) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(margins_top_left, "in"))
top_left

# b. SIM conjugation estimate plot ------------------------------------------------------------
specific_estimate = "SIM"
estimate_color = c(rgb(205, 102, 0, maxColorValue = 255)) #darkorange3
estimate_fill = c(rgb(255, 229, 204, maxColorValue = 255))
bottom_left <- ggplot(data = dat1 %>% filter(time < Incubation_time & time > T_1), aes(x = time, y = SIM)) +
  geom_hline(yintercept = Set_conjugation_rate, linetype = set_conjugation_rate_linetype, color = set_conjugation_rate_color, size = set_conjugation_rate_line_size)+
  geom_line(size = SIM_line_size, color = estimate_color, alpha = 0.5) +
  geom_boxplot(data = tau_TR_frame %>% filter(estimate == specific_estimate), 
               aes(x=time, y=conjugationestimate, group = time), 
               notch=FALSE, colour = estimate_color, fill = estimate_fill,
               size = box_line_size, outlier.shape = outlier_shape, outlier.size = outlier_point_size, outlier.alpha = 0.25)+
  geom_hpline(data = meantable %>% filter(estimate == specific_estimate), aes(x=time, y=conjugationestimate, group = time), 
              color = mean_line_color, size = mean_line_size, width = SIM_mean_line_length)+
  scale_x_continuous(limits = c(x_axis_min, x_axis_max), breaks = x_axis_breaks)+
  scale_y_log10(limits=c(SIM_y_axis_min, SIM_y_axis_max), breaks = scales::trans_breaks("log10", function(x) 10^x, n = SIM_y_axis_number_of_ticks), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  ylab(SIM_y_axis_label) + 
  xlab(x_axis_label)+
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text = element_text(size = text_size_axis_tick_labels))+
  theme(legend.position = "none") +
  theme(plot.margin = unit(margins_bottom_left, "in"))
bottom_left

# c. Density plot for LDM panel ------------------------------------------------------------
top_right <- ggplot(data = dat1 %>% filter(time < Incubation_time, Cell.type == 'Transconjugant'), aes(x = time, y = Density, color = Cell.type)) +
  geom_line(size = density_line_size, linetype = density_linetype) +
  geom_line(data = dat1 %>% filter(time < Incubation_time, Cell.type != 'Transconjugant'), position=position_dodge(width=density_line_dodge), size = density_line_size, linetype = density_linetype)+
  geom_line(data = stochastic_frame %>% filter(time < Incubation_time, Cell.type == 'Transconjugant'), aes(x = time, y = Density, group = interaction(Cell.type, replicate), color = Cell.type), size = stochastic_line_size) +
  scale_color_manual(values = density_cell_colors) +
  scale_x_continuous(limits = c(x_axis_min, x_axis_max), breaks = x_axis_breaks)+
  scale_y_log10(limits = c(density_y_axis_min, density_y_axis_max), breaks = scales::trans_breaks("log10", function(x) 10^x, n = density_y_axis_number_of_ticks), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  ylab(density_y_axis_label) + 
  theme(axis.title.x = element_blank()) +
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text.y = element_text(size = text_size_axis_tick_labels))+
  theme(axis.text.x = element_blank()) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(margins_top_left, "in"))
top_right

# d. LDM conjugation estimate plot ------------------------------------------------------------
specific_estimate = "LDM"
estimate_color = "tan4"
estimate_fill = c(rgb(226, 210, 195, maxColorValue = 255))
bottom_right <- ggplot(tau_TR_frame %>% filter(estimate == specific_estimate), aes(x=time, y=conjugationestimate, group = time)) +
  geom_hline(yintercept = Set_conjugation_rate, linetype = set_conjugation_rate_linetype, color = set_conjugation_rate_color, size = set_conjugation_rate_line_size)+
  geom_boxplot(notch=FALSE, colour = estimate_color, fill = estimate_fill,
               size = box_line_size, outlier.shape = outlier_shape, outlier.size = outlier_point_size, outlier.alpha = 0.25)+
  geom_hpline(data = meantable %>% filter(estimate == specific_estimate), aes(x=time, y=conjugationestimate, group = time), 
              color = mean_line_color, size = mean_line_size, width = SIM_mean_line_length)+
  scale_x_continuous(limits = c(x_axis_min, x_axis_max), breaks = x_axis_breaks)+
  scale_y_log10(limits=c(SIM_y_axis_min, SIM_y_axis_max), breaks = scales::trans_breaks("log10", function(x) 10^x, n = SIM_y_axis_number_of_ticks), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  ylab(SIM_y_axis_label) + 
  xlab(x_axis_label)+
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text = element_text(size = text_size_axis_tick_labels))+
  theme(legend.position = "none") +
  theme(plot.margin = unit(margins_bottom_left, "in"))
bottom_right

# Fig assemble ------------------------------------------------------------
Figleft <- (top_left / bottom_left) + 
  plot_layout(widths = unit(c(plot_width), units = c('in')), heights = unit(c(density_height, SIM_height), c('in', 'in')))

Figright <- (top_right / bottom_right) + 
  plot_layout(widths = unit(c(plot_width), units = c('in')), heights = unit(c(density_height, SIM_height), c('in', 'in')))

Fig <- (Figleft | Figright)  +
  plot_annotation(tag_levels = list(c('b', 'c', 'd', 'e'))) & 
  theme(plot.tag.position = c(-0.03, 1.1), plot.tag = element_text(size = Figure_label_size))

save_plot("Fig.pdf", plot = Fig, base_width = Figure_width, base_height = Figure_height)
