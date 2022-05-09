# Packages ----------------------------------------------------------
library(tidyverse)
library(cowplot)
library(patchwork)
library(egg)
library(grid)
library(ungeviz)
theme_set(theme_cowplot())

# Plot filtering ----------------------------------------------------------
Set_conjugation_rate = 1e-6 # The simulated conjugation rate
Figure_height = 6
Figure_width = 7
Figure_label_size = 10
plot_width = 6
text_size_axis_title = 8
text_size_axis_tick_labels = 8
axis_line_size = 0.75/(ggplot2::.pt*72.27/96) # pt size for the axis lines
axis_tick_size = 0.75/(ggplot2::.pt*72.27/96)
axis_tick_lengths = 0.03
x_axis_min = 0
x_axis_max = 4
x_axis_breaks = c(0,1,2,3,4)
x_axis_label = "Time (hr)"
margins_top_left <- c(0,0,0,0.0)
margins_bottom_left <- c(0.1,0,0,0.0)
margins_top_right <- c(0,0,0,0.0)
margins_bottom_right <- c(0.35,0,0,0.0)
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
SIM_line_size = 1.5/(ggplot2::.pt*72.27/96)
SIM_y_axis_number_of_ticks = 4
SIM_y_axis_min = 1e-10
SIM_y_axis_max = 1e-3
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
specific_estimate = "LDM"
estimate_color = "tan4"
estimate_fill = c(rgb(226, 210, 195, maxColorValue = 255))

# W LDM variance ----------
csvfiles.var <- list.files("Time_series/", full.names = TRUE)
csvfileslist.var <- lapply(csvfiles.var, read.csv, header = TRUE, stringsAsFactors = FALSE)
treatments <- strsplit(csvfiles.var,"/")
well_names <- rep(NA, length(treatments)) #This is to get the time points
treatment_names <- rep(NA, length(treatments)) #This is to get the time points
for (i in 1:length(treatments)){
  filename <- treatments[[i]][length(treatments[[i]])]
  well_names[i] <- strsplit(filename,'_')[[1]][1]
  treatment_names[i] <- strsplit(filename,'_')[[1]][2]}
for (l in 1:length(csvfileslist.var)) {
  csvfileslist.var[[l]] <- csvfileslist.var[[l]][2:101,] #grab the long incubation populations
  csvfileslist.var[[l]] <- subset(csvfileslist.var[[l]], select = c('LDM'))
  csvfileslist.var[[l]]['time'] <- rep(treatment_names[l], nrow(csvfileslist.var[[l]]))
  csvfileslist.var[[l]]['W'] <- rep(well_names[l], nrow(csvfileslist.var[[l]]))
}
var_frame <- bind_rows(csvfileslist.var)
var_frame$W <- factor(var_frame$W, 
                                 levels = c('T83', 'T82', 'T81', 'T9', 'T80'), 
                                 labels = c(5, 10, 50, 100, 500))
LDM.Var <- var_frame %>%
  group_by(time, W) %>%
  mutate(zerocount = ifelse(LDM == 0, 1, 0),
         infcount = ifelse(LDM == Inf, 1, 0), 
         positiveLDM = ifelse(LDM == 0 | LDM == Inf, 0, 1),
         LDM.filter = ifelse(LDM == 0 | LDM == Inf, NA, LDM)) %>%
  summarise(Var.LDM = var(LDM.filter, na.rm = T),
            Counts.0 = sum(zerocount),
            Counts.Inf = sum(infcount),
            Counts.pos = sum(positiveLDM))

LDM.Var.tile <- ggplot(data = LDM.Var %>% filter(time < 4), aes(x = as.numeric(time), y = W, fill = log10(Var.LDM))) +
  geom_tile()+
  scale_x_continuous(limits = c(x_axis_min, x_axis_max), breaks = x_axis_breaks)+
  ylab('W') + 
  xlab(x_axis_label)+
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text = element_text(size = text_size_axis_tick_labels))+
  theme(plot.margin = unit(margins_bottom_left, "in"))+
  labs(fill = paste("estimate\nvariance\n(log10)", '\n'))+
  theme(legend.title = element_text(size = text_size_axis_title))+
  theme(legend.text  = element_text(size = text_size_axis_title))+
  scale_fill_gradientn(colours = c("black", 'light grey'), na.value = 'white')
LDM.Var.tile

LDM.infinite.tile <- ggplot(data = LDM.Var, aes(x = as.numeric(time), y = W, fill = Counts.Inf)) +
  geom_tile()+
  scale_x_continuous(limits = c(x_axis_min, x_axis_max), breaks = x_axis_breaks)+
  ylab('W') + 
  xlab(x_axis_label)+
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text = element_text(size = text_size_axis_tick_labels))+
  theme(plot.margin = unit(margins_bottom_left, "in"))+
  labs(fill = paste("number of\ninfinite\nestimates",'\n'))+
  theme(legend.title = element_text(size = text_size_axis_title))+
  theme(legend.text  = element_text(size = text_size_axis_title))
LDM.infinite.tile

LDM.zero.tile <- ggplot(data = LDM.Var, aes(x = as.numeric(time), y = W, fill = Counts.0)) +
  geom_tile()+
  scale_x_continuous(limits = c(x_axis_min, x_axis_max), breaks = x_axis_breaks)+
  ylab('W') + 
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
  theme(legend.text  = element_text(size = text_size_axis_title))
LDM.zero.tile

# W averaging ----------
csvfiles <- list.files("T_crit_estimates/", full.names = TRUE)
csvfileslist <- lapply(csvfiles, read.csv, header = TRUE, stringsAsFactors = FALSE)
treatments <- strsplit(csvfiles,"/")
well_names <- rep(NA, length(treatments)) #This is to get the time points
for (i in 1:length(treatments)){
  filename <- treatments[[i]][length(treatments[[i]])]
  well_names[i] <- strsplit(filename,'_')[[1]][1]}
treatment_info <- read.csv(file = 'Treatments_reps_with_incubation.csv')
for (l in 1:length(csvfileslist)) {
  csvfileslist[[l]] <- csvfileslist[[l]][2:nrow(csvfileslist[[l]]),] #grab the long incubation populations
  csvfileslist[[l]] <- subset(csvfileslist[[l]], select = c('LDM'))
  add <- rep(seq(10), 500/(treatment_info$Num_sims_I2[treatment_info$Treatment_ID == well_names[l]]))
  csvfileslist[[l]]['rep'] <- c(add , rep(NA, nrow(csvfileslist[[l]]) - length(add)))
  csvfileslist[[l]] <- csvfileslist[[l]]  %>%
    mutate(LDM = ifelse(LDM == Inf, NA, LDM), 
           positiveLDM = ifelse(LDM == 0 | LDM == Inf, 0, 1), 
           zeroLDM = ifelse(LDM == 0, 1, 0),
           infiniteLDM = ifelse(LDM == Inf, 1, 0)) %>%
    group_by(rep) %>%
    summarise(N = n(), 
              positiveLDM = sum(positiveLDM, na.rm = T),
              LDM = mean(LDM, na.rm = T),
              zeroLDM = sum(zeroLDM, na.rm = T),
              infiniteLDM = sum(infiniteLDM, na.rm = T)) %>%
    filter(!(is.na(rep)))
  csvfileslist[[l]]['W'] <- rep(well_names[l], nrow(csvfileslist[[l]]))
}

Wframe <- bind_rows(csvfileslist)
Wframe$W <- factor(Wframe$W, 
                      levels = c('T83', 'T82', 'T81', 'T9', 'T80'), 
                      labels = c(5, 10, 50, 100, 500))
#write.csv(x = Wframe, file = 'rep_averages_per_W.csv')

W_samereps <- ggplot(Wframe, aes(x=W, y=LDM)) +
  geom_hline(yintercept = Set_conjugation_rate, linetype = set_conjugation_rate_linetype, color = set_conjugation_rate_color, size = set_conjugation_rate_line_size)+
  geom_boxplot(notch=FALSE, colour = estimate_color, fill =estimate_fill,
               size = box_line_size, outlier.shape = outlier_shape, outlier.size = outlier_point_size, outlier.alpha = 0.25)+
  geom_point(color = estimate_color, position = position_jitter(width = 0.15))+
  scale_y_log10(limits=c(10**-6.5, 10**-5.5), breaks = scales::trans_breaks("log10", function(x) 10^x, n = 3), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  ylab(SIM_y_axis_label) +
  xlab('W')+
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text = element_text(size = text_size_axis_tick_labels))+
  theme(legend.position = "none") +
  theme(plot.margin = unit(margins_bottom_left, "in"))
W_samereps

Fig_var <- (LDM.Var.tile / LDM.infinite.tile/LDM.zero.tile/ W_samereps) + 
  plot_layout(widths = unit(c(2), units = c('in')), heights = unit(c(1.5, 1.5, 1.5,1.5), c('in', 'in')))
Fig_var <- (Fig_var)  +
  plot_annotation(tag_levels = list(c('a', 'b','c','d'))) & 
  theme(plot.tag.position = c(0, 1.05), plot.tag = element_text(size = Figure_label_size))
save_plot("Fig_var.pdf", plot = Fig_var, base_width = Figure_width/1.75, base_height = 7.8)

