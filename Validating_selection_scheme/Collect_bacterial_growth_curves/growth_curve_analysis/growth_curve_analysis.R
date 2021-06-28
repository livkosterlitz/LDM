library(ggplot2)
library(Rmisc)
library(dplyr)
library(lattice)
library(plyr)
library(cowplot)
library(gridExtra)
library(grid)
library(egg)
theme_set(theme_cowplot())

dat <- read.csv("recorded_density_plating.csv")

dat <- dat %>%
  group_by(Strain, Time) %>%
  summarise(N = n(),
            CFUs = mean(Density),
            SD = sd(Density),
            SE = SD/sqrt(N))

dat$Strain <- factor(dat$Strain)
levels(dat$Strain)

colors_light <- c('firebrick1', 'blue2')

p1 <- ggplot(dat, aes(x=Time, y=CFUs, color=Strain)) + 
  geom_line(linetype = "dashed", size = 0.4668623442372146) +
  geom_errorbar(aes(ymin=CFUs-SE, ymax=CFUs+SE), width=0, size = 0.4668623442372146) +
  geom_point(aes(shape=Strain, color=Strain, fill=Strain, size=Strain)) +
  scale_color_manual(values=colors_light) +
  scale_fill_manual(values=colors_light)+
  scale_size_manual(values=c(2,1.5,1.5,3))+
  scale_shape_manual(values=c(16,25,24,18)) +
  scale_y_log10(limits=c(min(dat$CFUs), 1e9),
                breaks=10^seq(0, 9, by=1),
                labels=seq(0, 9, by=1)) +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(axis.line.y = element_line(size = 0.3734899)) +
  theme(axis.line.x = element_line(size = 0.3734899)) +
  theme(axis.ticks = element_line(size = 0.3734899))+
  theme(axis.text.x = element_text(margin=margin(1,0,0,0,"pt")),
        axis.text.y = element_text(margin=margin(0,1,0,0,"pt")))+
  theme(axis.ticks.length=unit(.025, "in"))+
  theme(plot.margin = margin(.2, 0, 0, .05, "in"))+
  expand_limits(x = 0, y = 0) +
  theme(axis.text = element_text(size = 8))

p1

Figp1_fixed <- set_panel_size(p1, width  = unit(2, "in"), height = unit(1.35, "in"))

Figure_main <- plot_grid(Figp1_fixed)

Figure_main

#create common x and y labels

y.grob <- textGrob("cell density \n [log [(CFUs/mL)+1)]]",
                   gp=gpar(fontsize=10), rot=90)

x.grob <- textGrob("Time (day)", 
                   gp=gpar(fontsize=10))

#add common axis to plot

Figure <- grid.arrange(arrangeGrob(Figure_main, left = y.grob, bottom = x.grob))

save_plot("Growth_curve_figure.pdf", plot = Figure, base_width = 2.7, base_height = 2)

