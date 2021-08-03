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
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input_file <- opt$file
dat <- read.csv(input_file)

dat <- dat %>%
  filter(Counts != 'x', Control != 'y') %>%
  select(-Control, -Predicted.counts, -Predicted.density) %>%
  mutate(Density = as.numeric(Density)) %>%  
  group_by(Strain, Cell.type, Time) %>%
  summarise(N = n(),
            CFUs = mean(Density),
            SD = sd(Density),
            SE = SD/sqrt(N))

dat$Cell.type <- factor(dat$Cell.type)
dat$Cell.type <- factor(dat$Cell.type, levels = c("D", "R"))
write.csv(dat, file = paste(substr(input_file, 1, nchar(input_file)-4), '_analyzed.csv'))

growth <- data.frame(matrix(ncol = 4, nrow = 0))
x <- c("Strain", "Cell.Type", "Time", "growth_rate")
colnames(growth) <- x

for (i in levels(dat$Cell.type)){
  temp <- dat %>%
    filter(Cell.type == i)
  
  growth_frame <- data.frame(matrix(ncol = 4, nrow = (nrow(temp)-1)))
  x <- c("Strain", "Cell.type", "Time", "growth_rate")
  colnames(growth_frame) <- x
  
  for(i in 1:(nrow(temp)-1)){
    growth_frame$Strain[i] <- as.character(temp$Strain[i])
    growth_frame$Cell.type[i] <- as.character(temp$Cell.type[i])
    growth_frame$Time[i] <- (temp$Time[i]+temp$Time[i+1])/2
    growth_frame$growth_rate[i] <- log((temp$CFUs[i+1]/temp$CFUs[i]))/(temp$Time[i+1]-temp$Time[i])
  }
  
  growth <- rbind(growth, growth_frame)
}
write.csv(growth, file = paste(substr(input_file, 1, nchar(input_file)-4), '_growthrates.csv'))

colors_light <- c('firebrick1', 'blue2')

p1 <- ggplot(dat, aes(x=Time, y=CFUs, color=Cell.type)) + 
  geom_line(linetype = "dashed", size = 0.4668623442372146) +
  geom_errorbar(aes(ymin=CFUs-SE, ymax=CFUs+SE), width=0, size = 0.4668623442372146) +
  geom_point(aes(shape=Cell.type, color=Cell.type, fill=Cell.type, size=Cell.type)) +
  scale_color_manual(values=colors_light) +
  scale_fill_manual(values=colors_light)+
  scale_size_manual(values=c(2,1.5,1.5,3))+
  scale_shape_manual(values=c(16,25,24,18)) +
  scale_y_log10(limits=c(min(dat$CFUs), 1e9),
                breaks=10^seq(0, 9, by=1),
                labels=seq(0, 9, by=1)) +
  theme(legend.position = "none")+
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
  theme(axis.text = element_text(size = 8))  +
  ylab("cell density \n [log [(CFUs/mL)+1)]]") +
  theme(axis.title.y = element_text(size = 10))


p1


colors_light <- c('firebrick1', 'blue2')

p2 <- ggplot(growth, aes(x=Time, y=growth_rate, color=Cell.type)) + 
  geom_line(linetype = "dashed", size = 0.4668623442372146) +
  geom_point(aes(shape=Cell.type, color=Cell.type, fill=Cell.type, size=Cell.type)) +
  scale_color_manual(values=colors_light) +
  scale_fill_manual(values=colors_light)+
  scale_size_manual(values=c(2,1.5,1.5,3))+
  scale_shape_manual(values=c(16,25,24,18)) +
  theme(legend.title=element_blank())+
  theme(legend.text=element_text(size=8))+
  theme(axis.title.x=element_blank()) +
  theme(axis.line.y = element_line(size = 0.3734899)) +
  theme(axis.line.x = element_line(size = 0.3734899)) +
  theme(axis.ticks = element_line(size = 0.3734899))+
  theme(axis.text.x = element_text(margin=margin(1,0,0,0,"pt")),
        axis.text.y = element_text(margin=margin(0,1,0,0,"pt")))+
  theme(axis.ticks.length=unit(.025, "in"))+
  theme(plot.margin = margin(0.2, 0.5, 0, .1, "in"))+
  expand_limits(x = 0, y = 0) +
  ylab("growth rate") +
  theme(axis.text = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 10))

p2



Figp1_fixed <- set_panel_size(p1, width  = unit(2, "in"), height = unit(1.35, "in"))
Figp2_fixed <- set_panel_size(p2, width  = unit(2, "in"), height = unit(1.35, "in"))

Figure_main <- plot_grid(Figp1_fixed, Figp2_fixed)

Figure_main

#create common x and y labels

x.grob <- textGrob("Time (day)", 
                   gp=gpar(fontsize=10))

#add common axis to plot

Figure <- grid.arrange(arrangeGrob(Figure_main, bottom = x.grob))

save_plot("Growth_curve_figure.pdf", plot = Figure, base_width = 6, base_height = 2.2)

