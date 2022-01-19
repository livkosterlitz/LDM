library(tidyverse)
library(readxl)
library(cowplot)
library(grid)
library(cowplot)
library(egg)
theme_set(theme_cowplot())

dat <- read_xlsx("Example_data/recorded_density_plating.xlsx")

dat_f <- dat %>% 
  select(-Predicted.counts, -Predicted.density, -Counting.day, -Volume.plated, 
         -Plate.label,-Plate.type, -Strain) %>%
  mutate(Counts = as.numeric(Counts)) %>% 
  filter(!is.na(Counts), Counts != 0) %>%
  group_by(Day, Set, Experiment, Time, Plate.cell.type) %>%
  summarise(N = n(), 
            CFUs = mean(Density),
            SD = sd(Density),
            SE = SD/sqrt(N))  %>%
  mutate(Time = as.numeric(Time))

dat_g <- dat_f %>%
  arrange(Time) %>%
  group_by(Day, Set, Experiment, Plate.cell.type) %>%
  mutate(Growthrate = log((CFUs/lag(CFUs, default = first(CFUs))))/(Time-lag(Time, default = 0)))

colors_light <- c('firebrick1', 'blue2')

p1 <- ggplot(data = dat_g, aes(x = Time, y = CFUs, color = Plate.cell.type))+
  geom_point() +
  geom_errorbar(aes(ymin=CFUs-SE, ymax=CFUs+SE, width = 0)) +
  geom_line() +
  theme(legend.position = "none")+
  scale_color_manual(values=colors_light) +
  ylab("cell density per ml") +
  scale_y_log10(limits=c(min(dat_g$CFUs), max(dat_g$CFUs)),
                breaks=10^seq(0, 9, by=1),
                labels=seq(0, 9, by=1)) +
  xlab("Time (hours)") 

p2 <- ggplot(data = dat_g, aes(x = Time, y = Growthrate, color = Plate.cell.type))+
  geom_point() + 
  geom_line() +
  theme(legend.title=element_blank()) +
  scale_color_manual(values=colors_light) +
  ylab("Growth rate") +
  xlab("Time (hours)") 

Figp1_fixed <- set_panel_size(p1, width  = unit(2, "in"), height = unit(1.35, "in"))
Figp2_fixed <- set_panel_size(p2, width  = unit(2, "in"), height = unit(1.35, "in"))

Figure_main <- plot_grid(Figp1_fixed, Figp2_fixed)

Figure <- grid.arrange(arrangeGrob(Figure_main))

save_plot("Plots/Growth_curve_figure.pdf", plot = Figure, base_width = 6.5, base_height = 2.5)
