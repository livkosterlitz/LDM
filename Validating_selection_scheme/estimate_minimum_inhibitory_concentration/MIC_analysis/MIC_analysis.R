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

dat <- read.csv("recorded_MIC_data.csv")

dat <- dat %>%
  group_by(Name, Plate, Drug) %>%
  summarise(N = n(),
            Conc = mean(Concentration))

dat$Name <- factor(dat$Name)
dat$Name <- factor(dat$Name, levels = c("Donor", "Recipient", "Mating"))
levels(dat$Name)

colors_light <- c('firebrick1', 'blue2', 'darkorchid2')

p1 <- ggplot(dat %>% filter(Drug=='Tet'),
  aes(x = Name, y=Conc, color = Name, fill = Name))+ 
  geom_bar(stat = "identity")+
  scale_color_manual(values=colors_light) +
  scale_fill_manual(values=colors_light) +
  scale_y_continuous(trans = "log2",
                     breaks=2^seq(0, (max(dat$Conc)/2), by=2)) +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  facet_wrap(~ Plate) +
  theme(axis.text = element_text(size = 8)) +
  theme(strip.text.x = element_text(size = 8)) + 
  geom_hline(yintercept=7.5, linetype="solid", color = "firebrick1") +
  geom_hline(yintercept=60, linetype="dashed", color = "firebrick1")
p1


p2 <- ggplot(dat %>% filter(Drug=='Str'),
             aes(x = Name, y=Conc, color = Name, fill = Name))+ 
  geom_bar(stat = "identity")+
  scale_color_manual(values=colors_light) +
  scale_fill_manual(values=colors_light) +
  scale_y_continuous(trans = "log2",
                     breaks=2^seq(0, (max(dat$Conc)/2), by=2)) +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  facet_wrap(~ Plate) +
  theme(axis.text = element_text(size = 8)) +
  theme(strip.text.x = element_text(size = 8)) + 
  geom_hline(yintercept=25, linetype="solid", color = "blue2") +
  geom_hline(yintercept=200, linetype="dashed", color = "blue2")
p2

Figp1_fixed <- set_panel_size(p1, width  = unit(2, "in"), height = unit(2, "in"))
Figp2_fixed <- set_panel_size(p2, width  = unit(2, "in"), height = unit(2, "in"))

tet.y.grob <- textGrob("Tet drug concentration (ug/ml)",
                   gp=gpar(fontsize=10), rot=90)

str.y.grob <- textGrob("Str drug concentration (ug/ml)",
                       gp=gpar(fontsize=10), rot=90)

Figp1_fixed <- grid.arrange(arrangeGrob(Figp1_fixed, left = tet.y.grob))
Figp2_fixed <- grid.arrange(arrangeGrob(Figp2_fixed, left = str.y.grob))

Figure_main <- plot_grid(Figp1_fixed, Figp2_fixed, ncol = 1, align = "v")

Figure_main


save_plot("MIC_figure.pdf", plot = Figure_main
          , base_width = 5, base_height = 6)

