# Time sweep functions -------------------------------------------------
## Load data function
Time_sweep_load_data <- function(treatments_filename,
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
  }
  Data_frame <- bind_rows(csvfileslist)
  Data_frame$estimate <- factor(Data_frame$estimate, levels=c("LDM", "SIM", "TDR",  "ASM"), ordered = TRUE, labels=c("LDM", "SIM", "TDR",  "ASM"))
  return(Data_frame)
}

## Mean line
Time_sweep_mean_calculation <- function(treatments_filename,
                                        zero_value) {
  csvfiles <- list.files(treatments_filename, full.names = TRUE)
  csvfileslist <- lapply(csvfiles, read.csv, header = TRUE, stringsAsFactors = FALSE)
  treatments <- strsplit(csvfiles,"/")
  treatment_names <- rep(NA, length(treatments)) #This is to get the time points 
  for (i in 1:length(treatments)){
    filename <- treatments[[i]][length(treatments[[i]])]
    treatment_names[i] <- as.numeric(strsplit(filename,'_')[[1]][2])
  }
  meantable = data.frame()
  new_zero_value = zero_value
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
  meantable$estimate <- factor(meantable$estimate, levels=c("LDM", "SIM", "TDR",  "ASM"), ordered = TRUE, labels=c("LDM", "SIM", "TDR",  "ASM"))
  return(meantable)
}

## Box plot
Time_box_plot <- function(sim_data_frame, 
                          Set_conjugation_rate, 
                          set_conjugation_rate_linetype = 'dashed', 
                          set_conjugation_rate_color = c(rgb(169, 169, 169, maxColorValue = 255)), 
                          set_conjugation_rate_line_size = 0.75/(ggplot2::.pt*72.27/96), 
                          basecase_color= c(rgb(169, 169, 169, maxColorValue = 255)), 
                          basecase_color_fill = c(rgb(221, 221, 221, maxColorValue = 255)), 
                          box_line_size = 0.75/(ggplot2::.pt*72.27/96) , 
                          outlier_shape = 19, 
                          outlier_point_size = 0, 
                          mean_dataframe, 
                          mean_line_color = 'black', 
                          mean_line_size = 0.75/(ggplot2::.pt*72.27/96), 
                          mean_line_length = 0.08, 
                          parameter_xaxis_x_axis_label, 
                          text_size_axis_title = 8, 
                          axis_line_size = 0.75/(ggplot2::.pt*72.27/96), 
                          axis_tick_size = 0.75/(ggplot2::.pt*72.27/96), 
                          axis_tick_lengths = 0.03, 
                          text_size_axis_tick_labels = 8, 
                          margins = c(0.1, 0.2, 0.1, 0.2),
                          x_axis_label,
                          y_axis_min, 
                          y_axis_max, 
                          x_axis_min = 0,
                          x_axis_max = 9,
                          x_axis_breaks = c(0,3,6,9),
                          y_axis_number_of_ticks = 4, 
                          color_fill = c(c(rgb(226, 210, 195, maxColorValue = 255)), c(rgb(255, 229, 204, maxColorValue = 255)), c(rgb(205, 249, 247, maxColorValue = 255)), c(rgb(211, 216, 175, maxColorValue = 255))), #light tan and light orange
                          color = c('tan4', c(rgb(205, 102, 0, maxColorValue = 255)), c(rgb(0, 139, 139, maxColorValue = 255)), c(rgb(85, 107, 47, maxColorValue = 255)))) #tan4, darkorange3, green, cyan)
{p <- ggplot(sim_data_frame, aes(x=time, y=conjugationestimate, group = interaction(time, estimate), color = estimate, fill = estimate)) +
  geom_hline(yintercept = Set_conjugation_rate, linetype = set_conjugation_rate_linetype, color = set_conjugation_rate_color, size = set_conjugation_rate_line_size)+
  geom_boxplot(notch=FALSE, size = box_line_size, outlier.shape = outlier_shape, outlier.size = outlier_point_size, outlier.alpha = 0.25)+
  geom_hpline(data = mean_dataframe, color = mean_line_color, size = mean_line_size, width = mean_line_length)+
  scale_x_continuous(limits = c(x_axis_min, x_axis_max), breaks = x_axis_breaks)+
  scale_y_log10(limits=c(y_axis_min, y_axis_max), breaks = scales::trans_breaks("log10", function(x) 10^x, n = y_axis_number_of_ticks), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())+
  facet_grid(~estimate)+
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text = element_text(size = text_size_axis_tick_labels))+
  theme(legend.position = "none") +
  theme(plot.margin = unit(margins, "in"))+
  scale_fill_manual(values = color_fill)+
  scale_color_manual(values = color, guide = "none")+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(size = text_size_axis_tick_labels))
return(p)
}

# Parameter sweep functions -------------------------------------------------
## Load data function
Parameter_sweep_load_data <-
  function(treatments_filename,
           simulation_estimates_folder,
           new_zero_value,
           parameter_for_xaxis_values) {
    treatment_description <-read.csv(treatments_filename, header = TRUE)
    csvfiles <- list.files(simulation_estimates_folder, full.names = TRUE)
    select_estimates <- c("LDM", "SIM", "TDR", "ASM")
    estimate_labels <- c("LDM", "SIM", "TDR", "ASM")
    csvfileslist <-lapply(csvfiles,read.csv,header = TRUE,stringsAsFactors = FALSE)
    treatments <- strsplit(csvfiles, "/")
    treatment_names <- rep(NA, length(treatments))
    for (i in 1:length(treatments)) {
      filename <- treatments[[i]][length(treatments[[i]])]
      treatment_names[i] <- strsplit(filename, '_')[[1]][1]}
    for (l in 1:length(csvfileslist)) {
      csvfileslist[[l]] <- csvfileslist[[l]][-1, ]
      csvfileslist[[l]] <-subset(csvfileslist[[l]], select = select_estimates)
      csvfileslist[[l]]['LDM'][csvfileslist[[l]]['LDM'] == 0] <-new_zero_value
      csvfileslist[[l]]['SIM'][csvfileslist[[l]]['SIM'] == 0] <-new_zero_value
      csvfileslist[[l]]['TDR'][csvfileslist[[l]]['TDR'] == 0] <-new_zero_value
      csvfileslist[[l]]['ASM'][csvfileslist[[l]]['ASM'] == 0] <-new_zero_value
      csvfileslist[[l]][parameter_for_xaxis_values] <- rep(treatment_description[parameter_for_xaxis_values][treatment_description['Treatment_ID'] == treatment_names[l]], nrow(csvfileslist[[l]]))
      csvfileslist[[l]] <-pivot_longer(csvfileslist[[l]],cols = select_estimates, names_to = "estimate",values_to = "conjugationestimate")
    }
    Data_frame <- bind_rows(csvfileslist)
    Data_frame$estimate <-factor(Data_frame$estimate , levels = select_estimates, ordered = TRUE, labels = estimate_labels)
    return(Data_frame)
  }
## Mean line
Parameter_sweep_mean_calculation <-
  function(treatments_filename,
           simulation_estimates_folder,
           new_zero_value,
           parameter_for_xaxis_values) {
    treatment_description <-read.csv(treatments_filename, header = TRUE)
    csvfiles <- list.files(simulation_estimates_folder, full.names = TRUE)
    select_estimates <- c("LDM", "SIM", "TDR", "ASM")
    estimate_labels <- c("LDM", "SIM", "TDR", "ASM")
    csvfileslist <-lapply(csvfiles,read.csv,header = TRUE,stringsAsFactors = FALSE)
    treatments <- strsplit(csvfiles, "/")
    treatment_names <- rep(NA, length(treatments))
    for (i in 1:length(treatments)){
      filename <- treatments[[i]][length(treatments[[i]])]
      treatment_names[i] <- strsplit(filename,'_')[[1]][1]}
    meantable = data.frame()
    for (l in 1:length(csvfileslist)) {
      csvfileslist[[l]] <- csvfileslist[[l]][-1,]
      csvfileslist[[l]] <- subset(csvfileslist[[l]], select = select_estimates)
      csvfileslist[[l]][parameter_for_xaxis_values] <- rep(treatment_description[parameter_for_xaxis_values][treatment_description['Treatment_ID'] == treatment_names[l]], nrow(csvfileslist[[l]]))
      csvfileslist[[l]] <- pivot_longer(csvfileslist[[l]], cols=select_estimates , names_to = "estimate", values_to = "conjugationestimate")
      meanlist <- aggregate(csvfileslist[[l]]$conjugationestimate ~ csvfileslist[[l]]$estimate, FUN = mean)
      colnames(meanlist) <- c("estimate", "conjugationestimate")
      meanlist[parameter_for_xaxis_values] <- rep(treatment_description[parameter_for_xaxis_values][treatment_description['Treatment_ID'] == treatment_names[l]], nrow(meanlist))
      meantable <- rbind(meantable, meanlist)
    }
    meantable$estimate <- factor(meantable$estimate, levels=select_estimates, ordered = TRUE, labels = estimate_labels)
    return(meantable)
  }
## Box plot
Box_plot <- function(sim_data_frame, 
                     Set_conjugation_rate, 
                     set_conjugation_rate_linetype = 'dashed', 
                     set_conjugation_rate_color = c(rgb(169, 169, 169, maxColorValue = 255)), 
                     set_conjugation_rate_line_size = 0.75/(ggplot2::.pt*72.27/96), 
                     basecase_color= c(rgb(169, 169, 169, maxColorValue = 255)), 
                     basecase_color_fill = c(rgb(221, 221, 221, maxColorValue = 255)), 
                     box_line_size = 0.75/(ggplot2::.pt*72.27/96) , 
                     outlier_shape = 19, 
                     outlier_point_size = 0, 
                     mean_dataframe, 
                     mean_line_color = 'black', 
                     mean_line_size = 0.75/(ggplot2::.pt*72.27/96), 
                     mean_line_length = 0.21, 
                     parameter_xaxis_x_axis_label, 
                     text_size_axis_title = 8, 
                     axis_line_size = 0.75/(ggplot2::.pt*72.27/96), 
                     axis_tick_size = 0.75/(ggplot2::.pt*72.27/96), 
                     axis_tick_lengths = 0.03, 
                     text_size_axis_tick_labels = 8, 
                     margins = c(0.1, 0.2, 0.1, 0.2),
                     x_axis_label,
                     y_axis_min, 
                     y_axis_max, 
                     y_axis_number_of_ticks = 4, 
                     color_fill = c(c(rgb(226, 210, 195, maxColorValue = 255)), c(rgb(255, 229, 204, maxColorValue = 255)), c(rgb(205, 249, 247, maxColorValue = 255)), c(rgb(211, 216, 175, maxColorValue = 255))), #light tan and light orange
                     color = c('tan4', c(rgb(205, 102, 0, maxColorValue = 255)), c(rgb(0, 139, 139, maxColorValue = 255)), c(rgb(85, 107, 47, maxColorValue = 255)))) #tan4, darkorange3, green, cyan) 
{ sim_data_frame <- sim_data_frame %>% rename('parameter_xaxis' = colnames(sim_data_frame)[1])
  mean_dataframe <- mean_dataframe %>% rename('parameter_xaxis' = colnames(mean_dataframe)[3])
  p<-ggplot(sim_data_frame %>% filter(parameter_xaxis != 1), aes(x=parameter_xaxis, y=conjugationestimate, group = interaction(parameter_xaxis, estimate), fill = estimate))+
  geom_hline(yintercept = Set_conjugation_rate, linetype = set_conjugation_rate_linetype, color = set_conjugation_rate_color, size = set_conjugation_rate_line_size)+
  geom_boxplot(data = sim_data_frame, aes(x=parameter_xaxis, y=conjugationestimate, group = interaction(parameter_xaxis, estimate)), color = basecase_color, fill = basecase_color_fill, notch=FALSE, size = box_line_size, outlier.shape = outlier_shape, outlier.size = outlier_point_size, outlier.alpha = 0.25)+
  geom_boxplot(notch=FALSE, aes(colour = estimate), size = box_line_size, outlier.shape = outlier_shape, outlier.size = outlier_point_size, outlier.alpha = 0.25)+
  geom_hpline(data = mean_dataframe, color = mean_line_color, size = mean_line_size, width = mean_line_length)+
  xlab(x_axis_label)+
  theme(axis.title.y = element_blank())+
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text = element_text(size = text_size_axis_tick_labels))+
  theme(legend.position = "none") +
  theme(plot.margin = unit(margins, "in"))+
  facet_grid(~estimate)+
  scale_x_log10()+
  scale_y_log10(limits = c(y_axis_min, y_axis_max), breaks = trans_breaks("log10", function(x) 10^x, n = y_axis_number_of_ticks), labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_manual(values = color_fill)+
  scale_color_manual(values = color, guide = "none")+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(size = text_size_axis_tick_labels))
return(p)
}

## Box plot
Gamma_box_plot <- function(sim_data_frame, 
                     Set_conjugation_rate, 
                     set_conjugation_rate_linetype = 'dashed', 
                     set_conjugation_rate_color = c(rgb(169, 169, 169, maxColorValue = 255)), 
                     set_conjugation_rate_line_size = 0.75/(ggplot2::.pt*72.27/96), 
                     basecase_color= c(rgb(169, 169, 169, maxColorValue = 255)), 
                     basecase_color_fill = c(rgb(221, 221, 221, maxColorValue = 255)), 
                     box_line_size = 0.75/(ggplot2::.pt*72.27/96) , 
                     outlier_shape = 19, 
                     outlier_point_size = 0, 
                     mean_dataframe, 
                     mean_line_color = 'black', 
                     mean_line_size = 0.75/(ggplot2::.pt*72.27/96), 
                     mean_line_length = 0.69, 
                     parameter_xaxis_x_axis_label, 
                     text_size_axis_title = 8, 
                     axis_line_size = 0.75/(ggplot2::.pt*72.27/96), 
                     axis_tick_size = 0.75/(ggplot2::.pt*72.27/96), 
                     axis_tick_lengths = 0.03, 
                     text_size_axis_tick_labels = 8, 
                     margins = c(0.1, 0.2, 0.1, 0.2), 
                     x_axis_label,
                     y_axis_min, 
                     y_axis_max, 
                     y_axis_number_of_ticks = 4, 
                     color_fill = c(c(rgb(226, 210, 195, maxColorValue = 255)), c(rgb(255, 229, 204, maxColorValue = 255)), c(rgb(205, 249, 247, maxColorValue = 255)), c(rgb(211, 216, 175, maxColorValue = 255))), #light tan and light orange
                     color = c('tan4', c(rgb(205, 102, 0, maxColorValue = 255)), c(rgb(0, 139, 139, maxColorValue = 255)), c(rgb(85, 107, 47, maxColorValue = 255)))) #tan4, darkorange3, green, cyan) 
{ sim_data_frame <- sim_data_frame %>% rename('parameter_xaxis' = colnames(sim_data_frame)[1])
mean_dataframe <- mean_dataframe %>% rename('parameter_xaxis' = colnames(mean_dataframe)[3])
p<-ggplot(sim_data_frame %>% filter(parameter_xaxis != 1e-6), aes(x=parameter_xaxis, y=conjugationestimate, group = interaction(parameter_xaxis, estimate), fill = estimate))+
  geom_hline(yintercept = Set_conjugation_rate, linetype = set_conjugation_rate_linetype, color = set_conjugation_rate_color, size = set_conjugation_rate_line_size)+
  geom_boxplot(data = sim_data_frame, aes(x=parameter_xaxis, y=conjugationestimate, group = interaction(parameter_xaxis, estimate)), color = basecase_color, fill = basecase_color_fill, notch=FALSE, size = box_line_size, outlier.shape = outlier_shape, outlier.size = outlier_point_size, outlier.alpha = 0.25)+
  geom_boxplot(notch=FALSE, aes(colour = estimate), size = box_line_size, outlier.shape = outlier_shape, outlier.size = outlier_point_size, outlier.alpha = 0.25)+
  geom_hpline(data = mean_dataframe, color = mean_line_color, size = mean_line_size, width = mean_line_length)+
  xlab(x_axis_label)+
  theme(axis.title.y = element_blank())+
  theme(axis.title = element_text(size = text_size_axis_title)) +
  theme(axis.line.y = element_line(size = axis_line_size)) +
  theme(axis.line.x = element_line(size = axis_line_size)) +
  theme(axis.ticks = element_line(size = axis_tick_size))+
  theme(axis.ticks.length = unit(axis_tick_lengths, 'in'))+
  theme(axis.text = element_text(size = text_size_axis_tick_labels))+
  theme(legend.position = "none") +
  theme(plot.margin = unit(margins, "in"))+
  facet_grid(~estimate)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = y_axis_number_of_ticks), labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(limits = c(y_axis_min, y_axis_max), breaks = trans_breaks("log10", function(x) 10^x, n = y_axis_number_of_ticks), labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_manual(values = color_fill)+
  scale_color_manual(values = color, guide = "none")+
  theme(strip.background = element_blank())+
  theme(strip.text = element_text(size = text_size_axis_tick_labels))
return(p)
}