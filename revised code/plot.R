library(openxlsx)
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(GGally)
library(patchwork)
library(posterior)
library(bayesplot)
library(odin2)
library(dust2)
library(monty)

rm(list = ls())

gen_sir <- odin('sir.R')

file <- 'mcmc/sensitivity/est-phi/'

samples1 <- readRDS(paste0(file, 'EV71-samples-2.rds'))  # 0
samples2 <- readRDS(paste0(file, 'CA16-samples-5.rds'))  # 0
samples3 <- readRDS(paste0(file, 'CA6-samples-5.rds'))  # 0

output_func0 <- function(serotype, samples0){
  
  samples <- samples0
  samples_thin <- if (serotype == 'EV71') {
    monty_samples_thin(samples, thinning_factor = 160, burnin = 2000)} else{
      monty_samples_thin(samples, thinning_factor = 80, burnin = 1000)
      }
    
    
  samples_df <- posterior::as_draws_df(samples_thin)
  fit_pars <- as.list(apply(samples_df, 2, mean))[1:(ncol(samples_df)-3)]
  
  time_len <- if (serotype == 'CA6') {14*52} else  {20*52}
  start_year <- if (serotype == 'CA6') {2011} else  {2005}

  data <- read.xlsx('data/hfmd derived ev by week.xlsx') %>% 
    filter(virus == serotype) %>% arrange(week.no) %>% 
    mutate(time = week.no, cases = value) %>% select(time, cases) %>% 
    filter(!is.na(cases)) %>% mutate(time = time - (start_year-2005) * 52)
  
  
  pop_data <- read.xlsx('data/Japan population.xlsx', sheet = 'Sheet2') %>% 
    mutate(birth_rate = Live.births/(Total*52), death_rate = Deaths/(Total*52),
           birth = Live.births/(52), death = Deaths/(52)) %>%      
    
    filter(Year >= start_year) %>% 
    
    mutate(time = 52 * (Year - start_year) + 26) %>% 
    dplyr::select(Year, time, Total, birth_rate, death_rate, birth, death)
  
  
  idx_data <- read.csv('data/OxCGRT_timeseries_StringencyIndex_v1.csv') %>% 
    filter(CountryName == 'Japan') %>% dplyr::select(-(2:7)) %>% 
    melt(id.vars = 'CountryCode', variable.name = "date") %>% 
    mutate(Date = seq(as.Date('2020-01-01'),as.Date('2023-02-28'), by = '1 day'),
           week = week(Date),
           week = ifelse(week == 53 , 52, week),
           week.no = (year(Date) - start_year) * 52 + week) %>%
    filter(! is.na(value)) %>% 
    group_by(week.no) %>% summarise(idx = mean(value))
  
  
  time_before_endidx <- (2022- start_year + 1) * 52  # 2005-2022 end
  t_values <- seq(1, 52*2, 1) + time_before_endidx    # 2023 and 2024
  Bt_value <- rowMeans(splines2::iSpline(t_values, degree = 1, df = 5)) 
  
  
  fixed_pars <- list(
    period = 52,
    gamma = 1,
    N_0 = pop_data$Total[1],
    birth_time = c(0, pop_data$time, time_len), 
    birth_value = c(pop_data$birth_rate[1], pop_data$birth_rate, pop_data$birth_rate[nrow(pop_data)]),
    death_time = c(0, pop_data$time, time_len), 
    death_value = c(pop_data$death_rate[1], pop_data$death_rate, pop_data$death_rate[nrow(pop_data)]),
    
    idx_start_time = (2020 - start_year) * 52 + 1,  #  from 2005 to 2020
    idx_end_time = time_before_endidx,
    idx_end_value = 33.33,
    
    idx_time = c(0, idx_data$week.no, time_len), 
    idx_value = c(0, idx_data$idx, idx_data$idx[nrow(idx_data)]),
    
    Bt_time = c(0, t_values),
    Bt_value = c(0, Bt_value))
  
  pars_list <- c(fixed_pars, fit_pars)
  
  filter <- dust_filter_create(gen_sir(), data = data, time_start = 0, n_particles = 1000, n_threads = 4)
  
  dust_likelihood_run(filter, pars_list)
  
  dim(samples_thin$observations$trajectories)
  
  # trajectory fitting output ---- 
  
  traj_fit <- dust_unpack_state(filter, samples_thin$observations$trajectories)
  
  inc <- t(array(traj_fit$incidence, c(nrow(data), 200)))
  sus <- t(array(traj_fit$S, c(nrow(data), 200)))
  total <- t(array(traj_fit$N, c(nrow(data), 200)))
  
  dim(inc)
  
  p <- c(samples_thin$pars['p', , ])
  cases_mod <- inc * p
  
  sus <- as.data.frame(array(traj_fit$S, c(nrow(data), 200))) %>% 
    mutate(week = seq(1, time_len, 1)) %>% 
    pivot_longer(cols = 1:200, names_to = 'sample.no', values_to = 'sus')
  
  total <- as.data.frame(array(traj_fit$N, c(nrow(data), 200))) %>% 
    mutate(week = seq(1, time_len, 1)) %>% 
    pivot_longer(cols = 1:200, names_to = 'sample.no', values_to = 'total')
  
  cases <- as.data.frame(t(cases_mod)) %>% 
    mutate(week = seq(1, time_len, 1)) %>% 
    pivot_longer(cols = 1:200, names_to = 'sample.no', values_to = 'case_mod')
  
  traj <- sus %>% left_join(total) %>% left_join(cases) %>% 
    left_join(data, by = c('week' = 'time'))
  
  return(list(traj = traj, samples = samples_thin, samples_df = samples_df, pars_list = pars_list))
  
}

data1 <- output_func0(serotype = 'EV71', samples0 = samples1)
data2 <- output_func0(serotype = 'CA16', samples0 = samples2)
data3 <- output_func0(serotype = 'CA6', samples0 = samples3)


# diag ----

jpeg(paste0(file, 'trace and density plot.tiff'), width = 9, height = 11, units = "in", res = 426)
((mcmc_trace(data1$samples_df) + labs(title = 'A. EV-A71')) + theme(axis.text.x = element_text(angle = 45)) + mcmc_dens(data1$samples_df) + theme(axis.text.x = element_text(angle = 45))) / 
  ((mcmc_trace(data2$samples_df) + labs(title = 'B. CVA16')) + theme(axis.text.x = element_text(angle = 45)) + mcmc_dens(data2$samples_df)+ theme(axis.text.x = element_text(angle = 45))) / 
  ((mcmc_trace(data3$samples_df) + labs(title = 'C. CVA6')) + theme(axis.text.x = element_text(angle = 45)) + mcmc_dens(data3$samples_df)+ theme(axis.text.x = element_text(angle = 45)))
dev.off()


jpeg(paste0(file, 'density.tiff'), width = 6, height = 4, units = "in", res = 426)
par(mfrow = c(3,1), oma = c(0,0,0,0), mar = c(4,4,2,1))
hist(drop(data1$samples$density), ylab = '', xlab = 'EV-A71 log-likelihood', main = NULL)
hist(drop(data2$samples$density), ylab = '', xlab = 'CVA16 log-likelihood', main = NULL)
hist(drop(data3$samples$density), ylab = '', xlab = 'CVA6 log-likelihood', main = NULL)
dev.off()


jpeg(paste0(file, "ev71-pairs.tiff"), width = 6, height = 5, units = "in", res = 426)
pairs(data1$samples_df[, c(1:10)], main = 'EV-A71', pch = 19, cex = 0.5, oma = c(2,2,5,2), gap = 0, cex.labels = 1.5)
dev.off()

jpeg(paste0(file, "ev71-acf.tiff"), width = 6, height = 3.5, units = "in", res = 426)
par(mfrow = c(2,5), mar = c(4,4,2,1))
coda::autocorr.plot(data1$samples_df[, c(1:10)], auto.layout = F, lag.max = 30)
dev.off()

jpeg(paste0(file, "ca16-pairs.tiff"), width = 6, height = 5, units = "in", res = 426)
pairs(data2$samples_df[, c(1:10)], main = 'CVA16', pch = 19, cex = 0.5, oma = c(2,2,5,2), gap = 0, cex.labels = 1.5)
dev.off()

jpeg(paste0(file, "ca16-acf.tiff"), width = 6, height = 3.5, units = "in", res = 426)
par(mfrow = c(2,5), mar = c(4,4,2,1))
coda::autocorr.plot(data2$samples_df[, c(1:10)], auto.layout = F, lag.max = 30)
dev.off()

jpeg(paste0(file, "ca6-pairs.tiff"), width = 6, height = 5, units = "in", res = 426)
pairs(data3$samples_df[, c(1:10)], main = 'CVA6', pch = 19, cex = 0.5, oma = c(2,2,5,2), gap = 0, cex.labels = 1.5)
dev.off()

jpeg(paste0(file, "ca6-acf.tiff"), width = 6, height = 3.5, units = "in", res = 426)
par(mfrow = c(2,5), mar = c(4,4,2,1))
coda::autocorr.plot(data3$samples_df[, c(1:10)], auto.layout = F, lag.max = 30)
dev.off()


pars_estimates <- function(data0, serotype){
  
  lik <- array(data0$samples$density, c(200))
  
  tables <- 
    posterior::summarise_draws(data0$samples_df) %>%
    select(1,2,6,7) %>% filter(variable != 'chain') %>%
    mutate(across(where(is.numeric), ~ ifelse(abs(.) >= 0.01 & abs(.) <= 100,
                                              sprintf('%.2f', .),
                                              sprintf('%.2e', .))
    )) %>% 
    
    rbind(c('log_lik', round(mean(lik),0), round(quantile(lik, 0.025),0), round(quantile(lik,0.975),0))) %>% 
    mutate(serotype = serotype, para = paste0(mean, " (", q5, ", ", q95, ")")) %>% 
    select(serotype, variable, para)
  
  return(tables)
}



pars_tab <- rbind.data.frame(pars_estimates(data1, 'EV71'),
                             pars_estimates(data2, 'CA16'),
                             pars_estimates(data3, 'CA6')) %>% 
  pivot_wider(names_from = serotype, values_from = para)


# write.xlsx(pars_tab, paste0(file, 'pars table.xlsx'))


# traj ----

plot.data <- function(data0){
   a <- data0$traj %>% 
     mutate(sus_re = (sus / 1e+6),
            add = floor(min(sus_re)),
            sus_re = sus_re - add,
            add_times = max(case_mod / sus_re))
   return(a)
}
 

pdata1 <- plot.data(data1)
pdata2 <- plot.data(data2)
pdata3 <- plot.data(data3)

  
col1 <- c('#ADCDEC',"gray65",'#ED0000FF')
labs <- c(2005,'', 2007,'', 2009,'', 2011,'', 2013,'', 2015,'',
          2017,'', 2019,'', 2021,'', 2023, '', 2025)
theme1 <- 
  theme(panel.background = element_rect(colour = "black", fill = 'white', size = 1),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        plot.margin = margin(0.2,0.3,0.2,0.2, "cm"),
        plot.title = element_text(colour = "black",size = 18, face = "bold"),
        axis.text = element_text(colour = "black",size = 18,hjust = 0.5,vjust = 0.5),
        axis.title = element_text(size = 17,vjust = 0, hjust = 0.5,face = "bold"))



p1 <- 
  ggplot(pdata1) +
  geom_line(aes(x = week, y = sus_re * max(add_times), group = sample.no), size = 0.8, color = col1[1])+
  geom_line(aes(x = week, y = case_mod, group = sample.no), size = 0.8, color = col1[2])+
  geom_line(aes(x = week, y = cases), size = 0.6, color = col1[3]) +
  scale_y_continuous(name = 'EV-A71 cases', 
                     sec.axis = sec_axis(trans = ~./unique(pdata1$add_times) + unique(pdata1$add),
                                         name='Susceptibles\n(million)', breaks = seq(7, 13, 2))) +
  scale_x_continuous(name = NULL, limits = c(0, 20*52), expand = c(0,10), 
                     breaks = seq(0, 20*52, 52), labels = labs) +
  theme1


p2 <- 
  ggplot(pdata2) +
  geom_line(aes(x = week, y = sus_re * max(add_times), group = sample.no), size = 0.8, color = col1[1])+
  geom_line(aes(x = week, y = case_mod, group = sample.no), size = 0.8, color = col1[2])+
  geom_line(aes(x = week, y = cases), size = 0.6, color = col1[3]) +
  scale_y_continuous(name = 'CVA16 cases', 
                     sec.axis = sec_axis(trans = ~./unique(pdata2$add_times) + unique(pdata2$add),
                                         name='Susceptibles\n(million)')) +
  scale_x_continuous(name = NULL, limits = c(0, 20*52), expand = c(0,10), 
                     breaks = seq(0, 20*52, 52), labels = labs) +
  theme1


p3 <- 
  ggplot() +
  geom_line(data = pdata3, aes(x = week + 312, y = sus_re * max(add_times), group = sample.no), size = 0.8, color = col1[1])+
  geom_line(data = pdata3, aes(x = week + 312, y = case_mod, group = sample.no), size = 0.8, color = col1[2])+
  geom_line(data = pdata3, aes(x = week + 312, y = cases), size = 0.6, color = col1[3]) +
  scale_y_continuous(name = 'CVA6 cases', 
                     sec.axis = sec_axis(trans = ~./unique(pdata3$add_times) + unique(pdata3$add),
                                         name='Susceptibles\n(million)')) +
  scale_x_continuous(name = NULL, limits = c(0, 20*52), expand = c(0,10), 
                     breaks = seq(0, 20*52, 52), labels = labs) +
  geom_point(aes(x = c(30, 30, 30), y = c(120,90,60)), color = col1, size = 4) +
  annotate('text', x = rep(50,3), y = c(120,90,60), 
           label = c('Susceptibles','Predicted cases', 'Observed cases'), size = 6, hjust = 0) + 
  theme1


jpeg(paste0(file, 'simulation.tiff'), width = 12, height = 10, units = "in", res = 426)
p1/p2/p3
dev.off()


# susceptible ----


s.71 <- pdata1 %>% 
  mutate(p = sus / total) %>% 
  group_by(week) %>% 
  summarise(p.m = mean(p), p.5 = quantile(p, 0.025), p.95 = quantile(p, 0.975),
            sus.m = mean(sus), sus.5 = quantile(sus, 0.025), sus.95 = quantile(sus, 0.975)) %>% 
  mutate(week.no = week + 104,
         year = 2003 + floor((week.no - 1)/52),
         circle = case_when(year >= 2019 ~ 5,
                            year >= 2015 ~ 4,
                            year >= 2011 ~ 3,
                            year >= 2007 ~ 2,
                            year >= 2005 ~ 1),
         circle.week = (ifelse(year >= 2007, week.no - 208 * (circle-1), week.no)),
         circle = as.character(circle)) 


s.16 <- pdata2 %>% 
  mutate(p = sus / total) %>% 
  group_by(week) %>% 
  summarise(p.m = mean(p), p.5 = quantile(p, 0.025), p.95 = quantile(p, 0.975),
            sus.m = mean(sus), sus.5 = quantile(sus, 0.025), sus.95 = quantile(sus, 0.975)) %>% 
  mutate(week.no = week,
         year = 2005 + floor((week.no - 1)/52),
         circle = case_when(year >= 2020 ~ 5,
                            year >= 2016 ~ 4,
                            year >= 2012 ~ 3,
                            year >= 2009 ~ 2,
                            year >= 2005 ~ 1),
         circle.week = week.no - 208 * (circle-1),
         circle.week = ifelse(year >= 2012, circle.week + 52, circle.week),
         circle = as.character(circle))

s.6 <- pdata3 %>% 
  mutate(p = sus / total) %>% 
  group_by(week) %>% 
  summarise(p.m = mean(p), p.5 = quantile(p, 0.025), p.95 = quantile(p, 0.975),
            sus.m = mean(sus), sus.5 = quantile(sus, 0.025), sus.95 = quantile(sus, 0.975)) %>% 
  mutate(week.no = week,
         year = 2011 + floor((week.no - 1)/52),
         circle = case_when(year >= 2023 ~ 7,
                            year >= 2020 ~ 6,
                            year >= 2018 ~ 5,
                            year >= 2016 ~ 4,
                            year >= 2014 ~ 3,
                            year >= 2012 ~ 2,
                            year >= 2011 ~ 1),
         circle.week = case_when(year == 2011 ~ week.no + 52,
                                 year >= 2012 ~ (week.no + 52) - 52*2*(circle - 1),
                                 year >= 2020 ~ (week.no + 52) - 52*3*(circle - 1),
                                 year >= 2023 ~ (week.no + 52) - 52*2*(circle - 1)),
         circle.week = ifelse(year >= 2023, circle.week - 52, circle.week),
         circle = as.character(circle))


s.71.s <- pdata1 %>% 
  mutate(week.no = week + 104,
         year = 2003 + floor((week.no - 1)/52),
         circle = case_when(year >= 2019 ~ 5,
                            year >= 2015 ~ 4,
                            year >= 2011 ~ 3,
                            year >= 2007 ~ 2,
                            year >= 2005 ~ 1),
         circle.week = (ifelse(year >= 2007, week.no - 208 * (circle-1), week.no)),
         circle = as.character(circle)) %>% 
  
  mutate(period = ifelse(circle < 5, 'pre', 'post')) %>% 
  group_by(sample.no, circle) %>% slice_max(sus) %>% 
  group_by(period) %>% 
  summarise(sus.m = mean(sus), p.m = mean(sus/total),
            sus.5 = quantile(sus, 0.025), p.5 = quantile(sus/total, 0.025),
            sus.95 = quantile(sus, 0.975), p.95 = quantile(sus/total, 0.975)) %>% 
  mutate(p = paste0(sprintf('%.2f', 100*p.m), ' (95% CrI: ', sprintf('%.2f', 100*p.5), ', ', sprintf('%.2f', 100*p.95), ')'),
         sus = paste0(sprintf('%.3g', sus.m/1e6), ' (95% CrI: ', sprintf('%.3g', sus.5/1e6), ', ', sprintf('%.3g', sus.95/1e6), ')')) %>% 
  mutate(sero = 'EV71') %>% select(10,1,8,9)


s.16.s <- pdata2 %>% 
  mutate(week.no = week,
         year = 2005 + floor((week.no - 1)/52),
         circle = case_when(year >= 2020 ~ 5,
                            year >= 2016 ~ 4,
                            year >= 2012 ~ 3,
                            year >= 2009 ~ 2,
                            year >= 2005 ~ 1),
         circle.week = week.no - 208 * (circle-1),
         circle.week = ifelse(year >= 2012, circle.week + 52, circle.week),
         circle = as.character(circle)) %>% 
  
  mutate(period = ifelse(circle < 5, 'pre', 'post')) %>% 
  group_by(sample.no, circle) %>% slice_max(sus) %>% 
  group_by(period) %>% 
  summarise(sus.m = mean(sus), p.m = mean(sus/total),
            sus.5 = quantile(sus, 0.025), p.5 = quantile(sus/total, 0.025),
            sus.95 = quantile(sus, 0.975), p.95 = quantile(sus/total, 0.975)) %>% 
  mutate(p = paste0(sprintf('%.2f', 100*p.m), ' (95% CrI: ', sprintf('%.2f', 100*p.5), ', ', sprintf('%.2f', 100*p.95), ')'),
         sus = paste0(sprintf('%.3g', sus.m/1e6), ' (95% CrI: ', sprintf('%.3g', sus.5/1e6), ', ', sprintf('%.3g', sus.95/1e6), ')')) %>% 
  mutate(sero = 'CA16') %>% select(10,1,8,9)



s.6.s <- pdata3 %>% 
  mutate(week.no = week,
         year = 2011 + floor((week.no - 1)/52),
         circle = case_when(year >= 2023 ~ 7,
                            year >= 2020 ~ 6,
                            year >= 2018 ~ 5,
                            year >= 2016 ~ 4,
                            year >= 2014 ~ 3,
                            year >= 2012 ~ 2,
                            year >= 2011 ~ 1),
         circle.week = case_when(year == 2011 ~ week.no + 52,
                                 year >= 2012 ~ (week.no + 52) - 52*2*(circle - 1),
                                 year >= 2020 ~ (week.no + 52) - 52*3*(circle - 1),
                                 year >= 2023 ~ (week.no + 52) - 52*2*(circle - 1)),
         circle.week = ifelse(year >= 2023, circle.week - 52, circle.week),
         circle = as.character(circle)) %>% 
  
  filter(circle != 7) %>% 
  
  mutate(period = ifelse(circle < 6, 'pre', 'post')) %>% 
  group_by(sample.no, circle) %>% slice_max(sus) %>% 
  group_by(period) %>% 
  summarise(sus.m = mean(sus), p.m = mean(sus/total),
            sus.5 = quantile(sus, 0.025), p.5 = quantile(sus/total, 0.025),
            sus.95 = quantile(sus, 0.975), p.95 = quantile(sus/total, 0.975)) %>% 
  mutate(p = paste0(sprintf('%.2f', 100*p.m), ' (95% CrI: ', sprintf('%.2f', 100*p.5), ', ', sprintf('%.2f', 100*p.95), ')'),
         sus = paste0(sprintf('%.3g', sus.m/1e6), ' (95% CrI: ', sprintf('%.3g', sus.5/1e6), ', ', sprintf('%.3g', sus.95/1e6), ')')) %>% 
  mutate(sero = 'CA6') %>% select(10,1,8,9)


sus_table <- rbind.data.frame(s.71.s, s.16.s, s.6.s)  
  

# write.xlsx(sus_table, paste0(file, 'sus.summary.xlsx'))

theme2 <- 
  theme(panel.background = element_rect(colour = "black", fill = 'white', size = 1),
        panel.grid.major.x = element_line(color = 'gray', linetype = "dashed"),
        axis.ticks.length = unit(0.4,"lines"),
        plot.margin = margin(0.2,0.4,0.2,0.2, "cm"),
        plot.title = element_text(colour = "black",size = 18, face = "bold"),
        axis.text = element_text(colour = "black",size = 18,hjust = 0,vjust = 0.5),
        axis.title = element_text(size = 17,vjust = 0, hjust = 0.5,face = "bold"),
        legend.position = 'bottom',
        legend.direction = "vertical",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.key = element_rect(fill = NA, color = NA),
        legend.key.height = unit(0.75, "cm"))

ps1 <-
  ggplot(data = s.71) +
  geom_line(aes(x = circle.week , y = 100*p.m, color = circle),size = 1 ) +
  geom_ribbon(aes(x = circle.week, ymin = 100*p.5, ymax = 100*p.95, fill = circle), alpha = 0.1) +
  scale_color_manual(values = brewer.pal(9,"Set1")[c(1,5,3,2,4)], name = 'Cycles',
                     labels = c('2005-2006', '2007-2010', '2011-2014', '2015-2018', '2019-2024')) +
  scale_fill_manual(values = brewer.pal(9,"Set1")[c(1,5,3,2,4)], name = 'Cycles',
                     labels = c('2005-2006', '2007-2010', '2011-2014', '2015-2018', '2019-2024')) +
  
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 10))) + # title.position = "left")) +
  scale_y_continuous(name = 'Susceptibles (%)') + 
  scale_x_continuous(name = NULL,  breaks = seq(0, 6*52, 52), labels = c(1:7)) +
  labs(title = 'A. EV-A71') +
  theme2

ps2 <- 
  ggplot(data = s.16) +
  geom_line(aes(x = circle.week , y = 100*p.m, color = circle),size = 1 ) +
  geom_ribbon(aes(x = circle.week, ymin = 100*p.5, ymax = 100*p.95, fill = circle), alpha = 0.1) +
  scale_color_manual(values = brewer.pal(9,"Set1")[c(1,5,3,2,4)], name = 'Cycles',
                     labels = c('2005-2008', '2009-2011', '2012-2015', '2016-2019', '2020-2024')) +
  scale_fill_manual(values = brewer.pal(9,"Set1")[c(1,5,3,2,4)], name = 'Cycles',
                     labels = c('2005-2008', '2009-2011', '2012-2015', '2016-2019', '2020-2024')) +
  
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 10))) + # title.position = "left")) +
  scale_y_continuous(name = NULL) + 
  scale_x_continuous(name = NULL,  breaks = seq(0, 6*52, 52), labels = c(1:7)) +
  labs(title = 'B. CVA16') +
  theme2

ps3 <- 
  ggplot(data = s.6) +
  geom_line(aes(x = circle.week , y = 100*p.m, color = circle),size = 1 ) +
  geom_ribbon(aes(x = circle.week, ymin = 100*p.5, ymax = 100*p.95, fill = circle), alpha = 0.1) +
  scale_color_manual(values = brewer.pal(9,"Set1")[c(1,5,3,2,4,7,8)], name = 'Cycles',
                     labels = c('2011', '2012-2013', '2014-2015', '2016-2017', '2018-2019', '2020-2022','2023-2024')) +
  scale_fill_manual(values = brewer.pal(9,"Set1")[c(1,5,3,2,4,7,8)], name = 'Cycles',
                     labels = c('2011', '2012-2013', '2014-2015', '2016-2017', '2018-2019', '2020-2022','2023-2024')) +
  
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 10))) + # title.position = "left")) +
  scale_y_continuous(name = NULL) + 
  scale_x_continuous(name = NULL,  breaks = seq(0, 6*52, 52), labels = c(1:7)) +
  labs(title = 'C. CVA6') +
  theme2


jpeg(paste0(file, "susceptible.tiff"), width = 12, height = 6, units = "in", res = 426)
ps1 + ps2 + ps3
dev.off()


# transmission -----


idx_data <- read.csv('data/OxCGRT_timeseries_StringencyIndex_v1.csv') %>% 
  filter(CountryName == 'Japan') %>% select(-(2:7)) %>% 
  melt(id.vars = 'CountryCode', variable.name = "date") %>% 
  mutate(Date = seq(as.Date('2020-01-01'),as.Date('2023-02-28'), by = '1 day'),
         week = lubridate::week(Date),
         week = ifelse(week == 53 , 52, week),
         week.no = (year(Date)-2005)*52 + week) %>%
  filter(! is.na(value)) %>% 
  group_by(week.no) %>% summarise(idx = mean(value))

mean_idx <- mean(idx_data$idx)


mcmc_cr <- function(data0, serotype){
  cr_data <- data0$samples_df %>% 
    mutate(cr = 1 / (1 + exp(-(c_0 - c_1 * mean_idx))),
           cr0 = 1 / (1 + exp(-(c_0))),
           redc = 100 * (1-cr),
           redc0 = 100 * (1-cr0),
           dif = redc - redc0,
           prop = 100 * (dif / redc),
           serotype = serotype) %>% 
    select(serotype, cr, cr0, redc, redc0, dif, prop) %>% 
    pivot_longer(cols = c('cr', 'cr0', 'redc', 'redc0', 'dif', 'prop') , names_to = 'cr')
  
  return(cr_data)
}


all_cr <- rbind.data.frame(mcmc_cr(data1, 'EV71'),
                           mcmc_cr(data2, 'CA16'),
                           mcmc_cr(data3, 'CA6')) %>% 
  mutate(serotype = factor(serotype, levels = c('EV71', 'CA16', 'CA6')))

# write.xlsx(all_cr, 'output/posterior cr.xlsx')


all_cr_sum <- all_cr %>% 
  group_by(serotype, cr) %>% 
  summarise(m = mean(value),
            low = quantile(value, 0.025), 
            up = quantile(value, 0.975)) %>% 
  mutate(across(where(is.numeric), ~ signif(., 3))) %>% 
  mutate(comb = paste0(m, ' (', low, ', ', up, ')')) 

sum_cr <- all_cr_sum %>% 
  filter(cr %in% c('cr', 'cr0')) %>% select(serotype, cr, comb) %>% 
  pivot_wider(names_from = 'serotype', values_from = 'comb') %>% 
  rename(variable = cr) %>% 
  rbind(pars_tab)

pars_tab_final <- sum_cr[c(3:11, 1,2,12), ]

# write.xlsx(pars_tab_final, paste0(file, 'all parameters.xlsx'))

sum_red <- all_cr_sum %>% 
  filter(! cr %in% c('cr', 'cr0')) %>% select(serotype, cr, comb) %>% 
  pivot_wider(names_from = 'cr', values_from = 'comb') %>% 
  select(1,4,5,2,3)

write.xlsx(sum_red, paste0(file, 'average redc summary.xlsx'))


all_cr_data <- all_cr %>% left_join(all_cr_sum %>% select(1:3)) %>% 
  filter(cr %in% c('cr', 'cr0'))

col <- brewer.pal(9,"Set1")[1:2]
labs <- c('COVID-19 Stringency Index dependent\naverage cr from 2020-2022', 
          'COVID-19 Stringency Index independent\naverage cr from 2020-2022')

post_cr <-
  ggplot(all_cr_data) +
  geom_density(aes(x = value, color = cr, fill = cr), alpha = 0.3) + 
  geom_vline(aes(xintercept = m, color = cr), linetype = "dashed", size = 1) +
  facet_wrap(~ serotype, scales = 'free') +
  labs(x = 'cr', y = 'Posterior distribution') + 
  scale_fill_manual(name = '', values = col, labels = labs) + 
  scale_color_manual(name = '', values = col, labels = labs) + 
  theme(panel.background = element_rect(colour = "black", fill = 'white', size = 1),
        panel.grid = element_blank(),
        panel.spacing = unit(0.5, "cm"),
        plot.margin = margin(0.1,0,0,0.1, "cm"),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x = unit(0.4,"lines"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black",size = 14,hjust = 0.5,vjust = 0.5, angle = 45),
        axis.title = element_text(size = 14,vjust = 0, hjust = 0.5,face = "bold"),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.key.size = unit(1, 'cm'),
        legend.spacing.x = unit(1, "cm"),
        strip.text = element_text(size = 13,face = 'bold', margin = margin(0, 0, 0, 0))
  )

jpeg(paste0(file, 'posterior cr.tiff'), width = 10, height = 4, units = "in", res = 426)
post_cr
dev.off()


# cr over time ----


idx_end_value <- 33.33


data.cr <- function(data0, sero){
  
  time_before_endidx <- (2022- 2005 + 1) * 52  # 2005-2022 end
  t_values <- seq(1, 52*2, 1) + time_before_endidx    # 2023 and 2024
  Bt_value <- rowMeans(splines2::iSpline(t_values, degree = 1, df = 5)) 
  Bt <- cbind.data.frame(week.no = t_values, Bt = Bt_value)
  
  mean_cr0 <- c()
  low_cr0 <- c()
  up_cr0 <- c()
  week.no0 <- c()
  
  a <- data0$samples_df
  
  for (i in (781:1040)) {
    
    week.no <- i
    
    if (i <= 936) {
      m_cr <- mean(1 / (1 + exp(-(a$c_0 - a$c_1 * idx_data$idx[idx_data$week.no == i]))))
      l_cr <- quantile(1 / (1 + exp(-(a$c_0 - a$c_1 * idx_data$idx[idx_data$week.no == i]))), 0.025)
      u_cr <- quantile(1 / (1 + exp(-(a$c_0 - a$c_1 * idx_data$idx[idx_data$week.no == i]))), 0.975)
      
    } else {
      m_cr <- mean(1 / (1 + exp(-(a$c_0 - a$c_1 * idx_end_value + a$c_2 * Bt$Bt[Bt$week.no == i]))))
      l_cr <- quantile(1 / (1 + exp(-(a$c_0 - a$c_1 * idx_end_value + a$c_2 * Bt$Bt[Bt$week.no == i]))), 0.025)
      u_cr <- quantile(1 / (1 + exp(-(a$c_0 - a$c_1 * idx_end_value + a$c_2 * Bt$Bt[Bt$week.no == i]))), 0.975)
    }
    
    mean_cr0 <- rbind(mean_cr0, m_cr)
    low_cr0 <- rbind(low_cr0, l_cr)
    up_cr0 <- rbind(up_cr0, u_cr)
    week.no0 <- rbind(week.no0, week.no)
    
    cr_data <- cbind.data.frame(week.no0, mean_cr0, low_cr0, up_cr0) %>% 
      mutate(serotype = sero)
    colnames(cr_data) <- c('week.no', 'mean_cr', 'low_cr', 'up_cr', 'serotype')
    
  }
  
  return(cr_data)
}

cr_data1 <- rbind.data.frame(data.cr(data1, 'EV-A71'), 
                             data.cr(data2, 'CVA16'))

data.cr.ca6 <- function(data0, sero){
  
  time_before_endidx <- (2022- 2011 + 1) * 52  # 2005-2022 end
  t_values <- seq(1, 52*2, 1) + time_before_endidx    # 2023 and 2024
  Bt_value <- rowMeans(splines2::iSpline(t_values, degree = 1, df = 5)) 
  Bt <- cbind.data.frame(week.no = t_values, Bt = Bt_value)
  
  idx_data_ca6 <- idx_data %>% mutate(week.no = week.no - 312)
  
  mean_cr0 <- c()
  low_cr0 <- c()
  up_cr0 <- c()
  week.no0 <- c()
  
  a <- data0$samples_df
  
  for (i in ((781:1040)-312)) {
    
    week.no <- i
    
    if (i <= (936-312)) {
      m_cr <- mean(1 / (1 + exp(-(a$c_0 - a$c_1 * idx_data_ca6$idx[idx_data_ca6$week.no == i]))))
      l_cr <- quantile(1 / (1 + exp(-(a$c_0 - a$c_1 * idx_data_ca6$idx[idx_data_ca6$week.no == i]))), 0.025)
      u_cr <- quantile(1 / (1 + exp(-(a$c_0 - a$c_1 * idx_data_ca6$idx[idx_data_ca6$week.no == i]))), 0.975)
      
    } else {
      m_cr <- mean(1 / (1 + exp(-(a$c_0 - a$c_1 * idx_end_value + a$c_2 * Bt$Bt[Bt$week.no == i]))))
      l_cr <- quantile(1 / (1 + exp(-(a$c_0 - a$c_1 * idx_end_value + a$c_2 * Bt$Bt[Bt$week.no == i]))), 0.025)
      u_cr <- quantile(1 / (1 + exp(-(a$c_0 - a$c_1 * idx_end_value + a$c_2 * Bt$Bt[Bt$week.no == i]))), 0.975)
    }
    
    mean_cr0 <- rbind(mean_cr0, m_cr)
    low_cr0 <- rbind(low_cr0, l_cr)
    up_cr0 <- rbind(up_cr0, u_cr)
    week.no0 <- rbind(week.no0, week.no)
    
    cr_data <- cbind.data.frame(week.no0, mean_cr0, low_cr0, up_cr0) %>% 
      mutate(serotype = sero)
    colnames(cr_data) <- c('week.no', 'mean_cr', 'low_cr', 'up_cr', 'serotype')
    
  }
  
  return(cr_data)
}


cr_data2 <- data.cr.ca6(data3, 'CVA6') %>% mutate(week.no = week.no + 312)
cr_data <- rbind.data.frame(cr_data1, cr_data2) %>% 
  mutate(serotype = factor(serotype, levels = c('EV-A71', 'CVA16', 'CVA6')))

date <- as.Date(c('2020/4/7','2020/5/25',      #soe1  14 21
                  '2021/1/8','2021/3/21',      #soe2  1  12
                  '2021/4/25','2021/9/30',
                  '2020/3/2', '2022/12/31'))    #soe3  17 39
week.value <- week(date)
week.no <- c(780+week.value[1], 780+week.value[2], 
             780+52+week.value[3], 780+52+week.value[4],
             780+52+week.value[5], 780+52+week.value[6],
             780+week.value[7], 780+52*2+week.value[8]-1)

fig.cr <-
  ggplot() + 
  annotate("rect",xmin=week.no[1],xmax = week.no[2],ymin=0,ymax=60,fill="grey95",lty=0.1)+
  annotate("rect",xmin=week.no[3],xmax = week.no[4],ymin=0,ymax=60,fill="grey95",lty=0.1)+
  annotate("rect",xmin=week.no[5],xmax = week.no[6],ymin=0,ymax=60,fill="grey95",lty=0.1)+
  
  geom_line(data = idx_data, aes(x = week.no, y = idx), color = 'gray30', size = 1.25) +
  geom_line(data = cr_data, aes(x = week.no, y = 100*(1-mean_cr), color = serotype), size = 1.25) +
  geom_ribbon(data = cr_data, aes(x = week.no, ymin = 100*(1-up_cr), ymax = 100*(1-low_cr), fill = serotype),
              alpha = 0.3) +
  
  annotate('text', x = week.no[8] + 2, y = 33.33, label = c('Stringency Index'),
           cex = 5, vjust = 0.5, hjust = 0) +
  
  annotate('text', x=(week.no[c(1,3,5)] + week.no[c(2,4,6)])/2, y=0,
           cex = 5, label=c('SoE1','SoE2','SoE3'), vjust = 0, hjust = 0.5) +
  
  scale_x_continuous(expand = c(0,6), breaks = seq(780, 1040, 52), # labels = c(2020:2025),
                     labels = c('Jan 2020',  'Jan 2021', 'Jan 2022', 'Jan 2023', 'Jan 2024',  'Jan 2025')) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~., name = 'COVID-19 Stringency Index')) +
  labs(x = '', y = "Transmission rate reduction %") +
  theme(panel.background = element_rect(colour = "black", fill = 'white', size = 1),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        plot.margin = margin(0.1,0,0,0.1, "cm"),
        plot.title = element_text(colour = "black",size = 14, face = "bold"),
        axis.text = element_text(colour = "black",size = 14,hjust = 0.5,vjust = 0.5),
        axis.title = element_text(size = 14,vjust = 0, hjust = 0.5,face = "bold"),
        legend.position = c(0.9,0.75),
        legend.title = element_blank(),
        legend.key = element_rect(colour = "white"),
        legend.text = element_text(size = 13),
        legend.key.size = unit(1, 'cm'))


jpeg(paste0(file, 'reduction.tiff'), width = 9, height = 4, units = "in", res = 426)
fig.cr
dev.off()


cr_data_sum <- cr_data %>% 
  mutate(stage = case_when(week.no %in% c(794:801) ~ 'soe1',
                           week.no %in% c(834:844) ~ 'soe2',
                           week.no %in% c(849:871) ~ 'soe3',
                           week.no > 936 ~ 'post-idx')) %>% 
  mutate(stage = ifelse(is.na(stage), 'non-soe', stage),
         mean_redx = 100*(1-mean_cr), low_redx = 100*(1-up_cr), up_redx = 100*(1-low_cr)) 


cr_data_sum.min <- cr_data_sum %>% 
  group_by(serotype, stage) %>% 
  slice_min(mean_redx) %>% 
  mutate(min_redx = paste0(signif(mean_redx,3), ' (', signif(low_redx,3), ', ', signif(up_redx,3), ')')) %>% 
  select(serotype, stage, min_redx)

cr_data_sum.max <- cr_data_sum %>% 
  group_by(serotype, stage) %>% 
  slice_max(mean_redx) %>% 
  mutate(max_redx = paste0(signif(mean_redx,3), ' (', signif(low_redx,3), ', ', signif(up_redx,3), ')')) %>% 
  select(serotype, stage,max_redx) %>% unique()


cr_data_sum.min.all <- cr_data_sum %>% 
  group_by(serotype) %>% 
  slice_min(mean_redx) %>% 
  mutate(min_redx = paste0(signif(mean_redx,3), ' (', signif(low_redx,3), ', ', signif(up_redx,3), ')'),
         stage = 'All') %>% 
  select(serotype, stage, min_redx)

cr_data_sum.max.all <- cr_data_sum %>% 
  group_by(serotype) %>% 
  slice_max(mean_redx) %>% 
  mutate(max_redx = paste0(signif(mean_redx,3), ' (', signif(low_redx,3), ', ', signif(up_redx,3), ')'),
         stage = 'All') %>% 
  select(serotype, stage, max_redx) %>% unique()

range.table <- rbind.data.frame(cr_data_sum.min.all, cr_data_sum.min) %>% 
  left_join(rbind.data.frame(cr_data_sum.max.all, cr_data_sum.max)) %>% 
  mutate(reduction = paste0(min_redx, ' - ', max_redx)) %>% 
  select(1,2,5) %>% 
  pivot_wider(names_from = serotype, values_from = reduction)

write.xlsx(range.table, paste0(file, 'reduction range by stages.xlsx'))


cr_data_sum.soe.mean <- cr_data_sum %>% 
  filter(stage != 'post-idx') %>% 
  mutate(stage = ifelse(stage == 'non-soe', stage, 'soe')) %>% 
  group_by(serotype, stage) %>% 
  summarise(redx_mean_value = mean(mean_redx),
            redx_5_value = quantile(mean_redx, 0.025),
            redx_95_value = quantile(mean_redx, 0.975)) %>% 
  mutate(reduced_trans = paste0(sprintf('%.3g', redx_mean_value), ' (95% CrI:', 
                                sprintf('%.3g', redx_5_value), ', ', 
                                sprintf('%.3g', redx_95_value), ')')) %>%
  select(serotype, stage, reduced_trans)

# write.xlsx(cr_data_sum.soe.mean, paste0(file, 'reduction soe-nonsoe.xlsx'))


