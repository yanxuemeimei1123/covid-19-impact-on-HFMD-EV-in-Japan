

time_len <- if (serotype == 'CA6') {14*52} else  {20*52}
start_year <- if (serotype == 'CA6') {2011} else  {2005}

initial_values <- if (serotype == 'CA6')
{c(50, 0.2, 0.5, 0.02, 2e-4, 5e-4, 100, 2, 0.01, 2)} else
{c(8.2, 0.01, 0.5, 0.08, 5e-5, 2e-4, 20, 2, 0.004, 0.02)}

vcv <- if (serotype == 'CA6')
{diag(c(10, 0.01, 0.01, 1e-4, 1e-8, 1e-8, 100, 0.1, 0.001, 0.1), 10)} else
{diag(c(1, 1e-4, 0.01, 1e-4, 1e-10, 1e-8, 10, 0.1, 1e-6, 1e-4), 10)}


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

