
pop_data <- read.xlsx('data/Japan population.xlsx', sheet = 'Sheet2') %>% 
  mutate(birth_rate = Live.births/(Total*52), death_rate = Deaths/(Total*52),
         birth = Live.births/(52), death = Deaths/(52)) %>%         
  mutate(time = case_when(Year == 2024 ~ (20*52),
                          Year > 2005 ~ ((Year-2005)*52+1),
                          Year == 2005 ~ 0)) %>% 
  dplyr::select(Year, time, Total, birth_rate, death_rate, birth, death)



idx_data <- read.csv('data/OxCGRT_timeseries_StringencyIndex_v1.csv') %>% 
  filter(CountryName == 'Japan') %>% dplyr::select(-(2:7)) %>% 
  melt(id.vars = 'CountryCode', variable.name = "date") %>% 
  mutate(Date = seq(as.Date('2020-01-01'),as.Date('2023-02-28'), by = '1 day'),
         week = week(Date),
         week = ifelse(week == 53 , 52, week),
         week.no = (year(Date)-2005)*52 + week) %>%
  filter(! is.na(value)) %>% 
  group_by(week.no) %>% summarise(idx = mean(value))
  
t_values <- seq(1, 52*2, 1)
Bt_value <- rowMeans(splines2::iSpline(t_values, degree = 1, df = 5)) 
# matplot(t_values, Bt_value, type = 'l')

sir_parameters <- function(serotype, initial_value) {
  
  if (serotype == 'CA6'){
    
    pop <- pop_data %>% filter(Year >= 2011) %>% mutate(time = time - 312)
    idx_start_time <- idx_data$week.no[1] - 312   #  from 2011
    idx_end_time <- idx_data$week.no[length(idx_data$week.no)] - 312   # from 2011
    idx0 <- rep(0, 468)       # idx from step 0-9*52  2011-2019
    Bt0 <- rep(0, 12*52)      # Bt from step 0-12*52  2011-2022
    
  } else {
    
    pop <- pop_data
    idx_start_time <- idx_data$week.no[1]  #  from 2005
    idx_end_time <- idx_data$week.no[length(idx_data$week.no)]
    idx0 <- rep(0, 780)       # idx from step 0-780
    Bt0 <- rep(0, 18*52)      # Bt from step 0-18*52  2005-2022
    
  }
  
  N_0 <- pop$Total[1]
  idx_end_value <- idx_data$idx[length(idx_data$idx)]
  
  # time-vary death and birth rate
  birth_rate_step <- as.data.frame(spline(pop$time, pop$birth_rate, n = (52*length(pop$Year))))$y
  death_rate_step <- as.data.frame(spline(pop$time, pop$death_rate, n = (52*length(pop$Year))))$y
  
  idx_step <- c(idx0, idx_data$idx)
  
  Bt_step <- c(Bt0, Bt_value)
  
  beta_0 <- initial_value[1]
  beta_s <- initial_value[2]
  s_0 <- initial_value[3]
  i_0 <- initial_value[4]
  p <- initial_value[5]
  w <- initial_value[6]
  c_0 <- initial_value[7]
  c_1 <- initial_value[8]
  c_2 <- initial_value[9]
  
  # Final list of all user-defined parameters, vectors and arrays
  list(beta_0 = beta_0,
       beta_s = beta_s,
       s_0 = s_0,
       i_0 = i_0,
       p = p,
       w = w, 
       c_0 = c_0,
       c_1 = c_1,
       c_2 = c_2,
       
       N_0 = N_0,
       idx_start_time = idx_start_time,
       idx_end_time = idx_end_time,
       idx_end_value = idx_end_value,
       idx_step = idx_step,
       Bt_step = Bt_step,
       birth_rate_step = birth_rate_step,
       death_rate_step = death_rate_step
       )
}

data_ob <- function(serotype){
  
  if (serotype == 'CA6'){
    reduce_length <- 312
    
  } else { 
    reduce_length <- 0
  }
  
  data_ob <- read.xlsx('data/hfmd derived ev by week.xlsx') %>% 
    filter(virus == serotype) %>% filter(!is.na(value)) %>% 
    arrange(week.no) %>% 
    mutate(week.no = week.no - reduce_length) %>% 
    dplyr::select(year, virus, week.no, value)
  
  return(data_ob)
}

