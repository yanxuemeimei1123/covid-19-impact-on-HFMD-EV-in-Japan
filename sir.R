# gen_sir <- odin.dust::odin_dust({
  
  # Parameters unit daily
  beta_0 <- user()     
  beta_s <- user()   
  phi <- user(29/52)       # phase of the transmission rate (the time of annual seasonal peak, AUG or SEP)
  period <- user(52)       # period of seasonality 52 weeks
  gamma <- user(1)         # recovery rate
  
  birth_rate_step[] <- user()
  dim(birth_rate_step) <- user()
  death_rate_step[] <- user()
  dim(death_rate_step) <- user()
  
  birth_rate <- if (step >= length(birth_rate_step)) 
    birth_rate_step[length(birth_rate_step)] else birth_rate_step[step + 1]
  death_rate <- if (step >= length(death_rate_step)) 
    death_rate_step[length(death_rate_step)] else death_rate_step[step + 1]
  
  # Definition of the time-step and output as "time"  unit = 1 day
  dt <- 1
  initial(time) <- 0
  update(time) <- (step + 1) * dt
  
  # Initial values
  s_0 <- user()              # proportions of S
  i_0 <- user()              # proportions of I
  N_0 <- user()              # total population in 2005 year
  
  initial(N) <- N_0    
  initial(S) <- floor(N_0 * s_0)
  initial(I) <- floor(N_0 * i_0)
  initial(R) <- N_0 - floor(N_0 * s_0) - floor(N_0 * i_0)
  
  idx_step[] <- user()
  dim(idx_step) <- user()
  idx <- if (step >= length(idx_step)) 
    idx_step[length(idx_step)] else idx_step[step + 1]
  
  idx_end_value <- user()
  idx_end_time <- user()
  idx_start_time <- user()   # school close

  # a vector at time t under the matrix of monotonic non-decreasing spline function
  Bt_step[] <- user()
  dim(Bt_step) <- user()
  Bt <- if (step >= length(Bt_step)) 
    Bt_step[length(Bt_step)] else Bt_step[step + 1]
  
  c_0 <- user()
  c_1 <- user()
  c_2 <- user()

  x1 <- c_0 - c_1 * idx
  x2 <- c_0 - c_1 * idx_end_value + c_2 * Bt
  cr1 <- 1 / (1 + exp(-x1))  # expit(x)
  cr2 <- 1 / (1 + exp(-x2))  # expit(x)


  beta <- if (time >= idx_start_time && time <= idx_end_time)
    cr1 * beta_0 * (1 + beta_s * cos(2 * 3.141593 * (time  - phi * period) / period)) else
      if (time > idx_end_time)
        cr2 * beta_0 * (1 + beta_s * cos(2 * 3.141593 * (time  - phi * period) / period)) else
          beta_0 * (1 + beta_s * cos(2 * 3.141593 * (time  - phi * period) / period))

  w <- user()    # imported case
  foi <- beta * (I + w) / N 

  
  # Two types of events for S, so competing hazards.
  p_SI_event <- 1 - exp(-(foi + death_rate) * dt)  # transition probabilities of (death and I) from S
  p_SI_death <- if (foi < 0) 1 else (death_rate / (foi + death_rate))
  
  n_events_S <- rbinom(S, p_SI_event)
  n_deaths_S <- if (n_events_S < 1) 0 else rbinom(n_events_S, p_SI_death)
  n_infections_S <- n_events_S - n_deaths_S
  
  # Two types of events for I, so competing hazards.
  p_IR_event <- 1 - exp(-(gamma + death_rate) * dt)  
  p_IR_death <- death_rate / (death_rate + gamma)
  
  n_events_I <- rbinom(I, p_IR_event)
  n_deaths_I <- if (n_events_I < 1) 0 else rbinom(n_events_I, p_IR_death)
  n_recoveries_I <- n_events_I - n_deaths_I
  
  # one type of events for R
  n_deaths_R <- rbinom(R, 1 - exp(-death_rate * dt)) 
  
  n_births <- rbinom(N, 1 - exp(-birth_rate * dt))
  
  # update for the next time step
  update(N) <- S + I + R
  update(S) <- S + n_births - n_events_S
  update(I) <- I + n_infections_S - n_events_I
  update(R) <- R + n_recoveries_I - n_deaths_R
  
  # Daily infections incidence
  update(I_inc) <- n_infections_S
  initial(I_inc) <- 0
  
  # reported case
  
  p <- user()        # reporting probability per infection
  update(I_case) <- n_infections_S * p
  initial(I_case) <- 0
  
  update(case_mol) <- rpois(n_infections_S * p)
  initial(case_mol) <- 0
  

# })
  
  