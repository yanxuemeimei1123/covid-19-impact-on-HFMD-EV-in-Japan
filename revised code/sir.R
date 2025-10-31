
# update for the next time step

update(N) <- S + I + R
update(S) <- S + n_birth - n_SI_event
update(I) <- I + n_SI - n_IR_event
update(R) <- R + n_IR - n_R_death
update(incidence) <- n_SI + incidence

# transition probabilities of (death and I) from S

beta <- if (time >= idx_start_time && time <= idx_end_time)
  cr1 * beta_0 * (1 + beta_s * cos(2 * pi * (time  - phi * period) / period)) else
    if (time > idx_end_time)
      cr2 * beta_0 * (1 + beta_s * cos(2 * pi * (time  - phi * period) / period)) else
        beta_0 * (1 + beta_s * cos(2 * pi * (time  - phi * period) / period))

foi <- beta * (I + w) / N

p_SI_event <- 1 - exp(-(foi + death_rate) * dt)  
p_SI_death <- death_rate / (foi + death_rate)

n_SI_event <- Binomial(S, p_SI_event)
n_SI_death <- if (n_SI_event < 1) 0 else Binomial(n_SI_event, p_SI_death)
n_SI <- n_SI_event - n_SI_death

p_IR_event <- 1 - exp(-(gamma + death_rate) * dt)  
p_IR_death <- death_rate / (death_rate + gamma)

n_IR_event <- Binomial(I, p_IR_event)
n_IR_death <- if (n_IR_event < 1) 0 else Binomial(n_IR_event, p_IR_death)
n_IR <- n_IR_event - n_IR_death

n_R_death <- Binomial(R, 1 - exp(-death_rate * dt)) 

n_birth <- Binomial(N, 1 - exp(-birth_rate * dt))

# Initial values

initial(N) <- N_0    
initial(S) <- floor(N_0 * s_0)
initial(I) <- floor(N_0 * i_0)
initial(R) <- N_0 - floor(N_0 * s_0) - floor(N_0 * i_0)
initial(incidence, zero_every = 1) <- 0

# Parameters unit daily

period <- parameter()     # period of seasonality 52 weeks
gamma <- parameter()      # recovery rate
phi <- parameter()

beta_0 <- parameter()
beta_s <- parameter()

N_0 <- parameter()        # total population in 2005 year
s_0 <- parameter()        # proportions of S
i_0 <- parameter()        # proportions of I
w <- parameter()          # imported case
idx_end_value <- parameter()
idx_end_time <- parameter()
idx_start_time <- parameter()   

birth_rate <- interpolate(birth_time, birth_value, "spline")
birth_time <- parameter(constant = TRUE)
birth_value <- parameter(constant = TRUE)
dim(birth_time, birth_value) <- parameter(rank = 1)

death_rate <- interpolate(death_time, death_value, "spline")
death_time <- parameter(constant = TRUE)
death_value <- parameter(constant = TRUE)
dim(death_time, death_value) <- parameter(rank = 1)

idx <- interpolate(idx_time, idx_value, "constant")
idx_time <- parameter(constant = TRUE)
idx_value <- parameter(constant = TRUE)
dim(idx_time, idx_value) <- parameter(rank = 1)

Bt_step <- interpolate(Bt_time, Bt_value, "constant")
Bt_time <- parameter(constant = TRUE)
Bt_value <- parameter(constant = TRUE)
dim(Bt_time, Bt_value) <- parameter(rank = 1)

c_0 <- parameter()
c_1 <- parameter()
c_2 <- parameter()

x1 <- c_0 - c_1 * idx
x2 <- c_0 - c_1 * idx_end_value + c_2 * Bt_step
cr1 <- 1 / (1 + exp(-x1))  # expit(x)
cr2 <- 1 / (1 + exp(-x2))  # expit(x)


# observation model
p <- parameter()
noise <- Exponential(rate = 1e6)
cases <- data()
cases ~ Poisson(incidence * p + noise)

