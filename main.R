
library(openxlsx)
library(readr)
library(dplyr)
library(coda)
library(lubridate)
library(reshape2)

gen_sir <- odin.dust::odin_dust('sir.R')

serotype <- 'EV71'
#                beta_0, beta_s,s_0,  i_0,   p, w, c0,c1,c2
initial_value <- c(8.5, 0.002, 0.08, 2e-4, 2e-4, 1, 2, 0, 0)
proposal_matrix <- diag(c(1, 1e-6, 1e-4, 1e-8, 1e-8, 10, 0.1, 0.01, 0.01), 9)

source('sir_parameters.R')

pars <- sir_parameters(serotype = serotype, initial_value = initial_value)
pars[c(1:13)]

data_observed <- data_ob(serotype = serotype)
data_particle <- mcstate::particle_filter_data(data_observed, time = "week.no", rate = 1, initial_time = 0)

source('sir_functions.R')

prior_list <- prior(initial_value)

mcmc_pars <- mcstate::pmcmc_parameters$new(prior_list, proposal_matrix, transform)

mcmc1 <- mcstate::pmcmc(mcmc_pars, filter = deterministic, control = control)
sum(mcmc(mcmc1))

mcmc2 <- pmcmc(mcmc1, filter = deterministic, n_steps = 5000)
sum(mcmc(mcmc2))

mcmc3 <- pmcmc(mcmc2, filter = deterministic, n_steps = 5000)
sum(mcmc(mcmc3))

mcmc4 <- pmcmc(mcmc3, filter = deterministic, n_steps = 5000)
sum(mcmc(mcmc4))

mcmc5 <- pmcmc(mcmc4, filter = deterministic, n_steps = 5000)
sum(mcmc(mcmc5))

plot(mcmc(mcmc5)[, c(3:12)])
mcmc_trace(as.data.frame(mcmc(mcmc5))) / mcmc_dens(as.data.frame(mcmc(mcmc5)))

# write_rds(mcmc5, paste0('1mcmc_', serotype, '.rds'))

pmcmc1 <- pmcmc(mcmc0, filter = stochastic, n_steps = 5000)
sum(mcmc(pmcmc1))
summary(mcmc(pmcmc1))

pmcmc2 <- pmcmc(pmcmc1, filter = stochastic, n_steps = 5000)
sum(mcmc(pmcmc2))
summary(mcmc(pmcmc2))

pmcmc3 <- pmcmc(pmcmc2, filter = stochastic, n_steps = 5000)
sum(mcmc(pmcmc3))
summary(mcmc(pmcmc3))

pmcmc4 <- pmcmc(pmcmc3, filter = stochastic, n_steps = 5000)
sum(mcmc(pmcmc4))
summary(mcmc(pmcmc4))
write_rds(pmcmc4, paste0('data/pmcmc_', serotype, '.rds'))


