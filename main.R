library(openxlsx)
library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(odin2)
library(dust2)
library(monty)
library(lubridate)

rm(list = ls())

sero_list <- c('EV71', 'CA16', 'CA6')
serotype <- sero_list[1]

gen_sir <- odin('sir.R')

source('parameters.R')


packer <- monty_packer(c("beta_0", 'beta_s', 'phi', 's_0', 'i_0', 'p', 'w', 'c_0', 'c_1', "c_2"),
                       fixed = fixed_pars)

prior <- monty_dsl({
  
  beta_0 ~ Uniform(min = 0, max = 80)
  beta_s ~ Uniform(min = 0, max = 1)
  phi ~ Uniform(min = 0, max = 1)
  s_0 ~ Uniform(min = 0, max = 0.5)
  i_0 ~ Uniform(min = 0, max = 0.01)
  p ~ Uniform(min = 0, max = 0.01)
  # eta ~ Exponential(mean = 1000)
  w ~ Uniform(min = 0, max = 1000)
  c_0 ~ Uniform(min = -10, max = 10)
  c_1 ~ Uniform(min = 0, max = 10)
  c_2 ~ Uniform(min = 0, max = 10)
  
})

runner <- monty_runner_callr(4)


filter <- dust_filter_create(gen_sir, data = data, time_start = 0, n_particles = 1000, n_threads = 4)

likelihood <- dust_likelihood_monty(filter, packer, save_trajectories = TRUE)
posterior <- likelihood + prior

sampler <- monty_sampler_random_walk(vcv, rerun_every = 50, rerun_random = TRUE)
samples <- monty_sample(posterior, sampler, initial = samples, n_steps = 5000, n_chains = 4, runner = runner)



