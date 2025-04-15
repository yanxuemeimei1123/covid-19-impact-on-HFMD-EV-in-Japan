
index <- function(info) {
  list(run = c(I_case = info$index$I_case), 
       state = c(t = info$index$time,
                 N = info$index$N,
                 S = info$index$S,
                 I = info$index$I,
                 I_case = info$index$I_case))  
}


compare <- function(state, observed, pars = NULL) {
  exp_noise <- 1e6  # A small amount of noise is added to prevent zero expectations

  incidence_modelled <- state['I_case', , drop = TRUE]
  incidence_observed <- observed$value

  lambda <- incidence_modelled +
    rexp(n = length(incidence_modelled), rate = exp_noise)
  dpois(x = incidence_observed, lambda = lambda, log = TRUE)
}


prior <- function(initial_value){
  list(
    mcstate::pmcmc_parameter("beta_0", initial_value[1], prior = function(p)
      dunif(p, min = 0, max = 80, log = TRUE)),
    
    mcstate::pmcmc_parameter("beta_s", initial_value[2], min = 0, max = 1, prior = function(p)
      dbeta(p, 2, 10, log = TRUE)),
    
    mcstate::pmcmc_parameter("s_0", initial_value[3], min = 0, max = 1, prior = function(p)
      dbeta(p, 2, 10, log = TRUE)),
    
    mcstate::pmcmc_parameter("i_0", initial_value[4], min = 0, max = 1, prior = function(p)
      dbeta(p, 2, 100, log = TRUE)),
    
    mcstate::pmcmc_parameter("p", initial_value[5], min = 0, max = 1, prior = function(p)
      dbeta(p, 2, 100, log = TRUE)),
    
    mcstate::pmcmc_parameter("w", initial_value[6], prior = function(p)
      dunif(p, min = 0, max = 500, log = TRUE)),
    
    mcstate::pmcmc_parameter("c_0", initial_value[7], prior = function(p)
      dunif(p, min = -10, max = 10, log = TRUE)),
    
    mcstate::pmcmc_parameter("c_1", initial_value[8], prior = function(p)
      dunif(p, min = 0, max = 10, log = TRUE)),
    
    mcstate::pmcmc_parameter("c_2", initial_value[9], prior = function(p)
      dunif(p, min = 0, max = 10, log = TRUE))
    
  )
}
# transform

make_transform <- function(pars_list = pars_list, beta_0, beta_s, s_0, i_0, p, w, 
                           c_0, c_1, c_2, idx_step, Bt_step,
                           birth_rate_step, death_rate_step) {
  pars_list <- pars
  
  function(theta) {
    list(beta_0 = theta[["beta_0"]],
         beta_s = theta[["beta_s"]],
         s_0 = theta[["s_0"]],
         i_0 = theta[["i_0"]],
         p = theta[["p"]],
         w = theta[["w"]],
         c_0 = theta[["c_0"]],
         c_1 = theta[["c_1"]],
         c_2 = theta[["c_2"]],
         
         N_0 = pars_list$N_0,
         idx_end_value = pars_list$idx_end_value,
         idx_start_time = pars_list$idx_start_time,
         idx_end_time = pars_list$idx_end_time,
         idx_step = pars_list$idx_step,
         Bt_step = pars_list$Bt_step,
         birth_rate_step = pars_list$birth_rate_step,
         death_rate_step = pars_list$death_rate_step)
  }
}

transform <- make_transform(pars_list = pars_list, beta_0, beta_s, s_0, i_0, p, w, 
                            c_0, c_1, c_2, idx_step, Bt_step,
                            birth_rate_step, death_rate_step)

# filter

deterministic <- mcstate::particle_deterministic$new(data = data_particle, model = gen_sir,
                                                     compare = compare, index = index)

stochastic <- mcstate::particle_filter$new(data = data_particle, model = gen_sir,
                                           n_particles = 1000L, 
                                           compare = compare, index = index)

control <- mcstate::pmcmc_control(
  n_steps = 5000,
  save_state = F,
  save_trajectories = F,
  progress = TRUE,
  n_chains = 4,
  n_workers = 4,
  n_threads_total = 8,
  # use_parallel_seed = FALSE,
  rerun_every = 50,
  rerun_random = TRUE
)


pmcmc <- function(pmcmc0, filter, n_steps){
  
  i <- which.max(pmcmc0$probabilities[, "log_posterior"])
  new_initial <- pmcmc0$pars[i, ]

  prior <- prior(new_initial)
  
  proposal_matrix <- cov(pmcmc0$pars)
  
  mcmc_pars <- mcstate::pmcmc_parameters$new(prior, proposal_matrix, transform)
  
  control <- mcstate::pmcmc_control(
    n_steps = n_steps,
    save_state = F,
    save_trajectories = F,
    progress = TRUE,
    n_chains = 4,
    n_workers = 4,
    n_threads_total = 8,
    # use_parallel_seed = FALSE,
    rerun_every = n_steps/100,
    rerun_random = TRUE
  )
  
  pmcmc <- mcstate::pmcmc(mcmc_pars, filter = filter, control = control)
  return(pmcmc)
  
}


mcmc <- function(pmcmc){
  mcmc <- coda::as.mcmc(cbind(pmcmc$probabilities, pmcmc$pars))
  return(mcmc)
}

sum <- function(mcmc){
  es <- round(coda::effectiveSize(mcmc),4)
  ar <- round(1 - coda::rejectionRate(mcmc),4)
  s <- rbind(es,ar)
  return(s)
}


