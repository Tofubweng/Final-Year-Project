library(odin.dust)
library(mcstate)
library(coda)
gen_sir <- odin.dust::odin_dust("odindust_alpha.R")

incidence <- read.csv("alphainfectionsandreinfections100days.csv")

dt <- 1
sir_data <- mcstate::particle_filter_data(data = incidence,time = "days",rate = 1 / dt)
rmarkdown::paged_table(sir_data)

plot(incidence$days, incidence$delta.infections, type = "p", xlab = "Days", ylab = "New cases")


ll_pois <- function(obs, model) {
  exp_noise <- 1e6
  if (is.na(obs)) {
    # Creates vector of zeros in ll with same length, if no data
    ll_obs <- numeric(length(model))
  } else {
    lambda <- model +
      rexp(n = length(model), rate = exp_noise)
    ll_obs <- dpois(x = obs, lambda = lambda, log = TRUE)
  }
  ll_obs
}

combined_compare <- function(state, observed, pars = NULL) {
  ll_alpha_infections <- ll_pois(observed$alpha.infections, state[6, , drop = TRUE])
  ll_other_infections <- ll_pois(observed$other.infections, state[10, , drop = TRUE])
  ll_alpha_reinfections <- ll_pois(observed$alpha.reinfections, state[22, , drop = TRUE])
  ll_other_reinfections <- ll_pois(observed$other.reinfections, state[18, , drop = TRUE])
  ll_deaths <- ll_pois(observed$DEATH_COUNT, state[28,,drop = TRUE])
  ll_alpha_infections + ll_other_infections + ll_alpha_reinfections + ll_other_reinfections+ll_deaths
}

n_particles <- 100
filter <- mcstate::particle_filter$new(data = sir_data,
                                       model = gen_sir,
                                       n_particles = n_particles,
                                       compare = combined_compare,
                                       seed = 1L)
filter$run(save_history = TRUE, pars = list())
filter$history()

gen_sir$new(pars = list(), time = 0, n_particles = 1L)$info()

plot_particle_filter <- function(history, true_history, times, obs_end = NULL) {
  if (is.null(obs_end)) {
    obs_end <- max(times)
  }
  
  par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
  cols <- c(I1 = "#8c8cd9", I2 = "#cc0044")
  matplot(times, t(history[10, , -1]), type = "l", xlab = "Time", ylab = "Number of individuals", col = cols[["I2"]],ylim = c(0,50000))
  matlines(times, t(history[6, , -1]), type = "l", col = cols[["I1"]], lty = 1)
  matpoints(times[1:obs_end], true_history[,13 , -1] , pch = 19, col = cols[["I1"]])
  matpoints(times[1:obs_end], true_history[,15 , -1] , pch = 19, col = cols[["I2"]])
  legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")
}

true_history <- incidence
plot_particle_filter(filter$history(), true_history, incidence$days)

beta1 <- mcstate::pmcmc_parameter("beta1", 0.5, min = 0, max = 1)
beta2 <- mcstate::pmcmc_parameter("beta2", 0.08, min = 0, max = 1)
beta3 <- mcstate::pmcmc_parameter("beta3", 0.3, min = 0, max = 1)
beta4 <- mcstate::pmcmc_parameter("beta4", 0.0495, min = 0, max = 1)
beta5 <- mcstate::pmcmc_parameter("beta5", 0.5, min = 0, max = 1)
beta6 <- mcstate::pmcmc_parameter("beta6", 0.08, min = 0, max = 1)
beta7 <- mcstate::pmcmc_parameter("beta7", 0.3, min = 0, max = 1)
beta8 <- mcstate::pmcmc_parameter("beta8", 0.0495, min = 0, max = 1)
alpha <- mcstate::pmcmc_parameter("alpha", 0.0555, min = 0, max = 1)
eta <- mcstate::pmcmc_parameter("eta", 0.0027, min = 0, max = 1)
epsilon_alpha <-  mcstate::pmcmc_parameter("epsilon_alpha", 0.0862, min = 0, max = 1)
rho <- mcstate::pmcmc_parameter("rho", 0.1, min = 0, max = 1)
epsilon_LV1 <- mcstate::pmcmc_parameter("epsilon_LV1", 0.895, min = 0, max = 1)
epsilon_LVA1 <- mcstate::pmcmc_parameter("epsilon_LVA1", 0.895, min = 0, max = 1)
epsilon_LV2 <- mcstate::pmcmc_parameter("epsilon_LV2", 0.519, min = 0, max = 1)
epsilon_LVA2 <- mcstate::pmcmc_parameter("epsilon_LVA2", 0.519, min = 0, max = 1)

omega <- mcstate::pmcmc_parameter("omega", 0.25, min = 0, max = 1)
p <- mcstate::pmcmc_parameter("p", 0.15, min = 0.01, max = 1)
q <- mcstate::pmcmc_parameter("q", 0.15, min = 0.01, max = 1)
cfr1 <- mcstate::pmcmc_parameter("cfr1", 0.01, min = 0, max = 0.05)
cfr2 <- mcstate::pmcmc_parameter("cfr2", 0.01, min = 0, max = 0.05)
delta1 <- mcstate::pmcmc_parameter("delta1", 0.1, min = 0, max = 1)
gamma1 <- mcstate::pmcmc_parameter("gamma1", 0.9, min = 0, max = 1)
delta2 <- mcstate::pmcmc_parameter("delta2", 0.1, min = 0, max = 1)
gamma2 <- mcstate::pmcmc_parameter("gamma2", 0.9, min = 0, max = 1)
epsilon_LR1 <- mcstate::pmcmc_parameter("epsilon_LR1", 0.519, min = 0, max = 1)
epsilon_LR1A <- mcstate::pmcmc_parameter("epsilon_LR1A", 0.519, min = 0, max = 1)
epsilon_LR2 <- mcstate::pmcmc_parameter("epsilon_LR2", 0.519, min = 0, max = 1)
epsilon_LR2A <- mcstate::pmcmc_parameter("epsilon_LR2A", 0.519, min = 0, max = 1)
delta12 <- mcstate::pmcmc_parameter("delta12", 0.1, min = 0, max = 1)
gamma12 <- mcstate::pmcmc_parameter("gamma12", 0.9, min = 0, max = 1)
delta21 <- mcstate::pmcmc_parameter("delta21", 0.1, min = 0, max = 1)
gamma21 <- mcstate::pmcmc_parameter("gamma21", 0.9, min = 0, max = 1)

proposal_matrix <- diag(0.01,33)
mcmc_pars <- mcstate::pmcmc_parameters$new(list(beta1 = beta1, beta2 = beta2, beta3 = beta3,
                                                beta4 = beta4, beta5 = beta5, beta6 = beta6,
                                                beta7 = beta7, beta8 = beta8, alpha = alpha,
                                                eta = eta, epsilon_alpha = epsilon_alpha,
                                                rho = rho, epsilon_LV1 = epsilon_LV1,
                                                epsilon_LVA1 = epsilon_LVA1, epsilon_LV2 = epsilon_LV2,
                                                epsilon_LVA2 = epsilon_LVA2, omega = omega, 
                                                p = p, q = q, cfr1 = cfr1, cfr2 = cfr2, delta1 = delta1,
                                                gamma1 = gamma1, delta2 = delta2, gamma2 = gamma2,
                                                epsilon_LR1 = epsilon_LR1, epsilon_LR1A = epsilon_LR1A,
                                                epsilon_LR2 = epsilon_LR2, epsilon_LR2A = epsilon_LR2A,
                                                delta12 = delta12, gamma12 = gamma12, delta21 = delta21,
                                                gamma21 = gamma21), proposal_matrix)
n_steps <- 500
n_burnin <- 200
control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE,
  save_trajectories = TRUE,
  progress = TRUE)
pmcmc_run <- mcstate::pmcmc(mcmc_pars, filter, control = control)

processed_chains <- mcstate::pmcmc_thin(pmcmc_run, burnin = n_burnin, thin = 2)
parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
parameter_mean_hpd

mcmc1 <- coda::as.mcmc(cbind(pmcmc_run$probabilities, pmcmc_run$pars))

summary(mcmc1)

plot(mcmc1)

plot_particle_filter(pmcmc_run$trajectories$state, true_history, incidence$days)

#alpha infections plot
matplot(incidence$days, t(pmcmc_run$trajectories$state[6, , -1]), type = "l", xlab = "Time", ylab = "Number of individuals")
matpoints(incidence$days, true_history[,13 , -1] , pch = 19, col = "#cc0044")

#Other infections plot 
matplot(incidence$days, t(pmcmc_run$trajectories$state[10, , -1]), type = "l", xlab = "Time", ylab = "Number of individuals")
matpoints(incidence$days, true_history[,15 , -1] , pch = 19, col = "#cc0044")

# alpha reinfections plot

matplot(incidence$days, t(pmcmc_run$trajectories$state[22, , -1]), type = "l", xlab = "Time", ylab = "Number of individuals")
matpoints(incidence$days, true_history[,14 , -1] , pch = 19, col = "#cc0044")

#Other reinfections plot

matplot(incidence$days, t(pmcmc_run$trajectories$state[18, , -1]), type = "l", xlab = "Time", ylab = "Number of individuals", ylim = c(0,50))
matpoints(incidence$days, true_history[,16 , -1] , pch = 19, col = "#cc0044")


#Deaths plot

matplot(incidence$days, t(pmcmc_run$trajectories$state[28, , -1]), type = "l", xlab = "Time", ylab = "Number of individuals")
matpoints(incidence$days, true_history[,6 , -1] , pch = 19, col = "#cc0044")
