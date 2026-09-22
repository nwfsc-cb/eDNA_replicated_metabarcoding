n.chains <- 3
n.warm <- 10#500
n.samps <- 10#2000

stanMod = stan(file = here("Zenodo","Stan Simulation","1_marker",MOD_NAME),
               model_name = "Single-Marker",
               data = stan_data,
               verbose = FALSE,
               chains = n.chains,
               thin = 1,
               warmup = n.warm,
               iter = n.warm+n.samps,
               init_r = 1,
               control = list(max_treedepth=10,
                              adapt_init_buffer = 75),
               #                stepsize=0.01,
               #                adapt_delta=0.8,
               #                metric="diag_e"),
               #refresh = 100,
               boost_lib = NULL,
               sample_file="tmpW.csv"
)

pars <- rstan::extract(stanMod, permuted = TRUE)
names(pars)
# get_adaptation_info(stanMod)
samp_params <- get_sampler_params(stanMod)

#samp_params 
stanMod_summary <- NULL
stanMod_summary$log_D <- summary(stanMod,pars="log_D")$summary
stanMod_summary$log_F <- summary(stanMod,pars="log_F")$summary
stanMod_summary$log_F_mb1 <- summary(stanMod,pars="log_F_mb1")$summary
stanMod_summary$log_Lambda_copy <- summary(stanMod,pars="log_Lambda_copy")$summary
stanMod_summary$Lambda_copy <- summary(stanMod,pars="Lambda_copy")$summary
stanMod_summary$pi_samp_mb1 <- summary(stanMod,pars="pi_samp_mb1")$summary
#stanMod_summary$phi_param <- summary(stanMod,pars=c("beta0","beta1","log_gamma1","gamma1","gamma0","sigma_epsilon"))$summary
stanMod_summary$phi_param <- summary(stanMod,pars=c("beta0","beta1","log_gamma1","gamma1","gamma0","psi_epsilon"))$summary
stanMod_summary$phi_RE <- summary(stanMod,pars=c("epsilon"))$summary 
stanMod_summary$phi_val <- summary(stanMod,pars=c("log_phi_mb1","phi_mb1"))$summary
 
div <- 0
for(XX in 1:n.chains){
  count_divergences <- div + sum(samp_params[[XX]][(n.warm+1):nrow(samp_params[[XX]]),"divergent__"])
}

stanMod_sampling <- list(samp_params = samp_params,
                         count_divergences = count_divergences)