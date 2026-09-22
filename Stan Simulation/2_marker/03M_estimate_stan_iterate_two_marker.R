n.chains <- 3
n.warm <- 500
n.samps <- 1500

stanMod = stan(file = paste0("./Estimation/",MOD_NAME),
               model_name = "Multi-Marker",
               data = stan_data,
               verbose = FALSE,
               chains = n.chains,
               thin = 1,
               warmup = n.warm,
               iter = n.warm+n.samps,
               init_r = 1,
               control = list(max_treedepth=10,
                              adapt_init_buffer = 75,
                              adapt_delta=0.80),
               #                stepsize=0.01,
               
               #                metric="diag_e"),
               #refresh = 100,
               boost_lib = NULL,
               sample_file="tmpV.csv"
)
# IF YOU ARE LOOKING FOR LINES TO DEBUG:
# stanmod <- readLines('Estimation/Multi_marker_amp_fix_v7.stan')
# stanmod <- stanmod[stanmod != ""]
# ## what rstan calls each line 25:
# stanmod



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
stanMod_summary$phi_param <- summary(stanMod,pars=c("beta0","beta1","gamma0","gamma1","log_gamma1","sigma_epsilon"))$summary
stanMod_summary$phi_RE <- summary(stanMod,pars=c("epsilon"))$summary 
stanMod_summary$phi_val <- summary(stanMod,pars=c("log_phi_mb1","phi_mb1"))$summary
stanMod_summary$log_phi_mb1 <- summary(stanMod,pars="log_phi_mb1")$summary
stanMod_summary$log_phi_mb2 <- summary(stanMod,pars="log_phi_mb2")$summary
stanMod_summary$phi_mb1 <- summary(stanMod,pars="phi_mb1")$summary
stanMod_summary$phi_mb2 <- summary(stanMod,pars="phi_mb2")$summary

 
# traceplot(stanMod,pars=c("beta0","beta1"),inc_warmup=FALSE)
# traceplot(stanMod,pars=c("log_F[5,10]"))
# pairs(stanMod,pars=c("beta0","beta1","gamma0","gamma1"))


div <- 0
for(XX in 1:n.chains){
  div <- div + sum(samp_params[[XX]][(n.warm+1):nrow(samp_params[[XX]]),"divergent__"])
}
count_divergences <- div
  
stanMod_sampling <- list(samp_params = samp_params,
                         count_divergences = count_divergences)

