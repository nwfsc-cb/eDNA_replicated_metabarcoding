
library(parallel)
cl<-makeCluster(3,outfile="") # Register cores

# initialization function 
nimble_inits <- function() {
  N_site <- nimble_constants$N_site
  N_sp <- nimble_constants$N_sp
  N_obs_mb1 <- nimble_constants$N_obs_mb1
  N_obs_mb2 <- nimble_constants$N_obs_mb2
  N_site_mock <- nimble_constants$N_site_mock
  N_run <- nimble_constants$N_run  # Number of sequencing runs
  
  # Initialize d_mb based on count_mb
  d_mb1_init <- as.matrix((nimble_data$count_mb1 > 0) * 1)
  d_mb2_init <- as.matrix((nimble_data$count_mb2 > 0) * 1)

  # Initialize mock detection states
  d_mock_mb1_init <- as.matrix((nimble_data$count_mock_mb1 > 0) * 1)
  d_mock_mb2_init <- as.matrix((nimble_data$count_mock_mb2 > 0) * 1)
  
  list(
    # log_D needs to be initialized with a fairly narrow range of values
    log_D = matrix(rnorm(N_site * N_sp, mean = 0, sd = 1), nrow = N_site, ncol = N_sp),
    
    # Run-specific Phi parameters (now shared across mock and unknown within runs)
    log_phi1_mb1 = rnorm(N_run, 1, 0.1),
    phi0_mb1 = runif(N_run, 0.8, 2),
    log_phi1_mb2 = rnorm(N_run, 1, 0.1),
    phi0_mb2 = runif(N_run, 0.8, 2),
    
    # Run-specific beta0s
    beta0_mb1 = runif(N_run, -0.5, 0),
    beta0_mb2 = runif(N_run, -0.5, 0),
    beta1_mb1 = runif(N_run, -10,-5),
    beta1_mb2 = runif(N_run, -10,-5),
    
    # Site-level standardized random effects
    epsilon_site_mb1 = matrix(rnorm(N_site*N_sp, 0, 0.0001),nrow=N_site,ncol=N_sp),
    epsilon_site_mb2 = matrix(rnorm(N_site*N_sp, 0, 0.0001),nrow=N_site,ncol=N_sp),

    # Variance parameters
    sigma_site_mb1 = runif(1, .2,1.3),
    sigma_site_mb2 = runif(1, .2,1.3),

    # Detection states
    d_mb1 = d_mb1_init,
    d_mb2 = d_mb2_init,
    d_mock_mb1 = d_mock_mb1_init,
    d_mock_mb2 = d_mock_mb2_init,
    
    # Amplification efficiency parameters
    amp_alpha_mb1 = rnorm(N_sp,0,0.01),
    amp_alpha_mb2 = rnorm(N_sp,0,0.01),
    amp_tau_mb1 = runif(1,0.005,0.015),
    amp_tau_mb2 = runif(1,0.015,0.025)
  )
}

run_Nimble_MCMC <- function(chain_info, nimble_data){
  library(nimble)
  
  # Custom beta-binomial distribution
  dbetabin <- nimbleFunction(
    run = function(x = integer(0), 
                   size = integer(0),
                   alpha = double(0), 
                   beta = double(0),
                   log_p_zero = double(0),   
                   log_p_pos = double(0),    
                   log_n_choose_k = double(0),
                   d_state = integer(0),      
                   log = integer(0, default = 0)) {
      returnType(double(0))
      
      if (x == 0) {
        if (d_state == 0) {
          if (log) return(0)
          else return(1)
        } else {
          log_bb_zero <- log_n_choose_k + 
            lgamma(size + beta) + lgamma(alpha + beta) -
            lgamma(alpha + beta + size) - lgamma(beta)
          log_prob <- log_bb_zero
        }
      } else {
        if (d_state == 0) {
          if (log) return(-Inf)
          else return(0)
        } else {
          log_bb_pos <- log_n_choose_k + 
            lgamma(x + alpha) + lgamma(size - x + beta) + lgamma(alpha + beta) - 
            lgamma(alpha + beta + size) - lgamma(alpha) - lgamma(beta)
          log_prob <- log_bb_pos
        }
      }
      
      if (log) return(log_prob)
      else return(exp(log_prob))
    }
  )
  # Random beta-binomial generator function
  rbetabin <- nimbleFunction(
    run = function(n = integer(0),
                   size = integer(0),
                   alpha = double(0), 
                   beta = double(0),
                   log_p_zero = double(0),
                   log_p_pos = double(0),
                   log_n_choose_k = double(0),
                   d_state = integer(0)) {
      returnType(integer(0))
      if (d_state == 0) {
        return(0)
      }
      if (runif(1) < exp(log_p_zero)) {
        return(0)
      } else {
        p <- rbeta(1, alpha, beta)
        return(rbinom(1, size, p))
      }
    }
  )
  
  registerDistributions(list(
    dbetabin = list(
      BUGSdist = "dbetabin(size, alpha, beta, log_p_zero, log_p_pos, log_n_choose_k, d_state)",
      types = c('value = integer(0)', 
                'size = integer(0)', 
                'alpha = double(0)', 
                'beta = double(0)', 
                'log_p_zero = double(0)', 
                'log_p_pos = double(0)', 
                'log_n_choose_k = double(0)',
                'd_state = integer(0)'),
      discrete = TRUE,
      pqAvail = FALSE
    )
  ))
  # Centered log-ratio only for detectable species
  clr_alpha <- nimbleFunction(
    run = function(amp_alpha = double(1),     # amplification efficiencies for all species (on log scale)
                   mb_ind = double(1)) {     # sample-level indicator for this sample
      returnType(double(1))
      
      n_sp <- length(amp_alpha)
      clr_result <- numeric(n_sp)
      
      present_alphas <- numeric(0)
      
      for(j in 1:n_sp) {
        if(mb_ind[j] == 1) {  # detectable 
          present_alphas <- c(present_alphas, amp_alpha[j])  # amp_alpha is already log scale
        }else{present_alphas <- present_alphas}
      }
      mean_log_alpha <- mean(present_alphas)
      
      # Calculate CLR for all species, but only meaningful for present ones
      for(j in 1:n_sp) {
        clr_result[j] <- amp_alpha[j] - mean_log_alpha  
      }
      return(clr_result)
    }
  )
  
  assign('dbetabin', dbetabin, envir = .GlobalEnv)
  assign('rbetabin', rbetabin, envir = .GlobalEnv)
  assign('clr_alpha', clr_alpha, envir = .GlobalEnv)
  
  inits <- chain_info$inits
  
  modelCode <- nimbleCode({
    
    # Priors for log_D
    for(m in 1:N_site) {
      for(j in 1:N_sp) {
        log_D[m,j] ~ dnorm(0, sd = 10)
      }
    }
    
    # Run-specific priors 
    for(r in 1:N_run) {
      # Parameters to estimate phi's
      beta0_mb1[r] ~ dnorm(0, sd = 1)
      beta0_mb2[r] ~ dnorm(0, sd = 1)
      beta1_mb1[r] ~ T(dnorm(-10, sd = 5), -Inf, 0)
      beta1_mb2[r] ~ T(dnorm(-10, sd = 5), -Inf, 0)
      phi0_mb1[r] ~ dnorm(1.5, sd = .25)
      phi0_mb2[r] ~ dnorm(1.5, sd = .25)
      log_phi1_mb1[r] ~ dnorm(1.25, sd = .25)
      log_phi1_mb2[r] ~ dnorm(1.25, sd = .25)
      phi1_mb1[r] <- exp(log_phi1_mb1[r])
      phi1_mb2[r] <- exp(log_phi1_mb2[r])
      
    }
    # Standardized site-level random effects (mapped to appropriate run)
    for(m in 1:N_site) {
      for(j in 1:N_sp){
        epsilon_site_mb1[m,j] ~ dnorm(0, sigma_site_mb1)
        epsilon_site_mb2[m,j] ~ dnorm(0, sigma_site_mb2)
        # Derived site-level effects
        # epsilon_site_mb1[m,j] <- epsilon_site_mb1_z[m,j] * sigma_site_mb1
        # epsilon_site_mb2[m,j] <- epsilon_site_mb2_z[m,j] * sigma_site_mb2
      }
    }
    
    sigma_site_mb1 ~ T(dnorm(0, sd = 1), 0, Inf)
    sigma_site_mb2 ~ T(dnorm(0, sd = 1), 0, Inf)
    amp_tau_mb1 ~ T(dnorm(0, sd = 0.01), 0, Inf)
    amp_tau_mb2 ~ T(dnorm(0, sd = 0.01), 0, Inf)
    
    # Amplification efficiency priors 
    for(j in 1:(N_sp-1)) {
      amp_alpha_mb1[j] ~ dnorm(0, amp_tau_mb1)
      amp_alpha_mb2[j] ~ dnorm(0, amp_tau_mb2)
      # amp_alpha_mb1_mock[j] <- amp_alpha_mb1[j] 
      # amp_alpha_mb2_mock[j] <- amp_alpha_mb2[j] 
    }
    amp_alpha_mb1[N_sp] <- 0
    amp_alpha_mb2[N_sp] <- 0
    # amp_alpha_mb1_mock[N_sp] <- 0
    # amp_alpha_mb2_mock[N_sp] <- 0
    # 
    # Deterministic assignment for unknown samples
    # amp_alpha_mb1[1:N_sp] <- amp_alpha_mb1_mock[1:N_sp]
    # amp_alpha_mb2[1:N_sp] <- amp_alpha_mb2_mock[1:N_sp]
    

    ########### "Known" Mock Communities (for estimating alphas) ###########
    # Site-level calculations 
    for(m in 1:N_site_mock) {
      clr_alpha_mock_mb1[m,1:N_sp] <- clr_alpha(amp_alpha_mb1[1:N_sp], MB1_mock_all[1:N_sp])
      clr_alpha_mock_mb2[m,1:N_sp] <- clr_alpha(amp_alpha_mb2[1:N_sp], MB2_mock_all[1:N_sp])
      for(j in 1:N_sp) {
        # Marker 1
        log_phi_mock_mb1[m,j] <- beta0_mb1[seq_run_mock_mb1[m]] + 
          beta1_mb1[seq_run_mock_mb1[m]] * clr_alpha_mock_mb1[m,j] +
          log_Lambda_site_mock_mb1[m] +
          phi0_mb1[seq_run_mock_mb1[m]] * exp(-phi1_mb1[seq_run_mock_mb1[m]] * log_K_site_mock_mb1[m, j])
        phi_mock_mb1[m,j] <- exp(min(20, log_phi_mock_mb1[m,j]))
        # Marker 2
        log_phi_mock_mb2[m,j] <- beta0_mb2[seq_run_mock_mb2[m]] +
          beta1_mb2[seq_run_mock_mb1[m]] * clr_alpha_mock_mb2[m,j] +
          log_Lambda_site_mock_mb2[m] +
          phi0_mb2[seq_run_mock_mb2[m]] * exp(-phi1_mb2[seq_run_mock_mb2[m]] * log_K_site_mock_mb2[m, j])
        phi_mock_mb2[m,j] <- exp(min(20, log_phi_mock_mb2[m,j]))
      }
    }
    # Known Mock community likelihood - Marker 1
    for(ii in 1:N_obs_mock_mb1) {
      amp_alpha_mb1_ref_mock[ii] <- amp_alpha_mb1[ref_sp_mock_mb1[ii]]
      for(j in 1:N_sp) {
        log_val_mock_mb1[ii, j] <- alr_mock_mb1[ii, j] + N_pcr * (amp_alpha_mb1[j] - amp_alpha_mb1_ref_mock[ii])
        d_mock_mb1[ii,j] ~ dbern(exp(log_1mexp_mock_mb1[ii,j]))
      }
      log_max_val_mock_mb1[ii] <- max(log_val_mock_mb1[ii,1:N_sp])
      for(j in 1:N_sp) {
        exp_shifted_mock_mb1[ii,j] <- exp(log_val_mock_mb1[ii,j] - log_max_val_mock_mb1[ii]) * d_mock_mb1[ii,j]
      }
      sum_exp_shifted_mock_mb1[ii] <- sum(exp_shifted_mock_mb1[ii,1:N_sp])
      for(j in 1:N_sp) {
        pi_mock_mb1_cond[ii,j] <- exp_shifted_mock_mb1[ii,j] / sum_exp_shifted_mock_mb1[ii]
        alpha_mock_mb1[ii,j] <- pi_mock_mb1_cond[ii,j] * phi_mock_mb1[mock_map[ii],j]
        beta_mock_mb1[ii,j] <- phi_mock_mb1[mock_map[ii],j] - alpha_mock_mb1[ii,j]
        count_mock_mb1[ii,j] ~ dbetabin(
          reads_mock_mb1[ii], alpha_mock_mb1[ii,j], beta_mock_mb1[ii,j],
          log_p_zero_mock_mb1[ii,j], log_p_pos_mock_mb1[ii,j],
          log_mock_mb1_n_k[ii,j], d_mock_mb1[ii,j]
        )
      }
    }
    
    # Known Mock community likelihood - Marker 2
    for(ii in 1:N_obs_mock_mb2) {
      amp_alpha_mb2_ref_mock[ii] <- amp_alpha_mb2[ref_sp_mock_mb2[ii]]
      for(j in 1:N_sp) {
        log_val_mock_mb2[ii, j] <- alr_mock_mb2[ii, j] + N_pcr * (amp_alpha_mb2[j] - amp_alpha_mb2_ref_mock[ii])
        d_mock_mb2[ii,j] ~ dbern(exp(log_1mexp_mock_mb2[ii,j]))
      }
      log_max_val_mock_mb2[ii] <- max(log_val_mock_mb2[ii,1:N_sp])
      for(j in 1:N_sp) {
        exp_shifted_mock_mb2[ii,j] <- exp(log_val_mock_mb2[ii,j] - log_max_val_mock_mb2[ii]) * d_mock_mb2[ii,j]
      }
      sum_exp_shifted_mock_mb2[ii] <- sum(exp_shifted_mock_mb2[ii,1:N_sp])
      for(j in 1:N_sp) {
        pi_mock_mb2_cond[ii,j] <- exp_shifted_mock_mb2[ii,j] / sum_exp_shifted_mock_mb2[ii]
        alpha_mock_mb2[ii,j] <- pi_mock_mb2_cond[ii,j] * phi_mock_mb2[mock_map[ii],j]
        beta_mock_mb2[ii,j] <- phi_mock_mb2[mock_map[ii],j] - alpha_mock_mb2[ii,j]
        count_mock_mb2[ii,j] ~ dbetabin(
          reads_mock_mb2[ii], alpha_mock_mb2[ii,j], beta_mock_mb2[ii,j],
          log_p_zero_mock_mb2[ii,j], log_p_pos_mock_mb2[ii,j],
          log_mock_mb2_n_k[ii,j], d_mock_mb2[ii,j]
        )
      }
    }
    ############  Unknown Communities ##############
    # Matrix operations for unknown samples
    log_F[1:N_site_samp, 1:N_sp] <- X_D_to_F[1:N_site_samp, 1:N_site] %*% log_D[1:N_site, 1:N_sp]
    # log_F_mb1[1:N_obs_mb1, 1:N_sp] <- X_MB1_ss[1:N_obs_mb1, 1:N_site_samp] %*% log_F[1:N_site_samp, 1:N_sp]
    # log_F_mb2[1:N_obs_mb2, 1:N_sp] <- X_MB2_ss[1:N_obs_mb2, 1:N_site_samp] %*% log_F[1:N_site_samp, 1:N_sp]
    
    # MB1 calculations 
    for(m in 1:N_site_samp) {
      # Calculate log_lambda_K values first
      for(j in 1:N_sp) {
        # Cache all exponentials and logs at once
        temp_exp_mb1[m,j] <- exp(log_F[m,j])
        log_1mexp_mb1[m,j] <- log(1 - exp(-temp_exp_mb1[m,j]))
        log_p_pos_mb1[m,j] <- log_1mexp_mb1[m,j] # log probability of presence in aliquot
        log_p_zero_mb1[m,j] <- -temp_exp_mb1[m,j] # log probability of absence in aliquot
        lambda_mb1[m,j] <- MB1_all[j] * temp_exp_mb1[m,j] # log copies
        log_lambda_K_mb1[m,j] <- log_F[m,j] - log_1mexp_mb1[m,j] # log expected copies in aliquot (from the zero-truncated Poisson)
      }
      
    }
    for(i in 1:N_obs_mb1) {
      # extract reference values
      log_lambda_K_mb1_ref[i] <- log_lambda_K_mb1[samp_map1[i],ref_sp_mb1[i]]
      amp_alpha_mb1_ref[i] <- amp_alpha_mb1[ref_sp_mb1[i]]
      # log-ratios
      for(j in 1:N_sp) {
        log_val_mb1[i,j] <- (log_lambda_K_mb1[samp_map1[i],j] - log_lambda_K_mb1_ref[i]) + 
          N_pcr * (amp_alpha_mb1[j] - amp_alpha_mb1_ref[i])
        # Discrete latent presence in aliquot
        d_mb1[i,j] ~ dbern(exp(log_1mexp_mb1[samp_map1[i],j]))
        # Copies in aliquot (adjusted for presence)
        lambda_K_d_mb1[i,j] <- d_mb1[i,j] * MB1_all[j] * exp(log_lambda_K_mb1[samp_map1[i],j])
      }
      log_max_val_mb1[i] <- max(log_val_mb1[i,1:N_sp])
    }
    
    for(m in 1:N_site_samp) {
      # Total copies in sample
      log_Lambda_mb1[m] <- log(sum(lambda_mb1[m,1:N_sp])) 
      clr_alpha_mb1[m,1:N_sp] <- clr_alpha(amp_alpha_mb1[1:N_sp], 
                                           MB1_all[1:N_sp])
      for(j in 1:N_sp) {
        # Calculate phi_mb
        log_phi_mb1[m,j] <- beta0_mb1[seq_run_mb1_sites[m]] + 
          beta1_mb1[seq_run_mb1_sites[m]] * clr_alpha_mb1[m,j] +
          log_Lambda_mb1[m] +
          phi0_mb1[seq_run_mb1_sites[m]] * exp(-phi1_mb1[seq_run_mb1_sites[m]] * log_lambda_K_mb1[m,j]) +
          epsilon_site_mb1[m,j] 
        phi_mb1[m,j] <- exp(min(20,log_phi_mb1[m,j]))
      }
    }
    for(i in 1:N_obs_mb1) {
      for(j in 1:N_sp) {
        # Apply discrete state masking to ALR values
        exp_shifted_mb1[i,j] <- exp(log_val_mb1[i,j] - log_max_val_mb1[i]) * d_mb1[i,j]
      }
      sum_exp_shifted_mb1[i] <- sum(exp_shifted_mb1[i,1:N_sp])
      for(j in 1:N_sp) {
        # Softmax to calculate pi
        pi_samp_mb1_cond[i,j] <- exp_shifted_mb1[i,j] / sum_exp_shifted_mb1[i]
        # Calculate alpha and beta params of beta-binomial from phi and pi
        alpha_mb1[i,j] <- pi_samp_mb1_cond[i,j] * phi_mb1[samp_map1[i],j]
        beta_mb1[i,j] <- phi_mb1[samp_map1[i],j] - alpha_mb1[i,j]
        # The Likelihood
        count_mb1[i,j] ~ dbetabin(
          reads_mb1[i], alpha_mb1[i,j], beta_mb1[i,j],
          log_p_zero_mb1[samp_map1[i],j], log_p_pos_mb1[samp_map1[i],j],
          log_mb1_n_k[i,j],d_mb1[i,j]
        )
      }
    }
    
    # Marker 2 calculations 
    for(m in 1:N_site_samp) {
      # Calculate log_lambda_K values first
      for(j in 1:N_sp) {
        # Cache all exponentials and logs at once
        temp_exp_mb2[m,j] <- exp(log_F[m,j])
        log_1mexp_mb2[m,j] <- log(1 - exp(-temp_exp_mb2[m,j]))
        log_p_pos_mb2[m,j] <- log_1mexp_mb2[m,j] # log probability of presence in aliquot
        log_p_zero_mb2[m,j] <- -temp_exp_mb2[m,j] # log probability of absence in aliquot
        lambda_mb2[m,j] <- MB2_all[j] * temp_exp_mb2[m,j] # log copies
        log_lambda_K_mb2[m,j] <- log_F[m,j] - log_1mexp_mb2[m,j] # log expected copies in aliquot (from the zero-truncated Poisson)
      }
    }
    for(i in 1:N_obs_mb2) {
      # extract reference values
      log_lambda_K_mb2_ref[i] <- log_lambda_K_mb2[samp_map2[i],ref_sp_mb2[i]]
      amp_alpha_mb2_ref[i] <- amp_alpha_mb2[ref_sp_mb2[i]]
      for(j in 1:N_sp) {
        log_val_mb2[i,j] <- (log_lambda_K_mb2[samp_map2[i],j] - log_lambda_K_mb2_ref[i]) + 
          N_pcr * (amp_alpha_mb2[j] - amp_alpha_mb2_ref[i])
        # Discrete latent presence in aliquot
        d_mb2[i,j] ~ dbern(exp(log_1mexp_mb2[samp_map2[i],j]))
        # Copies in aliquot (adjusted for presence)
        lambda_K_d_mb2[i,j] <- d_mb2[i,j] * MB2_all[j] * exp(log_lambda_K_mb2[samp_map2[i],j])
      }
      log_max_val_mb2[i] <- max(log_val_mb2[i,1:N_sp])
    }
    for(m in 1:N_site_samp) {
      # Total copies in sample
      log_Lambda_mb2[m] <- log(sum(lambda_mb2[m,1:N_sp])) 
      clr_alpha_mb2[m,1:N_sp] <- clr_alpha(amp_alpha_mb2[1:N_sp], 
                                           MB2_all[1:N_sp])
      # Apply discrete state masking to log values
      
      for(j in 1:N_sp) {
        # Calculate phi_mb
        log_phi_mb2[m,j] <- beta0_mb2[seq_run_mb2_sites[m]] + 
          beta1_mb2[seq_run_mb2_sites[m]] * clr_alpha_mb2[m,j] +
          log_Lambda_mb2[m] +
          phi0_mb2[seq_run_mb2_sites[m]] * exp(-phi1_mb2[seq_run_mb2_sites[m]] * log_lambda_K_mb2[m,j]) +
          epsilon_site_mb2[m,j] 
        phi_mb2[m,j] <- exp(min(20,log_phi_mb2[m,j]))
      }
    }
    for(i in 1:N_obs_mb2) {
      for(j in 1:N_sp) {
        # Apply discrete state masking to ALR values
        exp_shifted_mb2[i,j] <- exp(log_val_mb2[i,j] - log_max_val_mb2[i]) * d_mb2[i,j]
      }
      sum_exp_shifted_mb2[i] <- sum(exp_shifted_mb2[i,1:N_sp])
      for(j in 1:N_sp) {
        # Softmax to calculate pi
        pi_samp_mb2_cond[i,j] <- exp_shifted_mb2[i,j] / sum_exp_shifted_mb2[i]
        # Calculate alpha and beta params of beta-binomial from phi and pi
        alpha_mb2[i,j] <- pi_samp_mb2_cond[i,j] * phi_mb2[samp_map2[i],j]
        beta_mb2[i,j] <- phi_mb2[samp_map2[i],j] - alpha_mb2[i,j]
        # The Likelihood
        count_mb2[i,j] ~ dbetabin(
          reads_mb2[i], alpha_mb2[i,j], beta_mb2[i,j],
          log_p_zero_mb2[samp_map2[i],j], log_p_pos_mb2[samp_map2[i],j],
          log_mb2_n_k[i,j],d_mb2[i,j]
        )
      }
    }
  })
  
  model <- nimbleModel(code = modelCode,
                       data = nimble_data,
                       constants = nimble_constants,
                       inits = inits,
                       check = TRUE,
                       calculate=FALSE)
  
  # Clean up any leftovers from previous runs
  if (exists("cmodel")) {try(nimble:::clearCompiled(cmodel), silent = TRUE)}
  if (exists("cmcmc")) {try(nimble:::clearCompiled(cmcmc), silent = TRUE)}
  gc()  
  
  # Compile
  cmodel <- compileNimble(model)
  
  # Configure MCMC 
  mcmc_conf <- configureMCMC(model, useConjugacy=FALSE)
  
  # Remove default samplers for log_D
  mcmc_conf$removeSamplers("log_D")
  # mcmc_conf$removeSamplers("epsilon_site_mb1")
  # mcmc_conf$removeSamplers("epsilon_site_mb2")
  # mcmc_conf$removeSamplers("sigma_site_mb1")
  # mcmc_conf$removeSamplers("sigma_site_mb2")
  mcmc_conf$removeSamplers("amp_tau_mb1")
  mcmc_conf$removeSamplers("amp_tau_mb2")
  # mcmc_conf$removeSamplers(c(paste0("amp_alpha_mb1_raw[1:", nimble_constants$N_sp-1, "]")))
  # mcmc_conf$removeSamplers(c(paste0("amp_alpha_mb2_raw[1:", nimble_constants$N_sp-1, "]")))
  
  # Add block samplers for log_D by site
  for (i in 1:nimble_constants$N_site){
    targets <- paste0("log_D[", i, ", 1:", nimble_constants$N_sp, "]")
    mcmc_conf$addSampler(target=targets, "AF_slice")
  }
  
  # Add back multivariate samplers 
  # mcmc_conf$addSampler(target=c("epsilon_site_mb1"), "AF_slice")
  # mcmc_conf$addSampler(target=c("epsilon_site_mb2"), "AF_slice")
  # mcmc_conf$addSampler(target=c(paste0("epsilon_site_mb2_z[1:", nimble_constants$N_site, "]")), "AF_slice")
  # mcmc_conf$addSampler(target=c("sigma_site_mb1[2:3]", paste0("epsilon_site_mb1_z[1:", nimble_constants$N_site, "]")), "AF_slice")
  # mcmc_conf$addSampler(target=c("sigma_site_mb2[2:3]", paste0("epsilon_site_mb2_z[1:", nimble_constants$N_site, "]")), "AF_slice")
  
  # Change samplers for run-specific parameters
  for (s in 1:nimble_constants$N_site_samp) {
    run_targets_mb1 <- c(paste0("epsilon_site_mb1[", s, ",]"))
    run_targets_mb2 <- c(paste0("epsilon_site_mb2[", s, ",]"))
    mcmc_conf$addSampler(target=run_targets_mb1, "AF_slice")
    mcmc_conf$addSampler(target=run_targets_mb2, "AF_slice")
  }
  
  mcmc_conf$addSampler(target=c("amp_tau_mb1",paste0("amp_alpha_mb1[1:", nimble_constants$N_sp-1, "]")), "AF_slice")
  mcmc_conf$addSampler(target=c("amp_tau_mb2",paste0("amp_alpha_mb2[1:", nimble_constants$N_sp-1, "]")), "AF_slice")
  
  # Change samplers for run-specific parameters
  for (r in 1:nimble_constants$N_run) {
    run_targets_mb1 <- c(paste0("beta0_mb1[", r, "]"), paste0("beta1_mb1[", r, "]"),paste0("phi0_mb1[", r, "]"), paste0("log_phi1_mb1[", r, "]"))
    run_targets_mb2 <- c(paste0("beta0_mb2[", r, "]"), paste0("beta1_mb2[", r, "]"),paste0("phi0_mb2[", r, "]"), paste0("log_phi1_mb2[", r, "]"))
    mcmc_conf$removeSamplers(run_targets_mb1)
    mcmc_conf$removeSamplers(run_targets_mb2)
    mcmc_conf$addSampler(target=run_targets_mb1, "AF_slice")
    mcmc_conf$addSampler(target=run_targets_mb2, "AF_slice")
  }
  
  # Remove unneeded samplers (i.e., where d_mb==1)
  fixed_idx1 <- which(nimble_data$count_mb1 > 0, arr.ind = TRUE)
  for(k in 1:nrow(fixed_idx1)) {
    i <- fixed_idx1[k,1]
    j <- fixed_idx1[k,2]
    mcmc_conf$removeSamplers(paste0("d_mb1[",i,",",j,"]"))
  }
  fixed_idx2 <- which(nimble_data$count_mb2 > 0, arr.ind = TRUE)
  for(k in 1:nrow(fixed_idx2)) {
    i <- fixed_idx2[k,1]
    j <- fixed_idx2[k,2]
    mcmc_conf$removeSamplers(paste0("d_mb2[",i,",",j,"]"))
  }
  fixed_idx1 <- which(nimble_data$count_mock_mb1 > 0, arr.ind = TRUE)
  for(k in 1:nrow(fixed_idx1)) {
    i <- fixed_idx1[k,1]
    j <- fixed_idx1[k,2]
    mcmc_conf$removeSamplers(paste0("d_mock_mb1[",i,",",j,"]"))
  }
  fixed_idx2 <- which(nimble_data$count_mock_mb2 > 0, arr.ind = TRUE)
  for(k in 1:nrow(fixed_idx2)) {
    i <- fixed_idx2[k,1]
    j <- fixed_idx2[k,2]
    mcmc_conf$removeSamplers(paste0("d_mock_mb2[",i,",",j,"]"))
  }
  
  # Add monitors
  # mcmc_conf$addMonitors("log_phi_mb1")
  # mcmc_conf$addMonitors("log_phi_mb2")
  # mcmc_conf$addMonitors("log_phi_mock_mb1")
  # mcmc_conf$addMonitors("log_phi_mock_mb2")
  mcmc_conf$addMonitors("amp_alpha_mb1")
  mcmc_conf$addMonitors("amp_alpha_mb2")
  mcmc_conf$addMonitors("epsilon_site_mb1")
  mcmc_conf$addMonitors("epsilon_site_mb2")

  # Build and compile MCMC
  mcmc <- buildMCMC(mcmc_conf)
  cmcmc <- compileNimble(mcmc, project = model)
 
  # Run MCMC
  system.time(samples <- runMCMC(cmcmc, niter = 24000 , nburnin = 4000, thin = 1))
  return(samples)
}

chain_info <- list(
  list(inits = nimble_inits()),
  list(inits = nimble_inits()),
  list(inits = nimble_inits())
)

clusterExport(cl, c("nimble_data", "nimble_constants", "chain_info"))

chain_output <- parLapply(
  cl = cl, X = chain_info,
  fun = run_Nimble_MCMC, 
  nimble_data = nimble_data)
stopCluster(cl)


nimble_samples <- do.call(rbind, chain_output)

View(MCMCvis::MCMCsummary(chain_output))


MCMCvis::MCMCtrace(chain_output,ind=TRUE,params= c("sigma_site_mb1","sigma_site_mb2","amp_tau_mb1","amp_tau_mb2"))

MCMCvis::MCMCtrace(chain_output,ind=TRUE,params= c("epsilon_site_mb1"))
