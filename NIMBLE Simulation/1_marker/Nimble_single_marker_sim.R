# Script to 
# 1. Simulate single marker replicated metabarcoding data
# 2. Call the single marker discrete latent parameter Nimble estimation model 
# 3. Summarize results
library(nimble)
library(tidyverse)
library(data.table)
library(here)
source(here("Zenodo","shared_functions","process_mb.R"))
source(here("Zenodo","shared_functions","sim_functions.R"))

ind_sim <- 1  # How many datasets to simulate?
SIMSET=c(1:ind_sim) # Vector for looping
TOT.CONC <- c(10, 30, 100, 300, 1000, 3000, 10000) # Total true concentrations of each "site" (biological sample) in each dataset
n_site = length(TOT.CONC)
N_sp <- 30 # How many species
N_Rep <- 3 # How many technical replicates per "site"

# Pull a vector of seeds to be used in the simulation  
SEED <- sample(1:1e5,ind_sim,replace=FALSE)

THEME = "Nimble_"
NOM <- paste0(THEME,ind_sim,"_Single_Marker_",n_site,"site_",N_sp,"sp_conc.Rdata")

sim_scen <- data.frame(N_SP = rep(N_sp,n_site),
                       TOT.CONC = TOT.CONC,
                       CONC.DIST = rep("very_high_skew",n_site),
                       N_READS = rep(100000,n_site),
                       amp_sd = rep(0.01),
                       replicates = rep(N_Rep,n_site))

post_tot_conc_all <- NULL
conc_dat_all <- NULL
nimble_data_list = NULL
nimbleMod_summary_list = NULL

for(j in SIMSET){# loop over number of independent sims
  set.seed(SEED[j])
  print(paste(j,"of",ind_sim))
  OUT <- list()
  library(tidyverse)
  
 
  #######################################################################
  ######################## Simulation 1
  ###### N_sp sim... Reads fixed. 
  ######################################################################
  #####################################################################
  A1  <- sim_amp(amp_mean = 0.9,
                 amp_sd = sim_scen$amp_sd[1] ,
                 N_species = sim_scen$N_SP[1])
  
  for( i in 1:nrow(sim_scen)){
    
    N_SP <- sim_scen$N_SP[i]
    TOT.CONC <- sim_scen$TOT.CONC[i]
    CONC.DIST <- sim_scen$CONC.DIST[i]
    N_READS <- sim_scen$N_READS[i]
    amp_sd <- sim_scen$amp_sd[i]
    replicates <- sim_scen$replicates[i]
    
    conc_dat <- sim_conc(N_species=N_SP,
                         conc_dist = "very_high_skew",
                         tot_conc = TOT.CONC)
    
    # Ok Assign a amplification rate to each species for 2 marker
    
    conc_dat$MB1_alpha <- A1$a_val
    
    # OK.  Simulate 3 replicates for each marker
    MB1 <-sim_MB(out = conc_dat,
                 a_val=conc_dat$MB1_alpha,
                 replicates = replicates,
                 N_tot_reads_min = N_READS/2,
                 N_tot_reads_max = N_READS,
                 N_pcr = 40)
    
    ##############################################################
    # Write to file:
    ##############################################################
    out <- list(conc_dat=conc_dat,MB1=MB1)
    
    OUT[[i]] <- process_mb1(out,site_name=paste0("X",i))
  }
  
  # These are the true simulation values.
  CONC <- list()
  
  for(i in 1:n_site){
    nommy <- as.name(paste0("X",i))
    CONC[[nommy]] <- OUT[[i]]$conc_dat
  }               
  
  
  # These are the observations
  N_pcr <- out$MB1$N_pcr
  
  MB1 <- OUT[[1]]$MB1
  for(i in 2:n_site){
    MB1 <- rbind(MB1,OUT[[i]]$MB1)
  }  
  
  MB1 <- MB1 %>% rename(site_id =site_idx)
  
  # Make long-form 
  MB1_long <- data.table::melt(data.table(MB1),id.vars =c("sp_id","samp_idx", "site_id"),
                               variable.name = "rep", value.name = "count")
  
  MB1_long <- MB1_long %>% as.data.frame() %>% 
    mutate(sp_id = as.numeric(sp_id),
           samp_idx = as.numeric(samp_idx),
           count = as.numeric(count))
  
  MB1_long <- MB1_long %>% mutate(marker = "MB1")
  
  MB_all_long <- MB1_long %>% as.data.frame()
  
  ##############################################
  ###### UNKNOWN SAMPLES
  ###### Make Indices for species, sample, site
  ##############################################
  # Find the unique sites
  site <- MB1_long %>% 
    distinct(site_id) %>% mutate(site_idx = 1:nrow(.))
  N_site <- nrow(site)
  
  # Marker
  marker <- MB_all_long %>% distinct(marker) %>% mutate(marker_idx = 1:nrow(.))
  N_marker <- nrow(marker)
  
  # Find the unique sites in any marker
  site_marker <- MB1_long %>% distinct(marker,site_id)
  site_marker <- site_marker %>% mutate(site_marker_idx = 1:nrow(.)) %>% 
    left_join(.,site) %>% left_join(.,marker)
  
  species <- MB_all_long %>% distinct(sp_id) %>% mutate(sp_idx = 1: nrow(.))
  N_sp <- nrow(species)
  
  #############################################
  # add relevant indices to data.
  #############################################
  MB_all_long <- MB_all_long %>% left_join(.,species) %>% left_join(.,site_marker) 
  
  site_samp <- MB_all_long %>% distinct(site_idx, samp_idx)
  site_samp <- site_samp %>% mutate(site_samp_idx = 1:nrow(.)) %>% as.data.frame()
  N_site_samp <- nrow(site_samp)
  
  MB_all_long <- MB_all_long %>% left_join(.,site_samp) 
  
  #############################################
  # Make MARKER-SPECIFIC OBSERVATION MATRICES 
  #############################################
  # Make wide form matrices for each marker (columns= species, rows =replicates + site-sample pairs)
  MB1_dat <-pivot_wider(MB_all_long %>% filter(marker_idx==1),id_cols=c("site_samp_idx","rep"),
                        values_from = count,
                        names_from = sp_idx)
  N_obs_mb1 <- nrow(MB1_dat)
  
  #############################################
  # make design matrices for mapping site_level concentrations 
  # to site_sample concentrations
  #############################################
  form <- "samp_idx ~ 0 + factor(site_idx)"
  model_frame   <- model.frame(form, site_samp)  
  X_D_to_F <- model.matrix(as.formula(form), model_frame)
  
  #############################################
  # make design matrices for mapping site_sample pairs to observations
  #############################################
  form <- "rep ~ 0 + factor(site_samp_idx)"
  model_frame   <- model.frame(form, MB1_dat)  
  X_MB1_ss <- model.matrix(as.formula(form), model_frame)
  
  # Summarise which species were observed in each sample.
  TMP <- MB_all_long %>% group_by(marker_idx,site_samp_idx,sp_idx) %>% summarise(Sum = sum(count)) %>% 
    mutate(IND = ifelse(Sum>0,1,0))
  MB1_ind <- pivot_wider(TMP %>% filter(marker_idx==1),id_cols=c("site_samp_idx"),
                         values_from = IND,
                         names_from = sp_idx)
  
  #### PULL IN AMPLIFICATION EFFICIENCIES AND PRETEND THEY ARE KNOWN FOR NOW.
  # Define the reference species for each sample based on the most commonly observed species in each sample and marker
  ref_dat <-  MB_all_long %>% group_by(marker_idx,site_samp_idx,sp_idx) %>% summarise(Sum = sum(count)) %>% 
    ungroup() %>% group_by(marker_idx,site_samp_idx) %>% mutate(Max=max(Sum)) %>% filter(Sum==Max)
  
  ref_dat_all <- MB_all_long %>% group_by(marker_idx,sp_idx) %>% summarise(Sum = sum(count)) %>% 
    ungroup() %>% group_by(marker_idx) %>% mutate(Max=max(Sum)) %>% filter(Sum==Max)
  
  ref_sp_mb1 <- ref_dat %>% filter(marker_idx==1) %>% left_join(MB1_dat,.) %>% 
    dplyr::select(marker_idx,site_samp_idx,sp_idx)
  
  # Choose species 33 to be alpha == 0 for marker 1 and 27 to be alpha==0 for marker 2
  # These are treated as known and so the choice is somewhat arbitrary.
  amp_alpha_mb1 <- log(out$conc_dat$MB1_alpha + 1) - log(out$conc_dat$MB1_alpha[ref_dat_all$sp_idx[1]] + 1)
  
  # Here are the vector of alphas for use in the analysis.
  amp_alpha_mb1
  
  # Summarize observe reads for each replicate.
  Reads <- MB_all_long %>% group_by(site_idx,marker_idx,samp_idx,rep) %>% summarise(N_read = sum(count))
  MB_all_long <- MB_all_long %>% left_join(.,Reads)
  
  reads_mb1 <- MB1_dat %>% dplyr::select(-site_samp_idx,-rep) %>% rowSums()
  
  # Calculate combinatorial coefficient.
  log_mb1_n_k <- MB1_dat %>% dplyr::select(-site_samp_idx,-rep) %>% as.matrix()
  
  for(i in 1:nrow(MB1_dat)){
    tmp <-(MB1_dat %>% dplyr::select(-site_samp_idx,-rep))[i,] %>% unlist() %>% c()
    log_mb1_n_k[i,] <- as.numeric(lchoose(reads_mb1[i],tmp))
  }
  
  ############################################################
  #### Get files ready for nimble
  ############################################################
  
  # Format for NIMBLE -------------------------------------------------------
  
  nimble_data <- list(
    # Counters
    # Model parameters
    N_pcr = N_pcr,
    # Observations of read counts  
    count_mb1 = as.matrix(MB1_dat %>% dplyr::select(-site_samp_idx,-rep)),
    reads_mb1 = reads_mb1,
    # Combinatorial coefficients
    log_mb1_n_k = log_mb1_n_k,
    # Amplification efficiencies for each marker
    amp_alpha_mb1 = amp_alpha_mb1,
    # Indicators for species detection by marker
    MB1_ind = as.matrix(MB1_ind %>% ungroup() %>% dplyr::select(-site_samp_idx)),
    MB1_all =  as.matrix(MB1_ind %>% ungroup() %>% dplyr::select(-site_samp_idx)) %>% apply(.,2, function(col) as.numeric(any(col > 0))) ,
    # Design matrices for mapping samples to their replicates
    X_MB1_ss = matrix(X_MB1_ss,nrow=N_obs_mb1,ncol=N_site_samp),
    X_D_to_F = matrix(X_D_to_F,nrow=N_site_samp,ncol=N_site)
  )
  
  # Constants for NIMBLE (design matrices and fixed parameters)
  nimble_constants <- list(
    # Index vectors to avoid dynamic indexing
    # N_marker = as.numeric(N_marker),
    N_site = as.numeric(N_site),
    N_sp = as.numeric(N_sp),
    N_site_samp = as.numeric(N_site_samp),
    N_obs_mb1 = as.numeric(N_obs_mb1),
    # Mapping indices
    # site_map = site_samp$site_idx,
    samp_map1 = as.numeric(ref_sp_mb1$site_samp_idx),
    ref_sp_mb1 = as.numeric(ref_sp_mb1$sp_idx)
  )
  # Packages
  library(nimble)
  library(MCMCvis)
  library(parallel)
  cl<-makeCluster(2,outfile="")
  
  # initialization function
  nimble_inits <- function() {
    N_site <- nimble_constants$N_site
    N_sp <- nimble_constants$N_sp
    N_obs_mb1 <- nimble_constants$N_obs_mb1
    N_site_mock <- nimble_constants$N_site_mock
    # Initialize d_mb based on count_mb 
    d_mb1_init <- as.matrix((nimble_data$count_mb1 > 0) * 1)
    
    list(
      log_D = matrix(rnorm(N_site * N_sp, mean = 0, sd = 1), 
                     nrow = N_site, ncol = N_sp),
      beta0 = runif(1,-.2, .20),
      beta1 = runif(1,-.3, 0),
      log_phi1 = rnorm(1,1.5,.1),
      phi0 = abs(rnorm(1.5,.25)),
      #epsilon_site_mb1_z = rnorm(N_site*N_sp, 0, 1),
      epsilon_species_mb1_z = matrix(rnorm(N_site * N_sp, 0, 1), nrow = N_site, ncol = N_sp),
      #sigma_site_mb1 = runif(1, 0.1, 0.5),
       sigma_species_mb1  = runif(1, 0.1, 0.2),
      d_mb1 = (as.matrix(MB1_dat %>% dplyr::select(-site_samp_idx,-rep)) > 0) * 1
    )
  }
  
  chain_info <- list(
    list(inits = nimble_inits()),
    list(inits = nimble_inits()))
  
  run_Nimble_MCMC <- function(chain_info, nimble_data){
    #  zero-inflated beta-binomial with discrete states and sample-level zeroes
    library(nimble)
    cleanup_nimble <- function() {
      if (exists("cmodel")) {
        try(nimble:::clearCompiled(cmodel), silent = TRUE)
      }
      if (exists("cmcmc")) {
        try(nimble:::clearCompiled(cmcmc), silent = TRUE)
      }
      gc()
    }
    
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
    clr_alpha <- nimbleFunction(
      run = function(amp_alpha = double(1),     # amplification efficiencies for all species (on log scale)
                     mb_ind = double(1)) {     # sample-level indicator for this sample
        returnType(double(1))
        
        n_sp <- length(amp_alpha)
        clr_result <- numeric(n_sp)
        
        present_alphas <- numeric(0)
        
        for(j in 1:n_sp) {
          if(mb_ind[j] == 1) {  # detectable AND present
            present_alphas <- c(present_alphas, amp_alpha[j])  # amp_alpha is already log scale
          }else{present_alphas <- present_alphas}
        }
        # Calculate mean of log amplification efficiencies (geometric mean in original scale)
        mean_log_alpha <- mean(present_alphas)
        
        # Calculate CLR for all species, but only meaningful for present ones
        for(j in 1:n_sp) {
          clr_result[j] <- amp_alpha[j] - mean_log_alpha  # Simple difference since already on log scale
        }
        return(clr_result)
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
    
    assign('dbetabin', dbetabin, envir = .GlobalEnv)
    assign('rbetabin', rbetabin, envir = .GlobalEnv)
    assign('clr_alpha', clr_alpha, envir = .GlobalEnv)
    
    inits <- chain_info$inits
    
    # Enhanced NIMBLE model with proper two-level zero-inflation
    modelCode <- nimbleCode({
      
      # Parameters - priors
      for(i in 1:N_site) {
        for(j in 1:N_sp) {
          log_D[i,j] ~ dnorm(0, sd = 10)
        }
      }
      
      # Phi parameters
      beta0 ~ dnorm(0, sd = 1)
      beta1 ~ T(dnorm(-10, sd = 5),-Inf,0)
      phi0 ~ dnorm(2, sd = 1)
      log_phi1 ~ dnorm(1, sd = .75)
      phi1 <- exp(log_phi1)
      
      
      # Random Effects Sigmas
      #sigma_site_mb1 ~ T(dnorm(0, sd = 1), 0, Inf)
      sigma_species_mb1 ~ T(dnorm(0, sd = 1), 0, Inf)
      # for(m in 1:N_site) {
      #   sigma_species_mb1[m] ~ T(dnorm(0, sd = 1), 0, Inf)
      # }
      # Random Effects
      for(m in 1:N_site) {
        # epsilon_site_mb1_z[m] ~ dnorm(0, 1)
        # epsilon_site_mb1[m] <- epsilon_site_mb1_z[m] * sigma_site_mb1
        for(j in 1:N_sp) {
          epsilon_species_mb1_z[m,j] ~ dnorm(0, 1)
          epsilon_species_mb1[m,j] <- epsilon_species_mb1_z[m,j] * sigma_species_mb1
        }
      }
      
      # Matrix operations 
      log_F[1:N_site_samp, 1:N_sp] <- X_D_to_F[1:N_site_samp, 1:N_site] %*% log_D[1:N_site, 1:N_sp]
      log_F_mb1[1:N_obs_mb1, 1:N_sp] <- X_MB1_ss[1:N_obs_mb1, 1:N_site_samp] %*% log_F[1:N_site_samp, 1:N_sp]
      
      
      # MB1 calculations 
      for(i in 1:N_obs_mb1) {
        # Calculate ALL log_lambda_K values first
        for(j in 1:N_sp) {
          # Cache all exponentials and logs at once
          temp_exp_mb1[i,j] <- exp(log_F_mb1[i,j])
          log_1mexp_mb1[i,j] <- log(1 - exp(-temp_exp_mb1[i,j]))
          log_lambda_K_mb1[i,j] <- log_F_mb1[i,j] - log_1mexp_mb1[i,j]
          lambda_mb1[i,j] <- MB1_all[j] * temp_exp_mb1[i,j] # log copies
          # Cache sample-level probabilities
          log_p_pos_mb1[i,j] <- log_1mexp_mb1[i,j]
          log_p_zero_mb1[i,j] <- -temp_exp_mb1[i,j]
        }
        
        # extract reference values
        log_lambda_K_mb1_ref[i] <- log_lambda_K_mb1[i,ref_sp_mb1[i]]
        amp_alpha_mb1_ref[i] <- amp_alpha_mb1[ref_sp_mb1[i]]
        
        for(j in 1:N_sp) {
          # Additive log ratios
          log_val_mb1[i,j] <- (log_lambda_K_mb1[i,j] - log_lambda_K_mb1_ref[i]) + 
            N_pcr * (amp_alpha_mb1[j] - amp_alpha_mb1_ref[i])
          # Aliquot-level discrete DNA-presence state 
          d_mb1[i,j] ~ dbern(exp(log_1mexp_mb1[i,j]))
        }
        
        # Apply discrete state masking to log values
        log_Lambda_mb1[i] <- log(sum(lambda_mb1[i,1:N_sp])) 
        log_max_val_m1[i] <- max(log_val_mb1[i,1:N_sp])
        clr_alpha_mb1[i,1:N_sp] <- clr_alpha(amp_alpha_mb1[1:N_sp], 
                                             MB1_all[1:N_sp])
        for(j in 1:N_sp) {
          # Discrete state masking
          exp_shifted_m1[i,j] <- exp(log_val_mb1[i,j] - log_max_val_m1[i]) * d_mb1[i,j]
          # Phi calculations
          log_phi_mb1[i,j] <- beta0 +
            beta1 * clr_alpha_mb1[i,j] +
            log_Lambda_mb1[i] +
            phi0 * exp(-phi1 * log_lambda_K_mb1[i, j]) +
            #epsilon_site_mb1[samp_map1[i]] #+
           epsilon_species_mb1[samp_map1[i], j]
          phi_mb1[i,j] <- exp(min(20,log_phi_mb1[i,j]))
        }
        
        sum_exp_shifted_m1[i] <- sum(exp_shifted_m1[i,1:N_sp])
        
        for(j in 1:N_sp) {
          # Softmax
          pi_samp_mb1_cond[i,j] <- exp_shifted_m1[i,j] / sum_exp_shifted_m1[i]
          
          # Beta-binomial parameters
          alpha_mb1[i,j] <- pi_samp_mb1_cond[i,j] * phi_mb1[i,j]
          beta_mb1[i,j] <- phi_mb1[i,j] - alpha_mb1[i,j]
          
          # Likelihood 
          count_mb1[i,j] ~ dbetabin(
            reads_mb1[i], alpha_mb1[i,j], beta_mb1[i,j],
            log_p_zero_mb1[i,j], log_p_pos_mb1[i,j],
            log_mb1_n_k[i,j],d_mb1[i,j]
          )
        }
      }
    })
    
    
    # Model building
    model <- nimbleModel(code = modelCode,
                         data = nimble_data,
                         constants = nimble_constants,
                         inits = inits,
                         check = TRUE)
    
    
    cleanup_nimble()
    
    cmodel <- compileNimble(model)
    
    # Configure MCMC
    mcmc_conf <- configureMCMC(model, useConjugacy=FALSE)
    # Remove default samplers for log_D
    mcmc_conf$removeSamplers("log_D")
    
    # Add block samplers for log_D by site
    for (i in 1:nimble_constants$N_site){
      targets <- paste0("log_D[", i, ", 1:", nimble_constants$N_sp, "]")
      mcmc_conf$addSampler(target=targets, "AF_slice")
    }
    fixed_idx1 <- which(nimble_data$count_mb1 > 0, arr.ind = TRUE)
    for(k in 1:nrow(fixed_idx1)) {
      i <- fixed_idx1[k,1]
      j <- fixed_idx1[k,2]
      mcmc_conf$removeSamplers(paste0("d_mb1[",i,",",j,"]"))
    }
    # Remove default samplers for the hierarchical epsilon structure
    #mcmc_conf$removeSamplers("epsilon_site_mb1_z")
     mcmc_conf$removeSamplers("epsilon_species_mb1_z")
    
    # Remove default samplers for variance parameters
    # mcmc_conf$removeSamplers("sigma_species_mb1")
    
    
    # Add site-level epsilon samplers (these are univariate)
    #mcmc_conf$addSampler(target=c("sigma_site_mb1",paste0("epsilon_site_mb1_z[1:", nimble_constants$N_site, "]")), "AF_slice")
    
    
    # Add block samplers for species-level epsilon by site, including site-specific variance
    for (i in 1:nimble_constants$N_site) {
    # Block sample species effects and their variance together for each site
    epsilon_targets_mb1 <- paste0("epsilon_species_mb1_z[", i, ", 1:", nimble_constants$N_sp, "]")
    #sigma_target_mb1 <- paste0("sigma_species_mb1[", i, "]")
    mcmc_conf$addSampler(target=c(epsilon_targets_mb1), "AF_slice")
    }
    
    # Add block samplers for phi parameters
    run_targets_mb1 <- c("beta0","beta1",
                         "phi0",
                         "log_phi1")
    mcmc_conf$removeSamplers(run_targets_mb1)
    mcmc_conf$addSampler(target=run_targets_mb1, "AF_slice")
    
    mcmc_conf$addMonitors("phi_mb1[]")
    # Build and compile MCMC
    mcmc <- buildMCMC(mcmc_conf)
    cmcmc <- compileNimble(mcmc, project = model)
    # Run MCMC
    #samples <- runMCMC(cmcmc, niter = 3000, nburnin = 1500, thin = 1)
    samples <- runMCMC(cmcmc, niter = 2, nburnin = 1, thin = 1)
    
    return(samples)
  }
  
  ### THIS RUNS THE MODEL AND SPITS OUT MCMC CHAINS.
  clusterExport(cl, c("nimble_data","nimble_constants","chain_info"))
  
  chain_output <- parLapply(cl = cl, X = chain_info,
                            fun = run_Nimble_MCMC, 
                            nimble_data = nimble_data);stopCluster(cl)
  nimbleMod_summary = MCMCvis::MCMCsummary(chain_output)
  
  nimble_samples <- do.call(rbind, chain_output)
  
  ## OK.  Make some simple comparisons with the concentration values.
  ### GO GET THE SIMULATION DATA.
  conc_dat <- NULL
  for(i in 1:n_site){
    tmp <- data.frame(sp_id = CONC[[i]]$sp_id,
                      conc_val = CONC[[i]]$conc_val,
                      site_id = CONC[[i]]$site_idx,
                      conc_prop = CONC[[i]]$conc_prop,
                      conc_total = CONC[[i]]$tot_conc)
    conc_dat <- rbind(conc_dat,tmp)
  }
  # Merge in correct indexes
  conc_dat <- left_join(conc_dat,site)
  
  # Summarize log_D
  
  logD_columns <- colnames(nimble_samples[,grep("log_D",colnames(nimble_samples))])
  
  sp_Conc_summary <- function(data, funct,... ) {
    # Initialize result matrix
    n_samples <- nrow(data)
    site_sums <- matrix(0, nrow = n_site, ncol = N_sp)
    
    # For each site, sum across all species
    for(site in 1:n_site) {
      site_columns <- paste0("log_D[", site, ", ", 1:N_sp, "]")
      
      # Check which columns exist (in case some are missing)
      existing_cols <- site_columns[site_columns %in% colnames(data)]
      
      if(length(existing_cols) > 0) {
        site_sums[site, ] <- apply(exp(data[, existing_cols, drop = FALSE]),2,FUN=funct,...)
      }
    }
    
    # Add column names
    # colnames(site_sums) <- paste0("site_", 1:n_site, "_sum")
    
    return(site_sums)
  } 
  
  A <- sp_Conc_summary(nimble_samples[,logD_columns],mean) %>% as_tibble() %>%
    mutate(site_idx = 1:nrow(.)) %>%
    pivot_longer(.,-site_idx,names_to="sp_id",values_to ="Mean")
  B <- sp_Conc_summary(nimble_samples[,logD_columns],sd) %>% as_tibble() %>%
    mutate(site_idx = 1:nrow(.)) %>%
    pivot_longer(.,-site_idx,names_to="sp_id",values_to="SD")
  C <- sp_Conc_summary(nimble_samples[,logD_columns],median) %>% as_tibble() %>%
    mutate(site_idx = 1:nrow(.)) %>%
    pivot_longer(.,-site_idx,names_to="sp_id",values_to="Median")
  D <- sp_Conc_summary(nimble_samples[,logD_columns],quantile, probs = .05) %>% as_tibble() %>%
    mutate(site_idx = 1:nrow(.)) %>%
    pivot_longer(.,-site_idx,names_to="sp_id",values_to="q_05")
  E <- sp_Conc_summary(nimble_samples[,logD_columns],quantile, probs = .95)  %>% as_tibble() %>%
    mutate(site_idx = 1:nrow(.)) %>%
    pivot_longer(.,-site_idx,names_to="sp_id",values_to="q_95")
  
  post <- left_join(A,B) %>% left_join(.,C) %>% left_join(.,D) %>% left_join(.,E)
  post <- post %>% mutate(sp_idx = sp_id, sp_idx= gsub('V','',sp_idx))
  
  conc_dat_est <- left_join(conc_dat %>% mutate(sp_id=as.character(sp_id)),
                            post %>% dplyr::select(-sp_id),by = join_by(sp_id == sp_idx,site_idx) )
  
  
  site_Conc_sum <- function(data, n_site, N_sp) {
    # Initialize result matrix
    n_samples <- nrow(data)
    site_sums <- matrix(0, nrow = n_samples, ncol = n_site)
    
    # For each site, sum across all species
    for(site in 1:n_site) {
      site_columns <- paste0("log_D[", site, ", ", 1:N_sp, "]")
      
      # Check which columns exist (in case some are missing)
      existing_cols <- site_columns[site_columns %in% colnames(data)]
      
      if(length(existing_cols) > 0) {
        site_sums[, site] <- rowSums(exp(data[, existing_cols, drop = FALSE]))
      }
    }
    
    # Add column names
    colnames(site_sums) <- paste0("site_", 1:n_site, "_sum")
    
    return(site_sums)
  }
  
  
  # Calculate the total Concentration
  post_tot_conc <- site_Conc_sum(nimble_samples,n_site,N_sp) %>% as_tibble()
  colnames(post_tot_conc) <- c(paste0("X",1:n_site))
  post_tot_conc_summ <- post_tot_conc %>% pivot_longer(.,everything(),names_to = "site_id",values_to="tot_conc") %>% 
    group_by(site_id) %>% 
    summarise(Mean=mean(tot_conc),
              Median=median(tot_conc),
              SD = sd(tot_conc),
              q_05=quantile(tot_conc,probs=0.05),
              q_95=quantile(tot_conc,probs=0.95))
  
  post_tot_conc_summ <- data.frame(tot_conc_true=conc_dat$conc_total,site_id=conc_dat$site_id) %>% 
    distinct(site_id,tot_conc_true) %>% left_join(.,post_tot_conc_summ)
  
  
  calculate_phi_mb1_post <- function(data, N_obs_mb1, N_sp) {
    # Get phi_mb1 columns
    phi_mb1_cols <- grep("phi_mb1", colnames(data))
    phi_mb1_data <- data[, phi_mb1_cols]
    
    # Calculate column means
    phi_mb1_means <- colMeans(phi_mb1_data)
    
    # Create a mapping from column names to matrix positions
    col_names <- names(phi_mb1_means)
    
    # Extract observation and species indices from column names
    obs_indices <- as.numeric(gsub("phi_mb1\\[(\\d+), \\d+\\]", "\\1", col_names))
    sp_indices <- as.numeric(gsub("phi_mb1\\[\\d+, (\\d+)\\]", "\\1", col_names))
    
    # Initialize result matrix
    phi_mb1_matrix <- matrix(NA, nrow = N_obs_mb1, ncol = N_sp)
    
    # Fill matrix using vectorized indexing
    for(i in seq_along(phi_mb1_means)) {
      obs_idx <- obs_indices[i]
      sp_idx <- sp_indices[i]
      phi_mb1_matrix[obs_idx, sp_idx] <- phi_mb1_means[i]
    }
    
    # Add row and column names
    rownames(phi_mb1_matrix) <- paste0("obs_", 1:N_obs_mb1)
    colnames(phi_mb1_matrix) <- paste0("V", 1:N_sp)
    
    return(phi_mb1_matrix)
  }
  
  ######### Summarise the phi values
  A <- calculate_phi_mb1_post(nimble_samples, N_obs_mb1, N_sp) %>% as_tibble() %>% 
    mutate(site_idx = 1:nrow(.)) %>% 
    pivot_longer(.,-site_idx,names_to="sp_id",values_to ="Mean_phi")
  
  phi_mb1 <- A %>% mutate(sp_idx = sp_id, sp_idx= gsub('V','',sp_idx)) 
  
  conc_dat_est <- left_join(conc_dat_est,phi_mb1 %>% dplyr::select(-sp_id),
                            by=join_by(site_idx==site_idx,sp_id==sp_idx))
  
  
  BRKS = c(0.001,0.01,0.1,1,10,30,100,300,1000,3000)
  
  ggplot(conc_dat_est) +
    geom_point(aes(x=Mean,y=Mean_phi,color=as_factor(site_idx))) +
    scale_y_continuous(trans="log",breaks=BRKS) +
    scale_x_continuous(trans="log",breaks=BRKS) 

  ## Real simple plots
  
  p_true_est_conc <- ggplot(conc_dat_est) +
    geom_point(aes(x=conc_val,y=Mean,color=site_id))+
    geom_errorbar(aes(x=conc_val,ymin=q_05,ymax=q_95,color=site_id),width=0,alpha=0.8)+
    scale_x_continuous(trans="log","True conc.",breaks=BRKS) +
    scale_y_continuous(trans="log","Est conc.",breaks=BRKS) +
    coord_cartesian(xlim=c(min(BRKS),max(BRKS)),ylim=c(min(BRKS),max(BRKS))) +
    geom_abline(intercept=0,slope=1,linetype="dashed", color="red")+
    theme_bw()
  
  BRKS = c(10,30,100,300,1000,3000,10000)
  p_true_est_tot_conc <- ggplot(post_tot_conc_summ) +
    geom_point(aes(x=tot_conc_true,y=Mean,color=site_id))+
    geom_errorbar(aes(x=tot_conc_true,ymin=q_05,ymax=q_95,color=site_id),width=0,alpha=0.8)+
    scale_x_continuous(trans="log","True total conc.",breaks=BRKS) +
    scale_y_continuous(trans="log","Est total conc.",breaks=BRKS) +
    #coord_cartesian(xlim=c(min(BRKS),max(BRKS)),ylim=c(min(BRKS),max(BRKS))) +
    geom_abline(intercept=0,slope=1,linetype="dashed", color="red")+
    theme_bw()
  
  ### Save some output cross estimation loop.
  
  conc_dat_all <- rbind(conc_dat_all,conc_dat_est %>% mutate(rep_id = j))
  post_tot_conc_all <- rbind(post_tot_conc_all,post_tot_conc_summ %>% mutate(rep_id = j))
  
  
  ### Retain the input data for each replicate
  nimble_data_list[[j]] <- list(nimble_data=nimble_data,rep_id = j)
  nimbleMod_summary_list[[j]]  <- list(nimbleMod_summary,rep_id=j)}

# Summarize match between simulated and estimated.
post_tot_conc_all1 <- post_tot_conc_all %>% mutate(prop_error = (Median - tot_conc_true) / tot_conc_true)
conc_dat_all1 <- conc_dat_all %>% mutate(prop_error = (Median - conc_val) / conc_val)

post_tot_conc_all_summ <- post_tot_conc_all1 %>% group_by(site_id,tot_conc_true) %>% 
  summarise(grand_mean = mean(Mean),
            q_05= quantile(Mean,probs=0.05),
            q_25= quantile(Mean,probs=0.25),
            grand_median = median(Mean),
            q_75= quantile(Mean,probs=0.75),
            q_95= quantile(Mean,probs=0.95),
            grand_prop_error = mean(prop_error),
            q_05_error= quantile(prop_error,probs=0.05),
            q_25_error= quantile(prop_error,probs=0.25),
            grand_median_error = median(prop_error),
            q_75_error= quantile(prop_error,probs=0.75),
            q_95_error= quantile(prop_error,probs=0.95)) 

BRKS = c(0.001,0.01,0.1,1,10,30,100,300,1000,3000,10000)
BRKS2 = c(0.01,0.1,1,10,30,100,300,1000,3000,10000)
BRKS.prop = c(-2,-1,-0.5,-0.25,0,0.25,0.5,1,2)

# conc_dat_all <- conc_dat_all %>% filter(!rep_id %in% c(35,28,15))
# post_tot_conc_all <- post_tot_conc_all %>% filter(!rep_id %in% c(35,28,15))

p_true_est_conc_rep <- ggplot(conc_dat_all1) +
  geom_point(aes(x=conc_val,y=Median,color=site_id))+
  geom_errorbar(aes(x=conc_val,ymin=q_05,ymax=q_95,color=site_id),width=0,alpha=0.1)+
  scale_x_continuous(trans="log","True conc.",breaks=BRKS) +
  scale_y_continuous(trans="log","Est conc.",breaks=BRKS) +
  coord_cartesian(xlim=c(min(BRKS),max(BRKS)),ylim=c(min(BRKS),max(BRKS))) +
  geom_abline(intercept=0,slope=1,linetype="dashed", color="red")+
  theme_bw()

p_true_est_tot_conc_rep <- ggplot(post_tot_conc_all1) +
  geom_point(aes(x=tot_conc_true,y=Median,color=site_id),alpha=0.5)+
  geom_errorbar(aes(x=tot_conc_true,ymin=q_05,ymax=q_95,color=site_id),width=0,alpha=0.3)+
  scale_x_continuous(trans="log", "True total conc.",breaks=BRKS) +
  scale_y_continuous(trans="log","Est total conc.",breaks=BRKS) +
  #coord_cartesian(xlim=c(min(BRKS),max(BRKS)),ylim=c(min(BRKS),max(BRKS))) +
  geom_abline(intercept=0,slope=1,linetype="dashed", color="red")+
  theme_bw()

p_prop_error_conc_rep <- ggplot(conc_dat_all1) +
  geom_point(aes(x=conc_val,y=prop_error),alpha=0.5,shape=".")+
  #geom_jitter(aes(x=tot_conc_true,y=prop_error),height=0,width=0.5,alpha=0.5)+#geom_errorbar(aes(x=tot_conc_true,ymin=q_05,ymax=q_95,color=site_id),width=0,alpha=0.3)+
  scale_x_continuous(trans="log","True total conc.",breaks=BRKS) +
  scale_y_continuous("Proportional Error",breaks=BRKS.prop) +
  coord_cartesian(xlim=c(min(BRKS2),NA),ylim=c(-1,max(BRKS.prop))) +
  geom_hline(yintercept =0, linetype="dashed", color="red")+
  facet_wrap(~as.factor(conc_total)) +
  theme_bw()

p_prop_error_tot_conc_rep <- ggplot(post_tot_conc_all1) +
  geom_boxplot(aes(x=tot_conc_true,y=prop_error,group=tot_conc_true),alpha=0.5)+
  geom_jitter(aes(x=tot_conc_true,y=prop_error),height=0,width=0.25,alpha=0.5)+#geom_errorbar(aes(x=tot_conc_true,ymin=q_05,ymax=q_95,color=site_id),width=0,alpha=0.3)+
  scale_x_continuous(trans="log", "True total conc.",breaks=BRKS) +
  #scale_y_continuous(trans="log","Proportional Error",breaks=BRKS) +
  #coord_cartesian(xlim=c(min(BRKS),max(BRKS)),ylim=c(min(BRKS),max(BRKS))) +
  geom_hline(yintercept =0, linetype="dashed", color="red")+
  theme_bw()

plots = list( p_true_est_conc_rep=p_true_est_conc_rep,
              p_true_est_tot_conc_rep = p_true_est_tot_conc_rep,
              p_prop_error_conc_rep =p_prop_error_conc_rep,
              p_prop_error_tot_conc_rep =p_prop_error_tot_conc_rep)

OUT <- list(nimbleMod_summary_list = nimbleMod_summary_list,
            conc_dat_all = conc_dat_all , 
            post_tot_conc_all= post_tot_conc_all)

#save(OUT,file="FILE_NAME.RData")
