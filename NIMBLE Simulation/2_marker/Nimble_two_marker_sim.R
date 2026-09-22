# Script to 
# 1. Simulate single marker replicated metabarcoding data
# 2. Call the single marker discrete latent parameter Nimble estimation model 
# 3. Summarize results
library(nimble)
library(tidyverse)
library(data.table)
library(parallel)
library(MCMCvis)
library(here)

### Functions for simulations of metabarcoding
source(here("Zenodo","shared_functions","process_mb.R"))
source(here("Zenodo","shared_functions","sim_functions.R"))

set.seed(101)
ind_sim <- 1  # How many datasets to simulate?
SIMSET=c(1:ind_sim) # Vector for looping
TOT.CONC <- c(10, 30, 100, 300, 1000, 3000, 10000) # Total true concentrations of each "site" (biological sample) in each dataset
n_site = length(TOT.CONC)
N_sp <- 50 # How many species
N_Rep <- 3 # How many technical replicates per "sample"

# Pull a vector of seeds to be used in the simulation  
SEED <- sample(1:1e5,ind_sim,replace=FALSE)

# Get ready...
sim_scen <- data.frame(N_SP = rep(N_sp,n_site),
                       TOT.CONC = TOT.CONC,
                       CONC.DIST = rep("very_high_skew",n_site),
                       N_READS = rep(100000,n_site),
                       amp_sd = rep(0.01,n_site),
                       replicates = rep(N_Rep,n_site))

post_tot_conc_all <- NULL
conc_dat_all <- NULL
nimble_data_list = NULL
nimbleMod_summary_list = NULL

# Go! Within this block, 
for(j in SIMSET){# loop over number of independent sim
  set.seed(SEED[j])
  print(paste(j,"of",ind_sim))
  OUT <- list()
 
  # Simulate a single set of concentrations 
  # Ok Assign a amplification rate to each species for 2 marker
  A1  <- sim_amp(amp_mean = 0.9,
                 amp_sd = sim_scen$amp_sd[1] ,
                 N_species = sim_scen$N_SP[1])
  A2  <- sim_amp(amp_mean = 0.85,
                 amp_sd = sim_scen$amp_sd[1],
                 N_species = sim_scen$N_SP[1])
  # Modify A1 so 1:15 are 0
  A1$a_val[1:15] <- 0
  # Modify A2 so 36:50 are 0
  A2$a_val[36:sim_scen$N_SP[1]] <- 0
  #######################################################################
  ######################## Simulation 1
  ###### 50 sp. Reads fixed. 
  ######################################################################
  #####################################################################
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
    
    
    
    #This means that species 16 to 35  are shared between the two species.
    conc_dat$MB1_alpha <- A1$a_val
    conc_dat$MB2_alpha <- A2$a_val
    
    # OK.  Simulate 3 replicates for each marker
    MB1 <-sim_MB(out = conc_dat,
                 a_val=conc_dat$MB1_alpha,
                 replicates = replicates,
                 N_tot_reads_min = N_READS/2,
                 N_tot_reads_max = N_READS,
                 N_pcr = 40)
    
    MB2 <-sim_MB(out = conc_dat,
                 a_val=conc_dat$MB2_alpha,
                 replicates = replicates,
                 N_tot_reads_min = N_READS/2,
                 N_tot_reads_max = N_READS,
                 N_pcr = 40)
    
    ##############################################################
    # Write to file:
    ##############################################################
    out <- list(conc_dat=conc_dat,MB1=MB1,MB2=MB2)
    OUT[[i]] <- process_mb2(out,site_name=paste0("X",i))
    
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
  MB2 <- OUT[[1]]$MB2
  for(i in 2:n_site){
    MB1 <- rbind(MB1,OUT[[i]]$MB1)
    MB2 <- rbind(MB2,OUT[[i]]$MB2)
  }  
  
  MB1 <- MB1 %>% rename(site_id =site_idx)
  MB2 <- MB2 %>% rename(site_id =site_idx)
  
  # Make long-form 
  MB1_long <- data.table::melt(data.table(MB1),id.vars =c("sp_id","samp_idx", "site_id"),
                               variable.name = "rep", value.name = "count")
  MB2_long <- data.table::melt(data.table(MB2),id.vars =c("sp_id","samp_idx", "site_id"),
                               variable.name = "rep", value.name = "count")
  
  MB1_long <- MB1_long %>% as.data.frame() %>% 
    mutate(sp_id = as.numeric(sp_id),
           samp_idx = as.numeric(samp_idx),
           count = as.numeric(count))
  MB2_long <- MB2_long %>% as.data.frame() %>% 
    mutate(sp_id = as.numeric(sp_id),
           samp_idx = as.numeric(samp_idx),
           count = as.numeric(count))
  
  MB1_long <- MB1_long %>% mutate(marker = "MB1")
  MB2_long <- MB2_long %>% mutate(marker = "MB2")
  
  MB_all_long <- rbind(MB1_long, MB2_long) %>% as.data.frame()
  
  ##############################################
  ###### UNKNOWN SAMPLES
  ###### Make Indices for species, sample, site
  ##############################################
  # Find the unique sites
  site <- rbind(MB1_long , MB2_long) %>% 
    distinct(site_id) %>% mutate(site_idx = 1:nrow(.))
  N_site <- nrow(site)
  
  # Marker
  marker <- MB_all_long %>% distinct(marker) %>% mutate(marker_idx = 1:nrow(.))
  N_marker <- nrow(marker)
  
  # Find the unique sites in any marker
  site_marker <- rbind(MB1_long %>% distinct(marker,site_id),
                       MB2_long %>%  distinct(marker,site_id))
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
  
  MB2_dat <-pivot_wider(MB_all_long %>% filter(marker_idx==2),id_cols=c("site_samp_idx","rep"),
                        values_from = count,
                        names_from = sp_idx)
  N_obs_mb2 <- nrow(MB2_dat)
  
  
  #
  
  
  
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
  
  # make design matrices
  form <- "rep ~ 0 + factor(site_samp_idx)"
  model_frame   <- model.frame(form, MB2_dat)  
  X_MB2_ss <- model.matrix(as.formula(form), model_frame)
  
  # Summarise which species were observed in each sample.
  TMP <- MB_all_long %>% group_by(marker_idx,site_samp_idx,sp_idx) %>% summarise(Sum = sum(count)) %>% 
    mutate(IND = ifelse(Sum>0,1,0))
  MB1_ind <- pivot_wider(TMP %>% filter(marker_idx==1),id_cols=c("site_samp_idx"),
                         values_from = IND,
                         names_from = sp_idx)
  MB2_ind <- pivot_wider(TMP %>% filter(marker_idx==2),id_cols=c("site_samp_idx"),
                         values_from = IND,
                         names_from = sp_idx)
  
  #### PULL IN AMPLIFICATION EFFICIENCIES AND PRETEND THEY ARE KNOWN FOR NOW.
  # Define the reference species for each sample based on the most commonly observed species in each sample and marker
  ref_dat <-  MB_all_long %>% group_by(marker_idx,site_samp_idx,sp_idx) %>% summarise(Sum = sum(count)) %>% 
    ungroup() %>% group_by(marker_idx,site_samp_idx) %>% mutate(Max=max(Sum)) %>% filter(Sum==Max)
  
  ref_sp_mb1 <- ref_dat %>% filter(marker_idx==1) %>% left_join(MB1_dat,.) %>% 
    dplyr::select(marker_idx,site_samp_idx,sp_idx)
  ref_sp_mb2 <- ref_dat %>% filter(marker_idx==2) %>% left_join(MB2_dat,.) %>% 
    dplyr::select(marker_idx,site_samp_idx,sp_idx)
  
  
  
  # Choose species 23 to be alpha == 0 for marker 1 and 27 to be alpha==0 for marker 2
  # These are treated as known and so the choice is somewhat arbitrary.
  amp_alpha_mb1 <- log(out$conc_dat$MB1_alpha + 1) - log(out$conc_dat$MB1_alpha[23] + 1)
  amp_alpha_mb2 <- log(out$conc_dat$MB2_alpha + 1) - log(out$conc_dat$MB2_alpha[27] + 1)
  
  # Here are the vector of alphas for use in the analysis.
  amp_alpha_mb1
  amp_alpha_mb2
  
  
  # Summarize observe reads for each replicate.
  Reads <- MB_all_long %>% group_by(site_idx,marker_idx,samp_idx,rep) %>% summarise(N_read = sum(count))
  MB_all_long <- MB_all_long %>% left_join(.,Reads)
  
  reads_mb1 <- MB1_dat %>% dplyr::select(-site_samp_idx,-rep) %>% rowSums()
  reads_mb2 <- MB2_dat %>% dplyr::select(-site_samp_idx,-rep) %>% rowSums()
  
  # Calculate combinatorial coefficient.
  log_mb1_n_k <- MB1_dat %>% dplyr::select(-site_samp_idx,-rep) %>% as.matrix()
  log_mb2_n_k <- MB2_dat %>% dplyr::select(-site_samp_idx,-rep) %>% as.matrix()
  for(i in 1:nrow(MB1_dat)){
    tmp <-(MB1_dat %>% dplyr::select(-site_samp_idx,-rep))[i,] %>% unlist() %>% c()
    log_mb1_n_k[i,] <- as.numeric(lchoose(reads_mb1[i],tmp))
  }
  for(i in 1:nrow(MB2_dat)){
    tmp <-(MB2_dat %>% dplyr::select(-site_samp_idx,-rep))[i,] %>% unlist() %>% c()
    log_mb2_n_k[i,] <- as.numeric(lchoose(reads_mb2[i],tmp))
  }
  
  ### Make Indicator (0 if count>0, 1 if count==0)
  MB1_bias_ind <- MB1_dat %>% dplyr::select(-site_samp_idx,-rep)
  MB1_bias_ind[MB1_bias_ind>0] <- 99
  MB1_bias_ind[MB1_bias_ind==0] <- 1
  MB1_bias_ind[MB1_bias_ind==99] <- 0
  
  MB2_bias_ind <- MB2_dat %>% dplyr::select(-site_samp_idx,-rep)
  MB2_bias_ind[MB2_bias_ind>0] <- 99
  MB2_bias_ind[MB2_bias_ind==0] <- 1
  MB2_bias_ind[MB2_bias_ind==99] <- 0
  
  # Format for NIMBLE -------------------------------------------------------
  # Set concentration thresholds (these were missing from your Stan data)
  nimble_data <- list(
    # Counters
    # Model parameters
    N_pcr = N_pcr,
    
    # Observations of read counts  
    count_mb1 = as.matrix(MB1_dat %>% dplyr::select(-site_samp_idx,-rep)),
    count_mb2 = as.matrix(MB2_dat %>% dplyr::select(-site_samp_idx,-rep)),
    
    reads_mb1 = reads_mb1,
    reads_mb2 = reads_mb2,
    
    # Combinatorial coefficients for the betabinomial
    log_mb1_n_k = log_mb1_n_k,
    log_mb2_n_k = log_mb2_n_k,
    
    # Amplification efficiencies for each marker
    amp_alpha_mb1 = amp_alpha_mb1,
    amp_alpha_mb2 = amp_alpha_mb2,
    
    
    # Indicators for species detection by marker
    MB1_ind = as.matrix(MB1_ind %>% ungroup() %>% dplyr::select(-site_samp_idx)),
    MB1_all =  as.matrix(MB1_ind %>% ungroup() %>% dplyr::select(-site_samp_idx)) %>% apply(.,2, function(col) as.numeric(any(col > 0))) ,
    MB2_ind = as.matrix(MB2_ind %>% ungroup() %>% dplyr::select(-site_samp_idx)),
    MB2_all =  as.matrix(MB2_ind %>% ungroup() %>% dplyr::select(-site_samp_idx)) %>% apply(.,2, function(col) as.numeric(any(col > 0))),
    # Design matrices for mapping samples to their replicates
    X_MB1_ss = matrix(X_MB1_ss,nrow=N_obs_mb1,ncol=N_site_samp),
    X_MB2_ss = matrix(X_MB2_ss,nrow=N_obs_mb2,ncol=N_site_samp),
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
    N_obs_mb2 = as.numeric(N_obs_mb2),
    # Mapping indices
    samp_map1 = as.numeric(ref_sp_mb1$site_samp_idx),
    samp_map2 = as.numeric(ref_sp_mb2$site_samp_idx),
    ref_sp_mb1 = as.numeric(ref_sp_mb1$sp_idx),
    ref_sp_mb2 = as.numeric(ref_sp_mb2$sp_idx)
  )
  
  
  cl<-makeCluster(2,outfile="") # Register cores
  
  # Initialization function
  nimble_inits <- function() {
    N_site <- nimble_constants$N_site
    N_sp <- nimble_constants$N_sp
    N_obs_mb1 <- nimble_constants$N_obs_mb1
    N_obs_mb2 <- nimble_constants$N_obs_mb2
    
    # Initialize d_mb based on count_mb 
    d_mb1_init <- (as.matrix(MB1_dat %>% dplyr::select(-site_samp_idx,-rep)) > 0) * 1
    d_mb2_init <- (as.matrix(MB2_dat %>% dplyr::select(-site_samp_idx,-rep)) > 0) * 1
    
    list(
      # log_D needs to be initialized with a fairly narrow range of values (or possibly initialized conditional and d_mb_init)
      log_D = matrix(rnorm(N_site * N_sp, mean = 0, sd = 1), nrow = N_site, ncol = N_sp),
      beta0_mb1 = -abs(rnorm(1,0,0.2)),
      beta0_mb2 = -abs(rnorm(1,0,0.2)),
      beta1_mb1 = -abs(rnorm(1,0,0.3)),
      beta1_mb2 = -abs(rnorm(1,0,0.3)),
      log_phi1_mb1 = runif(1,0.8,2),
      log_phi1_mb2 = runif(1,0.8,2),
      phi0_mb1 = abs(rnorm(1,0.8,2)),
      phi0_mb2 = abs(rnorm(1,0.8,2)),
      d_mb1 = d_mb1_init,
      d_mb2 = d_mb2_init,
      epsilon_site_mb1_z = matrix(rnorm(N_site*N_sp, 0, 1),nrow=N_site,ncol=N_sp),
      epsilon_site_mb2_z = matrix(rnorm(N_site*N_sp, 0, 1),nrow=N_site,ncol=N_sp),
      sigma_site_mb1 = runif(1, 0.1, 0.5),
      sigma_site_mb2 = runif(1, 0.1, 0.5)
    )
  }
  # Nimble Parallel Run Function -------------------------------------------------------
  run_Nimble_MCMC <- function(chain_info, nimble_data){
    library(nimble)
    ### Custom distributions and functions    
    ### Beta-binomial density
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
    ### Beta-binomial random 
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
    ### Centered log-ratio for detectable species (i.e., species potentially detectable by marker)
    clr_alpha <- nimbleFunction(
      run = function(amp_alpha = double(1),     # amplification efficiencies for all species (on log scale)
                     mb_ind = double(1)) {     # sample-level indicator for this sample
        returnType(double(1))
        n_sp <- length(amp_alpha)
        clr_result <- numeric(n_sp)
        present_alphas <- numeric(0)
        for(j in 1:n_sp) {
          if(mb_ind[j] == 1) {  # detectable species
            present_alphas <- c(present_alphas, amp_alpha[j])  # amp_alpha is already log scale
          }else{present_alphas <- present_alphas}
        }
        # Calculate mean of log amplification efficiencies
        mean_log_alpha <- mean(present_alphas)
        # Calculate CLR for all species, but only meaningful for detectable ones
        for(j in 1:n_sp) {
          clr_result[j] <- amp_alpha[j] - mean_log_alpha  
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
    
    # ModelCode --------------------------------------------------------------
    modelCode <- nimbleCode({
      
      ############# Indexing and suffixes #####
      ## i = observations (i.e., unique technical replicates)
      ## j = species
      ## m = sites (or, equivalently, samples/communities); 
      ## *NB:Here N_site and N_site_samp are equivalent. But this need not be 
      ## *the case (i.e., within a site/community, we may have biological replicates or 
      ## *"site-samples" within a site, each with its own technical replication)
      ## _mb1 = metabarcoding marker 1
      ## _mb2 = metabarcoding marker 2
      
      ###### Priors ###################
      for(m in 1:N_site) {
        for(j in 1:N_sp) {
          log_D[m,j] ~ dnorm(0, sd = 10)
        }
      }
      
      # Phi parameters Marker 1
      beta0_mb1 ~ dnorm(0, sd = 1)
      beta1_mb1 ~ T(dnorm(-10, sd = 5),-Inf,0)
      phi0_mb1 ~ dnorm(2, sd = 1)
      log_phi1_mb1 ~ dnorm(1, sd = .75)
      phi1_mb1 <- exp(log_phi1_mb1)
      
      # Phi parameters Marker 2
      beta0_mb2 ~ dnorm(0, sd = 1) 
      beta1_mb2 ~ T(dnorm(-5, sd = 10),-Inf,0)
      phi0_mb2 ~ dnorm(2, sd = 1)
      log_phi1_mb2 ~ dnorm(1, sd = .75)
      phi1_mb2 <- exp(log_phi1_mb2)
      
      # Random Effects Sigmas
      sigma_site_mb1 ~ T(dnorm(0, sd = 1), 0, Inf)
      sigma_site_mb2 ~ T(dnorm(0, sd = 1), 0, Inf)
      
      # Random Effects
      for(m in 1:N_site) {
        for(j in 1:N_sp){
          epsilon_site_mb1_z[m,j] ~ dnorm(0, 1)
          epsilon_site_mb1[m,j] <- epsilon_site_mb1_z[m,j] * sigma_site_mb1
          epsilon_site_mb2_z[m,j] ~ dnorm(0, 1)
          epsilon_site_mb2[m,j] <- epsilon_site_mb2_z[m,j] * sigma_site_mb2
        }
      }
      
      ######## Matrix operations ########
      log_F[1:N_site, 1:N_sp] <- X_D_to_F[1:N_site, 1:N_site] %*% log_D[1:N_site, 1:N_sp]
      # log_F_mb1[1:N_obs_mb1, 1:N_sp] <- X_MB1_ss[1:N_obs_mb1, 1:N_site_samp] %*% log_F[1:N_site_samp, 1:N_sp]
      # log_F_mb2[1:N_obs_mb2, 1:N_sp] <- X_MB2_ss[1:N_obs_mb2, 1:N_site_samp] %*% log_F[1:N_site_samp, 1:N_sp]
      
      ######## MARKER 1 #################
      # MB1 calculations 
      # Site level
      for(m in 1:N_site) {
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
      # Observation-level Additive Log-Ratios,, Latent Presence Variables, and Copies
      for(i in 1:N_obs_mb1) {
        # extract reference values
        log_lambda_K_mb1_ref[i] <- log_lambda_K_mb1[samp_map1[i],ref_sp_mb1[i]]
        amp_alpha_mb1_ref[i] <- amp_alpha_mb1[ref_sp_mb1[i]]
        for(j in 1:N_sp) {
          log_val_mb1[i,j] <- (log_lambda_K_mb1[samp_map1[i],j] - log_lambda_K_mb1_ref[i]) + 
            N_pcr * (amp_alpha_mb1[j] - amp_alpha_mb1_ref[i])
          # Discrete binary for presence in aliquot
          d_mb1[i,j] ~ dbern(exp(log_1mexp_mb1[samp_map1[i],j]))
          # Copies in aliquot (adjusted for presence)
          lambda_K_d_mb1[i,j] <- d_mb1[i,j] * MB1_all[j] * exp(log_lambda_K_mb1[samp_map1[i],j])
        }
        log_max_val_mb1[i] <- max(log_val_mb1[i,1:N_sp])
      }
      # Site-level phi calculations
      for(m in 1:N_site) {
        # Total copies in sample
        log_Lambda_mb1[m] <- log(sum(lambda_mb1[m,1:N_sp])) 
        clr_alpha_mb1[m,1:N_sp] <- clr_alpha(amp_alpha_mb1[1:N_sp], 
                                             MB1_all[1:N_sp])
        for(j in 1:N_sp) {
          # Calculate phi_mb
          log_phi_mb1[m,j] <- beta0_mb1 + 
            beta1_mb1 * clr_alpha_mb1[m,j] +
            log_Lambda_mb1[m] +
            phi0_mb1 * exp(-phi1_mb1 * log_lambda_K_mb1[m,j])+
            epsilon_site_mb1[samp_map1[m],j] 
          phi_mb1[m,j] <- exp(min(20,log_phi_mb1[m,j]))
        }
      }
      # Another observation level loop
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
          ####### Marker 1 Likelihood #########
          count_mb1[i,j] ~ dbetabin(
            reads_mb1[i], alpha_mb1[i,j], beta_mb1[i,j],
            log_p_zero_mb1[samp_map1[i],j], log_p_pos_mb1[samp_map1[i],j],
            log_mb1_n_k[i,j],d_mb1[i,j]
          )
        }
      }
      
      ######## MARKER 2 #################
      # MB1 calculations 
      # Site level
      for(m in 1:N_site) {
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
      # Observation-level Additive Log-Ratios, Latent Presence Variables, and Copies
      for(i in 1:N_obs_mb2) {
        # extract reference values
        log_lambda_K_mb2_ref[i] <- log_lambda_K_mb2[samp_map2[i],ref_sp_mb2[i]]
        amp_alpha_mb2_ref[i] <- amp_alpha_mb2[ref_sp_mb2[i]]
        for(j in 1:N_sp) {
          log_val_mb2[i,j] <- (log_lambda_K_mb2[samp_map2[i],j] - log_lambda_K_mb2_ref[i]) + 
            N_pcr * (amp_alpha_mb2[j] - amp_alpha_mb2_ref[i])
          # Discrete binary for presence in aliquot
          d_mb2[i,j] ~ dbern(exp(log_1mexp_mb2[samp_map2[i],j]))
          # Copies in aliquot (adjusted for presence)
          lambda_K_d_mb2[i,j] <- d_mb2[i,j] * MB2_all[j] * exp(log_lambda_K_mb2[samp_map2[i],j])
        }
        log_max_val_mb2[i] <- max(log_val_mb2[i,1:N_sp])
      }
      # Site-level phi calculations
      for(m in 1:N_site) {
        # Total copies in sample
        log_Lambda_mb2[m] <- log(sum(lambda_mb2[m,1:N_sp])) 
        # Centered log-ratio
        clr_alpha_mb2[m,1:N_sp] <- clr_alpha(amp_alpha_mb2[1:N_sp], 
                                             MB2_all[1:N_sp])
        for(j in 1:N_sp) {
          # Calculate phi_mb
          log_phi_mb2[m,j] <- beta0_mb2 + 
            beta1_mb2 * clr_alpha_mb2[m,j] +
            log_Lambda_mb2[m] +
            phi0_mb2 * exp(-phi1_mb2 * log_lambda_K_mb2[m,j])+
            epsilon_site_mb2[samp_map2[m],j] 
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
          ####### Marker 2 Likelihood #########
          count_mb2[i,j] ~ dbetabin(
            reads_mb2[i], alpha_mb2[i,j], beta_mb2[i,j],
            log_p_zero_mb2[samp_map2[i],j], log_p_pos_mb2[samp_map2[i],j],
            log_mb2_n_k[i,j],d_mb2[i,j]
          )
        }
      }
    })
    
    
    # Nimble Model Building ---------------------------------------------------
    model <- nimbleModel(code = modelCode,
                         data = nimble_data,
                         constants = nimble_constants,
                         inits = inits,
                         check = TRUE)
    
    # Clean up any leftovers that might be clogging things up
    if (exists("cmodel")) {try(nimble:::clearCompiled(cmodel), silent = TRUE)}
    if (exists("cmcmc")) {try(nimble:::clearCompiled(cmcmc), silent = TRUE)}
    gc()  
    
    # Compile
    cmodel <- compileNimble(model)
    
    # Configure MCMC 
    mcmc_conf <- configureMCMC(model,useCongucacy = FALSE)
    
    # Replace default samplers with multivariate options
    mcmc_conf$removeSamplers(c("log_D[]"))
    for (i in 1:nimble_constants$N_site){
      targets = c()
      for (j in 1:nimble_constants$N_sp){
        targets[j]=paste0("log_D[",i,",",j,"]")
      }
      mcmc_conf$addSampler(target=targets,"AF_slice")
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
    # Add block samplers for phi parameters
    run_targets_mb1 <- c("beta0_mb1","beta1_mb1", "phi0_mb1","log_phi1_mb1")
    mcmc_conf$removeSamplers(run_targets_mb1)
    mcmc_conf$addSampler(target=run_targets_mb1, "AF_slice")    
    run_targets_mb2 <- c("beta0_mb2","beta1_mb2", "phi0_mb2","log_phi1_mb2")
    mcmc_conf$removeSamplers(run_targets_mb2)
    mcmc_conf$addSampler(target=run_targets_mb2, "AF_slice")
    
    # Remove default samplers for the hierarchical epsilon structure
    mcmc_conf$removeSamplers("epsilon_site_mb1_z")
    mcmc_conf$removeSamplers("epsilon_site_mb2_z")
    
    # Add site-level epsilon samplers (these are univariate)
    # mcmc_conf$addSampler(target=c("sigma_site_mb1",paste0("epsilon_site_mb1_z[1:", nimble_constants$N_site, "]")), "AF_slice")
    # mcmc_conf$addSampler(target=c("sigma_site_mb2",paste0("epsilon_site_mb2_z[1:", nimble_constants$N_site, "]")), "AF_slice")
    # 
    # Build and compile MCMC
    mcmc <- buildMCMC(mcmc_conf)
    cmcmc <- compileNimble(mcmc, project = model)
    
    # Run MCMC
    samples <- runMCMC(cmcmc, niter = 4000, nburnin = 2000, thin = 1)
    return(samples)
  }
  
  chain_info <- list(
    list(inits = nimble_inits()),
    list(inits = nimble_inits()))
  
  clusterExport(cl, c("nimble_data", "nimble_constants", "chain_info"))
  
  chain_output <- parLapply(
    cl = cl, X = chain_info,
    fun = run_Nimble_MCMC, 
    nimble_data = nimble_data)
  stopCluster(cl)
  
  nimbleMod_summary = MCMCvis::MCMCsummary(chain_output)
  
  nimble_samples <- do.call(rbind, chain_output)
  
  ## OK.  Make some simple comparisons with the concentration values.
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
  
  ### Save some output cross estimation loop.
  conc_dat_all <- rbind(conc_dat_all,conc_dat_est %>% mutate(rep_id = j))
  post_tot_conc_all <- rbind(post_tot_conc_all,post_tot_conc_summ %>% mutate(rep_id = j))
  
  ### Retain the input data for each replicate
  nimble_data_list[[j]] <- list(nimble_data=nimble_data,rep_id = j)
  nimbleMod_summary_list[[j]]  <- list(nimbleMod_summary,rep_id=j)}

# Just in case there are any open clusters lingering 
# Summarize match between simulated and estimated.
post_tot_conc_all <- post_tot_conc_all %>% mutate(prop_error = (Median - tot_conc_true) / tot_conc_true)
conc_dat_all <- conc_dat_all %>% mutate(prop_error = (Median - conc_val) / conc_val)

post_tot_conc_all <- post_tot_conc_all %>% select(-prop_error)
conc_dat_all <- conc_dat_all %>% select(-prop_error)


OUT <- list(nimbleMod_summary_list = nimbleMod_summary_list,
            conc_dat_all = conc_dat_all , 
            post_tot_conc_all= post_tot_conc_all)

#save(OUT,file="TwoMarker_11-25b.RData")



