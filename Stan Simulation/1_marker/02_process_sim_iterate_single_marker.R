
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


### Priors
prior_log_D_mu <- 0
prior_log_D_sig <- 10
prior_beta0_mu <- 0
prior_beta0_sig <- 1
prior_beta1_mu <- -10
prior_beta1_sig <- 5
prior_gamma0_mu <- 2
prior_gamma0_sig <- 1
prior_log_gamma1_mu <- 1
prior_log_gamma1_sig <- 0.75
prior_amp_sigma <- 0.01
prior_sigma_epsilon <- 1

############################################################
#### Get files ready for stan.
############################################################
stan_data <- list(
  # Counters
  N_sp = N_sp,
  N_marker = N_marker,
  N_site = N_site,
  N_site_samp = N_site_samp,
  N_pcr = N_pcr,
  N_obs_mb1 = N_obs_mb1,
  
  # Observations of read counts  
  count_mb1 = MB1_dat %>% dplyr::select(-site_samp_idx,-rep),
  reads_mb1 = reads_mb1,
  log_mb1_n_k = log_mb1_n_k, # n choose k coefficient mb1

   # design matrices for mapping samples to their replicates.
  X_MB1_ss = X_MB1_ss,

  X_D_to_F = X_D_to_F,
  # Indicator files for finding non-zero observations for each species
  # and making a total concentration estimate 
  MB1_ind = MB1_ind %>% ungroup() %>%  dplyr::select(-site_samp_idx),
  vec_one_sp = rep(1,N_sp), # a row_vector useful for calculating row-sums.
  site_samp_mb1_idx = ref_sp_mb1$site_samp_idx,

  # Amplification efficiency for each marker
  amp_alpha_mb1 = amp_alpha_mb1,

  N_pcr = N_pcr,
  ref_sp_mb1 = ref_sp_mb1$sp_idx,
  
  #Priors
  prior_log_D_mu = prior_log_D_mu,
  prior_log_D_sig = prior_log_D_sig,
  prior_beta0_mu = prior_beta0_mu,
  prior_beta0_sig = prior_beta0_sig,
  prior_beta1_mu = prior_beta1_mu,
  prior_beta1_sig = prior_beta1_sig,
  prior_gamma0_mu = prior_gamma0_mu,
  prior_gamma0_sig = prior_gamma0_sig,
  prior_log_gamma1_mu = prior_log_gamma1_mu,
  prior_log_gamma1_sig = prior_log_gamma1_sig,
  prior_amp_sigma = prior_amp_sigma, # SD of variability among amp efficiencies
  prior_sigma_epsilon = prior_sigma_epsilon
)

