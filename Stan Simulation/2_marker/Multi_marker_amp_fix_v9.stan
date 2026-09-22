
functions { // User-Defined Functions.
  real zero_like_func(real log_p_occur_zero, real log_p_occur_pos, real log_n_k, 
                        int count, int reads,
                        real alpha, real beta) {
                    real val ;
                    val = log_sum_exp(log_p_occur_zero,
                                    log_p_occur_pos +
                                    log_n_k +
                                    lbeta(count + alpha, reads - count + beta) -
                                    lbeta(alpha, beta)) ;
                    return val ;
  }

  real pos_like_func(real log_p_occur_pos, real log_n_k, 
                        int count, int reads,
                        real alpha, real beta) {
                    real val ;
                    val = log_p_occur_pos + 
                            log_n_k +
                            lbeta(count + alpha,reads - count + beta) -
                            lbeta(alpha, beta) ;
                    return val ;
  }
  
  vector alpha_func(vector Pi, vector log_phi){
                  vector[rows(Pi)] val;
                  val = Pi .* exp(log_phi);
                  return val ;
  }

  vector beta_func(vector alpha, vector Pi){
                  vector[rows(Pi)] val ;
                  val = alpha .* (1-Pi) ./  Pi ;
                  return val ;
  }
}

data { // 
  int N_sp; // Number of species in data
  int N_marker;  // Number of markers (only works for 2 as of now)
  int N_site;  // Number of observed samples for individual sites.
  int N_site_samp; // Number of site and sample combinations across all markers.
  int N_obs_mb1; // Number of observed rows in MB1 data
  int N_obs_mb2; // Number of observed rows in MB2 data
  real N_pcr; // number of PCR cycles used for each metabarcoding
  
  // Observed data of community matrices
    array[N_obs_mb1,N_sp] int count_mb1; // observed reads
    array[N_obs_mb2,N_sp] int count_mb2; // observed reads
  
    array[N_obs_mb1] int reads_mb1 ; // total reads
    array[N_obs_mb2] int reads_mb2 ; // total reads
    
    matrix[N_obs_mb1,N_sp] log_mb1_n_k ; // n choose k coefficient mb1
    matrix[N_obs_mb2,N_sp] log_mb2_n_k ; // n choose k coefficient mb2

   // Design matrices: field samples
    // map site-sample to observations.
    matrix[N_obs_mb1, N_site_samp] X_MB1_ss ; 
    matrix[N_obs_mb2, N_site_samp] X_MB2_ss ;
    // map site-level copy number to site-sample level copies.
    matrix[N_site_samp, N_site] X_D_to_F ;
  
  // These are matrices of selecting values of the F matrix that should be included
  // in calculating the total copies in each sample.
    matrix[N_site_samp,N_sp] MB1_ind ; 
    matrix[N_site_samp,N_sp] MB2_ind ; 
  
  // Reference species for each site-samp and marker
    array[N_obs_mb1] int ref_sp_mb1;
    array[N_obs_mb2] int ref_sp_mb2;
 
  /// These are unused indices for the moment.
    array[N_obs_mb1] int site_samp_mb1_idx; 
    array[N_obs_mb2] int site_samp_mb2_idx;
  
  // Read in amplification Estimates.
    vector[N_sp] amp_alpha_mb1;
    vector[N_sp] amp_alpha_mb2;
  
  // Priors
  real prior_log_D_mu;
  real prior_log_D_sig;

  real prior_beta1_mu;
  real prior_beta1_sig;
  
  real prior_beta0_mu;
  real prior_beta0_sig;
  
  real prior_gamma0_mu; 
  real prior_gamma0_sig; 
  
  real prior_log_gamma1_mu ;
  real prior_log_gamma1_sig ;
  
  #real prior_amp_tau ; # SD of variability among amp efficiencies
  real prior_sigma_epsilon; 
}

transformed data{
  array[N_obs_mb1] real log_reads_mb1 = log(reads_mb1);
  array[N_obs_mb2] real log_reads_mb2 = log(reads_mb2);

  matrix[N_site_samp,N_sp] kappa_tmp ;
  vector[N_site_samp] kappa_sum ;
  vector[N_site_samp] kappa_N ;
  matrix[N_site_samp,N_sp] kappa_mb1 ;
  matrix[N_site_samp,N_sp] kappa_mb2 ;
  
  /// This calculates the centered log-ratio for onluy species that were detected at least once
  /// in a sample.
  for(j in 1:N_site_samp){
    kappa_tmp[j,] = to_row_vector(amp_alpha_mb1) .* MB1_ind[j,] ;
    kappa_sum[j] = sum(kappa_tmp[j,]) ;
    kappa_N[j] = sum(MB1_ind[j,]) ;
    kappa_mb1[j,] = to_row_vector(amp_alpha_mb1) - (kappa_sum[j]/kappa_N[j]) ;
  }

  for(j in 1:N_site_samp){
    kappa_tmp[j,] = to_row_vector(amp_alpha_mb2) .* MB2_ind[j,] ;
    kappa_sum[j] = sum(kappa_tmp[j,]) ;
    kappa_N[j] = sum(MB2_ind[j,]) ;
    kappa_mb2[j,] = to_row_vector(amp_alpha_mb2) - (kappa_sum[j]/kappa_N[j]) ;
  }

}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  // Estimated amplification efficiencies.
  vector[N_marker] beta0 ;   // intercept of parameter determining overdispersion
                              // (log over-disp variance at copies = reads)
  vector<upper=0>[N_marker]  beta1 ;   // slope of parameter determining overdispersion
  
  vector[N_marker] log_gamma1 ;    // slope of parameter determining overdispersion
  vector<lower=0>[N_marker] gamma0 ;
  
  vector<lower=0>[N_marker] sigma_epsilon; // random effect SD
  //vector<lower=0>[N_marker] sigma_epsilon; // random effect SD
  array[N_marker] matrix[N_site,N_sp] epsilon_raw ; // random effect for each site-marker combination.
  
  matrix[N_site,N_sp] log_D; // True log-copies of each species at each site.
}

transformed parameters{
  matrix[N_site_samp,N_sp] log_F; // Long-form matrix for all combinations of sites and samples 
                                  // Allows for non-matching replicates in each marker.
                                  // This is the base value through which all markers are connect.
  
  // These are the set of observations (including replicates for each MB marker)
  matrix[N_obs_mb1,N_sp]  log_F_mb1 ; // Log-copies observed by marker 1 across all species.
  matrix[N_obs_mb2,N_sp]  log_F_mb2 ; // Log-copies observed by marker 2 across all species.

  # slope for phi parameter.
  vector[N_marker] gamma1 ;

  // Lambda (copies observed by each marker)
  matrix[N_site_samp,N_marker] log_Lambda_copy; 
  matrix[N_site_samp,N_marker] log_Lambda_K_copy; 

  array[N_marker] matrix[N_site_samp,N_sp] epsilon ; // random effect for each site-marker combination.

  matrix[N_obs_mb1,N_sp] pi_samp_mb1; // Proportion Mean mb1
  matrix[N_obs_mb2,N_sp] pi_samp_mb2; // Proportion Mean mb2 

  matrix[N_obs_mb1,N_sp] log_phi_mb1; // Overdispersion parameter for each observation mb1
  matrix[N_obs_mb2,N_sp] log_phi_mb2; // Overdispersion parameter for each observation mb2
  
  // Define the copies at each site as a function of covariates and spatial variables here:
  
  // PLACEHOLDER log_D = COV + SPACE;

 // Define the relationship between the site-level DNA copies and
  //  the sample within a site level.  
  // For now there are no replicate biological samples taken at each site, so this is a 1:1 relationship.
  
  log_F = X_D_to_F * log_D ; // You could add various random effects here if you have multiple biological samples.
    
 // Define the relationship between the copies for each marker and 
    // the shared, site-sample level copies.  Allow for a multiplicative offset if of interest. 
    // for now, there is no offset for any marker.
  log_F_mb1 = X_MB1_ss * log_F; 
  log_F_mb2 = X_MB2_ss * log_F;
  
  /// BEGIN CONVERSION OF COPIES TO PROPORTIONS.
  {// local variables for making reference species vectors
      vector[N_obs_mb1] log_lambda_K_mb1_ref ;
      vector[N_obs_mb2] log_lambda_K_mb2_ref ;
      
      vector[N_obs_mb1] amp_alpha_mb1_ref ;
      vector[N_obs_mb2] amp_alpha_mb2_ref ;
      
      matrix[N_obs_mb1,N_sp] log_lambda_K_mb1 ;
      matrix[N_obs_mb2,N_sp] log_lambda_K_mb2 ;

      matrix[N_obs_mb1,N_sp] log_val_m1 ;
      matrix[N_obs_mb2,N_sp] log_val_m2 ;
      
      // local variables
      real tmp_A;
      gamma1 = exp(log_gamma1);
      
    // This calculates the mean of the zero-truncated Poisson (zero observations are gone)
      log_lambda_K_mb1 = log_F_mb1  - log1m_exp(-exp(log_F_mb1)) ;
      log_lambda_K_mb2 = log_F_mb2  - log1m_exp(-exp(log_F_mb2)) ;
      
     // Make helper vectors observations
    for(i in 1:N_obs_mb1){
       log_lambda_K_mb1_ref[i] = log_lambda_K_mb1[i,ref_sp_mb1[i]];
       amp_alpha_mb1_ref[i] = amp_alpha_mb1[ref_sp_mb1[i]];
    }
    for(i in 1:N_obs_mb2){
       log_lambda_K_mb2_ref[i] = log_lambda_K_mb2[i,ref_sp_mb2[i]];
       amp_alpha_mb2_ref[i] = amp_alpha_mb2[ref_sp_mb2[i]];
    }

    // Calculate log-scale mean predictions for each species in each sample.
    for (n in 1:N_sp) {
      # field obs
      log_val_m1[,n] = (log_lambda_K_mb1[,n] - log_lambda_K_mb1_ref) + N_pcr*(amp_alpha_mb1[n] - amp_alpha_mb1_ref);
      log_val_m2[,n] = (log_lambda_K_mb2[,n] - log_lambda_K_mb2_ref) + N_pcr*(amp_alpha_mb2[n] - amp_alpha_mb2_ref);
    }
    // Convert to log-scale mean predictions to proportions for each species in each marker.
    for(m in 1:N_obs_mb1){
      pi_samp_mb1[m,] = to_row_vector(softmax(to_vector(log_val_m1[m,]))); // proportion of each taxon in field samples
    }
    for(m in 1:N_obs_mb2){
      pi_samp_mb2[m,] = to_row_vector(softmax(to_vector(log_val_m2[m,]))); // proportion of each taxon in field samples
    }

  //print(log_F);
  // calculate the total copies observed for each marker. (only include species that were observed non-zero times among replicates)
  for(i in 1:N_site_samp){
    log_Lambda_copy[i,1] = log(sum(exp(log_F[i,]) .* MB1_ind[i,])) ; 
    log_Lambda_copy[i,2] = log(sum(exp(log_F[i,]) .* MB2_ind[i,])) ; 
    
    log_Lambda_K_copy[i,1] = log(sum(exp(log_F[i,] - log1m_exp(-exp(log_F[i,]))) .* MB1_ind[i,])) ; 
    log_Lambda_K_copy[i,2] = log(sum(exp(log_F[i,] - log1m_exp(-exp(log_F[i,]))) .* MB2_ind[i,])) ;
  }

  for(i in 1:N_marker){
    for(j in 1:N_site){
      epsilon[i,j,] =  epsilon_raw[i,j,] * sigma_epsilon[i]; // -0.5*pow(tau_epsilon[i],2) +
    }
  }

  // calculate the phi values for each site-sample-species in each marker.
  // phi is constant for larger copy numbers

    for(i in 1: N_obs_mb1){
     log_phi_mb1[i,] = beta0[1] + rep_row_vector(log_Lambda_copy[site_samp_mb1_idx[i],1],N_sp) +  
                              fmax(gamma0[1] - gamma1[1] * log_lambda_K_mb1[i,],0) +
                              beta1[1] * kappa_mb1[site_samp_mb1_idx[i],] +
                              epsilon[1,site_samp_mb1_idx[i],] ;
    }
  #print("log_phi1",log_phi_mb1);
   for(i in 1: N_obs_mb2){
      log_phi_mb2[i,] = beta0[2] + rep_row_vector(log_Lambda_copy[site_samp_mb2_idx[i],2],N_sp) +  
                              fmax(gamma0[2] - gamma1[2] * log_lambda_K_mb2[i,],0) +
                              beta1[2] * kappa_mb2[site_samp_mb2_idx[i],] +
                              epsilon[2,site_samp_mb2_idx[i],] ;
    }
  #print("log_phi2",log_phi_mb1);
  }// end local variables.
}

model{
  {// LOCAL VARIABLES USEFUL FOR EVALUATING THE LIKELIHOODS
   /// Matrices used in calculating likelihoods. 
  matrix[N_obs_mb1,N_sp] log_p_occur_zero_mb1 ;
  matrix[N_obs_mb2,N_sp] log_p_occur_zero_mb2 ;
    
  matrix[N_obs_mb1,N_sp] log_p_occur_pos_mb1  ;
  matrix[N_obs_mb2,N_sp] log_p_occur_pos_mb2  ;

  matrix[N_obs_mb1,N_sp] alpha_mb1;
  matrix[N_obs_mb1,N_sp] beta_mb1;

  matrix[N_obs_mb2,N_sp] alpha_mb2;
  matrix[N_obs_mb2,N_sp] beta_mb2;

  // Calcuate appropriate p_occur for zero obs and for pos obs
  // These are in terms of the log-copies
    log_p_occur_zero_mb1 = -exp(log_F_mb1) ;
    log_p_occur_zero_mb2 = -exp(log_F_mb2) ;
    
    log_p_occur_pos_mb1 = log1m_exp(log_p_occur_zero_mb1) ;
    log_p_occur_pos_mb2 = log1m_exp(log_p_occur_zero_mb2) ;
  
  // Convert pi and phi to alpha and beta for likelihoood calculations...
    for(m in 1:N_sp){
      alpha_mb1[,m] = alpha_func(pi_samp_mb1[,m], log_phi_mb1[,m]) ;
      beta_mb1[,m]  = beta_func(alpha_mb1[,m], pi_samp_mb1[,m]) ;
    }
    for(m in 1:N_sp){
      alpha_mb2[,m] = alpha_func(pi_samp_mb2[,m],log_phi_mb2[,m]) ;
      beta_mb2[,m]  = beta_func(alpha_mb2[,m], pi_samp_mb2[,m]) ;
    }
    // Loop over likelihoods for metabarcoding marker 1
    for(i in 1:N_obs_mb1){
      for(j in 1:N_sp){
        if(count_mb1[i,j]==0){
          target += zero_like_func(log_p_occur_zero_mb1[i,j], log_p_occur_pos_mb1[i,j], log_mb1_n_k[i,j], 
                                    count_mb1[i,j], reads_mb1[i],
                                    alpha_mb1[i,j], beta_mb1[i,j]) ;
        }else{
          target += pos_like_func(log_p_occur_pos_mb1[i,j], log_mb1_n_k[i,j], 
                                    count_mb1[i,j], reads_mb1[i],
                                    alpha_mb1[i,j], beta_mb1[i,j]) ;
        }
      }
    }// end mb1 obs loop
    // Loop over likelihoods for metabarcoding marker 2
    for(i in 1:N_obs_mb2){
      for(j in 1:N_sp){
        if(count_mb2[i,j]==0){
          target += zero_like_func(log_p_occur_zero_mb2[i,j], log_p_occur_pos_mb2[i,j], log_mb2_n_k[i,j], 
                                    count_mb2[i,j], reads_mb2[i],
                                    alpha_mb2[i,j], beta_mb2[i,j]) ;
        }else{
          target += pos_like_func(log_p_occur_pos_mb2[i,j], log_mb2_n_k[i,j], 
                                    count_mb2[i,j], reads_mb2[i],
                                    alpha_mb2[i,j], beta_mb2[i,j]) ;
        }
      }
    } // end mb2 obs loop
  } // END LOCAL VARIABLES.
  
  /// Priors
    for(i in 1:N_sp){
      log_D[,i] ~ normal(prior_log_D_mu, prior_log_D_sig) ;
    }
    
    for(i in 1:N_marker){
      for(j in 1:N_site){
        epsilon_raw[i,j,] ~ std_normal();
      }
    }
    
    // PRIORS
    for(i in 1:N_marker){
        sigma_epsilon[i] ~ normal(0,prior_sigma_epsilon) ;
    }
    
    // log_beta1 ~ normal(prior_log_beta1_mu ,prior_log_beta1_sig);
    beta0 ~ normal(prior_beta0_mu , prior_beta0_sig);
    beta1 ~ normal(prior_beta1_mu , prior_beta1_sig);
    gamma0 ~ normal(prior_gamma0_mu, prior_gamma0_sig);
    log_gamma1 ~ normal(prior_log_gamma1_mu, prior_log_gamma1_sig);
}

generated quantities{
  matrix[N_site_samp,N_marker] Lambda_copy = exp(log_Lambda_copy);
  matrix[N_obs_mb1,N_sp] phi_mb1 = exp(log_phi_mb1);
  matrix[N_obs_mb2,N_sp] phi_mb2 = exp(log_phi_mb2);
}
