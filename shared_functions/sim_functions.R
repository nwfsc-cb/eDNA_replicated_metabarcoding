### Functions for simulations of metabarcoding
require(tidyverse)

sim_conc <- function(N_species = 20, 
                     conc_dist = "even", # options = mild_skew, mod_skew, high_skew, very_high_skew
                     tot_conc # total concentration in the sample.
                     ){
    
  # First, generate distribution concentration in a couple of discrete scenarios
  if(conc_dist!="even"){
    if(conc_dist == "mild_skew"){
       alpha_par <- 10
    }else if(conc_dist == "mod_skew"){
       alpha_par <- 5
    }else if(conc_dist == "high_skew"){
       alpha_par <- 1
    }else if(conc_dist == "very_high_skew"){
      alpha_par <- 0.5
    }
    stick_param <- rgamma(N_species,alpha_par,1)
    conc_prop <- stick_param / sum(stick_param)
  }
  
    if(conc_dist == "even"){
       conc_prop <- rep(1 / N_species,N_species) 
    }
  # conc_val is the thing that gets pushed forward.
    conc_val <- conc_prop * tot_conc

    out <- list(N_species= N_species,
                tot_conc = tot_conc,
                conc_prop = conc_prop,
                conc_val = conc_val)
    return(out)
}

# Function to simulated amplification
sim_amp <- function(amp_mean = 0.85,
                    amp_sd = 0,
                    N_species){
    if(amp_sd>0){
      # use method of moments to calculate beta distribution parameters.
      alpha = amp_mean* (amp_mean*(1-amp_mean)/amp_sd^2 - 1)
      beta  = alpha * (1-amp_mean) / amp_mean
      a_val <- rbeta(N_species,alpha,beta)
    }else{ 
      a_val <- rep(amp_mean,N_species)
    }
  
  clr_alpha <- a_val - exp(mean(log(a_val)))
  
  out <- list(a_val=a_val,
              clr_alpha=clr_alpha)
  return(out)
}


# Simulate Amplification and Metabarcoding sampling.
sim_MB <- function(out,
                   a_val,
                   replicates = 10,
                   N_tot_reads_min = 50000,
                   N_tot_reads_max = 100000,
                   N_pcr = 40){
    # First, simulate the number of copies in each replicate sample.
      # rows are species, columns are replicates.
      G <- rpois(replicates*length(out$conc_val),out$conc_val) %>% matrix(.,ncol=replicates)
  
      G_prop = G*0
      amp <- G * 0
      prop_true <- G * 0
      Y <- G * 0
      Y_prop <- G * 0
      # Amplification, calculate proportions, Calculate observed reads.
      
      for(i in 1:ncol(G)){
        if(sum(G[,i])>0){
          G_prop = G[,i] / sum(G[,i])
          amp[,i] <- exp(log(G[,i]) + N_pcr * log(1+a_val))
          prop_true[,i] <- amp[,i] / sum(amp[,i])
          Y[,i] <- rmultinom(1, runif(1,N_tot_reads_min,N_tot_reads_max),prop_true[,i])
          Y_prop[,i] <- Y[,i] / sum(Y[,i])
        }
      }
      # Calculate 
      
      out <- list(N_pcr =N_pcr,
                  G = G, # Realized DNA copies in a tube.
                  G_prop = G_prop, #  Proportion of copies in each tube.
                  prop_true = prop_true, # Proportion true after amplification
                  Y = Y, # Observed after sampling of sequences.
                  Y_prop = Y_prop
                  )

}

sim_scenarios <- function(N_species,
                          tot_conc_vec, 
                          conc_dist,
                          amp_mean, 
                          amp_sd,
                          n_tech_rep = 10,
                          n_comm_rep = 100,
                          N_tot_reads = c(50000,100000),
                          N_pcr = 40
                          ) {

  out_sim <- NULL
  
  for(i in 1:length(tot_conc_vec)){
    for(j in 1:n_comm_rep){
      # This represents the true DNA concentration in a tube.
      X <- sim_conc(N_species = N_species, 
                    conc_dist = conc_dist, # options = even mild_skew, mod_skew, high_skew, very_high_skew
                    tot_conc= tot_conc_vec[i])
      
      # This represents the amplification rate for each species.
      amp_X <- sim_amp(amp_mean = amp_mean, # among species
                       amp_sd=amp_sd, # among species
                       N_species=X$N_species) 
      
      # Make metabarcoding observations for an arbitrary number of replicate observations
      obs <- sim_MB(X,
                    a_val = amp_X$a_val,
                    replicates = n_tech_rep,
                    N_tot_reads_min = N_tot_reads[1],
                    N_tot_reads_max = N_tot_reads[2],
                    N_pcr = N_pcr)
      
      # summaries
      p_zero_G = rowSums(obs$G==0) / n_tech_rep
      p_zero_Y = rowSums(obs$Y==0) / n_tech_rep
      prop_mean = rowMeans(obs$Y_prop)
      prop_sd = apply(obs$Y_prop,1,sd)
      
      cond_prop_mean <- rep(0,X$N_species)
      cond_prop_sd <- rep(0,X$N_species)
      for(k in 1:X$N_species){
        G <- obs$G[k,] 
        Z <- obs$Y_prop[k,] 
        # conditional probabilities based on at least 1 copy being present in the aliquot.
        cond_prop_mean[k] <- mean(Z[G>0])
        cond_prop_sd[k] <- sd(Z[G>0])
      }
      
      tmp <- data.frame(sim_id = paste0(tot_conc_vec[i],"_",j),
                        conc = X$conc_val,
                        prop_conc = X$conc_prop,
                        tot_conc=tot_conc_vec[i],
                        p_zero_G,
                        p_zero_Y,
                        prop_mean,
                        prop_sd,
                        cond_prop_mean,
                        cond_prop_sd,
                        n_tech_rep,
                        a_val = amp_X$a_val,
                        clr_alpha = amp_X$clr_alpha)
      out_sim <- rbind(out_sim,tmp)
      
    } # end j loop
  } # end tot_conc loop

return(out_sim)
}
