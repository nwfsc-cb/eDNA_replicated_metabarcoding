# Convert NIMBLE output to Stan format for two-marker workflow
# This function converts NIMBLE MCMC samples to match the structure expected by the Stan post-processing pipeline


  # Extract dimensions from constants
  N_site <- nimble_constants$N_site
  N_sp <- nimble_constants$N_sp
  N_obs_mb1 <- nimble_constants$N_obs_mb1
  N_obs_mb2 <- nimble_constants$N_obs_mb2
  N_site_samp <- nimble_constants$N_site_samp
  
  # Get number of MCMC iterations
  n_iter <- nrow(nimble_samples)
  
  # Initialize pars list in Stan format
  pars <- list()
  
  # 1. Convert log_D (most important for concentration estimates)
  # NIMBLE: log_D[1:N_site, 1:N_sp] - each element stored as separate column
  # Stan format: array[n_iter, N_site, N_sp]
  
  pars$log_D <- array(NA, dim = c(n_iter, N_site, N_sp))
  
  # Fill array element by element to ensure correct mapping
  log_D_found <- 0
  for(i in 1:N_site) {
    for(j in 1:N_sp) {
      col_name <- paste0("log_D[", i, ", ", j, "]")
      if(col_name %in% colnames(nimble_samples)) {
        pars$log_D[, i, j] <- nimble_samples[, col_name]
        log_D_found <- log_D_found + 1
      }
    }
  }
  
  if(log_D_found == 0) {
    warning("No log_D parameters found in NIMBLE samples")
    pars$log_D <- NULL
  } else if(log_D_found < N_site * N_sp) {
    warning(paste("Only found", log_D_found, "of", N_site * N_sp, "log_D parameters"))
  }
  
  # 2. Convert phi parameters for both markers
  # phi_mb1: iterations x observations x species
  pars$phi_mb1 <- array(NA, dim = c(n_iter, N_obs_mb1, N_sp))
  phi_mb1_found <- 0
  
  for(i in 1:N_obs_mb1) {
    for(j in 1:N_sp) {
      col_name <- paste0("phi_mb1[", i, ", ", j, "]")
      if(col_name %in% colnames(nimble_samples)) {
        pars$phi_mb1[, i, j] <- nimble_samples[, col_name]
        phi_mb1_found <- phi_mb1_found + 1
      }
    }
  }
  
  if(phi_mb1_found == 0) {
    warning("No phi_mb1 parameters found in NIMBLE samples")
    pars$phi_mb1 <- NULL
  }
  
  pars$phi_mb2 <- array(NA, dim = c(n_iter, N_obs_mb2, N_sp))
  phi_mb2_found <- 0
  
  for(i in 1:N_obs_mb2) {
    for(j in 1:N_sp) {
      col_name <- paste0("phi_mb2[", i, ", ", j, "]")
      if(col_name %in% colnames(nimble_samples)) {
        pars$phi_mb2[, i, j] <- nimble_samples[, col_name]
        phi_mb2_found <- phi_mb2_found + 1
      }
    }
  }
  
  if(phi_mb2_found == 0) {
    warning("No phi_mb2 parameters found in NIMBLE samples") 
    pars$phi_mb2 <- NULL
  }
  
  # 3. Reconstruct log_Lambda_copy (total copies)
  if(!is.null(pars$log_D)) {
    X_D_to_F <- nimble_data$X_D_to_F
    MB1_ind <- nimble_data$MB1_ind
    MB2_ind <- nimble_data$MB2_ind
    log_Lambda_copy_array <- array(NA, dim = c(n_iter, N_site_samp, 2))
    
    for(iter in 1:n_iter) {
      # Get log_D for this iteration
      log_D_iter <- pars$log_D[iter, , ]
      
      # Calculate log_F: site-sample level copies
      log_F_iter <- X_D_to_F %*% log_D_iter
      
      # Calculate total copies for each marker at each site-sample
      for(i in 1:N_site_samp) {
        # Marker 1 (MB1)
        log_Lambda_copy_array[iter, i, 1] <- log(sum(exp(log_F_iter[i, ]) * MB1_ind[i, ]))
        # Marker 2 (MB2) 
        log_Lambda_copy_array[iter, i, 2] <- log(sum(exp(log_F_iter[i, ]) * MB2_ind[i, ]))
      }
    }
    
    pars$log_Lambda_copy <- log_Lambda_copy_array
  }
  
  param_names <- c("log_beta0_mb1", "log_beta0_mb2", "log_beta0_mock_mb1", "log_beta0_mock_mb2",
                   "log_phi0_mb1", "log_phi0_mb2", "phi1_mb1", "phi1_mb2", 
                   "amp_alpha_mb1", "amp_alpha_mb2")
  
  for(param in param_names) {
    if(param %in% colnames(nimble_samples)) {
      pars[[param]] <- nimble_samples[, param]
    }
  }
  
  for(marker in 1:2) {
    marker_name <- ifelse(marker == 1, "mb1", "mb2")
    amp_alpha_array <- array(NA, dim = c(n_iter, N_sp))
    amp_found <- 0
    for(j in 1:N_sp) {
      col_name <- paste0("amp_alpha_", marker_name, "[", j, "]")
      if(col_name %in% colnames(nimble_samples)) {
        amp_alpha_array[, j] <- nimble_samples[, col_name]
        amp_found <- amp_found + 1
      }
    }
    if(amp_found > 0) {
      # Stan format expects [species, iterations] for vectors
      pars[[paste0("amp_alpha_", marker_name)]] <- t(amp_alpha_array)
    }
  }

    # Set up variables needed for post-processing (mimicking Stan workflow)
  N_marker <- 2
  MARKER <- 2
  NAME1 <- "MFU"
  NAME2 <- "MV1"
  # Create stan_data equivalent for post-processing
  stan_data <- list(
    site_samp_mb1_idx = as.numeric(nimble_constants$samp_map1),
    site_samp_mb2_idx = as.numeric(nimble_constants$samp_map2),
    ref_sp_mb1 = nimble_constants$ref_sp_mb1,
    ref_sp_mb2 = nimble_constants$ref_sp_mb2
  )
  
  # Create temporary variables to match Stan workflow expectations
  stan_data_tmp <- stan_data
  
  ### Summarise log D (from Mock_analysis_03_post_summary.R)
  A <- apply(exp(pars$log_D),c(2,3),mean) %>% as_tibble() %>% 
    mutate(site_idx = 1:nrow(.)) %>% 
    pivot_longer(.,-site_idx,names_to="sp_id",values_to ="Mean")
  B <- apply(exp(pars$log_D),c(2,3),sd) %>% as_tibble() %>% 
    mutate(site_idx = 1:nrow(.)) %>% 
    pivot_longer(.,-site_idx,names_to="sp_id",values_to="SD")
  C <- apply(exp(pars$log_D),c(2,3),quantile,probs=0.5) %>% as_tibble() %>% 
    mutate(site_idx = 1:nrow(.)) %>% 
    pivot_longer(.,-site_idx,names_to="sp_id",values_to="Median")
  D <- apply(exp(pars$log_D),c(2,3),quantile,probs=0.025) %>% as_tibble() %>% 
    mutate(site_idx = 1:nrow(.)) %>% 
    pivot_longer(.,-site_idx,names_to="sp_id",values_to="q_2.5")
  E <- apply(exp(pars$log_D),c(2,3),quantile,probs=0.975) %>% as_tibble() %>% 
    mutate(site_idx = 1:nrow(.)) %>% 
    pivot_longer(.,-site_idx,names_to="sp_id",values_to="q_97.5")
  
  post <- left_join(A,B) %>% left_join(.,C) %>% left_join(.,D) %>% left_join(.,E)
  post <- post %>% mutate(sp_idx = sp_id, sp_idx= gsub('V','',sp_idx))
  
  ### Summarize concentrations
  conc.true <- dat.long.h %>% dplyr::select(Species,sp_idx,site_idx,dilution,mock,skew,mock.skew,species_copies) %>% 
    group_by(Species,sp_idx,site_idx,dilution,mock,skew,mock.skew) %>% 
    summarise(true.copies = mean(species_copies),sd.true.copies = sd(species_copies)) %>%
    mutate(true.copies.dil = true.copies /dilution) %>% 
    mutate(sp_idx = as.character(sp_idx)) %>% 
    left_join(.,post %>% dplyr::select(-sp_id))
  
  ### Summarize total copies across all species.
  if(!is.null(pars$log_Lambda_copy)) {
    log_Lambda_copy_post <- pars$log_Lambda_copy[,,1] %>% as.data.frame() %>%  
      mutate(marker_idx=1) %>% 
      pivot_longer(.,cols= !matches("marker_idx")) %>% 
      mutate(site_idx= substr(name,2,nchar(name)))
    
    log_Lambda_copy_post <- pars$log_Lambda_copy[,,2] %>% as.data.frame() %>%  
      mutate( marker_idx=2) %>% 
      pivot_longer(.,cols= !matches("marker_idx")) %>% 
      mutate(site_idx= substr(name,2,nchar(name))) %>% 
      rbind(log_Lambda_copy_post,.)
    
    # Fix species_names
    Lambda_copy_post <- log_Lambda_copy_post %>%  
      mutate(Lambda = exp(value)) %>% 
      group_by(site_idx,marker_idx) %>% 
      summarise(Mean = mean(Lambda),
                Median =median(Lambda),
                q2.5 = quantile(Lambda,probs=c(0.025)),
                q25 = quantile(Lambda,probs=c(0.25)),
                q75 = quantile(Lambda,probs=c(0.75)),
                q97.5 = quantile(Lambda,probs=c(0.975))) %>% 
      left_join(.,SITE %>% mutate(site_idx = as.character(site_idx)))
    
    # Merge in the true copy
    Lambda_copy_post <- left_join(Lambda_copy_post, dat.conc.summ) %>% mutate(tot_copy_dil = tot_copy / dilution)
  } else {
    Lambda_copy_post <- NULL
    warning("Could not calculate Lambda_copy_post - log_Lambda_copy not available")
  }
  
  ######### Summarise the phi values (if available)
  if(!is.null(pars$phi_mb1)) {
    A <- apply(pars$phi_mb1,c(2,3),mean) %>% as_tibble() %>% 
      mutate(site_samp_mb1_idx = stan_data_tmp$site_samp_mb1_idx) %>% 
      pivot_longer(.,-site_samp_mb1_idx,names_to="sp_id",values_to ="mean_phi_mb1") %>% 
      group_by(site_samp_mb1_idx,sp_id) %>% 
      summarize(Mean_phi_mb1 = mean(mean_phi_mb1))
    B <- apply(pars$phi_mb1,c(2,3),sd) %>% as_tibble() %>% 
      mutate(site_samp_mb1_idx = stan_data_tmp$site_samp_mb1_idx) %>% 
      pivot_longer(.,-site_samp_mb1_idx,names_to="sp_id",values_to ="sd_phi_mb1") %>% 
      group_by(site_samp_mb1_idx,sp_id) %>% 
      summarize(SD_phi_mb1 = mean(sd_phi_mb1))
    C <- apply(pars$phi_mb1,c(2,3),median) %>% as_tibble() %>% 
      mutate(site_samp_mb1_idx = stan_data_tmp$site_samp_mb1_idx) %>% 
      pivot_longer(.,-site_samp_mb1_idx,names_to="sp_id",values_to ="median_phi_mb1") %>% 
      group_by(site_samp_mb1_idx,sp_id) %>% 
      summarize(Median_phi_mb1 = mean(median_phi_mb1))
    
    phi_mb1 <- left_join(A,B) %>% left_join(.,C) %>% mutate(sp_idx = sp_id, sp_idx= gsub('V','',sp_idx))
    
    conc.true <- left_join(conc.true, phi_mb1 %>% dplyr::select(-sp_id) %>% rename(site_idx =site_samp_mb1_idx))
  }
  
  if(!is.null(pars$phi_mb2)) {
    A <- apply(pars$phi_mb2,c(2,3),mean) %>% as_tibble() %>% 
      mutate(site_samp_mb2_idx = stan_data_tmp$site_samp_mb2_idx) %>% 
      pivot_longer(.,-site_samp_mb2_idx,names_to="sp_id",values_to ="mean_phi_mb2") %>% 
      group_by(site_samp_mb2_idx,sp_id) %>% 
      summarize(Mean_phi_mb2 = mean(mean_phi_mb2))
    B <- apply(pars$phi_mb2,c(2,3),sd) %>% as_tibble() %>% 
      mutate(site_samp_mb2_idx = stan_data_tmp$site_samp_mb2_idx) %>% 
      pivot_longer(.,-site_samp_mb2_idx,names_to="sp_id",values_to ="sd_phi_mb2") %>% 
      group_by(site_samp_mb2_idx,sp_id) %>% 
      summarize(SD_phi_mb2 = mean(sd_phi_mb2))
    C <- apply(pars$phi_mb2,c(2,3),median) %>% as_tibble() %>% 
      mutate(site_samp_mb2_idx = stan_data_tmp$site_samp_mb1_idx) %>% 
      pivot_longer(.,-site_samp_mb2_idx,names_to="sp_id",values_to ="median_phi_mb1") %>% 
      group_by(site_samp_mb2_idx,sp_id) %>% 
      summarize(Median_phi_mb2 = mean(median_phi_mb1))
    
    phi_mb2 <- left_join(A,B) %>% left_join(.,C) %>% mutate(sp_idx = sp_id, sp_idx= gsub('V','',sp_idx))
    
    conc.true <- left_join(conc.true, phi_mb2 %>% dplyr::select(-sp_id)%>% rename(site_idx =site_samp_mb2_idx))
  }
  
  # Finalize conc.true
  conc.true <- conc.true %>% 
    mutate(est_conc_val_K = Mean / (1-exp(-Mean))) %>%
    group_by(site_idx) %>% 
    mutate(tot_copy = sum(Mean)) %>% 
    ungroup()
  
  # Plots of analysis output.
  # p_amp_MFU <- ggplot(amp_post) + 
  #               geom_point(aes(x=Species,y=mean,color=marker)) +
  #               geom_errorbar(aes(x=Species,ymin=X2.5.,ymax=X97.5.,color=marker),width=0) +
  #               geom_hline(yintercept = 0,color="red",linetype="dashed")+
  #               scale_y_continuous(expression(alpha)) +
  #               theme_bw() +            
  #               theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  # 
  BRKS <- c(0.01,0.1,1,10,100,1000,10000,100000)
  BRKS_char <- BRKS %>% as.character()
  
  ### Pred-Obs Concentration
  p_pred_v_obs1 <- ggplot(conc.true) +
    geom_point(aes(x = true.copies.dil, y = Mean, color = as.factor(dilution))) +
    geom_abline(slope=1,intercept=0) +
    facet_grid(skew~mock) +
    scale_color_viridis_d("Dilution",option="turbo")+
    scale_y_continuous("Predicted copies",trans="log",breaks=BRKS,labels=BRKS_char) +
    scale_x_continuous("True copies", trans="log",breaks=BRKS,labels=BRKS_char) +
    theme_bw()
  
  p_pred_v_obs2 <- ggplot(conc.true) +
    geom_point(aes(x = true.copies.dil, y = Mean, color = mock.skew)) +
    geom_abline(slope=1,intercept=0) +
    facet_wrap(~dilution,ncol=2) +
    scale_color_viridis_d("Community",option="turbo")+
    scale_y_continuous("Predicted copies",trans="log",breaks=BRKS,labels=BRKS_char) +
    scale_x_continuous("True copies", trans="log",breaks=BRKS,labels=BRKS_char) +
    theme_bw()
  
  ### Pred-Obs Total Conc
  p_pred_v_obs_tot_conc1 <- ggplot(Lambda_copy_post %>% filter(marker_idx == 1)) +
    geom_point(aes(x = tot_copy_dil, y = Mean, color = as.factor(dilution))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_errorbar(aes(x=tot_copy_dil, ymin=q2.5, ymax=q97.5, color=as.factor(dilution)), width=0) +
    scale_color_viridis_d("Dilution", option="turbo")+
    scale_y_continuous("Predicted copies", trans="log", breaks=BRKS, labels=BRKS_char) +
    scale_x_continuous("True copies", trans="log", breaks=BRKS, labels=BRKS_char) +
    ggtitle("Total Copies across all species") +
    theme_bw()
  
  #######
  # Extract phi shape parameters from each iteration
  #######
  p_phi_func_mb1 <- ggplot(conc.true) +
    #geom_point(aes(x=est_conc_val_K,y=Mean_phi,color=as.factor(conc_total),group=rep_id),alpha=0.3)+
    geom_point(aes(x=est_conc_val_K,y=Mean_phi_mb1,color=tot_copy),alpha=0.3)+
    geom_line(aes(x=est_conc_val_K,y=Mean_phi_mb1,color=tot_copy,group=site_idx),alpha=0.3)+
    scale_y_continuous(expression(phi),trans="log",breaks=BRKS) +
    scale_x_continuous(expression(lambda[K]), trans="log",breaks=BRKS) +
    scale_color_viridis_c("Total DNA \nCopies",option="turbo",trans="log",breaks=BRKS) +
    facet_wrap(~dilution,ncol=2) +
    ggtitle(paste("MB1=",NAME1)) +
    theme_bw()
  if(MARKER ==2){
    p_phi_func_mb2 <- ggplot(conc.true) +
      #geom_point(aes(x=est_conc_val_K,y=Mean_phi,color=as.factor(conc_total),group=rep_id),alpha=0.3)+
      geom_point(aes(x=est_conc_val_K,y=Mean_phi_mb2,color=tot_copy),alpha=0.3)+
      geom_line(aes(x=est_conc_val_K,y=Mean_phi_mb2,color=tot_copy,group=site_idx),alpha=0.3)+
      scale_y_continuous(expression(phi),trans="log",breaks=BRKS) +
      scale_x_continuous(expression(lambda[K]), trans="log",breaks=BRKS) +
      scale_color_viridis_c("Total DNA \nCopies",option="turbo",trans="log",breaks=BRKS) +
      facet_wrap(~dilution,ncol=2) +
      ggtitle(paste("MB2=",NAME2)) +
      theme_bw()  
  }
  
  plots <- list(
    p_pred_v_obs1 = p_pred_v_obs1,
    p_pred_v_obs2 = p_pred_v_obs2,
    p_pred_v_obs_tot_conc1 = p_pred_v_obs_tot_conc1,
    p_phi_func_mb1 = p_phi_func_mb1)
  if(MARKER==2){
    plots <- c(plots,
               list(p_phi_func_mb2 = p_phi_func_mb2))
  }
  
  
  
  # Write some things to file.
  pdf(here("figures","Mock Output",paste("MOCK OUTPUT TEMP NIMBLE;","EPS_site & EPS_sp",".pdf")),onefile=TRUE,width=8,height=8)
  print(p_pred_v_obs1)
  print(p_pred_v_obs2)
  print(p_pred_v_obs_tot_conc1)
  print(p_phi_func_mb1)
  if(MARKER==2){
    print(p_phi_func_mb2)
  }
  dev.off()
  
  Two_Marker <- list(
    stanMod_summary = NULL,  # NIMBLE doesn't have this
    conc.true = conc.true,
    Lambda_copy_post = Lambda_copy_post,
    stan_data = stan_data,
    pars = pars,
    plots = NULL,  # Generate these using Mock_analysis_04_plots.R
    samp_params = NULL  # NIMBLE doesn't have this
  )
  
