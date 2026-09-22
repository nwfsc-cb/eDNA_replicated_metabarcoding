
process_mb2 <- function(out, site_name){
  
  MB1 <- out$MB1
  MB2 <- out$MB2
  conc_dat <- out$conc_dat
  
  N_sp <- conc_dat$N_species
  dat_mb1 <- MB1$Y 
  dat_mb2 <- MB2$Y 
  colnames(dat_mb1) <- paste0("r",1:ncol(dat_mb1))
  colnames(dat_mb2) <- paste0("r",1:ncol(dat_mb2))
  
  conc_dat$site_idx = site_name
  conc_dat$sp_id = 1:N_sp
  dat_mb1 <- data.frame(cbind(sp_id = 1:N_sp,samp_idx = 1, site_idx = site_name, dat_mb1))
  dat_mb2 <- data.frame(cbind(sp_id = 1:N_sp,samp_idx = 1, site_idx = site_name, dat_mb2))
 
  ret <- list(conc_dat = conc_dat,
              MB1 = dat_mb1, 
              MB2 = dat_mb2, 
              N_pcr=out$N_pcr)
  return(ret) 
}


process_mb1 <- function(out, site_name){
  
  MB1 <- out$MB1
  conc_dat <- out$conc_dat
  
  N_sp <- conc_dat$N_species
  dat_mb1 <- MB1$Y 
  colnames(dat_mb1) <- paste0("r",1:ncol(dat_mb1))
  
  conc_dat$site_idx = site_name
  conc_dat$sp_id = 1:N_sp
  dat_mb1 <- data.frame(cbind(sp_id = 1:N_sp,samp_idx = 1, site_idx = site_name, dat_mb1))
  
  ret <- list(conc_dat = conc_dat,
              MB1 = dat_mb1, 
              N_pcr=out$N_pcr)
  return(ret) 
}