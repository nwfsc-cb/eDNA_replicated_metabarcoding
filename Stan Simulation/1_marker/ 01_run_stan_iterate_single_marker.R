# This is a script to call the STAN estimation model
library(here)
library(rstan)
library(tidyverse)
library(data.table)
options(mc.cores = parallel::detectCores()) # rstan options.
source(here("Zenodo","shared_functions","process_mb.R"))
source(here("Zenodo","shared_functions","sim_functions.R"))

MOD_NAME <- "Single_marker_amp_fix_v10.stan"

set.seed(101)
ind_sim <- 1
TOT.CONC <- c(10, 30,100, 300, 1000, 3000, 10000)
n_site = length(TOT.CONC)
N_sp <- 30  
AMP_SD <- 0.01

THEME = paste0("V10.fin_3_rep_amp-sd=",AMP_SD,"_")
NOM <- paste0(THEME,"_SINGLE-SP_",ind_sim,"=sim_2025-DATE",n_site,"site_",N_sp,"sp_conc.Rdata")

sim_scen <- data.frame(N_SP = rep(N_sp,n_site),
                           TOT.CONC = TOT.CONC,
                           CONC.DIST = rep("very_high_skew",n_site),
                           N_READS = rep(100000,n_site),
                           amp_sd = AMP_SD,
                           replicates = rep(3,n_site))

post_tot_conc_all <- NULL
conc_dat_all <- NULL
stan_data_list = NULL
divergences = NULL
stanMod_sampling_list = NULL
stanMod_summary_list = NULL

for(j in 1:ind_sim){# loop over number of independent sim
  print(paste(j,"of",ind_sim))
  OUT <- list()
  
  # Ok Assign a amplification rate to each species for 1 marker
  A1  <- sim_amp(amp_mean = 0.9,
                 amp_sd = sim_scen$amp_sd[1] ,
                 N_species = N_sp)
  
  for( i in 1:nrow(sim_scen)){
    source(here("Zenodo","Stan Simulation","1_marker","simulate_1marker_metabar_iterate.R"))
    OUT[[i]] <- process_mb1(out,site_name=paste0("X",i))
  }
  source(here("Zenodo","Stan Simulation","1_marker","02_process_sim_iterate_single_marker.R"))
  source(here("Zenodo","Stan Simulation","1_marker","03_estimate_stan_iterate_single_marker.R"))
  source(here("Zenodo","Stan Simulation","1_marker","04_post_process_stan_iterate_single_marker.R"))
}

## Pull out any replicates that had divergent transitions and throw them away.

THESE <- divergences %>% as.data.frame() %>% filter(count>0) %>% pull(rep_id)

post_tot_conc_all <- post_tot_conc_all %>% filter(!rep_id %in% c(THESE) )
conc_dat_all <- conc_dat_all %>% filter(!rep_id %in% c(THESE) )

# Summarize match between simulated and estimated.
post_tot_conc_all <- post_tot_conc_all %>% mutate(prop_error = (Median - tot_conc_true) / tot_conc_true)
conc_dat_all <- conc_dat_all %>% mutate(prop_error = (Median - conc_val) / conc_val)

post_tot_conc_all_summ <- post_tot_conc_all %>% group_by(site_id,tot_conc_true) %>% 
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
                      
BRKS = c(0.001,0.01,0.1,1,10,30,100,300,1000,3000,10000,30000)
BRKS2 = c(0.01,0.1,1,10,30,100,300,1000,3000,10000,30000)
BRKS.prop = c(-2,-1,-0.5,-0.25,0,0.25,0.5,1,2)

# conc_dat_all <- conc_dat_all %>% filter(!rep_id %in% c(35,28,15))
# post_tot_conc_all <- post_tot_conc_all %>% filter(!rep_id %in% c(35,28,15))

p_true_est_conc_rep  <- ggplot(conc_dat_all) +
  geom_point(aes(x=conc_val,y=Median,color=site_id))+
  geom_errorbar(aes(x=conc_val,ymin=q_05,ymax=q_95,color=site_id),width=0,alpha=0.1)+
  scale_x_continuous(trans="log","True conc.",breaks=BRKS) +
  scale_y_continuous(trans="log","Est conc.",breaks=BRKS) +
  coord_cartesian(xlim=c(min(BRKS),max(BRKS)),ylim=c(min(BRKS),max(BRKS))) +
  geom_abline(intercept=0,slope=1,linetype="dashed", color="red")+
  theme_bw()

p_true_est_tot_conc_rep <- ggplot(post_tot_conc_all) +
  geom_point(aes(x=tot_conc_true,y=Median,color=site_id),alpha=0.5)+
  geom_errorbar(aes(x=tot_conc_true,ymin=q_05,ymax=q_95,color=site_id),width=0,alpha=0.3)+
  scale_x_continuous(trans="log", "True total conc.",breaks=BRKS) +
  scale_y_continuous(trans="log","Est total conc.",breaks=BRKS) +
  #coord_cartesian(xlim=c(min(BRKS),max(BRKS)),ylim=c(min(BRKS),max(BRKS))) +
  geom_abline(intercept=0,slope=1,linetype="dashed", color="red")+
  theme_bw()

p_prop_error_conc_rep <- ggplot(conc_dat_all) +
  geom_point(aes(x=conc_val,y=prop_error),alpha=0.5,shape=".")+
  #geom_jitter(aes(x=tot_conc_true,y=prop_error),height=0,width=0.5,alpha=0.5)+#geom_errorbar(aes(x=tot_conc_true,ymin=q_05,ymax=q_95,color=site_id),width=0,alpha=0.3)+
  scale_x_continuous(trans="log","True total conc.",breaks=BRKS) +
  scale_y_continuous("Proportional Error",breaks=BRKS.prop) +
  coord_cartesian(xlim=c(min(BRKS2),NA),ylim=c(-1,max(BRKS.prop))) +
  geom_hline(yintercept =0, linetype="dashed", color="red")+
  facet_wrap(~as.factor(conc_total)) +
  theme_bw()

p_prop_error_tot_conc_rep <- ggplot(post_tot_conc_all) +
  geom_boxplot(aes(x=tot_conc_true,y=prop_error,group=tot_conc_true),alpha=0.5)+
  geom_jitter(aes(x=tot_conc_true,y=prop_error),height=0,width=0.25,alpha=0.5)+#geom_errorbar(aes(x=tot_conc_true,ymin=q_05,ymax=q_95,color=site_id),width=0,alpha=0.3)+
  scale_x_continuous(trans="log", "True total conc.",breaks=BRKS) +
  #scale_y_continuous(trans="log","Proportional Error",breaks=BRKS) +
  #coord_cartesian(xlim=c(min(BRKS),max(BRKS)),ylim=c(min(BRKS),max(BRKS))) +
  geom_hline(yintercept =0, linetype="dashed", color="red")+
  theme_bw()

 p_prop_error_tot_conc_rep2 <- ggplot(post_tot_conc_all) +
  geom_boxplot(aes(x=tot_conc_true,y=prop_error,group=tot_conc_true),alpha=0.5)+
  geom_jitter(aes(x=tot_conc_true,y=prop_error),height=0,width=0.25,alpha=0.5)+#geom_errorbar(aes(x=tot_conc_true,ymin=q_05,ymax=q_95,color=site_id),width=0,alpha=0.3)+
  geom_line(aes(x=tot_conc_true,y=prop_error,group=rep_id),alpha=0.5)+
  scale_x_continuous(trans="log", "True total conc.",breaks=BRKS) +
  #scale_y_continuous(trans="log","Proportional Error",breaks=BRKS) +
  #coord_cartesian(xlim=c(min(BRKS),max(BRKS)),ylim=c(min(BRKS),max(BRKS))) +
  geom_hline(yintercept =0, linetype="dashed", color="red")+
  theme_bw()

plots = list( p_true_est_conc_rep=p_true_est_conc_rep,
              p_true_est_tot_conc_rep = p_true_est_tot_conc_rep,
              p_prop_error_conc_rep =p_prop_error_conc_rep,
              p_prop_error_tot_conc_rep =p_prop_error_tot_conc_rep)

OUT <- list(stan_data_list = stan_data_list,
            stanMod= stanMod,
            divergences = divergences,
            stanMod_summary_list = stanMod_summary_list,
            stanMod_sampling_list = stanMod_sampling_list,
            conc_dat_all = conc_dat_all , 
            post_tot_conc_all= post_tot_conc_all,
            plots=plots)

#save(OUT,file=here("data","Model_est_1marker_sim",NOM))
 
