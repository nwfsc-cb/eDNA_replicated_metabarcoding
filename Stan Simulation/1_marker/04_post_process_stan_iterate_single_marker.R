
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

A <- apply(exp(pars$log_D),c(2,3),mean) %>% as_tibble() %>% 
  mutate(site_idx = 1:nrow(.)) %>% 
  pivot_longer(.,-site_idx,names_to="sp_id",values_to ="Mean")
B <- apply(exp(pars$log_D),c(2,3),sd) %>% as_tibble() %>% 
  mutate(site_idx = 1:nrow(.)) %>% 
  pivot_longer(.,-site_idx,names_to="sp_id",values_to="SD")
C <- apply(exp(pars$log_D),c(2,3),quantile,probs=0.5) %>% as_tibble() %>% 
  mutate(site_idx = 1:nrow(.)) %>% 
  pivot_longer(.,-site_idx,names_to="sp_id",values_to="Median")
D <- apply(exp(pars$log_D),c(2,3),quantile,probs=0.05) %>% as_tibble() %>% 
  mutate(site_idx = 1:nrow(.)) %>% 
  pivot_longer(.,-site_idx,names_to="sp_id",values_to="q_05")
E <- apply(exp(pars$log_D),c(2,3),quantile,probs=0.95) %>% as_tibble() %>% 
  mutate(site_idx = 1:nrow(.)) %>% 
  pivot_longer(.,-site_idx,names_to="sp_id",values_to="q_95")

post <- left_join(A,B) %>% left_join(.,C) %>% left_join(.,D) %>% left_join(.,E)
post <- post %>% mutate(sp_idx = sp_id, sp_idx= gsub('V','',sp_idx))

conc_dat_est <- left_join(conc_dat %>% mutate(sp_id=as.character(sp_id)),
                          post %>% dplyr::select(-sp_id),by = join_by(sp_id == sp_idx,site_idx) )

# Calculate the total Concentration
post_tot_conc <- apply(exp(pars$log_D),c(2),rowSums) %>% as_tibble()
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

conc_dat_tot_conc <- stanMod_summary$Lambda_conc

######### Summarise the phi values
A <- apply(pars$phi_mb1,c(2,3),mean) %>% as_tibble() %>% 
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

BRKS = c(10,30,100,300,1000,3000,10000,30000)
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
stan_data_list[[j]] <- list(stan_data=stan_data,rep_id = j)
stanMod_sampling_list[[j]] <- list(stanMod_sampling=stanMod_sampling,rep_id =j)
stanMod_summary_list[[j]]  <- list(stanMod_summary=stanMod_summary,rep_id=j)

if(j==1){
  divergences <- c(rep_id = j, count = count_divergences)
}else{
  divergences <- rbind(divergences, c(rep_id = j, count = count_divergences))
}
