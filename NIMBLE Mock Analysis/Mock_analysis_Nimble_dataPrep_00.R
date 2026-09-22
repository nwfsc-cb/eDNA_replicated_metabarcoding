### 
### Analyze the mock communities 
###

library(here)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(rstan)

## Read in the data.
# Pull in metabarcoding reads first.

#Species to omit
dat.mfu <- read.csv(here("Zenodo","Data","Mocks_combined_curated_MFU_taxon_table_wide.csv"))
dat.mv1 <- read.csv(here("Zenodo","Data","Mocks_combined_curated_MV1_taxon_table_wide.csv"))

dat.conc  <- read.csv(here("Zenodo","Data","mock_construction_concs.csv"))

dat.conc <- dat.conc %>% mutate(mock = substr(Mock_name,5,nchar(Mock_name))) %>% 
  filter(!Species %in% NO_SP) %>% # get rid of one messed up species.
  rename(skew=Even_skew) %>% 
  mutate(Species = ifelse(Species == "Artedius fenetralis","Artedius fenestralis",Species)) %>%  # Fix species_name typo
  group_by(Mock_name, skew, Species) %>% ## Combine Artedius fenestralis into one line
  summarise(Exp_mtDNA_prop = sum(Expected_mtDNA_prop),
            mock_conc_sp_copies_ul = sum(Mock_conc_species_copies_ul),
            mock_total_conc_copies_ul = mean(Mock_total_conc_copies_ul)) %>% 
  mutate(species_copies = mock_conc_sp_copies_ul * 2) %>% 
  mutate(mock = substr(Mock_name,5,9)) %>% 
  ungroup()

dat.conc.summ <- dat.conc %>% group_by(mock,skew) %>% 
  summarise(tot_conc = mean(mock_total_conc_copies_ul)) %>% 
  mutate(tot_copy = tot_conc*2)

dat.conc2 <- read.csv(here("Zenodo","Data","mock_ddPCR_quants.csv")) 
dat.conc2 <- dat.conc2 %>% filter(!mock=="NTC") %>% 
  mutate(skew = substr(mock,10,13),
         mock = substr(mock,4,8)) %>% 
  group_by(mock,skew,dilution) %>% 
  summarise(tot_conc_Mean_dd = mean(conc_sample_copies_ul),tot_conc_SD_dd = sd(conc_sample_copies_ul))

# Pull out labels from 
dat.mfu.long <- melt(data.table(dat.mfu), id.vars = c("BestTaxon", "Class")) %>% 
  filter(!BestTaxon %in% NO_SP) %>% 
  rename(ID=variable) %>% 
  mutate(ID = as.character(ID), 
         ID2 = gsub("_.*", "", ID),
         marker= substr(ID,1,3),
         group = substr(ID,5,7),
         mock = substr(ID,8,13), mock = gsub("\\.","",mock),
         skew = substr(ID,14,18), skew = gsub("\\.","",skew),
         mock.skew = paste0(mock,".",skew),
         dilution = sub(".*d","",ID2),dilution = gsub("\\..*","",dilution),
         dilution = as.numeric(dilution),
         tech_rep1 = substr(ID,nchar(ID2),nchar(ID2)),
         tech_rep2 = substr(ID,nchar(ID2)-2,nchar(ID2)-2),
         tech_rep = ifelse(tech_rep1 %in% c("h","l"),tech_rep2,tech_rep1),
         seq_depth = ifelse(grepl("h",ID),"h","h"),
         seq_depth = ifelse(grepl("l",ID),"l",seq_depth)) %>% 
  dplyr::select(-tech_rep1,-tech_rep2)

dat.mv1.long <- melt(data.table(dat.mv1), id.vars = c("BestTaxon", "Class")) %>% 
  filter(!BestTaxon %in% NO_SP) %>% 
  rename(ID=variable) %>% 
  mutate(ID = as.character(ID), 
         ID2 = gsub("_.*", "", ID),
         marker= substr(ID,1,3),
         group = substr(ID,5,7),
         mock = substr(ID,8,13), mock = gsub("\\.","",mock),
         skew = substr(ID,14,18), skew = gsub("\\.","",skew),
         mock.skew = paste0(mock,".",skew),
         dilution = sub(".*d","",ID2),dilution = gsub("\\..*","",dilution),
         dilution = as.numeric(dilution),
         tech_rep1 = substr(ID,nchar(ID2),nchar(ID2)),
         tech_rep2 = substr(ID,nchar(ID2)-2,nchar(ID2)-2),
         tech_rep = ifelse(tech_rep1 %in% c("h","l"),tech_rep2,tech_rep1),
         seq_depth = ifelse(grepl("h",ID),"h","h"),
         seq_depth = ifelse(grepl("l",ID),"l",seq_depth))%>% 
  dplyr::select(-tech_rep1,-tech_rep2)


## summarize to make sure things look right.
dat.mfu.long %>% distinct(marker,group,mock,dilution)
dat.mv1.long %>% distinct(marker,group,mock,dilution)

# Combine 
dat.long <- rbind(dat.mfu.long,dat.mv1.long) %>% dplyr::select(-ID2)
# Merge in the concentrations
dat.long <- dat.long %>% left_join(.,dat.conc.summ %>% dplyr::select(mock,skew,tot_conc))

### OK. calculate some summaries 
dat.summ <- dat.long %>% group_by(ID,marker,group,mock,skew,mock.skew,dilution,seq_depth) %>% 
  mutate(total.reads = sum(value)) %>% 
  ungroup() %>% 
  mutate(prop = value / total.reads)  
dat.summ <- dat.summ %>% 
  left_join(., dat.conc %>% ungroup() %>%  dplyr::select(mock,skew, BestTaxon = Species,species_copies) ) %>% 
  mutate(species_copies_dil = species_copies / dilution) %>% 
  mutate(tot_conc_dil = tot_conc / dilution) %>% 
  mutate(value_bin = ifelse(value>0,1,0),
         prop_copies = species_copies_dil / tot_conc_dil) %>% 
  left_join(.,dat.conc2) %>% 
  # derive dd_concentraion
  mutate(species_copies_dd = prop_copies * tot_conc_Mean_dd)

dat.conc.compare <- dat.summ %>% distinct(mock,skew,dilution,tot_conc_dil) %>% 
  left_join(.,dat.conc2)

# # Make a plot of dd concentration vs. just dilution series.
# ggplot(dat.conc.compare) +
#   geom_point(aes(x=tot_conc_dil,y=tot_conc_Mean_dd)) +
#   geom_errorbar(aes(x=tot_conc_dil,
#                     ymax=tot_conc_Mean_dd+tot_conc_SD_dd,
#                     ymin=tot_conc_Mean_dd-tot_conc_SD_dd),width=0) +
#   scale_x_continuous(trans="log",breaks=c(10,30,100,300,1000,3000,10000,30000,100000)) +
#   scale_y_continuous(trans="log",breaks=c(10,30,100,300,1000,3000,10000,30000,100000)) +
#   geom_abline(intercept=0,slope=1,linetype="dashed") +
#   theme_bw()

####
dat.reads <- dat.summ %>% 
  distinct(marker,group,mock,skew,dilution,seq_depth,total.reads) 
dat.reads.summ <- dat.reads %>% 
  group_by(marker,group,mock,skew,dilution,seq_depth) %>% 
  summarise(read.mean = mean(total.reads), 
            read.sd = sd(total.reads))

dat.reads <- left_join(dat.reads,dat.reads.summ)

dat.prop.summ <- dat.summ %>% group_by(BestTaxon,marker,group,mock,skew,mock.skew,dilution,seq_depth) %>% 
  summarise(mean_prop = mean(prop),sd_prop = sd(prop)) %>% 
  left_join(.,dat.conc.summ) %>% 
  mutate(tot_conc_dil = tot_conc / dilution)

############################################################
### Identify somes sample to use to estimate amplification efficiency.
############################################################

## DECLARE A REFERENCE SPECIES FOR AMP BIAS AND RENAME FOR EFFICIENCY.
REF.SP <- "Engraulis mordax" # Choose Anchovy

dat.long <- dat.long %>% mutate(Species = ifelse(BestTaxon == REF.SP,paste0("zzz_",REF.SP),BestTaxon))

# Choose the even, undiluted samples for estimating amp efficiency.
dat.mock <- dat.long %>% filter(seq_depth=="h",skew=="even",dilution == 1) %>% 
  left_join(.,dat.conc %>% dplyr::select(mock,skew,BestTaxon=Species,species_copies)) %>% 
  mutate(species_copies = ifelse(is.na(species_copies),0,species_copies))
# convert species that were observed with species_copies == 0.
# This is mostly for Chinook in mock2.
dat.mock <- dat.mock %>% mutate(value = ifelse(species_copies==0,0,value)) %>% arrange(Species)

# cull the samples from the bigger data.frame
dat.long <- dat.long %>% filter(!ID %in% dat.mock$ID) %>% arrange(Species)

# Find the unique sites for mocks

site_mock <- dat.mock %>% 
  distinct(mock,skew,dilution) %>% mutate(site_idx = 1:nrow(.))
N_site_mock <- nrow(site_mock)

# Marker
marker_mock <- dat.mock %>% distinct(marker) %>% mutate(marker_idx = 1:nrow(.))
N_marker_mock <- nrow(marker_mock)

# Find the unique sites in any marker
site_marker_mock <- dat.mock %>% distinct(marker,mock,skew,dilution)

site_marker_mock <- site_marker_mock %>% mutate(site_marker_idx = 1:nrow(.)) %>% 
  left_join(.,site_mock) %>% left_join(.,marker_mock)

N_site_samp_mock <- nrow(site_marker_mock)

species_mock <- dat.mock %>% distinct(Species) %>% mutate(Species = as.character(Species)) %>% 
  arrange(Species) %>% mutate(sp_idx = 1: nrow(.))
N_sp_mock <- nrow(species_mock)

SPECIES <- species_mock

#############################################
# add relevant indices to data.
#############################################
dat.mock <- dat.mock %>% left_join(.,species_mock) %>% left_join(.,site_marker_mock) 

# Make a site_samp index... in this case we only have one sample so site_samp will be identical to site_idx
site_samp_mock <- dat.mock %>% distinct(site_idx)
site_samp_mock <- site_samp_mock %>% mutate(site_samp_idx = 1:nrow(.)) %>% as.data.frame()
N_site_samp_mock <- nrow(site_samp_mock)

dat.mock <- dat.mock %>% left_join(.,site_samp_mock) 

# Arrange in species order so that the below matrices will be ordered correctly.
dat.mock <- dat.mock %>% arrange(sp_idx)

#############################################
# Make MARKER-SPECIFIC OBSERVATION MATRICES 
#############################################
# Make wide form matrices for each marker (columns= species, rows =replicates + site-sample pairs)
MFU_mock <-pivot_wider(dat.mock %>% filter(marker=="MFU"),id_cols=c("marker","site_samp_idx","tech_rep"),
                       values_from = value,
                       names_from = sp_idx)
N_obs_mb1_mock <- nrow(MFU_mock)

MV1_mock <-pivot_wider(dat.mock %>% filter(marker=="MV1"),id_cols=c("marker","site_samp_idx","tech_rep"),
                       values_from = value,
                       names_from = sp_idx)
N_obs_mb2_mock <- nrow(MV1_mock)

# Make vectors of total reads for mock communities.
Reads <- dat.mock %>% group_by(site_idx,marker_idx,site_samp_idx,tech_rep) %>% summarise(N_read = sum(value))
dat.mock <- dat.mock %>% left_join(.,Reads)

reads_mb1_mock <- MFU_mock %>% dplyr::select(-marker,-site_samp_idx,-tech_rep) %>% rowSums()
reads_mb2_mock <- MV1_mock %>% dplyr::select(-marker,-site_samp_idx,-tech_rep) %>% rowSums()


#############################################
# make design matrices for mapping site_sample pairs to observations
#############################################
form <- "site_samp_idx ~ 0 + factor(site_idx)"
model_frame   <- model.frame(form, site_samp_mock)  
X_D_to_F_mock <- model.matrix(as.formula(form), model_frame)


form <- "tech_rep ~ 0 + factor(site_samp_idx)"
model_frame   <- model.frame(form, MFU_mock)  
X_MFU_ss_mock <- model.matrix(as.formula(form), model_frame)

# make design matrices
form <- "tech_rep ~ 0 + factor(site_samp_idx)"
model_frame   <- model.frame(form, MV1_mock)  
X_MV1_ss_mock <- model.matrix(as.formula(form), model_frame)



###########################################################
###########################################################
###########################################################
###########################################################
####  MOve on to unknown samples.
####  Get the datasets together to feed a statistical model.
###########################################################
###########################################################
###########################################################
###########################################################

# only keep species in the mock community
dat.long <- dat.long %>% filter(Species %in% SPECIES$Species )

dat.mfu.wide <- dat.long %>% filter(marker == "MFU",seq_depth=="h") %>% 
  pivot_wider(id_cols = c("mock","skew","dilution","tech_rep","seq_depth"),
              names_from = BestTaxon, values_from=value) 

dat.mv1.wide <- dat.long %>% filter(marker == "MV1",seq_depth=="h") %>% 
  pivot_wider(id_cols = c("mock","skew","dilution","tech_rep","seq_depth"),
              names_from = BestTaxon, values_from=value)

dat.mfu.wide$N_read <- rowSums(dat.mfu.wide[,6:ncol(dat.mfu.wide)])
dat.mv1.wide$N_read <- rowSums(dat.mv1.wide[,6:ncol(dat.mv1.wide)])

dat.mv1.wide %>% as.data.frame() %>% head()

##############################################
# Find the unique sites
# filter for just high read replicates
dat.long.h <- dat.long %>% filter(seq_depth=="h") %>% 
  left_join(.,dat.conc %>% dplyr::select(mock,skew,BestTaxon=Species,species_copies)) %>% 
  mutate(species_copies = ifelse(is.na(species_copies),0,species_copies)) 

# get rid of Chinook from mock2 because of contamination from Diaphus.
dat.long.h <- dat.long.h %>% 
  mutate(value= ifelse(mock=="mock2" & BestTaxon == "Oncorhynchus tshawytscha",0,value))

site <- dat.long.h %>% 
  distinct(mock,skew,dilution) %>% mutate(site_idx = 1:nrow(.))
N_site <- nrow(site)
SITE <- site

# Marker
marker <- dat.long.h %>% distinct(marker) %>% mutate(marker_idx = 1:nrow(.))
N_marker <- nrow(marker)

# Find the unique sites in any marker
site_marker <- dat.long.h %>% distinct(marker,mock,skew,dilution)

site_marker <- site_marker %>% mutate(site_marker_idx = 1:nrow(.)) %>% 
  left_join(.,site) %>% left_join(.,marker)

species <- dat.long.h %>% distinct(Species) %>% arrange(Species) %>% mutate(sp_idx = 1: nrow(.))
N_sp <- nrow(species)

#############################################
# add relevant indices to data.
#############################################
dat.long.h <- dat.long.h %>% left_join(.,species) %>% left_join(.,site_marker) 

# Make a site_samp index... in this case we only have one sample so site_samp will be identical to site_idx
site_samp <- dat.long.h %>% distinct(site_idx)
site_samp <- site_samp %>% mutate(site_samp_idx = 1:nrow(.)) %>% as.data.frame()
N_site_samp <- nrow(site_samp)

dat.long.h <- dat.long.h %>% left_join(.,site_samp) %>% arrange(sp_idx)

#############################################
# Make MARKER-SPECIFIC OBSERVATION MATRICES 
#############################################
# Make wide form matrices for each marker (columns= species, rows =replicates + site-sample pairs)
MFU_dat <-pivot_wider(dat.long.h %>% filter(marker=="MFU"),id_cols=c("marker","site_samp_idx","tech_rep"),
                      values_from = value,
                      names_from = sp_idx)
N_obs_mb1 <- nrow(MFU_dat)

MV1_dat <-pivot_wider(dat.long.h %>% filter(marker=="MV1"),id_cols=c("marker","site_samp_idx","tech_rep"),
                      values_from = value,
                      names_from = sp_idx)
N_obs_mb2 <- nrow(MV1_dat)

#############################################
# make design matrices for mapping site_level copies 
# to site_sample copies
#############################################
form <- "site_samp_idx ~ 0 + factor(site_idx)"
model_frame   <- model.frame(form, site_samp)  
X_D_to_F <- model.matrix(as.formula(form), model_frame)

#############################################
# make design matrices for mapping site_sample pairs to observations
#############################################
form <- "tech_rep ~ 0 + factor(site_samp_idx)"
model_frame   <- model.frame(form, MFU_dat)  
X_MB1_ss <- model.matrix(as.formula(form), model_frame)

# make design matrices
form <- "tech_rep ~ 0 + factor(site_samp_idx)"
model_frame   <- model.frame(form, MV1_dat)  
X_MB2_ss <- model.matrix(as.formula(form), model_frame)

# Summarise which species were observed in each sample.
TMP <- dat.long.h %>% group_by(marker_idx,site_samp_idx,sp_idx) %>% summarise(Sum = sum(value)) %>% 
  mutate(IND = ifelse(Sum>0,1,0)) %>% arrange(sp_idx)
MB1_ind <- pivot_wider(TMP %>% filter(marker_idx==1),id_cols=c("site_samp_idx"),
                       values_from = IND,
                       names_from = sp_idx)
MB2_ind <- pivot_wider(TMP %>% filter(marker_idx==2),id_cols=c("site_samp_idx"),
                       values_from = IND,
                       names_from = sp_idx)

TMP.MOCK <- dat.mock %>% group_by(marker_idx,site_samp_idx,sp_idx) %>% summarise(Sum = sum(value)) %>% 
  mutate(IND = ifelse(Sum>0,1,0)) %>% arrange(sp_idx)
MB1_mock <- pivot_wider(TMP.MOCK %>% filter(marker_idx==1),id_cols=c("site_samp_idx"),
                       values_from = IND,
                       names_from = sp_idx)
MB2_mock <- pivot_wider(TMP.MOCK %>% filter(marker_idx==2),id_cols=c("site_samp_idx"),
                       values_from = IND,
                       names_from = sp_idx)

##### OK. Make Mock Matrix of known COPIES (not Conc) AND Mock matrix of observations
# use Engraulis mordax (anchovy) as the reference species for amplification efficiency (i.e. alpha = 0)
MFU_dat_mock_copy <-pivot_wider(dat.mock %>% filter(marker=="MFU"),id_cols=c("marker","site_samp_idx","tech_rep"),
                                values_from = species_copies,
                                names_from = sp_idx) 

MFU_dat_mock_copy_log <- MFU_dat_mock_copy
MFU_dat_mock_copy_log[,4:ncol(MFU_dat_mock_copy_log)] <- log(MFU_dat_mock_copy_log[,4:ncol(MFU_dat_mock_copy_log)] + 1e-10)

MV1_dat_mock_copy <-pivot_wider(dat.mock %>% filter(marker=="MV1"),id_cols=c("marker","site_samp_idx","tech_rep"),
                                values_from = species_copies,
                                names_from = sp_idx)
MV1_dat_mock_copy_log <- MV1_dat_mock_copy
MV1_dat_mock_copy_log[,4:ncol(MV1_dat_mock_copy_log)] <- log(MV1_dat_mock_copy_log[,4:ncol(MV1_dat_mock_copy_log)] + 1e-10)

#### PULL IN AMPLIFICATION EFFICIENCIES AND PRETEND THEY ARE KNOWN FOR NOW.
# Define the reference species for each sample based on the most commonly observed species in each sample and marker

# UNKNOWN SAMPLES.
ref_dat <-  dat.long.h %>% group_by(marker_idx,site_samp_idx,sp_idx) %>% summarise(Sum = sum(value)) %>% 
  ungroup() %>% group_by(marker_idx,site_samp_idx) %>% mutate(Max=max(Sum)) %>% filter(Sum==Max)

ref_sp_mb1 <- ref_dat %>% filter(marker_idx==1) %>% left_join(MFU_dat,.) %>% 
  dplyr::select(marker_idx,site_samp_idx,sp_idx)
ref_sp_mb2 <- ref_dat %>% filter(marker_idx==2) %>% left_join(MV1_dat,.) %>% 
  dplyr::select(marker_idx,site_samp_idx,sp_idx)

# MOCKS SAMPLES.
ref_dat <-  dat.mock %>% group_by(marker_idx,site_samp_idx,sp_idx) %>% summarise(Sum = sum(value)) %>% 
  ungroup() %>% group_by(marker_idx,site_samp_idx) %>% mutate(Max=max(Sum)) %>% filter(Sum==Max)

ref_sp_mb1_mock <- ref_dat %>% filter(marker_idx==1) %>% left_join(MFU_mock,.) %>% 
  dplyr::select(marker_idx,site_samp_idx,sp_idx)
ref_sp_mb2_mock <- ref_dat %>% filter(marker_idx==2) %>% left_join(MV1_mock,.) %>% 
  dplyr::select(marker_idx,site_samp_idx,sp_idx)

# Choose the final, ref species (Anchovy) to be alpha == 0 for marker 1 and marker 2 
# These are treated as known and so the choice is somewhat arbitrary.

# Summarize observe reads for each replicate.
Reads <- dat.long.h %>% group_by(site_idx,marker_idx,site_samp_idx,tech_rep) %>% summarise(N_read = sum(value))
dat.long.h <- dat.long.h %>% left_join(.,Reads)

reads_mb1 <- MFU_dat %>% dplyr::select(-marker,-site_samp_idx,-tech_rep) %>% rowSums()
reads_mb2 <- MV1_dat %>% dplyr::select(-marker,-site_samp_idx,-tech_rep) %>% rowSums()

# Calculate combinatorial coefficient. for unknowns and mocks
log_mb1_n_k <- MFU_dat %>% dplyr::select(-marker,-site_samp_idx,-tech_rep) %>% as.matrix()
log_mb2_n_k <- MV1_dat %>% dplyr::select(-marker,-site_samp_idx,-tech_rep) %>% as.matrix()
for(i in 1:nrow(MFU_dat)){
  tmp <-(MFU_dat %>% dplyr::select(-marker,-site_samp_idx,-tech_rep))[i,] %>% unlist() %>% c()
  log_mb1_n_k[i,] <- as.numeric(lchoose(reads_mb1[i],tmp))
}
for(i in 1:nrow(MV1_dat)){
  tmp <-(MV1_dat %>% dplyr::select(-marker,-site_samp_idx,-tech_rep))[i,] %>% unlist() %>% c()
  log_mb2_n_k[i,] <- as.numeric(lchoose(reads_mb2[i],tmp))
}
# MOCKS
log_mb1_n_k_mock <- MFU_mock %>% dplyr::select(-marker,-site_samp_idx,-tech_rep) %>% as.matrix()
log_mb2_n_k_mock <- MV1_mock %>% dplyr::select(-marker,-site_samp_idx,-tech_rep) %>% as.matrix()
for(i in 1:nrow(MFU_mock)){
  tmp <-(MFU_mock %>% dplyr::select(-marker,-site_samp_idx,-tech_rep))[i,] %>% unlist() %>% c()
  log_mb1_n_k_mock[i,] <- as.numeric(lchoose(reads_mb1_mock[i],tmp))
}
for(i in 1:nrow(MV1_mock)){
  tmp <-(MV1_mock %>% dplyr::select(-marker,-site_samp_idx,-tech_rep))[i,] %>% unlist() %>% c()
  log_mb2_n_k_mock[i,] <- as.numeric(lchoose(reads_mb2_mock[i],tmp))
}

### Make Indicator (0 if count>0, 1 if count==0)
MB1_bias_ind <- MFU_dat %>% dplyr::select(-marker,-site_samp_idx,-tech_rep)
MB1_bias_ind[MB1_bias_ind>0] <- 99
MB1_bias_ind[MB1_bias_ind==0] <- 1
MB1_bias_ind[MB1_bias_ind==99] <- 0

MB2_bias_ind <- MV1_dat %>% dplyr::select(-marker,-site_samp_idx,-tech_rep)
MB2_bias_ind[MB2_bias_ind>0] <- 99
MB2_bias_ind[MB2_bias_ind==0] <- 1
MB2_bias_ind[MB2_bias_ind==99] <- 0

############################################################
# N_pcr
############################################################
N_pcr = 43 # 

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
  N_obs_mb2 = N_obs_mb2,
  
  # Counters for Mocks
  N_site_mock = N_site_mock,
  N_site_samp_mock = N_site_samp_mock,
  N_obs_mb1_mock = N_obs_mb1_mock,
  N_obs_mb2_mock = N_obs_mb2_mock,
  
  # Mock communities
  count_mb1_mock = MFU_mock %>% dplyr::select(-marker,-site_samp_idx,-tech_rep),
  count_mb2_mock = MV1_mock %>% dplyr::select(-marker,-site_samp_idx,-tech_rep),
  
  log_copy_mock_mb1 = MFU_dat_mock_copy_log %>% dplyr::select(-marker,-site_samp_idx,-tech_rep), #known copies in mock samples
  log_copy_mock_mb2 = MV1_dat_mock_copy_log %>% dplyr::select(-marker,-site_samp_idx,-tech_rep), #known copies in mock samples
  
  log_Lambda_mock = MFU_dat_mock_copy %>% dplyr::select(-marker,-site_samp_idx,-tech_rep) %>% rowSums() %>% log(),
  
  reads_mb1_mock = reads_mb1_mock,
  reads_mb2_mock = reads_mb2_mock,
  
  log_mb1_n_k_mock = log_mb1_n_k_mock, # n choose k coefficient mb1
  log_mb2_n_k_mock = log_mb2_n_k_mock, # n choose k coefficient mb2
  
  # Observations of read counts  (Unknown Samples)
  count_mb1 = MFU_dat %>% dplyr::select(-marker,-site_samp_idx,-tech_rep),
  count_mb2 = MV1_dat %>% dplyr::select(-marker,-site_samp_idx,-tech_rep),
  
  reads_mb1 = reads_mb1,
  reads_mb2 = reads_mb2,
  
  log_mb1_n_k = log_mb1_n_k, # n choose k coefficient mb1
  log_mb2_n_k = log_mb2_n_k, # n choose k coefficient mb2
  
  # design matrices for mapping samples to their replicates.
  X_MB1_ss = X_MB1_ss,
  X_MB2_ss = X_MB2_ss,
  
  X_D_to_F = X_D_to_F,
  
  # Indicator files for finding non-zero observations for each species
  # and making a total concentration estimate 
  MB1_ind = MB1_ind %>% ungroup() %>% dplyr::select(-site_samp_idx),
  MB2_ind = MB2_ind %>% ungroup() %>% dplyr::select(-site_samp_idx),
  
  MB1_bias_ind = MB1_bias_ind,
  MB2_bias_ind = MB2_bias_ind,
  
  site_samp_mb1_idx = ref_sp_mb1$site_samp_idx,
  site_samp_mb2_idx = ref_sp_mb2$site_samp_idx,
  
  # Amplification efficiency for each marker
  # amp_alpha_mb1 = amp_alpha_mb1,
  # amp_alpha_mb2 = amp_alpha_mb2,
  
  ref_sp_mb1 = ref_sp_mb1$sp_idx,
  ref_sp_mb2 = ref_sp_mb2$sp_idx,
  
  ref_sp_mb1_mock = ref_sp_mb1_mock$sp_idx,
  ref_sp_mb2_mock = ref_sp_mb2_mock$sp_idx
  
  # Priors
  # prior_log_D_mu = prior_log_D_mu,
  # prior_log_D_sig = prior_log_D_sig,
  # prior_log_beta_phi0_mu = prior_log_beta_phi0_mu,
  # prior_log_beta_phi0_sig = prior_log_beta_phi0_sig,
  # prior_log_phi_int_mu = prior_log_phi_int_mu,
  # prior_log_phi_int_sig = prior_log_phi_int_sig,
  # prior_log_beta_phi1_mu = prior_log_beta_phi1_mu,
  # prior_log_beta_phi1_sig = prior_log_beta_phi1_sig,
  # prior_amp_tau = prior_amp_tau # SD of variability among amp efficiencies
)
#################

# data preparation for NIMBLE model,
# Assumes you've already run through prepping a "stan_data" list as in the 
# `01_run_stan.R` script.
# Format for NIMBLE -------------------------------------------------------
attach(stan_data)
# Set concentration thresholds (these were missing from your Stan data)
nimble_data <- list(
  # Counters
  # Model parameters
  N_pcr = 43,
  
  # Observations of read counts  
  count_mb1 = MFU_dat %>% dplyr::select(-marker,-site_samp_idx,-tech_rep),
  count_mb2 = MV1_dat %>% dplyr::select(-marker,-site_samp_idx,-tech_rep),
  # Mock communities
  count_mock_mb1 = MFU_mock %>% dplyr::select(-marker,-site_samp_idx,-tech_rep),
  count_mock_mb2 = MV1_mock %>% dplyr::select(-marker,-site_samp_idx,-tech_rep),
  
  reads_mb1 = reads_mb1,
  reads_mb2 = reads_mb2,
  
  # Combinatorial coefficients
  log_mb1_n_k = log_mb1_n_k,
  log_mb2_n_k = log_mb2_n_k,
  
  reads_mock_mb1 = reads_mb1_mock,
  reads_mock_mb2 = reads_mb2_mock,
  
  log_mock_mb1_n_k = log_mb1_n_k_mock, # n choose k coefficient mb1
  log_mock_mb2_n_k = log_mb2_n_k_mock, # n choose k coefficient mb2
  
  # Indicators for species detection by marker
  MB1_ind = as.matrix(MB1_ind %>% ungroup() %>% dplyr::select(-site_samp_idx)),
  MB2_ind = as.matrix(MB2_ind %>% ungroup() %>% dplyr::select(-site_samp_idx)),
  # Design matrices for mapping samples to their replicates
  X_MB1_ss = X_MB1_ss,
  X_MB2_ss = X_MB2_ss,
  X_D_to_F = matrix(X_D_to_F,nrow=N_site_samp,ncol=N_site),
  # X_MB1_ss_mock = X_MFU_ss_mock,
  # X_MB2_ss_mock = X_MV1_ss_mock,
  # X_D_to_F_mock = X_D_to_F_mock,
  log_D_mock = (MFU_dat_mock_copy_log %>% dplyr::select(-marker,-site_samp_idx,-tech_rep))[c(1,4,7),]
);any(is.null(nimble_data))
# Constants for NIMBLE (design matrices and fixed parameters)
nimble_constants <- list(
  
  # Index vectors to avoid dynamic indexing
  # N_marker = as.numeric(N_marker),
  N_site = as.numeric(N_site),
  N_sp = as.numeric(N_sp),
  N_site_samp = as.numeric(N_site_samp),
  N_obs_mb1 = as.numeric(N_obs_mb1),
  N_obs_mb2 = as.numeric(N_obs_mb2),
  # Counters for Mocks
  N_obs_mock_mb1 = N_obs_mb1_mock,
  N_obs_mock_mb2 = N_obs_mb2_mock,
  N_site_mock = N_site_mock,
  # N_site_samp_mock = N_site_samp_mock,
  # Mapping indices
  # site_map = site_samp$site_idx,
  samp_map1 = as.numeric(ref_sp_mb1$site_samp_idx),
  samp_map2 = as.numeric(ref_sp_mb2$site_samp_idx),
  ref_sp_mb1 = as.numeric(ref_sp_mb1$sp_idx),
  ref_sp_mb2 = as.numeric(ref_sp_mb2$sp_idx),
  ref_sp_mock_mb1 = as.numeric(ref_sp_mb1_mock$sp_idx),
  ref_sp_mock_mb2 = as.numeric(ref_sp_mb2_mock$sp_idx)
  # Vector of ones (utility)
  # vec_one_sp = rep(1, N_sp)
)
any(is.null(nimble_constants))


# Pre-compute mock community deterministic calculations
# These were previously calculated inside the NIMBLE model but can be pre-computed
# since log_D_mock is known data

# Extract log_D_mock for easier reference
log_D_mock <- nimble_data$log_D_mock # Add small constant to avoid log(0)

# Mock community matrix operations (pre-computed)
log_F_mock <- as.matrix(X_D_to_F_mock %*% as.matrix(log_D_mock))
log_F_mock_mb1 <- as.matrix(X_MFU_ss_mock) %*% log_F_mock
log_F_mock_mb2 <- as.matrix(X_MV1_ss_mock) %*% log_F_mock

# Pre-compute exponential transformations
temp_exp_mock_mb1 <- exp(log_F_mock_mb1)
temp_exp_mock_mb2 <- exp(log_F_mock_mb2)

# Pre-compute log(1 - exp(-x)) terms
log_1mexp_mock_mb1 <- log1p(- exp(-temp_exp_mock_mb1))
log_1mexp_mock_mb2 <- log1p(- exp(-temp_exp_mock_mb2))

# Pre-compute detection probabilities
log_p_pos_mock_mb1 <- log_1mexp_mock_mb1
log_p_pos_mock_mb2 <- log_1mexp_mock_mb2
log_p_zero_mock_mb1 <- -temp_exp_mock_mb1
log_p_zero_mock_mb2 <- -temp_exp_mock_mb2

# Pre-compute K concentrations
log_K_mock_mb1 <- log_F_mock_mb1 - ifelse(log_1mexp_mock_mb1>-Inf,log_1mexp_mock_mb1,-100)
log_K_mock_mb2 <- log_F_mock_mb2 - ifelse(log_1mexp_mock_mb2>-Inf,log_1mexp_mock_mb1,-100)

# Pre-compute reference species log_K values for each observation
ref_log_K_mock_mb1 <- numeric(nrow(log_K_mock_mb1))
ref_log_K_mock_mb2 <- numeric(nrow(log_K_mock_mb2))

for(ii in 1:nrow(log_K_mock_mb1)) {
  ref_log_K_mock_mb1[ii] <- log_K_mock_mb1[ii, nimble_constants$ref_sp_mock_mb1[ii]]
}

for(ii in 1:nrow(log_K_mock_mb2)) {
  ref_log_K_mock_mb2[ii] <- log_K_mock_mb2[ii, nimble_constants$ref_sp_mock_mb2[ii]]
}

# Pre-compute ALR (additive log-ratio) transformations
alr_mock_mb1 <- matrix(0, nrow = nrow(log_K_mock_mb1), ncol = ncol(log_K_mock_mb1))
alr_mock_mb2 <- matrix(0, nrow = nrow(log_K_mock_mb2), ncol = ncol(log_K_mock_mb2))

for(ii in 1:nrow(log_K_mock_mb1)) {
  alr_mock_mb1[ii, ] <- log_K_mock_mb1[ii, ] - ref_log_K_mock_mb1[ii]
}

for(ii in 1:nrow(log_K_mock_mb2)) {
  alr_mock_mb2[ii, ] <- log_K_mock_mb2[ii, ] - ref_log_K_mock_mb2[ii]
}

# Add these pre-computed values to nimble_data
nimble_data$log_F_mock_mb1 <- log_F_mock_mb1
nimble_data$log_F_mock_mb2 <- log_F_mock_mb2
# nimble_data$temp_exp_mock_mb1 <- temp_exp_mock_mb1
# nimble_data$temp_exp_mock_mb2 <- temp_exp_mock_mb2
nimble_data$log_1mexp_mock_mb1 <- log_1mexp_mock_mb1
nimble_data$log_1mexp_mock_mb2 <- log_1mexp_mock_mb2
nimble_data$log_p_pos_mock_mb1 <- log_p_pos_mock_mb1
nimble_data$log_p_pos_mock_mb2 <- log_p_pos_mock_mb2
nimble_data$log_p_zero_mock_mb1 <- log_p_zero_mock_mb1
nimble_data$log_p_zero_mock_mb2 <- log_p_zero_mock_mb2
nimble_data$log_K_mock_mb1 <- log_K_mock_mb1
nimble_data$log_K_mock_mb2 <- log_K_mock_mb2
# nimble_data$ref_log_K_mock_mb1 <- ref_log_K_mock_mb1
# nimble_data$ref_log_K_mock_mb2 <- ref_log_K_mock_mb2
nimble_data$alr_mock_mb1 <- alr_mock_mb1
nimble_data$alr_mock_mb2 <- alr_mock_mb2
nimble_constants$mock_map <- MFU_mock$site_samp_idx


# Remove log_D_mock from nimble_data since it's no longer needed in the model
# nimble_data$log_D_mock <- NULL
any(is.null(nimble_data))


# Pre-compute site-level values for mock communities from existing data
log_Lambda_site_mock_mb1 <- numeric(N_site_mock)
log_Lambda_site_mock_mb2 <- numeric(N_site_mock) 
log_K_site_mock_mb1 <- matrix(0, nrow = N_site_mock, ncol = N_sp)
log_K_site_mock_mb2 <- matrix(0, nrow = N_site_mock, ncol = N_sp)

# Extract unique site values for mocks from already computed values
for(m in 1:N_site_mock) {
  # Find first observation for this mock site
  obs_idx_mb1 <- which(MFU_mock$site_samp_idx == m)[1]
  obs_idx_mb2 <- which(MV1_mock$site_samp_idx == m)[1]
  
  # Use the already computed values from your existing preprocessing
  if (!is.na(obs_idx_mb1)) {
    # Compute lambda values for this observation
    lambda_vals_mb1 <- exp(log_F_mock_mb1[obs_idx_mb1, ])
    log_Lambda_site_mock_mb1[m] <- log(sum(lambda_vals_mb1))
    log_K_site_mock_mb1[m, ] <- log_K_mock_mb1[obs_idx_mb1, ]
  }
  
  if (!is.na(obs_idx_mb2)) {
    # Compute lambda values for this observation  
    lambda_vals_mb2 <- exp(log_F_mock_mb2[obs_idx_mb2, ])
    log_Lambda_site_mock_mb2[m] <- log(sum(lambda_vals_mb2))
    log_K_site_mock_mb2[m, ] <- log_K_mock_mb2[obs_idx_mb2, ]
  }
}


# Pre-compute mock community deterministic calculations
# These were previously calculated inside the NIMBLE model but can be pre-computed
# since log_D_mock is known data

# Extract log_D_mock for easier reference
log_D_mock <- nimble_data$log_D_mock # Add small constant to avoid log(0)

# Mock community matrix operations (pre-computed)
log_F_mock <- as.matrix(X_D_to_F_mock %*% as.matrix(log_D_mock))
log_F_mock_mb1 <- as.matrix(X_MFU_ss_mock) %*% log_F_mock
log_F_mock_mb2 <- as.matrix(X_MV1_ss_mock) %*% log_F_mock

# Pre-compute exponential transformations
temp_exp_mock_mb1 <- exp(log_F_mock_mb1)
temp_exp_mock_mb2 <- exp(log_F_mock_mb2)

# Pre-compute log(1 - exp(-x)) terms
log_1mexp_mock_mb1 <- log1p(- exp(-temp_exp_mock_mb1))
log_1mexp_mock_mb2 <- log1p(- exp(-temp_exp_mock_mb2))

# Pre-compute detection probabilities
log_p_pos_mock_mb1 <- log_1mexp_mock_mb1
log_p_pos_mock_mb2 <- log_1mexp_mock_mb2
log_p_zero_mock_mb1 <- -temp_exp_mock_mb1
log_p_zero_mock_mb2 <- -temp_exp_mock_mb2

# Pre-compute K concentrations
log_K_mock_mb1 <- log_F_mock_mb1 - ifelse(log_1mexp_mock_mb1>-Inf,log_1mexp_mock_mb1,-100)
log_K_mock_mb2 <- log_F_mock_mb2 - ifelse(log_1mexp_mock_mb2>-Inf,log_1mexp_mock_mb1,-100)

# Pre-compute reference species log_K values for each observation
ref_log_K_mock_mb1 <- numeric(nrow(log_K_mock_mb1))
ref_log_K_mock_mb2 <- numeric(nrow(log_K_mock_mb2))

for(ii in 1:nrow(log_K_mock_mb1)) {
  ref_log_K_mock_mb1[ii] <- log_K_mock_mb1[ii, nimble_constants$ref_sp_mock_mb1[ii]]
}

for(ii in 1:nrow(log_K_mock_mb2)) {
  ref_log_K_mock_mb2[ii] <- log_K_mock_mb2[ii, nimble_constants$ref_sp_mock_mb2[ii]]
}

# Pre-compute ALR (additive log-ratio) transformations
alr_mock_mb1 <- matrix(0, nrow = nrow(log_K_mock_mb1), ncol = ncol(log_K_mock_mb1))
alr_mock_mb2 <- matrix(0, nrow = nrow(log_K_mock_mb2), ncol = ncol(log_K_mock_mb2))

for(ii in 1:nrow(log_K_mock_mb1)) {
  alr_mock_mb1[ii, ] <- log_K_mock_mb1[ii, ] - ref_log_K_mock_mb1[ii]
}

for(ii in 1:nrow(log_K_mock_mb2)) {
  alr_mock_mb2[ii, ] <- log_K_mock_mb2[ii, ] - ref_log_K_mock_mb2[ii]
}

# Add these pre-computed values to nimble_data
nimble_data$log_F_mock_mb1 <- log_F_mock_mb1
nimble_data$log_F_mock_mb2 <- log_F_mock_mb2
# nimble_data$temp_exp_mock_mb1 <- temp_exp_mock_mb1
# nimble_data$temp_exp_mock_mb2 <- temp_exp_mock_mb2
nimble_data$log_1mexp_mock_mb1 <- log_1mexp_mock_mb1
nimble_data$log_1mexp_mock_mb2 <- log_1mexp_mock_mb2
nimble_data$log_p_pos_mock_mb1 <- log_p_pos_mock_mb1
nimble_data$log_p_pos_mock_mb2 <- log_p_pos_mock_mb2
nimble_data$log_p_zero_mock_mb1 <- log_p_zero_mock_mb1
nimble_data$log_p_zero_mock_mb2 <- log_p_zero_mock_mb2
nimble_data$log_K_mock_mb1 <- log_K_mock_mb1
nimble_data$log_K_mock_mb2 <- log_K_mock_mb2
# nimble_data$ref_log_K_mock_mb1 <- ref_log_K_mock_mb1
# nimble_data$ref_log_K_mock_mb2 <- ref_log_K_mock_mb2
nimble_data$alr_mock_mb1 <- alr_mock_mb1
nimble_data$alr_mock_mb2 <- alr_mock_mb2
nimble_constants$mock_map <- MFU_mock$site_samp_idx
# ADD TO nimble_data:
nimble_data$log_Lambda_site_mock_mb1 <- log_Lambda_site_mock_mb1
nimble_data$log_Lambda_site_mock_mb2 <- log_Lambda_site_mock_mb2  
nimble_data$log_K_site_mock_mb1 <- log_K_site_mock_mb1
nimble_data$log_K_site_mock_mb2 <- log_K_site_mock_mb2

nimble_data$MB1_all =  as.matrix(MB1_ind %>% ungroup() %>% dplyr::select(-site_samp_idx)) %>% apply(.,2, function(col) as.numeric(any(col > 0))) 
nimble_data$MB2_all =  as.matrix(MB2_ind %>% ungroup() %>% dplyr::select(-site_samp_idx)) %>% apply(.,2, function(col) as.numeric(any(col > 0)))
nimble_data$MB1_mock_all =  as.matrix(MB1_mock %>% ungroup() %>% dplyr::select(-site_samp_idx)) %>% apply(.,2, function(col) as.numeric(any(col > 0))) 
nimble_data$MB2_mock_all =  as.matrix(MB2_mock %>% ungroup() %>% dplyr::select(-site_samp_idx)) %>% apply(.,2, function(col) as.numeric(any(col > 0)))

# add sequencing run indicators

# Define number of sequencing runs
N_run <- 3  # Run 1: mock sites ; Run 2: sites 1:3; Run 3: sites 4+
# To use 1:10 dilutions as "known" communities (i.e., traditional mocks), set 1:N_site <= 6 below
# Create sequencing run indicators for unknown samples
seq_run_mb1_sites <- ifelse(1:N_site <= 3, 2, 3)
seq_run_mb2_sites <- ifelse(1:N_site <= 3, 2, 3)
# For mock samples - all mocks were in run 1
seq_run_mock_mb1 <- rep(1, N_site_mock)  # All mocks in run 1
seq_run_mock_mb2 <- rep(1, N_site_mock)  # All mocks in run 1

# Add to your nimble_constants list
nimble_constants$N_run <- N_run

# Add to your nimble_data list
# nimble_constants$seq_run_mb1 <- seq_run_mb1
# nimble_constants$seq_run_mb2 <- seq_run_mb2
nimble_constants$seq_run_mb1_sites <- seq_run_mb1_sites
nimble_constants$seq_run_mb2_sites <- seq_run_mb2_sites
nimble_constants$seq_run_mock_mb1 <- seq_run_mock_mb1
nimble_constants$seq_run_mock_mb2 <- seq_run_mock_mb2
