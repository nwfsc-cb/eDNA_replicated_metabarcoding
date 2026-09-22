library(here)
library(tidyverse)

# Simulate a single set of concentrations 
source(here("Simulation","sim_functions.R"))

#######################################################################
######################## Simulation 1
###### 50 sp. Reads fixed. 
######################################################################
#####################################################################

N_SP <- sim_scen$N_SP[i]
TOT.CONC <- sim_scen$TOT.CONC[i]
CONC.DIST <- sim_scen$CONC.DIST[i]
N_READS <- sim_scen$N_READS[i]
amp_sd <- sim_scen$amp_sd[i]
replicates <- sim_scen$replicates[i]

conc_dat <- sim_conc(N_species=N_SP,
              conc_dist = "very_high_skew",
              tot_conc = TOT.CONC)

conc_dat$MB1_alpha <- A1$a_val

# OK.  Simulate replicates for each marker
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
