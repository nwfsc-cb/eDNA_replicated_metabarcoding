Read_Me for code and data associated with the manuscript, “Deriving Quantitative Signals from Replicated DNA Metabarcoding”

This repository contains all of the necessary code and data to recreate the analyses in the main text and supplements.

The repository has the following main folders.  A read me within each subfolder provides more detailed description of the contents

/Data - Contains empirical data from metabarcoding analyses using the MiFish and MarVer1 DNA primers and information about constructing mock samples.

/NIMBLE Mock Analysis - Contains three scripts for processing the empirical data from mock samples and estimating a statistical model for DNA copies from metabarcoding observations

/NIMBLE Simulation - Contains scripts for generating simulations of metabarcoding data and estimating models for simulation output using Nimble (Model A as described in Supplement S1).  There are separate subfolders and scripts for analyzing data from a single marker (1_marker) and for combining data from two markers (2_marker)

/shared_functions - Contains a few function that are helpful for simulating and processing simulated metabarcoding data.

/Stan Simulation - Contains scripts for generating simulations of metabarcoding data and estimating models for simulation output using Stan (Model B as described in Supplement S1).  There are separate subfolders and scripts for analyzing data from a single marker (1_marker) and for combining data from two markers (2_marker)