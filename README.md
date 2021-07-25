# conservation_genomics
Documentation for analysis pipeline for conservation genomics workshop using Mimulus cardinalis

This repo outlines a genomics-based and a climate only strategy for selecting donor sites for assited migration.
Mimulus cardinalis (scarlet monkeyflower) is used a system to illustrate these two methods. We leaverage a whole genome sequening dataset for 402 popultions
across 55 sites throughout California and Oregon, as well as present and future climate data from Climate NA.


genomics_informed outlies genomic-based strategies that involve analysis of the output of BayPass associations between climate data and genomic SNP data.
If you are interested in just the most relevant illustrated results please see the R markdown file: Genomics_donor_selection.html_

If you would like to follow all the computational steps from snp tables and BayPass output through to all the final results see all five scirpts.
"1_manhattan_plots.R"  Shows manhhatan plots for three climate variables, and can be modified to show upto 9 climate variables
"2_filterSNPs.R"  Documents initital filtering of SNPs with a strong evidence of assocaition with climate (BayesFactor >10).
"3_SNP_env_table.R"  Shows how the SNP table for selected SNPS (BF>10) is joined to the the SNP table that has abudance information for both biallelic genotypes
"4_frequency_binary_tables.R"  Generates tables needed to carry out donor site selection. 
  This script  
  Tabulates proportion of alleles associated with each environmental variable
  Generate frequency tables and binary tables containing population level proportions and binary snp info for alleles associated with adaptation to climate change
"5_genomic_donor_selection.R" Provides maps and figures highlighting possible donor sites per SNP, for all SNPs and for SNPs missing from target site.


climate_only outlies climate only analyses that use climate layers and the coordinates of known locations of the target species.
If you are interested in just the final illustrated results please see: 3._find_donors.R

If you would like to follow all the computational steps from initial raw cliamte layer please see:
"1._climate_layer_processing.R"  Cuts down the size of cliamte layer to better fit a github repo.
"2._masking_climate.R"  Further constrians cliamte layers by cropping to approximate range extent of the species
"3._find_donors.R"  Constrains climate layers to match future climate change conditions for a target site illustrating an area of conservaiton concern. 
Known populations are overlaid to show which might be adapted to future cliamte at the target site.


"Data" provides all the nessary data. Links are given for downloads of large files. 

"Graphs" provides all the graphical outputs.

"Genomics_donor_selection.Rmd" is the code underlying the genomics R markdown file.
