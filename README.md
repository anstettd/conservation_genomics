# conservation_genomics
Documentation for analysis pipeline for conservation genomics workshop using Mimulus cardinalis

This repo outlines a genomics-based and a climate only strategy for selecting donor sites for assisted migration.
Mimulus cardinalis (scarlet monkeyflower) is used as a system to illustrate these two methods. We leverage a whole genome sequencing dataset for 402 populations
across 55 sites throughout California and Oregon, as well as present and future climate data from Climate NA.


The folder "genomics_informed" outlies genomic-based strategies that involve analysis of the output of BayPass associations (Genome-environment Association (GEA)) between climate data and genomic SNP data.
If you are interested in just the most relevant illustrated results please see the R markdown file: Genomics_donor_selection.html_

If you would like to follow all the computational steps from snp tables and BayPass output through to all the final results see all five scripts.
"1_manhattan_plots.R"  Shows manhhatan plots for three climate variables, and can be modified to show upto 9 climate variables
"2_filterSNPs.R"  Documents initial filtering of SNPs with strong evidence of association with climate (BayesFactor >10).
"3_SNP_env_table.R"  Shows how the SNP table for selected SNPS (BF>10) is joined to the the SNP table that has abundance information for both biallelic genotypes.
"4_frequency_binary_tables.R"  Generates tables needed to carry out donor site selection. 
  This script:  
  tabulates proportion of alleles associated with each environmental variable
  generates frequency tables and presence/absence tables containing population level proportions and presence/absence snp info for alleles associated with adaptation to climate change
"5_genomic_donor_selection.R" Provides maps and figures highlighting possible donor sites per SNP, for all SNPs and for SNPs missing from target sites.


climate_only outlies climate only analyses that use climate layers and the coordinates of known locations of the target species.
If you are interested in just the final illustrated results please see: 3._find_donors.R

If you would like to follow all the computational steps from initial raw climate layer please see:
"1._climate_layer_processing.R"  Cuts down the size of the climate layers to better fit a github repo.
"2._masking_climate.R"  Further constrians climate layers by cropping to the approximate range extent of the species.
"3._find_donors.R"  Constrains climate layers to match future climate change conditions for a target site illustrating an area of conservation concern. 
Known populations are overlaid to show which might be adapted to future climate at the target site.


"Data" provides all the necessary data. Links are given for downloads of large files. 

"Graphs" provides all the graphical outputs.

"Genomics_donor_selection.Rmd" is the code underlying the genomics R markdown file.
