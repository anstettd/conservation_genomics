##################################################################################
## Join SNP table with environment-snp associated data
## Author Daniel Anstett & Julia Anstett
## 
##
## Last Modified July 23, 2021

### OVERALL WORKFLOW:
## Assumes you have:
# A file with SNPS with a BF >10
# SNP Table and loci labels
## Produces: 
#SNP table (SNPs with BF >10) with all genotype abundance information
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)

# Import population, SNP, and loci data (BayPass input files)
# Files are large. Download files from : https://www.dropbox.com/sh/l6zm0av77fxsndg/AABYluDaUZE_OY1U32FqSgwUa?dl=0
# Place files in directory and substitute "Users/daniel_anstett/Dropbox/AM_Workshop/trim/" with your directory path
#pop_order<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.pop_order", 
#                      header=F, sep="\t")
snp<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table", 
                      header=F, sep=" ")
loci<-read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/baseline_filtered_variants.QUAL20_MQ40_AN80_MAF0.03_DP1SD.Baypass_table.loci", 
                 header=F, sep="\t")
colnames(loci) <- c("Chromosome","SNP") # lable columns

loci_united <- loci %>% unite(chr_snp,"Chromosome","SNP",sep="_")
loci_snp <-cbind(loci_united,snp) #add SNP labels to rows

#Import environment-SNP associations
env1 <- read_csv("Data/env1_adapt.csv")
env2 <- read_csv("Data/env2_adapt.csv")
env3 <- read_csv("Data/env3_adapt.csv")
env4 <- read_csv("Data/env4_adapt.csv")
env5 <- read_csv("Data/env5_adapt.csv")
env6 <- read_csv("Data/env6_adapt.csv")
env7 <- read_csv("Data/env7_adapt.csv")
env8 <- read_csv("Data/env8_adapt.csv")
env9 <- read_csv("Data/env9_adapt.csv")

#left_join snp table to each environment-snp dataframe to get associated loci data for each population
env1_loci <- left_join(env1,loci_snp, by="chr_snp")
env2_loci <- left_join(env2,loci_snp, by="chr_snp")
env3_loci <- left_join(env3,loci_snp, by="chr_snp")
env4_loci <- left_join(env4,loci_snp, by="chr_snp")
env5_loci <- left_join(env5,loci_snp, by="chr_snp")
env6_loci <- left_join(env6,loci_snp, by="chr_snp")
env7_loci <- left_join(env7,loci_snp, by="chr_snp")
env8_loci <- left_join(env8,loci_snp, by="chr_snp")
env9_loci <- left_join(env9,loci_snp, by="chr_snp")

#Write out snp table
write_csv(env1_loci, "Data/env1_loci.csv",col_names=TRUE)
write_csv(env2_loci, "Data/env2_loci.csv",col_names=TRUE)
write_csv(env3_loci, "Data/env3_loci.csv",col_names=TRUE)
write_csv(env4_loci, "Data/env4_loci.csv",col_names=TRUE)
write_csv(env5_loci, "Data/env5_loci.csv",col_names=TRUE)
write_csv(env6_loci, "Data/env6_loci.csv",col_names=TRUE)
write_csv(env7_loci, "Data/env7_loci.csv",col_names=TRUE)
write_csv(env8_loci, "Data/env8_loci.csv",col_names=TRUE)
write_csv(env9_loci, "Data/env9_loci.csv",col_names=TRUE)













