##################################################################################
## Make Manhattan plots for SNPS associated with each environmental variable
## Author Daniel Anstett & Julia Anstett
## 
##
## Last Modified July 23, 2021

### OVERALL WORKFLOW:
## Assumes you have: environment-SNP association tables with BayesFactors generated through BayPass
## Produces Manhattan plots 
###################################################################################


###################################################################################
#Import libraries
library(tidyverse)
library(qqman)

#Import chromosome size
chr_size <- read_csv("Data/chr_size.csv")
chr_size[,3] <- cumsum(chr_size$size) #get cumulative chromosome position
colnames(chr_size)[3] <- "poz"


#Import environment-SNP association tables
#Use files provided through external link: https://www.dropbox.com/sh/l6zm0av77fxsndg/AABYluDaUZE_OY1U32FqSgwUa?dl=0
#Set up in a directory and change "/Users/daniel_anstett/Dropbox/AM_Workshop/trim/" to match that directory
env1 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_1_trim.tsv",header=F, sep=" ")
env2 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_2_trim.tsv",header=F, sep=" ")
env3 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_3_trim.tsv",header=F, sep=" ")
env4 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_4_trim.tsv",header=F, sep=" ")
env5 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_5_trim.tsv",header=F, sep=" ")
env6 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_6_trim.tsv",header=F, sep=" ")
env7 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_7_trim.tsv",header=F, sep=" ")
env8 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_8_trim.tsv",header=F, sep=" ")
env9 <- read.table("/Users/daniel_anstett/Dropbox/AM_Workshop/trim/ENV_9_trim.tsv",header=F, sep=" ")

#Name Columns
colnames(env1) <- c("Chromosome","SNP","Env","BF")
colnames(env2) <- c("Chromosome","SNP","Env","BF")
colnames(env3) <- c("Chromosome","SNP","Env","BF")
colnames(env4) <- c("Chromosome","SNP","Env","BF")
colnames(env5) <- c("Chromosome","SNP","Env","BF")
colnames(env6) <- c("Chromosome","SNP","Env","BF")
colnames(env7) <- c("Chromosome","SNP","Env","BF")
colnames(env8) <- c("Chromosome","SNP","Env","BF")
colnames(env9) <- c("Chromosome","SNP","Env","BF")

# Set up data frame used to make the mahattan plot for each environmental variable
env1_united <- env1 %>% unite(chr_snp,Chromosome,SNP) #merge chomosome and snp ID into one combined ID
env1_united <- env1_united %>% select(chr_snp) #select only the combined ID
env1_united <- cbind(env1_united,env1) #bind combined ID back to the main dataframe for env1
#calculates the cumulatve basepair position based on chromsome sizes and SNP IDs
env1_united <- env1_united%>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env2_united <- env2 %>% unite(chr_snp,Chromosome,SNP)
env2_united <- env2_united %>% select(chr_snp)
env2_united <- cbind(env2_united,env2)
env2_united <- env2_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env3_united <- env3 %>% unite(chr_snp,Chromosome,SNP)
env3_united <- env3_united %>% select(chr_snp)
env3_united <- cbind(env3_united,env3)
env3_united <- env3_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env4_united <- env4 %>% unite(chr_snp,Chromosome,SNP)
env4_united <- env4_united %>% select(chr_snp)
env4_united <- cbind(env4_united,env4)
env4_united <- env4_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env5_united <- env5 %>% unite(chr_snp,Chromosome,SNP)
env5_united <- env5_united %>% select(chr_snp)
env5_united <- cbind(env5_united,env5)
env5_united <- env5_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env6_united <- env6 %>% unite(chr_snp,Chromosome,SNP)
env6_united <- env6_united %>% select(chr_snp)
env6_united <- cbind(env6_united,env6)
env6_united <- env6_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env7_united <- env7 %>% unite(chr_snp,Chromosome,SNP)
env7_united <- env7_united %>% select(chr_snp)
env7_united <- cbind(env7_united,env7)
env7_united <- env7_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env8_united <- env8 %>% unite(chr_snp,Chromosome,SNP)
env8_united <- env8_united %>% select(chr_snp)
env8_united <- cbind(env8_united,env8)
env8_united <- env8_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))

env9_united <- env9 %>% unite(chr_snp,Chromosome,SNP)
env9_united <- env9_united %>% select(chr_snp)
env9_united <- cbind(env9_united,env9)
env9_united <- env9_united %>% mutate(CHR=as.numeric(gsub("CE10_chr","",Chromosome))) %>%
  mutate(BP = ifelse(CHR == 1, SNP,
                     ifelse(CHR==2, SNP+chr_size$poz[1],
                            ifelse(CHR==3, SNP+chr_size$poz[2],
                                   ifelse(CHR==4, SNP+chr_size$poz[3],
                                          ifelse(CHR==5, SNP+chr_size$poz[4],
                                                 ifelse(CHR==6, SNP+chr_size$poz[5],
                                                        ifelse(CHR==7, SNP+chr_size$poz[6],
                                                               SNP+chr_size$poz[7]))))))))
###
#Make Manhattan Plots
#Done for example variables:
# env1 = Mean annual temperature (MAT)
# env2 = Mean annual precipitation (MAP)
# env5 = Hargreaves cumulative moisture deficit (CMD)
# !!! Warning will make a few minutes of computational time. The resulting image is rather large.
#ENV 1 - MAT 
axisdf_1 <- env1_united %>% group_by(CHR) %>% summarize(center=( max(BP) + min(BP) ) / 2 )

ggplot(env1_united, aes(x=BP, y=BF)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=0.5) +
  scale_color_manual(values = rep(c("black", "darkgoldenrod"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "red", size=0.9) +
  # custom X axis:
  scale_x_continuous(expand = c(0, 0), label = axisdf_1$CHR, breaks= axisdf_1$center) +
  scale_y_continuous(breaks = c(-20,-10,0,10,20,30,40)) +
  theme_classic() +
  labs(
    y = "Bayes Factor",
    x = "Position")+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=15, face="bold"),
    axis.text.y = element_text(size=15,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5)
  )



#ENV 2 - MAP 
axisdf_2 <- env2_united %>% group_by(CHR) %>% summarize(center=( max(BP) + min(BP) ) / 2 )

ggplot(env2_united, aes(x=BP, y=BF)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=0.5) +
  scale_color_manual(values = rep(c("black", "deepskyblue"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "red", size=0.9) +
  # custom X axis:
  scale_x_continuous(expand = c(0, 0), label = axisdf_2$CHR, breaks= axisdf_2$center) +
  scale_y_continuous(breaks = c(-20,-10,0,10,20,30,40)) +
  theme_classic() +
  labs(
    y = "Bayes Factor",
    x = "Position")+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=15, face="bold"),
    axis.text.y = element_text(size=15,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5)
  )




#ENV 5 - CMD 
axisdf_5 <- env5_united %>% group_by(CHR) %>% summarize(center=( max(BP) + min(BP) ) / 2 )

ggplot(env5_united, aes(x=BP, y=BF)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=0.5) +
  scale_color_manual(values = rep(c("black", "magenta3"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "deepskyblue", size=0.9) +
  # custom X axis:
  scale_x_continuous(expand = c(0, 0), label = axisdf_5$CHR, breaks= axisdf_5$center) +
  scale_y_continuous(breaks = c(-20,-10,0,10,20,30,40)) +
  theme_classic() +
  labs(
    y = "Bayes Factor",
    x = "Position")+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=15, face="bold"),
    axis.text.y = element_text(size=15,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5)
  )




################
#Test plot with only two chromosomes, for changing visualization
env1_united_chr1_2 <- env1_united %>% filter(CHR<3)

axisdf_1 <- env1_united_chr1_2 %>% group_by(CHR) %>% summarize(center=( max(BP) + min(BP) ) / 2 )

ggplot(env1_united_chr1_2, aes(x=BP, y=BF)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.8,size=0.5) +
  scale_color_manual(values = rep(c("black", "magenta2"), 22 )) +
  geom_hline(yintercept=10, linetype="dashed",color = "deepskyblue1", size=0.9) +
  # custom X axis:
  scale_x_continuous(expand = c(0, 0), label = axisdf_1$CHR, breaks= axisdf_1$center) +
  scale_y_continuous(breaks = c(-20,-10,0,10,20,30,40)) +
  theme_classic() +
    labs(
      y = "Bayes Factor",
      x = "Position")+
  theme(
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.8),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.x = element_text(size=15, face="bold"),
    axis.text.y = element_text(size=15,face="bold"),
    axis.title.x = element_text(color="black", size=20, vjust = 0.5, face="bold"),
    axis.title.y = element_text(color="black", size=20,vjust = 2, face="bold",hjust=0.5)
  )

  












