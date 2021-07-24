##################################################################################
## Select possible donor sites using genomics
## Make maps displaying genomic variation across landsacpe
## Author Daniel Anstett 
##
##
## Last Modified July 23, 2021

### OVERALL WORKFLOW:
## Assumes you have:
# Tables for the frequency and prescence/abscence of all loci association with climate change
# Coordinates of all populations
## Produces:
# Maps and figures highlighting possible donor sites per SNP, for all SNPs and for SNPs missing from target site

###################################################################################


###################################################################################
#Import libraries
library(tidyverse)
library(cowplot)
library(sf)
library(tmap)
library(rnaturalearth)
library(rnaturalearthdata)
#devtools::install_github("ropensci/rnaturalearthhires")
library(rnaturalearthhires)
library(rgeos)
library(geodist)
library(RColorBrewer)
library(VennDiagram)
library(viridis)

###################################################################################
#Import population data to calculate distance between sites
pop_var_raw <- read_csv("Data/paper_ID_site_select.csv")
deg_distance <- pop_var_raw %>% dplyr::select(Long,Lat)
geography <- geodist(deg_distance) 
geography_km <- as.data.frame(geography/1000)

#Import abundance results
freq_1 <- read_csv("Data/freq_1.csv")
freq_2 <- read_csv("Data/freq_2.csv")
freq_3 <- read_csv("Data/freq_3.csv")
freq_4 <- read_csv("Data/freq_4.csv")
freq_5 <- read_csv("Data/freq_5.csv")
freq_6 <- read_csv("Data/freq_6.csv")
freq_7 <- read_csv("Data/freq_7.csv")
freq_8 <- read_csv("Data/freq_8.csv")
freq_9 <- read_csv("Data/freq_9.csv")

#Import presence/abscence results
binary_1 <- read_csv("Data/freq_binary_1.csv")
binary_2 <- read_csv("Data/freq_binary_2.csv")
binary_3 <- read_csv("Data/freq_binary_3.csv")
binary_4 <- read_csv("Data/freq_binary_4.csv")
binary_5 <- read_csv("Data/freq_binary_5.csv")
binary_6 <- read_csv("Data/freq_binary_6.csv")
binary_7 <- read_csv("Data/freq_binary_7.csv")
binary_8 <- read_csv("Data/freq_binary_8.csv")
binary_9 <- read_csv("Data/freq_binary_9.csv")

# Get Proportion of SNPs Present
psp_1 <- as.data.frame(colMeans(binary_1[5:59],na.rm = TRUE))
psp_2 <- as.data.frame(colMeans(binary_2[5:59],na.rm = TRUE))
psp_3 <- as.data.frame(colMeans(binary_3[5:59],na.rm = TRUE))
psp_4 <- as.data.frame(colMeans(binary_4[5:59],na.rm = TRUE))
psp_5 <- as.data.frame(colMeans(binary_5[5:59],na.rm = TRUE))
psp_6 <- as.data.frame(colMeans(binary_6[5:59],na.rm = TRUE))
psp_7 <- as.data.frame(colMeans(binary_7[5:59],na.rm = TRUE))
psp_8 <- as.data.frame(colMeans(binary_8[5:59],na.rm = TRUE))
psp_9 <- as.data.frame(colMeans(binary_9[5:59],na.rm = TRUE))

#Put proportion of SNPs Present into population dataframe
pop_var <- cbind(pop_var_raw,psp_1,psp_2,psp_3,psp_4,psp_5,psp_6,psp_7,psp_8,psp_9)
colnames(pop_var)[6:14] <- c("MAT","MAP","PAS","EXT","CMD","Tave_wt","Tave_sm","PPT_wt","PPT_sm") 
###################################################################################
#Mapping setup

# California & Oregon Map Setup
states<-ne_states(country=c("canada","united states of america"),returnclass= "sf")
calo <- states %>%
  filter(name_en=="Oregon" | name_en=="California" | name_en=="Nevada")
#st_crs(calo) in WGS 1984

#Define CRS WGS 1984
EPSG4326<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" #setup WGS 1984 CRS

###################################################################################
###################################################################################
#Venn Diagram of SNP overlap between three climate variables

##### Annual #####
# MAT = Mean annual temperature (°C)
# MAP = Mean annual precipitation (mm)
# CMD = Hargreaves climatic moisture deficit (mm)

MAT_set <- pull(binary_1[,1])
MAP_set <- pull(binary_2[,1])
CMD_set <- pull(binary_5[,1])

ven_3 <- list(MAT_set,MAP_set,CMD_set)

VD_1<-venn.diagram(x=ven_3,
             category.names = c("MAT" , "MAP" , "CMD"),
             fill = c("yellow", "cyan","magenta"),
             #alpha = c(0.5, 0.5, 0.5),
             main.cex = 5,
             cat.cex = 0, cex=1.8,
             fontface = "bold",
             filename = NULL,
)

pdf("Graphs/venn_3.pdf")
grid.newpage()
grid.draw(VD_1)
dev.off()


###################################################################################
#Map distribution of the proportion of individual SNPs

# Map of P1, present in most localities
b1_snp2 <- freq_1 %>% filter(SNP==224717) #select abundance data for snp 2
#b1_snp2 <- binary_1 %>% filter(SNP==224717) #select binary data for snp 2, remove comment to see binary data
b1_snp2 <- b1_snp2[5:59] #take only snps values
b1_snp2 <- as.data.frame(t(b1_snp2)) #transpose and put in dataframe
pop_var_snp2 <- cbind(pop_var_raw,b1_snp2) #bind to lat/long population data
colnames(pop_var_snp2)[6] <- "MAT" 

MAT_points_snp2 <- pop_var_snp2 %>% dplyr::select(Long,Lat,MAT) #select relevant data
MAT_sf_snp2 <- st_as_sf(MAT_points_snp2,coords=c("Long","Lat"), crs=EPSG4326)

#Plot MAT associated allele snp2
tmap_mode("plot")
#tmap_mode("view")
snp2 <- tm_shape(calo)+
  tm_borders()+
  tm_shape(MAT_sf_snp2)+
  tm_bubbles(size = 0.18,col="MAT")+ 
  #tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.005)
snp2
tmap_save(snp2, filename = "Graphs/snp2_freq.pdf",width=5, height=6)

# Map of P25, missing from most localities
b1_snp25 <- freq_1 %>% filter(SNP==946612) #filter for snp
b1_snp25 <- b1_snp25[5:59] #take only snps values
b1_snp25 <- as.data.frame(t(b1_snp25)) #transpose and put in dataframe
pop_var_snp25 <- cbind(pop_var_raw,b1_snp25) #bind to lat/long population data
colnames(pop_var_snp25)[6] <- "MAT" 

MAT_points_snp25 <- pop_var_snp25 %>% dplyr::select(Long,Lat,MAT) #select relevant data
MAT_sf_snp25 <- st_as_sf(MAT_points_snp25,coords=c("Long","Lat"), crs=EPSG4326)

#Plot MAT associated allele snp25
tmap_mode("plot")
#tmap_mode("view")
snp25 <- tm_shape(calo)+
  tm_borders()+
  tm_shape(MAT_sf_snp25)+
  tm_bubbles(size = 0.18,col="MAT")+ 
  tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.005)
tmap_save(snp25, filename = "Graphs/snp25_freq.pdf",width=5, height=6)

###################################################################################

###################################################################################
#Map the presence or absence of individual SNPs 


# Map of P1, present in most localities
b1_snp2 <- binary_1 %>% filter(SNP==224717) #select abundance data for snp 2
#b1_snp2 <- binary_1 %>% filter(SNP==224717) #select binary data for snp 2, remove comment to see binary data
b1_snp2 <- b1_snp2[5:59] #take only snps values
b1_snp2 <- as.data.frame(t(b1_snp2)) #transpose and put in dataframe
pop_var_snp2 <- cbind(pop_var_raw,b1_snp2) #bind to lat/long population data
colnames(pop_var_snp2)[6] <- "MAT" 

MAT_points_snp2 <- pop_var_snp2 %>% dplyr::select(Long,Lat,MAT) #select relevant data
MAT_sf_snp2 <- st_as_sf(MAT_points_snp2,coords=c("Long","Lat"), crs=EPSG4326)

#Plot MAT associated allele snp2
tmap_mode("plot")
#tmap_mode("view")
snp2 <- tm_shape(calo)+
  tm_borders()+
  tm_shape(MAT_sf_snp2)+
  tm_bubbles(size = 0.18,col="MAT")+ 
  #tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.005)
snp2
tmap_save(snp2, filename = "Graphs/snp2_freq_p_a.pdf",width=5, height=6)

# Map of P25, missing from most localities
b1_snp25 <- binary_1 %>% filter(SNP==946612) #filter for snp
b1_snp25 <- b1_snp25[5:59] #take only snps values
b1_snp25 <- as.data.frame(t(b1_snp25)) #transpose and put in dataframe
pop_var_snp25 <- cbind(pop_var_raw,b1_snp25) #bind to lat/long population data
colnames(pop_var_snp25)[6] <- "MAT" 

MAT_points_snp25 <- pop_var_snp25 %>% dplyr::select(Long,Lat,MAT) #select relevant data
MAT_sf_snp25 <- st_as_sf(MAT_points_snp25,coords=c("Long","Lat"), crs=EPSG4326)

#Plot MAT associated allele snp25
tmap_mode("plot")
#tmap_mode("view")
snp25 <- tm_shape(calo)+
  tm_borders()+
  tm_shape(MAT_sf_snp25)+
  tm_bubbles(size = 0.18,col="MAT")+ 
  tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.005)
tmap_save(snp25, filename = "Graphs/snp25_freq_p_a.pdf",width=5, height=6)


###################################################################################
# Map Proportion of putatively adaptive SNPs Present

#Set up sf objects for snp proportions
#MAT
MAT_points <- pop_var %>% dplyr::select(Long,Lat,MAT)
MAT_sf <- st_as_sf(MAT_points,coords=c("Long","Lat"), crs=EPSG4326)
ggplot()+ geom_sf(data = MAT_sf)+ ggtitle("MAT points") #check data is set up properly

#MAP
MAP_points <- pop_var %>% dplyr::select(Long,Lat,MAP)
MAP_sf <- st_as_sf(MAP_points,coords=c("Long","Lat"), crs=EPSG4326)
ggplot()+ geom_sf(data = MAP_sf)+ ggtitle("MAP points") #check data is set up properly

#CMD
CMD_points <- pop_var %>% dplyr::select(Long,Lat,CMD)
CMD_sf <- st_as_sf(CMD_points,coords=c("Long","Lat"), crs=EPSG4326)
ggplot()+ geom_sf(data = CMD_sf)+ ggtitle("CMD points") #check data is set up properly

#Make Maps 
# Note: remove "#" from tmap_mode("view") to switch to the dynamic map viewer mode in R studio.
#MAT
tmap_mode("plot")
#tmap_mode("view")
MAT_all <- tm_shape(calo)+
  tm_borders()+
  tm_shape(MAT_sf)+
  tm_bubbles(size = 0.18,col="MAT")+ 
  #tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
MAT_all
tmap_save(MAT_all, filename = "Graphs/MAT_all.pdf",width=5, height=6)

#MAP
tmap_mode("plot")
#tmap_mode("view")
MAP_all <- tm_shape(calo)+
  tm_borders()+
  tm_shape(MAP_sf)+
  tm_bubbles(size = 0.15,col="MAP")+ 
  #tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
MAP_all
tmap_save(MAP_all, filename = "Graphs/MAP_all.pdf",width=5, height=6)

#CMD
tmap_mode("plot")
#tmap_mode("view")
CMD_all <- tm_shape(calo)+
  tm_borders()+
  tm_shape(CMD_sf)+
  tm_bubbles(size = 0.15,col="CMD")+ 
  #tm_dots(size=0.1,shape=1)+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
CMD_all
tmap_save(CMD_all, filename = "Graphs/CMD_all.pdf",width=5, height=6)



###################################################################################

##plot missing variation for P8 Deep Creek

#Setup sf object for P8
p8_only <- pop_var_raw %>% filter(Paper_ID==8) %>% dplyr::select(Long,Lat) #(Paper_ID can be changed to selection population of choice
p8_sf <- st_as_sf(p8_only,coords=c("Long","Lat"), crs=EPSG4326)

##Setup P8 MAT
binary1_miss_p8 <- binary_1 %>% filter(P8==0) #select only snps not in P8, can be changed to selection population of choice
psp_1_p8 <- as.data.frame(colMeans(binary1_miss_p8[5:59],na.rm = TRUE))
pop_var_p8 <- cbind(pop_var_raw,psp_1_p8)
colnames(pop_var_p8)[6] <- "MAT" 

#Population values for p8
MAT_points_p8 <- pop_var_p8 %>% dplyr::select(Long,Lat,MAT)
MAT_sf_p8 <- st_as_sf(MAT_points_p8,coords=c("Long","Lat"), crs=EPSG4326)

#Make Maps 
# Note: remove "#" from tmap_mode("view") to switch to the dynamic map viewer mode in R studio.
#MAT
tmap_mode("plot")
#tmap_mode("view")
miss_p8_MAT <- tm_shape(calo)+
  tm_borders()+
  tm_shape(MAT_sf_p8)+
  tm_bubbles(size = 0.15,col="MAT")+ 
  tm_shape(p8_sf)+
  tm_dots(size=0.3,shape=20,col= "#33FFFF")+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
miss_p8_MAT
tmap_save(miss_p8_MAT, filename = "Graphs/miss_p8_MAT.pdf",width=5, height=6)


##Setup P8 MAP
binary2_miss_p8 <- binary_2 %>% filter(P8==0) #select only snps not in P8
PSP_2_p8 <- as.data.frame(colMeans(binary2_miss_p8[5:59],na.rm = TRUE))
pop_2var_p8 <- cbind(pop_var_raw,PSP_2_p8)
colnames(pop_2var_p8)[6] <- "MAP" 

#Population values for p8
MAP_points_p8 <- pop_2var_p8 %>% dplyr::select(Long,Lat,MAP)
MAP_sf_p8 <- st_as_sf(MAP_points_p8,coords=c("Long","Lat"), crs=EPSG4326)

#MAP
tmap_mode("plot")
#tmap_mode("view")
miss_p8_MAP <- tm_shape(calo)+
  tm_borders()+
  tm_shape(MAP_sf_p8)+
  tm_bubbles(size = 0.15,col="MAP")+ 
  tm_shape(p8_sf)+
  tm_dots(size=0.3,shape=20,col= "#33FFFF")+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
miss_p8_MAP
tmap_save(miss_p8_MAP, filename = "Graphs/miss_p8_MAP.pdf",width=5, height=6)



##Setup P8 CMD
binary5_miss_p8 <- binary_5 %>% filter(P8==0) #select only snps not in P8
PSP_5_p8 <- as.data.frame(colMeans(binary5_miss_p8[5:59],na.rm = TRUE))
pop_5var_p8 <- cbind(pop_var_raw,PSP_5_p8)
colnames(pop_5var_p8)[6] <- "CMD" 

#Population values for p8
CMD_points_p8 <- pop_5var_p8 %>% dplyr::select(Long,Lat,CMD)
CMD_sf_p8 <- st_as_sf(CMD_points_p8,coords=c("Long","Lat"), crs=EPSG4326)

#CMD
tmap_mode("plot")
#tmap_mode("view")
miss_p8_CMD <- tm_shape(calo)+
  tm_borders()+
  tm_shape(CMD_sf_p8)+
  tm_bubbles(size = 0.15,col="CMD")+ 
  tm_shape(p8_sf)+
  tm_dots(size=0.3,shape=20,col= "#33FFFF")+
  tm_layout(legend.position = c(0.29, 0.73),legend.title.size = 0.001)
miss_p8_CMD
tmap_save(miss_p8_CMD, filename = "Graphs/miss_p8_CMD.pdf",width=5, height=6)



###################################################################################
#Distance to punitively adaptive snp for target site

##MAT
dist_p8_MAT <- data.frame() #Generate empty data frame

#select shortest distance to snp punitively adaptive to MAT 
for (i in 1:nrow(binary_1)){
  geography_p8 <- geography_km[,8] #get distance from p8, can be changed to select population of choice
  b1_t <- as.data.frame(t(binary_1[i,5:59])) #get line i from binary matrix
  geo_b1 <- as.data.frame(cbind(geography_p8,b1_t)) #bind columns together
  colnames(geo_b1) <- c("geo","snp")  #rename
  geo_b1 <-  geo_b1 %>% filter(snp==1) #select for presence of "adaptive" snps
  geo_b1 <-  geo_b1 %>% filter(geo==min(geo)) #select for minimum distance
  dist_p8_MAT[i,1] <- geo_b1[1,1] #write out in distance for i snp
}
colnames(dist_p8_MAT) <- "Distance"

#Make graphic
dist_MAT <- ggplot(data=dist_p8_MAT, aes(dist_p8_MAT$Distance))+
  geom_histogram(binwidth = 20)+
  scale_x_continuous(name="Distance (km)")+
  scale_y_continuous(name="Number of SNPs")+
  theme_classic()
dist_MAT + theme(
  axis.text.x = element_text(size=14,face="bold"),
  axis.text.y = element_text(size=14,face="bold"),
  axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=21,vjust = 2, face="bold",hjust=0.5))
ggsave(filename = "Graphs/dist_MAT.pdf",width=9, height=6)


##MAP
dist_p8_MAP <- data.frame() #Generate empty data frame

#select shortest distance to snp punitively adaptive to MAP 
for (i in 1:nrow(binary_2)){
  geography_p8 <- geography_km[,8] #get distance from p8
  b2_t <- as.data.frame(t(binary_2[i,5:59])) #get line i from binary matrix
  geo_b2 <- as.data.frame(cbind(geography_p8,b2_t)) #bind columns together
  colnames(geo_b2) <- c("geo","snp")  #rename
  geo_b2 <-  geo_b2 %>% filter(snp==1) #select for presence of "adaptive" snps
  geo_b2 <-  geo_b2 %>% filter(geo==min(geo)) #select for minimum distance
  dist_p8_MAP[i,1] <- geo_b2[1,1] #write out in distance for i snp
}
colnames(dist_p8_MAP) <- "Distance"

#Make graphic
dist_MAP <- ggplot(data=dist_p8_MAP, aes(dist_p8_MAP$Distance))+
  geom_histogram(binwidth = 20)+
  scale_x_continuous(name="Distance (km)")+
  scale_y_continuous(name="Number of SNPs")+
  theme_classic()
dist_MAP + theme(
  axis.text.x = element_text(size=14,face="bold"),
  axis.text.y = element_text(size=14,face="bold"),
  axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=21,vjust = 2, face="bold",hjust=0.5))
ggsave(filename = "Graphs/dist_MAP.pdf",width=9, height=6)


##CMD
dist_p8_CMD <- data.frame() #Generate empty data frame

#select shortest distance to snp punitively adaptive to CMD 
for (i in 1:nrow(binary_5)){
  geography_p8 <- geography_km[,8] #get distance from p8
  b5_t <- as.data.frame(t(binary_5[i,5:59])) #get line i from binary matrix
  geo_b5 <- as.data.frame(cbind(geography_p8,b5_t)) #bind columns together
  colnames(geo_b5) <- c("geo","snp")  #rename
  geo_b5 <-  geo_b5 %>% filter(snp==1) #select for presence of "adaptive" snps
  geo_b5 <-  geo_b5 %>% filter(geo==min(geo)) #select for minimum distance
  dist_p8_CMD[i,1] <- geo_b5[1,1] #write out in distance for i snp
}
colnames(dist_p8_CMD) <- "Distance"

#Make graphic
dist_CMD <- ggplot(data=dist_p8_CMD, aes(dist_p8_CMD$Distance))+
  geom_histogram(binwidth = 20)+
  scale_x_continuous(name="Distance (km)")+
  scale_y_continuous(name="Number of SNPs", limits=c(0,3000))+
  theme_classic()
dist_CMD + theme(
  axis.text.x = element_text(size=14,face="bold"),
  axis.text.y = element_text(size=14,face="bold"),
  axis.title.x = element_text(color="black", size=24, vjust = 0.5, face="bold"),
  axis.title.y = element_text(color="black", size=21,vjust = 2, face="bold",hjust=0.5))
ggsave(filename = "Graphs/dist_CMD.pdf",width=9, height=6)

###################################################################################







