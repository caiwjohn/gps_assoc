library(tidyverse)
library(useful)
library(ggfortify)
theme_set(theme_gray(base_size = 18))

####
# SUMMARY
#   1. Data is log_2(x+1) normalized
#   2. Standardized by samples
#   3. Inspect PC but no patterns evident, skeptical of this
####

setwd("~/programming/muchero_lab/gps_assoc/scripts")

#####
# LOADING AND CLEANING
#####
# Read in data
data<- read.csv("../local_data/Leaf_eQTN.csv", header=TRUE, sep=',')
rownames(data)<- data$SampleID
data<- as.data.frame(t(data[,-1]))
data<- cbind(SampleID= rownames(data), data)
data$SampleID<- as.character(data$SampleID)

# Convert to clean data format
dat<- gather(data, key="Transcript", value="Expr", -1,na.rm = TRUE)

# Inspect distri. of Expr values
##test<- sample_n(dat, size=100)
##test<- arrange(test, Expr)
##qplot(seq(1, nrow(test)), Expr, data=test)

# Log normalize expression values
dat$Expr<- log2(dat$Expr+1)

# Spread data and standardize
wide<- spread(dat, key= "Transcript", value = "Expr")
wide[,-1] <- t(scale(t(as.matrix(wide[,-1]))))
dat<- wide

# Match SampleID to GenoID
dat$SampleID<- gsub("\\.", "-", dat$SampleID)

# Load Phenotype data
## WARNING this may be introducing NA's
pheno<- readRDS("../local_data/clean_pheno.RDS") #%>%
        mutate(i = row_number()) %>% 
        spread(phenotype, value) %>%
        select(-i)

# Select garden, year and rep
pheno<- subset(pheno, garden=="CL" & rep==1 & year==2010)

# Remove duplicates
pheno<- pheno[which(!duplicated(pheno$Genotype_ID)),]

# Merge expression and covariates
annotated_dat<- merge(pheno, dat, by.x = "Genotype_ID", by.y= "SampleID")

#####
# Auto-partition groups
#####
# Define size of latitude groups
num_groups<- 4
lat_start<- range(annotated_dat$Latitude)[1]
quar_width<- (range(annotated_dat$Latitude)[2] - range(annotated_dat$Latitude)[1])/num_groups

# Partition
annotated_dat$lat_group[annotated_dat$Latitude<(lat_start+quar_width)]<-1
for(i in 2:(num_groups-1)){
  annotated_dat$lat_group[annotated_dat$Latitude>=(lat_start+((i-1)*quar_width)) & annotated_dat$Latitude<(lat_start+(i*quar_width))]<-i
  }
annotated_dat$lat_group[annotated_dat$Latitude>=(lat_start+((num_groups-1)*quar_width))]<-num_groups
annotated_dat$lat_group<- as.factor(annotated_dat$lat_group)

###
# Visually Partition
###
# Visualize latitudinal groups to manually separate
ggplot(annotated_dat, aes(Latitude, rep(1, nrow(annotated_dat)))) +
  geom_point()+
  ggtitle("Partitioning of Sample Latitudes into Clusters") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text = element_text(size=15)) +
  geom_vline(xintercept = 45.2, color= "green", size= 0.3) +
  geom_vline(xintercept = 47.5, color= "green", size= 0.3) +
  geom_vline(xintercept = 49, color= "green", size= 0.3) +
  geom_vline(xintercept = 50.4, color= "green", size= 0.3) +
  geom_vline(xintercept = 52.4, color= "green", size= 0.3)

# Assign partitions
annotated_dat$lat_group[annotated_dat$Latitude<45.2]<-1
annotated_dat$lat_group[annotated_dat$Latitude>=45.2 & annotated_dat$Latitude<47.5]<-2
annotated_dat$lat_group[annotated_dat$Latitude>=47.5 & annotated_dat$Latitude<49]<-3
annotated_dat$lat_group[annotated_dat$Latitude>=49 & annotated_dat$Latitude<50.4]<-4
annotated_dat$lat_group[annotated_dat$Latitude>=50.4 & annotated_dat$Latitude<52.4]<-5
annotated_dat$lat_group[annotated_dat$Latitude>=52.4]<-6
annotated_dat$lat_group<- as.factor(annotated_dat$lat_group)
  
# Rearrange columns to my liking
annotated_dat<- annotated_dat %>%
  select("Genotype_ID", "Latitude", "lat_group", everything()) %>%
  rename(bud_flush= CL2010_Bud_Flush_rep1, bud_set= CL2010_Bud_Set_score1_rep1)

saveRDS(annotated_dat, file="../data/cln_expr_covar.Rdata")
