library(tidyverse)

# Load covariates
pheno<- readRDS("../local_data/clean_pheno.RDS")

# Subset to select garden and rep and year
pheno<- subset(pheno, garden=="CL" & rep==1 & year==2010)

# Remove duplicate entries
pheno<- pheno[which(duplicated(pheno)==FALSE),]

# Select bud flush data
bud_flush<- subset(pheno, phenotype=="Bud.Flush") %>%
            select(Genotype_ID, value)

# Select latitude data
lat<- select(pheno, Genotype_ID, Latitude) %>%
      .[which(duplicated(.)==FALSE),]

# Select River data
river<- select(pheno, Genotype_ID, River) %>%
      .[which(duplicated(.)==FALSE),]

# Drop samples with no river label
river<- river[-which(river$River==""),]

# Make binary class out of River system
river$river_class<- 0
river <- within(river, river_class[River == 'Columbia'] <- 1)

# Drop original river column
river<- river[,-c(2)]

# Save data for use on cades
write.table(bud_flush, file="../cades/bud_flush.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(lat, file="../local_data//lat_2.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(river, file="../local_data/river.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
