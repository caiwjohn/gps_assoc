library(tidyverse)

# Load covariates
samples<- readRDS("../local_data/core_SNP.RDS")

# Keep just variables I need for GEMMA
river<- samples[,c("Genotype_ID", "numeric_class")]

# Build training and validation sets
river <- river%>%
  as_tibble() %>%
  group_by(numeric_class)

# Shuffle the order of samples
river <- river[sample(nrow(river)),]

# See how many of each class we have
length(which(river$numeric_class==1))
length(which(river$numeric_class==0))

# Select even number of each
river_balanced<- sample_n(river, 88)

# Save data for use on cades
write.table(river_balanced, file="../local_data/gemma_core_balanced.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
