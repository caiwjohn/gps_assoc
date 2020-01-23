library(tidyverse)

# Load budflush SNPs
bud_flush<- read.table("../local_data/budflush_sig_snps.csv", sep=",")
colnames(bud_flush)<- c("ID", "Chr", "bp", "ref", "alt", "info", "type", "geneID")

bud_flush$geneID <- as.character(bud_flush$geneID)

# Load lat SNPs
lat<- readRDS("../output/labeled_lat_snps.RDS")
lat$Gene.ID<- as.character(lat$Gene.ID)

# Compare SNPs
which(lat$Gene.ID %in% bud_flush$geneID)
