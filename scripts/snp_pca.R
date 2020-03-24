library(SNPRelate)
library(ggplot2)
library(tidyverse)
library(ggfortify)
library(ggdendro)
library(dendextend)
theme_set(theme_gray(base_size = 18))

######
# LOAD VCF, CONVERT AND RETAIN BIALLELIC
## SECTION TO BE RUN ON SERVER
######
# Define file path
vcf.fn <- "../cades/All_besc_filter.vcf"

# Convert file
snpgdsVCF2GDS(vcf.fn, "all.gds", method="biallelic.only")

# View created file
snpgdsSummary("../cades/all.gds")

######
# LD PRUNING AND PCA
## SECTION TO BE RUN ON SERVER
######
# Open file
genofile <- snpgdsOpen("../cades/all.gds")

# Try different LD thresholds for sensitivity analysis
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)

# Select IDs
snpset.id <- unlist(snpset)

# Run PCA
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)

# Save result
saveRDS(pca, "./pca_0.2.rds")

######
# INTERPRET PCA
######
# Load PCA obj
pca<- readRDS("../output/snp_analysis/pca/pca_core_results.RDS")

# Load covariates
covar<- readRDS("../local_data/core_SNP.RDS")

# Inspect PC %
pc.percent <- pca$varprop*100

# Create df for plotting pca
tab <- data.frame(sample.id = pca$sample.id,
                  PC1 = pca$eigenvect[,1],   
                  PC2 = pca$eigenvect[,2],    
                  stringsAsFactors = FALSE)

#####
# SINGLE PCA
#####
# Merge into one df for plotting
snp_pca<- merge(covar, tab, by.x = "Genotype_ID", by.y= "sample.id")

# Remove duplicates
snp_pca <- snp_pca[which(duplicated(snp_pca$Genotype_ID)==FALSE),]

# PC Plot using ggplot
ggplot(data= snp_pca, aes(PC1, PC2, colour=Latitude)) +
  geom_point() +
  ylab("PC2 (0.74% Var Explained)") +
  xlab("PC1 (1.30% Var Explained)") +
  ggtitle("Genotype PCA", subtitle="Mixed and Mislabeled Samples Removed") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle= element_text(hjust=0.5), 
        legend.text = element_text(size=12), legend.title = element_text(size = 12))# +
  #scale_colour_discrete(name="River System")

# Pairwise PC 1-4
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], labels=lbls)

# Determine mislabeled columbia
colum_mislab<- as.data.frame(snp_pca$Genotype_ID[which(snp_pca$PC1 < 0.025 & snp_pca$river_class=="Columbia")])
colnames(colum_mislab)<- c("samples")

# Get non-Columbia mislabeled
nonColum_mislab<- as.data.frame(snp_pca$Genotype_ID[which(snp_pca$PC1 > 0.025 & snp_pca$river_class=="non-Columbia")])
colnames(nonColum_mislab)<- c("samples")

# Combine
mislabeledSNP<- rbind(colum_mislab, nonColum_mislab)

# Save mislabeled names
write.table(mislabeledSNP, "../output/snp_analysis/snp_mislabeled_samples.txt", quote= F, row.names = F, col.names = F)

