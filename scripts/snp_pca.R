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
pca<- readRDS("../local_data/pca_above10_coreSamples_reducedSNP.RDS")

# Load covariates
covar<- readRDS("../local_data/clean_pheno.RDS")

# Inspect PC %
pc.percent <- pca$varprop*100

# Create df for plotting pca
tab <- data.frame(sample.id = pca$sample.id,
                  PC1 = pca$eigenvect[,1],   
                  PC2 = pca$eigenvect[,2],    
                  stringsAsFactors = FALSE)
#####
# HI THRU-PUT PCA
#####
# Iterate through all garden's, rep's and years
for (garden in unique(covar$garden)) {
  for (rep in unique(covar$rep)){
    for(year in unique(covar$year)){

      if(!is.na(rep)){
      # Subset to select garden and rep and year
      pheno1<- subset(covar, garden== garden & rep==rep & year==year)
      
      # Subset to select phenotypes
      pheno1<- subset(pheno1, phenotype=="Bud.Flush" )
      
      # Remove duplicates
      pheno1 <- pheno1[which(duplicated(pheno1)==FALSE),]
      
      # Spread phenotypes for plotting
      #pheno1<- spread(pheno1, phenotype, value)
      
      # Remove entries with empty River classifiation
      pheno1<- pheno1[which(pheno1$River!=""),]
      
      # Merge covar and tab
      pcplot<- merge(pheno1, tab, by.y= "sample.id", by.x = "Genotype_ID")
      
      # PC Plot using ggplot
      ggplot(data= pcplot, aes(PC2, PC1, colour=River)) +
        geom_point() +
        xlab("PC2 (0.47% Var Explained)") +
        ylab("PC1 (0.73% Var Explained)") +
        ggtitle("PCA of Poplar Genotypes") +
        theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.text = element_text(size=12),
              legend.title = element_text(size = 12),
              plot.subtitle = element_text( hjust = 0.5) ) +
        labs(colour="River System", subtitle = paste0(garden, " Garden, rep ", rep, " in ", year)) +
        ggsave(filename=paste0("../output/pca/", garden, "_", rep,"_", year), device= "png")
      }
    }
  }
}
#####
#####
# SINGLE PCA
#####
# Load river system data
river<- read.table(file="../local_data/river.txt", sep="\t")
colnames(river)<- c("Genotype_ID", "river_class")
river$river_class<- as.factor(river$river_class)

# Load structure population data
struct<- readRDS("../local_data/struct_pop.RDS")
struct<- struct[,c(1,4,5)]

## Perform this if sample ID's differ
tab$sample.id<- mgsub(tab$sample.id, "\\.", "-")

# Merge into one df for plotting
snp<- merge(river, tab, by.x = "Genotype_ID", by.y= "sample.id")

# Remove duplicates
snp <- snp[which(duplicated(snp$Genotype_ID)==FALSE),]

# Combine PCA data and structure pops
snp<- merge(snp, struct, by.x = "Genotype_ID", by.y = "GenotypeID")

# PC Plot using ggplot
ggplot(data= snp, aes(PC1, PC2, colour=River)) +
  geom_point() +
  ylab("PC2 (2.08% Var Explained)") +
  xlab("PC1 (91.61% Var Explained)") +
  ggtitle("Genotype PCA", subtitle="GEMMA SNPs above 10: ends removed") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle= element_text(hjust=0.5), 
        legend.text = element_text(size=12), legend.title = element_text(size = 12)) +
  scale_colour_discrete(name="River System")#,
                        #labels=c("Non-Columbia", "Columbia"))

# Pairwise PC 1-4
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], labels=lbls)

######
# RUN IBD
## SECTION TO BE RUN ON SERVER
######
genofile <- snpgdsOpen("../cades/all.gds")

# YRI samples
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
YRI.id <- sample.id[pop_code == "YRI"]

# Estimate IBD coefficients
ibd <- snpgdsIBDMoM(genofile, sample.id=YRI.id, snp.id=snpset.id,
                    maf=0.05, missing.rate=0.05, num.thread=2)

# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)

######
# VISUALIZE IBD
######
ibd.coeff<- readRDS("../local_data/ibd_coeff.RDS")

plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="YRI samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)

######
# HIERARCHICAL CLUS. IBD
######
# Isolate kinship
kin<- ibd.coeff[,-c(3,4)]

# Convert from kinship to distance proxy
kin$kinship<- (1-kin$kinship)

# Spread into matrix
kin_mat<- spread(kin, ID2, kinship, fill=1)

# Merge phen. into kinship
pheno_kin<- merge(covar, kin_mat, by.y = "ID1", by.x = "Genotype_ID", sort=FALSE)

# Make denodrogram
dend_kin<- as.dist(kin_mat[,-c(1)]) %>%          #pheno_kin[,-c(1:5)]) %>%
          hclust %>%
          #as.dendrogram() %>%
          plot(labels=FALSE, main= "Clustering of Kinship Coeff. for Poplar Genomes")

# Create labels df
tree_labels<- dendro_data(dend_kin, type = "rectangle")
tree_labels$labels<- merge(x= tree_labels$labels, y= kin_mat, by.x= "label", by.y= "ID1")

# Plot dendrogram with colored branches and legend
ggplot() +
  geom_segment(data=segment(tree_labels), aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_segment(data = tree_labels$segments %>%
                 filter(yend == 0) %>%
                 left_join(tree_labels$labels, by = "x"), aes(x=x, y=y.x, xend=xend, yend=yend, color = Diagnosis)) +
  geom_text(data = label(tree_labels), aes(x=x, y=y, label=label, colour = Diagnosis, hjust=0), size=1) +
  coord_flip() +
  scale_y_reverse(expand=c(0.2, 0)) +
  scale_colour_brewer(palette = "Spectral") + 
  theme_dendro()
