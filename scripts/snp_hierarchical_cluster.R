library(adegenet)
library(SNPRelate)
library(ggplot2)
library(tidyverse)
library(ggfortify)
library(ggdendro)
library(dendextend)
library(reshape2)
library(grid)
theme_set(theme_gray(base_size = 18))

######
# LOAD VCF, CONVERT AND RETAIN BIALLELIC
######
# Define file path
vcf.fn <- "../cades/data/filtered_vcfs/above10_SNPS.vcf"

# Convert file
snpgdsVCF2GDS(vcf.fn, "../cades/data/filtered_vcfs/above10.gds", method="biallelic.only")

######
# RUN IBD
######
genofile <- snpgdsOpen("../cades/data/filtered_vcfs/above10.gds")

# YRI samples
#sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
#YRI.id <- sample.id[pop_code == "YRI"]

# Estimate IBD coefficients
ibd <- snpgdsIBDMoM(genofile, num.thread=2)

# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)

######
# VISUALIZE IBD
######
plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="YRI samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)

######
# HIERARCHICAL CLUS. IBD
######
# Load phenotypes
## See alternative below
covar<- readRDS("../local_data/clean_pheno.RDS")

# Select only phenotypes I need
covar<- subset(covar, garden=='CL' & year==2010 & rep==1)

# Drop other phenotypes
covar<- covar[,1:6]
covar<- covar[which(!duplicated(covar)),]

# Only core sample set
covar<- covar[which(covar$River!=''),]

## If only interested in river system load that data directly
# Columbia river system encoded as '1' in this file
river<- read.table("../local_data/river.txt")
colnames(river)<- c("Genotype_ID", "River")
##

##
# BUILD DENDROGRAM
##
# Isolate kinship
kin<- ibd.coeff[,-c(3,4)]

# Load structure population data
struct<- readRDS("../local_data/struct_pop.RDS")
struct<- struct[,c(1,4,5)]

# Convert from kinship to distance proxy
kin$kinship<- (1-kin$kinship)

# Spread into matrix
kin_spread<- spread(kin, ID2, kinship, fill=1)

# Merge phen. into kinship
pheno_kin<- merge(struct, kin_spread, by.y = "ID1", by.x = "GenotypeID", sort=FALSE)

# Drop columns labeled as mixed
pheno_kin<- subset(pheno_kin, Struct_Pred != 'Mixed')

# Cluster
kin_mat<- as.matrix(pheno_kin[,-c(1,2,3)])
rownames(kin_mat)<- pheno_kin$Genotype_ID
kin_dendro <- as.dendrogram(hclust(d = dist(x = kin_mat)))

# Create dendro
dendro.plot <- ggdendrogram(data = kin_dendro, rotate = TRUE) +
               theme(axis.text.y = element_blank())

# Preview the plot
print(dendro.plot)

##
# BUILD HEATMAP
##
# Condense data to right format
kin_long <- melt(pheno_kin, id = c("GenotypeID", "River", "Struct_Pred"))

# Convert to factor
kin_long$River<- as.factor(kin_long$River)

# View heatmap
heatmap<- ggplot(data = kin_long, aes(x = variable, y = GenotypeID)) +
  geom_tile(aes(fill = River)) +
  ggtitle("Heatmap of Poplar Genotypes", subtitle="GEMMA Selected SNP Subset") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle= element_text(hjust=0.5), 
        legend.text = element_text(size=12), legend.title = element_text(size = 12),
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  scale_fill_discrete(name="River System",labels=c("Non-Columbia", "Columbia"))

print(heatmap)
##
# COMBINE HEATMAP AND DENDROGRAM
##
# Extract order from dendrogram to match
kin_order <- order.dendrogram(kin_dendro)

# Order the levels according to their position in the cluster
kin_long$GenotypeID <- factor(x = kin_long$GenotypeID,
                               levels = pheno_kin$GenotypeID[kin_order],
                               ordered = TRUE)

# Recreate heatmap in right order
heatmap.plot <- ggplot(data = kin_long, aes(x = variable, y = GenotypeID)) +
  geom_tile(aes(fill = River)) +
  ggtitle("Clustering Poplar Genotypes by Kinship", subtitle="GEMMA Selected SNP Subset") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle= element_text(hjust=0.5), 
        legend.text = element_text(size=12), legend.title = element_blank(),
        axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        legend.position = 'top')

# View full viz
grid.newpage()
print(heatmap.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.422, width = 0.2, height = 1))

