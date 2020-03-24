# Script to determine why my sample sizes are different between analyses
library(mgsub)

#####
# DIFFERENCES BETWEEN VCF SNP DATA AND PHENOTYPES
#####
# Load trait data
traits<- readRDS("../local_data/clean_pheno.RDS")

# Load river system data
river<- read.table(file="../local_data/river.txt", sep="\t")
colnames(river)<- c("Genotype_ID", "river_class")
river$river_class<- as.factor(river$river_class)

# Merge river and trait data to compare river class to system
pheno<- merge(river, traits[,c(1:4)])
pheno<- pheno[which(duplicated(pheno)==FALSE), ]
core<- pheno

# Load post-PCA samples (these exist in phenotype and VCF)
pca<- readRDS("../local_data/pca_above10_coreSamples.RDS")
SNPpca<- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[,1],   
                    PC2 = pca$eigenvect[,2],    
                    stringsAsFactors = FALSE)

# Standardize sample names
SNPpca$sample.id<- mgsub(SNPpca$sample.id, "\\.", "-")

# Check if all SNP samples exist in the phenotype data
all(SNPpca$sample.id %in% traits$Genotype_ID)

# Save the 80 sample names that exist in Phenotype but not VCF file
missing80<- core$Genotype_ID[which(!(core$Genotype_ID %in% SNPpca$sample.id))]
write.table(missing80, "../local_data/missing80_pheno2VCF.txt", quote = F, col.names= F, row.names= F, sep="\t")

# Reduce the core set to all those in both VCF and phenotype
core<- core[which(core$Genotype_ID %in% SNPpca$sample.id), ]

# Load Structure predictions
struct<- readRDS("../local_data/struct_pop.RDS")
struct<- struct[,c(1,5)]

# Merge structure pred with core set
gwas<- merge(core, struct, by.x = "Genotype_ID", by.y="GenotypeID")
gwas<- gwas[,c(1,3,2,6,4,5)]

# Save
saveRDS(gwas, "../local_data/sampleSet_GWAS.RDS")

#####
# DIFFERENCES BETWEEN CORE SET AND LEAF EXPRESSION
#####
# Load expr data
expr<- read.csv("../local_data/Leaf_eQTN.csv", header=TRUE, sep=',')
expr_samples<- colnames(expr)
expr_samples<- expr_samples[-c(1)]

# Fix name formatting
expr_samples<- mgsub(expr_samples, "\\.", "-")

# Extract and convert SNP samples
snp_samples<- as.character(core$Genotype_ID)

# Compare SNP and expr
## DEBUG if samples don't align with RNAseq file change the %in% order
overlap_snpExpr<- expr_samples[which(expr_samples %in% snp_samples)]

# Create sampleSet for RNA-seq analyses
rna<- gwas[which(gwas$Genotype_ID %in% overlap_snpExpr),]

# Save
saveRDS(rna, "../local_data/sampleSet_RNA.RDS")


