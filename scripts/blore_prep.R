library("qdapRegex")
library("topGO")

#####
# ADD LABELS TO VCF FILE
#####
# Run this in 'blore_data' dir while mounted locally to map SNPs to known IDs
#bedtools closest -k 2 -a all.vcf -b ../local_data/Ptrichocarpa_sorted_210_v3.0.gff >all_labeled.vcf

# Load data
snps<- read.table("../cades/data/blore_data/relaxed_SNPs.vcf", header = TRUE, comment.char = '')

# Fix first column name
colnames(snps)[1]<- "CHROM"

# Remove duplicates
snps<- snps[which(duplicated(snps)==F),]

# Add generic ID values to column
snps$ID <- paste0("SNP",seq.int(nrow(snps)))

# Move labels to front of df
#labeled_snp<- labeled_snp[,c(927:935, 1:926)]

# Assign meaningful column names
#colnames(labeled_snp)[1:9]<- c("Chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
#colnames(labeled_snp)[10:ncol(labeled_snp)]<- colnames(sigs)

# Extract POTRIIDs into their own column
#geneID<- unlist(ex_between(labeled_snp$attribute, c('Name='), c('((;)|($))'), fixed=FALSE))
#labeled_snp<- cbind(geneID, labeled_snp)

# Remove extra columns to create valid VCF
#labeled_snp$ID<- labeled_snp$geneID
#all_vcf<- labeled_snp[,-c(1:10)]

# Save unique gene list for functional analysis
#write.table(genes, "../output/columbia_cutoff10_mapped_genes.txt", quote=FALSE, row.names = FALSE, col.names = F)

#####
# Split selected SNPs into file for each loci
#####
# Split into df for each chr
sep<- split(snps, snps$CHROM)

# Inspect to determine number of chrom
sort(unique(snps$CHROM))

chr01<- as.data.frame(sep[1])
chr02<- as.data.frame(sep[2])
chr03<- as.data.frame(sep[3])
chr04<- as.data.frame(sep[4])
chr05<- as.data.frame(sep[5])
chr06<- as.data.frame(sep[6])
chr07<- as.data.frame(sep[7])
chr08<- as.data.frame(sep[8])
chr09<- as.data.frame(sep[9])
chr10<- as.data.frame(sep[10])
chr11<- as.data.frame(sep[11])
chr12<- as.data.frame(sep[12])
chr13<- as.data.frame(sep[13])
chr14<- as.data.frame(sep[14])
chr15<- as.data.frame(sep[15])
chr16<- as.data.frame(sep[16])
chr17<- as.data.frame(sep[17])
chr18<- as.data.frame(sep[18])
chr19<- as.data.frame(sep[19])

# Check assignment worked
corner(chr07)

# Assign colnames
colnames(chr01)<- colnames(snps)
colnames(chr03)<- colnames(snps)
colnames(chr04)<- colnames(snps)
colnames(chr05)<- colnames(snps)
colnames(chr06)<- colnames(snps)
colnames(chr07)<- colnames(snps)
colnames(chr08)<- colnames(snps)
colnames(chr09)<- colnames(snps)
colnames(chr10)<- colnames(snps)
colnames(chr11)<- colnames(snps)
colnames(chr12)<- colnames(snps)
colnames(chr13)<- colnames(snps)
colnames(chr14)<- colnames(snps)
colnames(chr15)<- colnames(snps)
colnames(chr16)<- colnames(snps)
colnames(chr17)<- colnames(snps)
colnames(chr18)<- colnames(snps)
colnames(chr19)<- colnames(snps)

# Save dataframes with over 100 SNPs for BLORE
write.table(chr01, "../cades/data/blore_data/Chr01.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr02, "../cades/data/blore_data/Chr02.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr03, "../cades/data/blore_data/Chr03.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr04, "../cades/data/blore_data/Chr04.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr05, "../cades/data/blore_data/Chr05.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr06, "../cades/data/blore_data/Chr06.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr07, "../cades/data/blore_data/Chr07.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr08, "../cades/data/blore_data/Chr08.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr09, "../cades/data/blore_data/Chr09.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr10, "../cades/data/blore_data/Chr10.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr11, "../cades/data/blore_data/Chr11.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr12, "../cades/data/blore_data/Chr12.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr13, "../cades/data/blore_data/Chr13.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr14, "../cades/data/blore_data/Chr14.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr15, "../cades/data/blore_data/Chr15.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr16, "../cades/data/blore_data/Chr16.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr17, "../cades/data/blore_data/Chr17.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr18, "../cades/data/blore_data/Chr18.vcf", quote=F, sep = "\t", row.names = F)
write.table(chr19, "../cades/data/blore_data/Chr19.vcf", quote=F, sep = "\t", row.names = F)

## !!
# open saved files manually and add '#' to start of each header line
## !!

#####
# Add phenotype to sample file
#####

# Load phenoype data
pheno<- read.table("../cades/blore_anal/river.txt")
colnames(pheno)<- c("GenotypeID", "pheno")

# Load sample file
sample<- read.table("../cades/blore_anal/finemap.sample", header = TRUE)

# Merge files
sample<- merge(sample, pheno, by.x = "ID_1", by.y= "GenotypeID", all.x = TRUE, sort=T)

# Change second header to binary format
sample$pheno[1]<- "B"

# Save sample file
write.table(sample, "../cades/blore_anal/finemap_pheno.sample", quote=F, row.names = F)
