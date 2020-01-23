library("topGO")
library("qdapRegex")

#####
# PREP DATA FOR BEDTOOLS
#####
# Load VCF format data for sigs
sigs<- read.table("../output/sigHeader.vcf", header = TRUE)

# Extract necessary columns and sort them
#bedt<- sigs[,c(1,2)]
sig_sort<- sigs[order(sigs[,1], sigs[,2]), ]

# Convert position to integer
sig_sort$POS<- as.integer(sig_sort$POS)

# Save
write.table(sig_sort, file="../output/bedt_in.vcf", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Load database file and sort
db<- read.table("../local_data/Ptrichocarpa_210_v3.0.gene.gff", header=FALSE)
db<- db[order(db[,1], db[,4]), ]

# Save db
write.table(db, file="../local_data/Ptrichocarpa_sorted_210_v3.0.gff", sep="\t", row.names = FALSE, 
            col.names = FALSE, quote = FALSE)


## Add a vcf header and a gff header to sorted files
## then run the subsequent CLI command:
## bedtools closest -k 2 -a bedt_in.vcf -b ../local_data/Ptrichocarpa_sorted_210_v3.0.gff

#####
# INTERPRET BEDTOOLS OUTPUT
#####
labeled_snp<- read.table("../output/annotated_sigs.vcf", header = FALSE)

# Move labels to front of df
labeled_snp<- labeled_snp[,c(927:935, 1:926)]

# Assign meaningful column names
colnames(labeled_snp)[1:9]<- c("Chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
colnames(labeled_snp)[10:ncol(labeled_snp)]<- colnames(sigs)

# Extract POTRIIDs into their own column
geneID<- unlist(rm_between(labeled_snp$attribute[2], c('Name='), c('((;)|($)'), extract=TRUE, fixed=FALSE))
labeled_snp<- cbind(geneID, labeled_snp)

saveRDS(labeled_snp, file="../output/labeled_lat_snps.RDS")

#####
# MAKE TABLE OF ANNOTATIONS
#####
# Load annotation file
role<- read.table("../local_data/Ptrichocarpa_210_v3.0.defline.txt", header = FALSE, sep = "\t")
colnames(role)<- c("geneID", "defLine", "function")

# Select genes of interest
sig_labels<- subset(role, role$geneID %in% labeled_snp$geneID)
sig_labels<- sig_labels[,c(1,3)]

# Make Table
library(dplyr)
library(knitr)
library(DT)
library(xtable)

# Fix column names
colnames(sig_labels)<- c("Gene", "Annotation")

# Create table
print(xtable(sig_labels, caption= "Latitude Gene Annotations",
             display = c("s", "s", "s"), include.rownames = F, 
             type = "latex"), file= "../output/lat_gene_annotations_table.tex")

write.table(sig_labels, file="../output/lat_gene_table.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")

#####
# PERFORM ENRICHMENT TESTS
## Most of the code below written by Tim Yates
#####
# set the output file
output_file<- "../output/enrichment_results.txt"
sink(output_file)

# read in the 'gene universe' file
geneID2GO <- readMappings(file = "../local_data/Ptrichocarpa_444_v3.1.annotation_info_gene_GO.txt")
geneUniverse <- names(geneID2GO)

# read in the genes of interest 
genesOfInterest <- as.character(labeled_snp$geneID)
genesOfInterest<- genesOfInterest[-which(is.na(genesOfInterest))]
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

# build the GOdata object in topGO
myGOdata <- new("topGOdata", description="Lat Association", ontology="BP", 
                allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata

# run the Fisher's exact tests
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
resultElim <- runTest(myGOdata, algorithm="elim", statistic="fisher")
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
resultParentchild <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")

# see how many results we get where weight01 gives a P-value <= 0.001:
mysummary <- summary(attributes(resultTopgo)$score <= 0.001)
numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.001







