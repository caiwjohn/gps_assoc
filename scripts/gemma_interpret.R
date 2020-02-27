library("tidyverse")
library(gdata)
library(splitstackshape)
library(rematch2)

## Make sure server folder is mounted!!

# Define path to results csv
fn<- "../cades/gps_assoc/gwas/columbia_gemma/columbia.csv"

# Define file-saving analysis prefix
save_tag<- "relaxed"

#####
# SELECT SIGSNPS AND SAVE FOR FURTHER USE
#####
# Load GWAS results
gwas<- read.table(fn, header=TRUE, sep=",")

# Determine Bonferroni Corrected p-val
pval<- 0.05/8300000
man_cutoff<- -log10(pval)

# Isolate sig SNPs
#sig<- gwas[which(gwas$P.value<pval),]

# Isolate unique gene IDs
id<- as.data.frame(unique(sig$Gene.ID))

# !
# Remove scaffold SNPs prior to saving
# !

# Save
saveRDS(sig, file=paste0("../output/", save_tag, "_SNPs.RDS"))

# Save GeneIDs
write.table(id, file=paste0("../output/", save_tag, "_geneID.txt"), quote=FALSE, row.names = FALSE)

#####
# COMPARE SELECTED SNPs AND LASSO TRANSCRIPTS
#####
# Load selected transcripts
transcripts<- read.table(file=paste0("../output/", save_tag, "_sel_trans.txt"), skip=1)
snps<- as.data.frame(sig$Gene.ID)

# Determine unique overlaps
overlap<- unique(merge(transcripts, snps, by.x = "V1", by.y = 'sig$Gene.ID'))

# Save as text for reference
write.table(overlap, file=paste0("../output/", save_tag, "_transcript_snp_overlap.txt"), quote=FALSE, row.names = FALSE)

#####
# BUILD SCRIPT TO EXTRACT SNPS FROM VCF TO HARVEST EFFECT TERM
#####
sig<- readRDS(paste0("../output/", save_tag, "_SNPs.RDS"))

# Remove duplicate positions
sig<- sig[-which(duplicated(sig[,1:2])),]

# Create df
extract_sig<- data.frame(command= rep(NA, nrow(sig)+1))

# Make first command
extract_sig$command[1]<- paste0("zgrep -P '#CHROM' Chr01_besc_filter.vcf > ", save_tag, "_SNPs.vcf")

# Iterate thru all
for (snp in 1:nrow(sig)) {
  extract_sig$command[snp+1]<- paste0("zgrep -P '", sig$Chr[snp],"\\t",sig$Pos[snp],"' ", sig$Chr[snp], "_besc_filter.vcf >> ", save_tag, "_SNPs.vcf")
}

# Save as bash script
write.table(extract_sig, file = paste0("../cades/",save_tag, "_extract.sh"), quote=FALSE, row.names = FALSE, col.names = FALSE)

#####
# LABEL EFFECT OF EXTRACTED SNPS
#####
sigs<- read.table("../output/sigHeader.vcf", header = TRUE)

# Isolate annotations
annot<- as.data.frame(sigs$INFO)

# Split labels and convert to char
tations<- concat.split(annot, split.col = 1, sep=";", fixed = FALSE, drop=TRUE) %>%
  mutate_all(as.character)

# Make df for storing effects
eff<- sigs[,1:5]
eff$EFF<- NA
  
#Iterate over each row and fill in effect term
for (row in 1:nrow(tations)){
  col<- grep("EFF=", tations[row,])
  if(length(col)==0){
    next()
  }
  else{
    eff$EFF[row]<- tations[row,col]
  }
}

# Select all SNPs with effect
regSNPs<- eff[which(!is.na(eff$EFF)),]





