library(useful)
library(WGCNA)
library(tidyverse)

#####
# PREP DATA
#####

# The following setting is important, do not omit
options(stringsAsFactors = FALSE)

# Load raw leaf expression
raw<- read.csv("../local_data/Leaf_eQTN.csv")
rownames(raw)<- raw$SampleID
raw<- raw[,-c(1)]
datExpr0<- as.data.frame(t(raw))

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

# If the above line is False, some samples/genes need to be removed
if (!gsg$allOK){
  
# Optionally, print the gene and sample names that were removed:
if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  
# Remove the offending genes and samples from the data:
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Cluster samples to look for outliers
sampleTree = hclust(dist(datExpr0), method = "average")

# Plot the sample tree
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)

# No outliers to remove so simply rename data for subsequent processing
datExpr<- datExpr0

#####
# THIS STEP NOT NECESSARY
#####
# Load trait data to incorporate into tree
pheno<- readRDS("../local_data/clean_pheno.RDS")

test<- pheno %>%
  group_by(phenotype) %>%
  mutate(row=row_number()) %>%
  pivot_wider(names_from = phenotype, values_from = value) %>%
  select(-row)
#####

#####
# FIT AND BUILD MODULES
#####

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

## Use these plots to pick a thresholding power, I chose 8

# Perform clustering
net = blockwiseModules(datExpr, power = 8,
TOMType = "unsigned", minModuleSize = 30,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = "femaleMouseTOM",
verbose = 3)

#####
# PLOT AND SAVE RESULTS
#####

# open a graphics window
sizeGrWindow(12, 9)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)

# Extract info and save
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
eigenGenes = net$MEs;
geneTree = net$dendrograms[[1]];
save(eigenGenes, moduleLabels, moduleColors, geneTree, file = "../output/leaf_expr_modules.RDS")













