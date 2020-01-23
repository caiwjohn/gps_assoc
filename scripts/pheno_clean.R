library(tidyverse)
library(useful)
library(textclean)

##
# NOTES
#   NA may be introduced while converting to numeric
##
# Read in phenotypes
orig<- read.csv("../local_data/phenotypes.csv", header=TRUE, sep=',')
colnames(covar)[1]<- "genotype_id"
covar$genotype_id<- as.character(covar$genotype_id)

# Preserve genotype and convert phenotypes to numeric
geno<- covar$genotype_id
covar<- as.data.frame(apply(covar, 2, as.numeric))
covar$genotype_id<- geno

# Define arrays for replacement
#####
patt=c("Bud_Flush", 
       "Bud_Set",
       "Cookie_Lignin",
       "Cookie_S.G.",
       "Cookie_Glucose_release",
       "Cookie_Xylose_release",
       "Cookie_combined_sugar_release",
       "Cookie_radial_shrinkage",
       "Cookie_longititudinal_shrinkage",
       "Cookie_specific_gravity",
       "Cookie_miosture_content_wood",
       "Cookie_miosture_content_bark",
       "Cookie_bark_thickness",
       "Core_Lignin",
       "Core_S.G",
       "Core_Glucose_release",
       "Core_Xylose_release",
       "Core_combined_sugar_release",
       "Branch_Angle",
       "Stem_Fresh_Weight",
       "Stem_Fresh_Volume",
       "Branch_Fresh_Weight",
       "Stem_Weight.to.Volume_ratio",
       "BLUPs_Lignin",
       "BLUPs_S.G",
       "BLUPs_Glucose_release",
       "BLUPs_Xylose_release",
       "BLUPs_Combined_Sugar_release",
       "NonCoppice_DBH",
       "NonCoppice_Ht",
       "NonCoppice_StemVolume",
       "NonCoppice_SYL.BR",
       "Leaf_Green_Weight",
       "Petiole_Length",
       "Petiole_Max_Diameter",
       "Petiole_Min_Diameter",
       "Coppice_StemCount",
       "Coppice_D20",
       "Coppice_Height",
       "Coppice_StemVolume",
       "NonCoppice_Height",
       "Coppice_Rust",
       "Coppice_DBH",
       "NonCoppice_Rust",
       "NonCoppice.Height_to_live_branch",
       "Wound_Healing",
       "Epicormic_branching",
       "chlorophyll_in_bark",
       "slug_browsing",
       "Cookies_Lignin",
       "Cookies_S.G.",
       "Cookies_Glucose_release",
       "Cookies_Xylose_release",
       "Cookies_combined_sugar_release",
       "Cookies_radial_shrinkage",
       "Cookies_longititudinal_shrinkage",
       "Cookies_specific_gravity",
       "Cookies_miosture_content_wood",
       "Cookies_miosture_content_bark",
       "Cookies_bark_thickness",
       "Apical_Dominance",
       "Stem_Volume",
       "Clear_Stem",
       "Bud_Break",
       "Coppice_No_of_stem_sprouts",
       "Stem_Vol",
       "No_of_stems",
       "Height_to_first_branch",
       "number_of_sylleptic_branches",
       "CW_SN",
       "CW_EW",
       "Crown_Area",
       "Stem.Volume_to_Crown.Area_Ratio",
       "Branch_Length",
       "Horizontal_Branch.Length",
       "Branch_Deflection",
       "Number_of_Side_Branches",
       "Branch_Color",
       "Tertiary_Branch",
       "Longest_Side_Branch.Length",
       "radial_shrinkage",
       "longititudinal_shrinkage",
       "specific_gravity",
       "miosture_content_wood",
       "miosture_content_bark",
       "bark_thickness",
       "Pretreated_Glucose_release",
       "Pretreated_Xylose_release",
       "Pretreated_Combined_sugar_release",
       "NoPretreated_Glucose_release",
       "NoPretreated_Xylose_release",
       "NoPretreated_Combined_sugar_release"
)

replce<- c("Bud.Flush", 
  "Bud.Set",
  "Cookie.Lignin",
  "Cookie.S.G_",
  "Cookie.Glucose.release",
  "Cookie.Xylose.release",
  "Cookie.combined.sugar.release",
  "Cookie.radial.shrinkage",
  "Cookie.longititudinal.shrinkage",
  "Cookie.specific.gravity",
  "Cookie.moisture.content.wood",
  "Cookie.moisture.content.bark",
  "Cookie.bark.thickness",
  "Core.Lignin",
  "Core.S.G",
  "Core.Glucose.release",
  "Core.Xylose.release",
  "Core.combined.sugar.release",
  "Branch.Angle",
  "Stem.Fresh.Weight",
  "Stem.Fresh.Volume",
  "Branch.Fresh.Weight",
  "Stem.Weight.to.Volume.ratio",
  "BLUPs.Lignin",
  "BLUPs.S.G",
  "BLUPs.Glucose.release",
  "BLUPs.Xylose.release",
  "BLUPs.Combined.Sugar.release",
  "NonCoppice.DBH",
  "NonCoppice.Ht",
  "NonCoppice.StemVolume",
  "NonCoppice.SYL.BR",
  "Leaf.Green.Weight",
  "Petiole.Length",
  "Petiole.Max.Diameter",
  "Petiole.Min.Diameter",
  "Coppice.StemCount",
  "Coppice.D20",
  "Coppice.Height",
  "Coppice.StemVolume",
  "NonCoppice.Height",
  "Coppice.Rust",
  "Coppice.DBH",
  "NonCoppice.Rust",
  "NonCoppice.Height.to.live.branch",
  "Wound.Healing",
  "Epicormic.branching",
  "chlorophyll.in.bark",
  "slug.browsing",
  "Cookie.Lignin",
  "Cookie.S.G_",
  "Cookie.Glucose.release",
  "Cookie.Xylose.release",
  "Cookie.combined.sugar.release",
  "Cookie.radial.shrinkage",
  "Cookie.longititudinal.shrinkage",
  "Cookie.specific.gravity",
  "Cookie.miosture.content.wood",
  "Cookie.miosture.content.bark",
  "Cookie.bark.thickness",
  "Apical.Dominance",
  "Stem.Volume",
  "Clear.Stem",
  "Bud.Break",
  "Coppice.No.of.stem.sprouts",
  "Stem.Vol",
  "No.of.stems",
  "Height.to.first.branch",
  "number.of.sylleptic.branches",
  "CW.SN",
  "CW.EW",
  "Crown.Area",
  "Stem.Volume.to.Crown.Area.Ratio",
  "Branch.Length",
  "Horizontal.Branch.Length",
  "Branch.Deflection",
  "Number.of.Side.Branches",
  "Branch.Color",
  "Tertiary.Branch",
  "Longest.Side.Branch.Length",
  "Cookie.radial.shrinkage",
  "Cookie.longititudinal.shrinkage",
  "Cookie.specific.gravity",
  "Cookie.moisture.content.wood",
  "Cookie.moisture.content.bark",
  "Cookie.bark.thickness",
  "Pretreated.Glucose.release",
  "Pretreated.Xylose.release",
  "Pretreated.Combined.sugar.release",
  "NoPretreated.Glucose.release",
  "NoPretreated.Xylose.release",
  "NoPretreated.Combined.sugar.release"
)

#####
# Replace messy phenotype names
for (string in 1:length(patt)) {
  cat("Replacing", patt[string], "\n")
  colnames(covar)<- str_replace_all(colnames(covar), patt[string], replce[string])
}

# Gather into cleaner format
covar<- gather(covar, key = "sample_info", value = "value", -genotype_id)

# Split the garden name into new column
covar<- separate(covar, sample_info, into = c("garden", "info"), sep = 2, extra= "merge")

# Split year and phenotype into new column
covar<- separate(covar, info, into=c("year", "phenotype", "info"), sep="_", extra = "merge" )

# Split info and rep into columns
covar<- separate(covar, info, into=c("info", "rep"), sep="_", fill ="left" )
covar$rep<- mgsub(covar$rep, c("rep1", "rep2", "rep3", "rep4"), c(1,2,3,4))
covar$rep<- as.numeric(covar$rep)

# Remove average rows
covar<- covar[-which(covar$info=="Average"),]
covar<- covar[-which(covar$info=="Aver"),]

# Drop info column
covar<- select(covar, -"info")

######

# Read GPS data
gps<- read.csv("../local_data/gps.csv", header=TRUE, sep=',')
gps$Genotype_ID<- as.character(gps$Genotype_ID)
colnames(gps)[8]<- "Altitude"
gps<- select(gps, -c(Comment, Collection, Genotype_No))

# Merge phenotypes and gps
pheno<- merge(gps, covar, by.x = "Genotype_ID", by.y = "genotype_id")
pheno<- pheno[,c(1:8, 10, 9, 11)]

# Save cleaned phenotypes
saveRDS(pheno, file="../local_data/clean_pheno.RDS")



