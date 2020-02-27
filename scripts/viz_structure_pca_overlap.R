library(tidyverse)
library(mgsub)
library(VennDiagram)
library(gridExtra)
library(ggpubr)
theme_set(theme_gray(base_size = 18))

# Load populations from fastSTRUCTURE runs
struc<- read.table("../local_data/structure/run_k.2.meanQ")

# Load bed file to determine sample order
samples<- read.table("../cades/gwas/lat_2_gemma/temp_bed/Chr01_besc_filter.fam")

# Combine sample labels and pops
struc<- cbind(samples[,1], struc)

# Name columns
colnames(struc)<- c("GenotypeID", "Pop1", "Pop2")

# Load phenotypes to get river system
pheno<- readRDS("../local_data/clean_pheno.RDS") %>%
          select(Genotype_ID, River) %>%
          unique(.)

# Merge river system in data
pops<- merge(struc, pheno, by.x = "GenotypeID", by.y= "Genotype_ID")

# Pred river system based on STRUCTURE population
pops$Struct_Pred[pops$Pop2 > 0.75]<- "Columbia"
pops$Struct_Pred[pops$Pop1 > 0.75]<- "non-Columbia"
pops$Struct_Pred[which(is.na(pops$Struct_Pred))] <- "Mixed"

# Replace River labels to match predictions
## Make sure textclean package isn't loaded, there are two mgsub functions
pops$River<- as.character(pops$River)
pops$River<- mgsub(pops$River, pattern=c("Puyallup", "Puyallup_Carbon", "Skyomish", "Skagit", "Skykomish"), 
                         replacement = c("non-Columbia"), recycle=T)

# Drop rows with no true river classification
pops<- pops[which(pops$River!=""),]

# Save data obj
saveRDS(pops, "../local_data/struct_pop.RDS")

# Determine accuracy of prediction method
acc<- length(which(pops$River==pops$Struct_Pred))/nrow(pops)
colnames(pops)[4:5]<- c("River System", "STRUCTURE Population")

# Create Venn Diagram of relationship
match_venn<- draw.pairwise.venn(area1 = 463, area2=7, cross.area=0 , category = c("Matching", "Not Matching"), euler.d = TRUE, sep.dist = 0.005, ext.text = FALSE, lty = "blank",
                          fill = c("lightskyblue", "brown1"), cat.pos = c(0,0), alpha = 0.8)

# Add titles
title=text_grob("Coherency Between River System and STRUCTURE Populations", face="bold", size= 18 )
grid.arrange(gTree(children=match_venn), top=title)

# Create plotting structure
freq<- data.frame(
              River= c("Columbia", "non-Columbia"),
              True= c(length(which(pops$True=="Columbia")), length(which(pops$True=="non-Columbia"))),
              Pred=c(length(which(pops$Pred=="Columbia")), length(which(pops$Pred=="non-Columbia")))
)

# Gather for plotting
barplot<- gather(freq, key="Status", value="Count", True, Pred)

# Plot data
ggplot(barplot, aes(River, Count))+
  geom_col(aes(fill=Status), position = "dodge")+
  labs(fill= "Population", title= "Populations from STRUCTURE vs River System", 
       x= "River System") + 
  scale_fill_discrete(name = "Population", labels = c("STRUCTURE", "River System")) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.text = element_text(size=12),
        legend.title = element_text(size = 12))
  
  





