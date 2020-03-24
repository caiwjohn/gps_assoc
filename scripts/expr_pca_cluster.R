library(tidyverse)
library(useful)
library(ggfortify)
library(mgsub)
theme_set(theme_gray(base_size = 18))

#####
# PCA
#####
# Load data
annotated_dat<- readRDS("../local_data/leaf_expr_lat_groups.RDS")
samples<- readRDS("../local_data/sampleSet_RNA.RDS")

# Change class label to a factor string
samples$numeric_class<- samples$river_class
samples$river_class<- as.character(samples$river_class)
samples$river_class<- mgsub(samples$river_class, c("1", "0"), c("Columbia", "non-Columbia"))

# Tidy samples
samples<- samples[,c(1:3, 7, 4:6)]

# Drop old annotations
annotated_dat<- annotated_dat[,-c(2,3,4,5)]

# Match genotype labels
annotated_dat$Genotype_ID<- gsub("\\.", "-", annotated_dat$Genotype_ID)

# Merge river class into expression data
expr<- merge(samples, annotated_dat)
expr$river_class<- as.factor(expr$river_class)

# Save as core expression data
saveRDS(expr, "../local_data/core_Expression.RDS")

# Drop samples with Mixed label
expr<- expr[which(expr$Struct_Pred!="Mixed"),]

#####
# Auto decomposition
#####
# Isolate expression data and compute
pca_results <- expr[,8:ncol(expr)] %>%
  prcomp(scale=F, center= T, tol= 0.01)

# Visualize PCA
autoplot(pca_results, data= expr, label= F, colour= "river_class") +
  labs(colour="River System", title= "Leaf Expression PCA", subtitle = "Mixed Samples Removed") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle= element_text(hjust=0.5), 
        legend.text = element_text(size=12), legend.title = element_text(size = 12)) +
  scale_colour_discrete(name="River System")

# Reduce transcripts to those selected from prediction
red_trans<- expr[,which(colnames(expr) %in% features$features)]

# Repeat with reduced set
pca_results <- red_trans %>%
  prcomp(scale=F, center= T, tol= 0.01)

# Combine pca with traits
plotDat<- cbind(pca_results$x, expr[,c(1, 3)])

# Visualize single dimensional PCA
ggplot(data= plotGene, aes(seq_along(Genotype_ID), red_trans, colour= river_class))+
  geom_point() +
  labs(x= "Samples", y= "Potri.010G079500" ,colour="River System", title= "Leaf Expression", subtitle= "LASSO Transcript, Core Samples") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.text = element_text(size=12),
        legend.title = element_text(size = 12), plot.subtitle = element_text(hjust = 0.5)) +
  scale_colour_discrete(name="River System")

# Visualize PCA
autoplot(pca_results, data= expr, label= F, colour= "river_class", loadings=F, loadings.label=F, scale=0) +
  labs( colour="River System", title= "Leaf Expression PCA", subtitle= "LASSO Transcripts, Core Samples") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.text = element_text(size=12),
        legend.title = element_text(size = 12), plot.subtitle = element_text(hjust = 0.5)) +
  scale_colour_discrete(name="River System")

# Extract samples that are mislabeled
pc<- as.data.frame(pca_results$x[,1:2])
pc$samples<- expr$Genotype_ID
pc$River<- expr$river_class

# Get mislabeled columbia sample names
colum_mislab<- as.data.frame(pc$samples[which(pc$PC1 < 1.35 & pc$River=="Columbia")])
colnames(colum_mislab)<- c("samples")

# Get non-Columbia mislabeled
nonColum_mislab<- as.data.frame(pc$samples[which(pc$PC1 > 1.36 & pc$River=="non-Columbia")])
colnames(nonColum_mislab)<- c("samples")

# Combine
mislabeled<- rbind(colum_mislab, nonColum_mislab)

# Save mislabeled names
write.table(mislabeled, "../output/rna/mislabeled_samples.txt", quote= F, row.names = F, col.names = F)

#####
# Hierarchical clustering
#####
library(cluster)
library(dendextend)
library(factoextra)

set.seed(123)

# Compute distances
dis<- dist(num_dat, method="euclidean")

# Cluster
hc<- hclust(dis, method= "complete")

# Use Elbow method to determine number of clusters
fviz_nbclust(num_dat, kmeans, method = "wss", k.max=20)

# Cut tree to optimal clusters
sub_grp <- cutree(hc, k = 12)

# Add cluster groups to data
dat<- dat %>%
  mutate(cluster = sub_grp)

# Visualize dendrogram with clusters labelled
plot(hc, cex = 0.6, hang=-1,main= "Dendrogram of Leaf Expression Samples")
rect.hclust(hc, k = 6, border = 2:5)

# Use GAP statistic method to determine optimal number of clusters
##gap_stat <- clusGap(num_dat, FUN = hcut, nstart = 25, K.max = 10, B = 50)
##fviz_gap_stat(gap_stat)



