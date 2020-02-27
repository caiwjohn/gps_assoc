library(tidyverse)
library(useful)
library(ggfortify)
theme_set(theme_gray(base_size = 18))

#####
# PCA
#####
# Load data
annotated_dat<- readRDS("../local_data/leaf_expr_lat_groups.RDS")
river<- read.table(file="../local_data/river.txt", sep="\t")
colnames(river)<- c("Genotype_ID", "river_class")

# Match genotype labels
annotated_dat$Genotype_ID<- gsub("\\.", "-", annotated_dat$Genotype_ID)

# Merge river class into expression data
expr<- merge(river, annotated_dat, by.x = "Genotype_ID", by.y = "Genotype_ID")
expr$river_class<- as.factor(expr$river_class)


# N.B.
## This expression data object has only samples that have valid river labels

###
# Manual Decomposition
###
# Isolate numeric data
num_dat<- expr[,7:ncol(expr)]
rownames(num_dat)<- expr$Genotype_ID

# Decompose and isolate components
decomp<- svd(num_dat)
d<- decomp[[1]]
u<- decomp[[2]]
v<- decomp[[3]]

# Compute total variance
total_var<- sum(d**2)

# Construct PCs
pc<- u %*% diag(d)
pc<- as.data.frame(pc)

# Join annotations to PCs
pc<- cbind(expr[,c("Genotype_ID", "river_class", "bud_flush", "bud_set")], pc)

# Visualize results
ggplot(data= pc, aes(V1, V2, colour= river_class))+
  geom_point() +
  labs(colour="River System", title= "Leaf Expression PCA", x= "PC1", y="PC2")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.text = element_text(size=12),
        legend.title = element_text(size = 12))

#####
# Auto decomposition
#####
# Isolate expression data and compute
pca_results <- expr[,7:ncol(annotated_dat)] %>%
  prcomp(scale=F, center= T, tol= 0.01)

# Visualize PCA
autoplot(pca_results, data= expr, label= F, colour= "river_class") +
  labs(colour="River System", title= "Leaf Expression PCA") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.text = element_text(size=12),
        legend.title = element_text(size = 12)) +
  scale_colour_discrete(name="River System",
                      labels=c("Non-Columbia", "Columbia"))

# Reduce transcripts to those selected from prediction
red_trans<- expr[,which(colnames(expr) %in% features$features)]

# Repeat with reduced set
pca_results <- red_trans %>%
  prcomp(scale=F, center= T, tol= 0.01)

# Visualize PCA
autoplot(pca_results, data= expr, label= F, colour= "river_class") +
  labs(colour="River System", title= "Leaf Expression PCA", subtitle= "Transcripts selected in LASSO Prediction") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.text = element_text(size=12),
        legend.title = element_text(size = 12), plot.subtitle = element_text(hjust = 0.5)) +
  scale_colour_discrete(name="River System",
                        labels=c("Non-Columbia", "Columbia"))
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

#####
# NMF
#####
library(NMF)

## Comments 
# 1. NMF scale samples by their mean value before log-transform, i.e. divide by its mean

# Shift all data up by 20 to avoid negative numbers
shifted<- num_dat+20

# Factorize
res<- nmf(shifted, rank=4)

# Visualize
basismap(res, subsetRow=TRUE)



