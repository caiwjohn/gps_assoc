library(tidyverse)
library(useful)
library(ggplot2)
library("glmnet")
library(lmtest)
library(textclean)
library(ROCR)
library(randomForest)

#####
# PREP DATA FOR PREDICTION
#####
# Load river system data
river<- read.table(file="../local_data/river.txt", sep="\t")
colnames(river)<- c("Genotype_ID", "river_class")
river$river_class<- as.factor(river$river_class)

# Load expr data
annotated_dat<- readRDS("../local_data/leaf_expr_lat_groups.RDS")

# Match genotype labels
annotated_dat$Genotype_ID<- gsub("\\.", "-", annotated_dat$Genotype_ID)

# Merge river class into expression data
expr<- merge(river, annotated_dat, by.x = "Genotype_ID", by.y = "Genotype_ID")

# Build training and validation sets
model_dat <- expr%>%
  as_tibble() %>%
  group_by(river_class)

# Shuffle the order of samples
model_dat <- model_dat[sample(nrow(model_dat)),]

# OPTIONAL: shuffle class labels
#shuffled_model_dat<- transform(model_dat, lat_group= sample(lat_group))

# Select training and validation
train<- sample_n(model_dat, 50)
test<- model_dat[which(!(model_dat$Genotype_ID %in% train$Genotype_ID)),]

# Isolate predictors and responses
train_resp<- as.factor(as.matrix(train[,c("river_class")]))
train_pred<- as.matrix(train[,7:ncol(train)])
test_resp<- as.factor(as.matrix(test[,c("river_class")]))
test_pred<- as.matrix(test[,7:ncol(test)])

#####
# LASSO PRED
#####
# Determine optimal lambda value
fit<- cv.glmnet(train_pred, train_resp, family = "binomial", 
                type.measure= "class", nfolds = 10)

# Plot lambda optimization
plot(fit)

# Find corresponding minimum misclassification
mse.min <- fit$cvm[fit$lambda == fit$lambda.min]

# View coefficients
myCoefs<- coef(fit, s = "lambda.min")

# Assemble into df
features <- data.frame(
  features = myCoefs@Dimnames[[1]][which(myCoefs != 0 )],
  coefs    = myCoefs[ which(myCoefs != 0 ) ]
)

# OPTIONAL: Save list of selected features
features<- features[-c(1),]
write.table(features$features, file ="../output/columbia_sel_trans.txt" , 
            sep="\t", row.names=FALSE, quote = FALSE)

# Predict classes for held-out data
model_pred<- as.numeric(predict(fit, newx = test_pred, s = "lambda.min", type="class"))

# Compute accuracy
acc<- length(which(model_pred==test_resp))/length(model_pred)






















