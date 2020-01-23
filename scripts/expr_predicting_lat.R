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
# Load data and remove lat_groups
annotated_dat<- readRDS("../local_data/leaf_expr_lat_groups.RDS")
annotated_dat<- select(annotated_dat, -one_of("lat_group"))

# Rearrange columns to my liking
annotated_dat<- annotated_dat %>%
  select("Genotype_ID", "Latitude", "lat_group", everything())

# Build training and validation sets
model_dat <- annotated_dat%>%
  as_tibble() %>%
  group_by(lat_group) %>%
  drop_na(lat_group)

# Shuffle the order of samples
model_dat <- model_dat[sample(nrow(model_dat)),]

# OPTIONAL: shuffle class labels
#shuffled_model_dat<- transform(model_dat, lat_group= sample(lat_group))

# Select training and validation
train<- sample_n(model_dat, 50)
test<- model_dat[which(!(model_dat$Genotype_ID %in% train$Genotype_ID)),]

# Isolate predictors and responses
train_resp<- as.factor(as.matrix(train[,c("lat_group")]))
train_pred<- as.matrix(train[,6:ncol(train)])
test_resp<- as.factor(as.matrix(test[,c("lat_group")]))
test_pred<- as.matrix(test[,6:ncol(test)])

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
write.table(features$features, file ="../output/shfled_sel_trans.txt" , 
            sep="\t", row.names=FALSE, quote = FALSE)


# Predict classes for held-out data
model_pred<- as.numeric(predict(fit, newx = test_pred, s = "lambda.min", type="class"))

# Compute accuracy
acc<- length(which(model_pred==test_resp))/length(model_pred)

#####
# DECISION TREE
#####
# Build random forest
rf<- randomForest(train_pred, train_resp, xtest = test_pred, ytest = test_resp)

# Compute test acc
length(which(rf$test$predicted==test_resp))/length(test_resp)

# Reduce to top 24 most selective transcripts
rank<- rf$importance
rank<- as.data.frame(rank)
rank$Transcripts<- rownames(rank)
rank<- arrange(rank, desc(MeanDecreaseGini))
train_pred_small<- as.data.frame(train_pred[,colnames(train_pred) %in% rank$Transcripts[1:24]])
test_pred_small<- as.data.frame(test_pred[,colnames(test_pred) %in% rank$Transcripts[1:24]])
dim(train_pred)

# Random forest with reduced transcripts
rf_reduced<- randomForest(train_pred_small, train_resp, xtest = test_pred_small, ytest = test_resp)

# Compute test acc
length(which(rf_reduced$test$predicted==test_resp))/length(test_resp)











