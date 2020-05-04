##########################################
############ GRO in isolation ############
##########################################
# Copyright 2020, Bérénice A. Benayoun
# 2020-05-06

#######        Week 7 code         #######

### Practical Machine Learning 101 ###

###################################################################################################################
### Set the working directory
setwd('/Users/berenice/Dropbox/USC_Work_Folder/Teaching/2019-2020/GRO-in-isolation/Github/Week7')
options(stringsAsFactors = F)

### load necessary libraries
library(caret)
library(MLmetrics)
library(kernlab)        # for SVM
library(randomForest)   # for Random Forest
library(beeswarm)       # overlaying points as beeswarm
library(ALL)            # data for example 2
library(biomaRt)        # to annotate results of example 2

########################################################################################
### Example #1
### Can we predict the species of iris based on 4 features:
### Sepal.Length, Sepal.Width, Petal.Length,Petal.Width

data(iris)

##############################################################
##### 1. partition data in training and testing set
##############################################################

set.seed(123456789)                                                      # setting the seed helps make your code doing "random" sampling reproducible

my.training.idx <- createDataPartition(iris$Species, p=0.67, list=FALSE) # 2/3 for training, 1/3 for testing
my.training.data <- iris[my.training.idx,]

summary(my.training.data)                                                # the classes are balanced

save(my.training.data, my.training.idx, file = paste0(Sys.Date(),"_Iris_Training_Data.RData"))

##############################################################
##### 2. train an SVM model
##############################################################

# use 10-fold cross-validation to build the model
my.ctrl.opt.svm <- trainControl(method          = "cv",
                                number          = 10,
                                allowParallel   = TRUE,
                                verbose         = F,
                                summaryFunction = multiClassSummary, # here, 3 species of iris, so multiclass. Check twoClassSummary for 2-class problems
                                classProbs      = TRUE)

# train model with caret train function
my.svm.fit       <- train( Species ~ .,
                           data       = my.training.data        ,
                           method     ="svmRadial"              ,  # SVM require kernels. Here, use of Gaussian kernel
                           trControl  = my.ctrl.opt.svm         ,
                           tuneLength = 10                      ,  # test 10 parameters combinations
                           metric     = "Mean_Balanced_Accuracy")
save(my.svm.fit, file = paste0(Sys.Date(),"_iris_SVM_model.RData"))

# Check the learned model
my.svm.fit

##############################################################
##### 3. train a RF model
##############################################################

set.seed(123456789)  # RF is based on random variable sampling

# use 10-fold cross-validation to build the model
my.ctrl.opt.rf <- trainControl(method          = "cv",
                               number          = 10,
                               allowParallel   = TRUE,
                               verbose         = F,
                               summaryFunction = multiClassSummary, # here, 3 species of iris, so multiclass. Check twoClassSummary for 2-class problems
                               classProbs      = TRUE)

# train model with caret train function
my.rf.fit <- train( Species ~ .,
                    data        =  my.training.data  ,
                    method      = "rf"               ,
                    importance  = TRUE               ,
                    trControl   =  my.ctrl.opt.rf    ,
                    tuneLength  = 10                 ,  # test 10 parameters combinations
                    metric      = "Mean_Balanced_Accuracy"               )
save(my.rf.fit, file = paste0(Sys.Date(),"_ALL_RF_model.RData"))

# Check the learned model
my.rf.fit

##############################################################
##### 4. Check model accuracy
##############################################################

###### Testing accuracy
# create function to extract info from both SVM and RF
get_acc_metrics <- function(my.mod.fit, my.testing) {
  my.mod.preds  <- predict(my.mod.fit, my.testing)
  my.confus.mat <- confusionMatrix(my.mod.preds,my.testing$Species) # prediction, then ref
  
  my.bal.acc <- mean(my.confus.mat$byClass[,"Balanced Accuracy"] )
  my.sens    <- mean(my.confus.mat$byClass[,"Sensitivity"]       )
  my.spe     <- mean(my.confus.mat$byClass[,"Specificity"]       )
  my.prec    <- mean(my.confus.mat$byClass[,"Precision"]         )
  
  my.results         <- c(my.bal.acc, my.sens, my.spe, my.prec)
  names(my.results)  <- c("Balanced Accuracy", "Sensitivity", "Specificity", "Precision")
  
  return(my.results)
}

# extrat withheld testing data
my.testing.data <- iris[-my.training.idx,]

get_acc_metrics(my.svm.fit, my.testing.data)
# Balanced Accuracy       Sensitivity       Specificity         Precision 
#       0.9375000         0.9166667         0.9583333         0.9191176 

get_acc_metrics(my.rf.fit, my.testing.data)
# Balanced Accuracy       Sensitivity       Specificity         Precision 
#        0.9531250         0.9375000         0.9687500         0.9385621 


##############################################################
##### 5. Variable Importance for RF model
##############################################################

# feature importance on tree-based models (caret method)
my.rf.Imp    <- varImp(my.rf.fit    , useModel = TRUE , scale = TRUE)
my.rf.Imp$importance

plot(my.rf.Imp)

# feature importance on for RF only (native to the model)
my.rf.varimps.native <- my.rf.fit$finalModel$importance[,4:5]
my.rf.varimps.native
########################################################################################



########################################################################################
### Example #2
### Can we predict the disease based on gene expression and other covariates

data(ALL)

ALL
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 12625 features, 128 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: 01005 01010 ... LAL4 (128 total)
# varLabels: cod diagnosis ... date last seen (21 total)
# varMetadata: labelDescription
# featureData: none
# experimentData: use 'experimentData(object)'
# pubMedIds: 14684422 16243790 
# Annotation: hgu95av2

# Let's see if we can learn whether the ALL is B or T derived

##############################################################
##### 1. Data preparation/feature engineering
##############################################################

# response/outcome variable
ALL$BT   # many subtypes as well, so we need to simplify

my.all.type <- rep(NA,length(ALL$BT))
my.all.type[grep("T",ALL$BT)] <- "T"
my.all.type[grep("B",ALL$BT)] <- "B"

my.all.type <- factor(my.all.type)

# features/covariate
ALL$sex
ALL$age

# gene expression information
all.exprs         <- exprs(ALL)

# remove control affy probes
all.exprs         <- all.exprs[-grep("AFF",rownames(all.exprs)),]

# select only probes with high variamnce
my.var.probes     <- apply(all.exprs,1,var) > quantile(apply(all.exprs,1,var), 0.90)   # top 10% most variable genes
all.exprs.select  <- all.exprs[my.var.probes,]
dim(all.exprs.select)
## 1256  128

# transpose so that patients are lines, and genes/features are columns
my.expres.feat <- data.frame(t(all.exprs.select))

# collate all features to be used
my.ALL.features <- cbind(my.all.type, ALL$sex,ALL$age,my.expres.feat)
colnames(my.ALL.features)[1:3] <- c("ALL_TYPE","Sex","Age")

# Check feature matrix
my.ALL.features[,1:10]
# there are some NAs!!

my.ALL.features.noNA <- my.ALL.features[apply(is.na(my.ALL.features),1,sum) == 0,] # get non NA examples

##############################################################
#### 2. partition data in training and testing set
##############################################################
set.seed(123456789)                                                      # setting the seed helps make your code doing "random" sampling reproducible

my.training.idx <- createDataPartition(my.ALL.features.noNA$ALL_TYPE, p=0.67, list=FALSE) # 2/3 for training, 1/3 for testing
my.training.data <- my.ALL.features.noNA[my.training.idx,]
save(my.training.data, my.training.idx, file = paste0(Sys.Date(),"_ALL_Training_Data.RData"))


##############################################################
##### 3. train a RF model
##############################################################

set.seed(123456789)  # RF is based on random variable sampling

# use 10-fold cross-validation to build the model
my.ctrl.opt.rf <- trainControl(method          = "cv",
                               number          = 10,
                               allowParallel   = TRUE,
                               verbose         = F,
                               summaryFunction = multiClassSummary, # use multiclasssummary to get balanced accuracyt twoClassSummary for 2-class problems (T-ALL vs B-ALL)
                               classProbs      = TRUE)

# train model with caret train function
my.rf.fit <- train( ALL_TYPE ~ .,
                    data        =  my.training.data,
                    method      = "rf",
                    importance  = TRUE,
                    trControl   =  my.ctrl.opt.rf,
                    tuneLength  = 10                  ,  # test 10 parameters combinations
                    metric      = "Balanced_Accuracy",
                    na.action   = na.omit)              # deal with NAs
save(my.rf.fit, file = paste0(Sys.Date(),"_RF_model.RData"))

# Check the learned model
my.rf.fit

###############################
##### 3. Check model accuracy
###############################

###### Testing accuracy
# extract withheld testing data
my.testing.data <- my.ALL.features.noNA[-my.training.idx,]

# create function to extract info from RF
my.mod.preds  <- predict(my.rf.fit, my.testing.data)
my.confus.mat <- confusionMatrix(my.mod.preds,my.testing.data$ALL_TYPE) # prediction, then ref

my.bal.acc <- my.confus.mat$byClass["Balanced Accuracy"]
my.sens    <- my.confus.mat$byClass["Sensitivity"]      
my.spe     <- my.confus.mat$byClass["Specificity"]      
my.prec    <- my.confus.mat$byClass["Precision"]        

my.results         <- c(my.bal.acc, my.sens, my.spe, my.prec)
names(my.results)  <- c("Balanced Accuracy", "Sensitivity", "Specificity", "Precision")


##############################################################
##### 4. Variable Importance for RF model
##############################################################

# feature importance on tree-based models (caret method)
my.rf.Imp    <- varImp(my.rf.fit    , useModel = TRUE , scale = TRUE)
my.rf.Imp$importance

varImpPlot(my.rf.fit$finalModel)

# feature importance on for RF only (native to the model)
my.rf.varimps.native <- my.rf.fit$finalModel$importance[,3:4]
my.rf.varimps.native <- data.frame(my.rf.varimps.native)

my.sort <- sort(my.rf.varimps.native$MeanDecreaseAccuracy, decreasing = T, index.return = T)

my.rf.top.vars <- my.rf.varimps.native[my.sort$ix[1:20],]


##################################################
### get predictor gene identity from BIOMART
# open a connection to biomart
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# see available data in biomart
listMarts(mart, host="www.biomart.org", path="/biomart/martservice", port=80, includeHosts = FALSE, archive=FALSE, verbose = FALSE)
listAttributes(mart)

# choose needed attributes
my.attributes <- c("ensembl_gene_id", "external_gene_name", "affy_hg_u95av2")

my.top.genes <- unlist(strsplit(rownames(my.rf.top.vars),"X"))
my.top.genes <- my.top.genes[my.top.genes != ""]

# get annotations
gene_annot.biomart <- getBM(my.attributes, filters = "affy_hg_u95av2", values = my.top.genes, mart)
unique(gene_annot.biomart$external_gene_name)
#  "GGA2"       "HLA-DRA"    "HLA-DMB"    "HLA-DMA"    "CD9"        "HLA-DPA1"   "CD3D"       "CD3G"      
#  "AL645941.2" "BLNK"       "TCL1A"      "TRAT1"      "BX248406.2" "HLA-DQB1"   "PRKCQ"      "LCK"       
#  "SH2D1A"     "CD74"       "CD247"      "CR753846.2" "CD79B"      "AL662796.1"
########################################################################################


