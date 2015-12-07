############################################
## Use H2O to create a random forest
##  against the entire data set in 
##  just a couple minutes
##
## This is a starter script, using defaults
##  so it can be improved. And RF may not be
##  the best algorithm for this problem. But
##  the script shows that it can be done in R
##  fairly quickly, in fact. And it will scale
##  well to adding many more columns.
##
## To better fit MAE, the log of the 
##  target has been used: log1p/expm1
##
## The final step blends with the Marshall-Palmer
##  benchmark 50/50.
##
## The API for h2o.randomForest is shown 
##  at the bottom of the script
###########################################

library(h2o)
library(data.table)
library(Metrics)
h2o.init(nthreads=-1)

calc_accelerate <- function(ref, minutes_past) {
    sort_min_index = order(minutes_past)
    diff <-
    c(ref[1:length(ref) - 1]) -
    c(ref[2:length(ref)])

    minutes_past <- minutes_past[sort_min_index]
    ref <- ref[sort_min_index]
  
  # calculate the length of time for which each reflectivity value is valid
    valid_time <-
        c(minutes_past[-length(minutes_past)], 60) -
        c(0, minutes_past[-length(minutes_past)])
    valid_time = valid_time / 60
    return(mean(diff / valid_time, na.rm=TRUE))
}

ref_diff <- function(ref) {
    diff <-
    c(ref[1:length(ref) - 1]) -
    c(ref[2:length(ref)])
    return(mean(diff, na.rm=TRUE))
}

mpalmer <- function(ref, minutes_past) {
  # order reflectivity values and minutes_past
  sort_min_index = order(minutes_past)
  minutes_past <- minutes_past[sort_min_index]
  ref <- ref[sort_min_index]
  
  # calculate the length of time for which each reflectivity value is valid
  valid_time <-
    c(minutes_past[-length(minutes_past)], 60) -
    c(0, minutes_past[-length(minutes_past)])
  valid_time = valid_time / 60

  # calculate hourly rain rates using marshall-palmer weighted by valid times
  return(sum(((10^(ref/10))/200)^0.625*valid_time, na.rm=TRUE))
}

rate_kdp <- function(Kdp, minutes_past) {
  
  # order Kdp values and minutes_past
  sort_min_index = order(minutes_past)
  minutes_past <- minutes_past[sort_min_index]
  Kdp <- Kdp[sort_min_index]
  
  # calculate the length of time for which each Kdp value is valid
  valid_time <-
    c(minutes_past[-length(minutes_past)], 60) -
    c(0, minutes_past[-length(minutes_past)])
  valid_time = valid_time / 60
  # calculate hourly rain rates using S and Z formula weighted by valid times
  return(sum((sign(Kdp)*(4.06)*(abs(Kdp) **.0866)*valid_time), na.rm=TRUE))
}


print(paste("reading training file:",Sys.time()))
train<-fread("../train.csv",select=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,24))

#Cut off outliers of Expected >= 69
train <- subset(train, Expected < 69)

summary(train)

# train$dt <- time_difference(train$minutes_past)
# train$mp5 <- marshall_palmer(train$Ref_5x5_50th)
# train$mp9 <- marshall_palmer(train$Ref_5x5_90th)

#Cut off Ref values < 0
train$Ref_5x5_50th[which(train$Ref_5x5_10th < 0)] <- NA
train$Ref_5x5_50th[which(train$Ref_5x5_50th < 0)] <- NA
train$Ref_5x5_90th[which(train$Ref_5x5_90th < 0)] <- NA
train$RefComposite[which(train$RefComposite < 0)] <- NA
train$RefComposite_5x5_10th[which(train$RefComposite_5x5_10th < 0)] <- NA
train$RefComposite_5x5_50th[which(train$RefComposite_5x5_50th < 0)] <- NA
train$RefComposite_5x5_90th[which(train$RefComposite_5x5_90th < 0)] <- NA
train$Ref[which(train$Ref < 0)] <- NA

cor(train, use = "pairwise.complete.obs")
summary(train)

trainHex<-as.h2o(train[,.(
    dist   = mean(radardist_km, na.rm = T),

    refArea1   = mean(Ref_5x5_10th, na.rm = T),
    varRefArea1   = var(Ref_5x5_10th, na.rm = T),
    refArea5   = mean(Ref_5x5_50th, na.rm = T),
    varRefArea5   = var(Ref_5x5_50th, na.rm = T),
    refArea9  = mean(Ref_5x5_90th, na.rm = T),
    varRefArea9  = var(Ref_5x5_90th, na.rm = T),

    meanRefcomp = mean(RefComposite,na.rm = T),
    varRefcomp = var(RefComposite,na.rm = T),
    meanRefcomp1 = mean(RefComposite_5x5_10th,na.rm = T),
    #varRefcomp1 = var(RefComposite_5x5_10th,na.rm = T),
    meanRefcomp5 = mean(RefComposite_5x5_50th,na.rm = T),
    varRefcomp5 = var(RefComposite_5x5_50th,na.rm = T),
    meanRefcomp9 = mean(RefComposite_5x5_90th,na.rm = T),
    varRefcomp9 = var(RefComposite_5x5_90th,na.rm = T),
    maxRefcomp = max(RefComposite, na.rm = T),
    minRefcomp = min(RefComposite, na.rm = T),

    rhoHV = mean(RhoHV, na.rm = T),
    rhoHV1 = mean(RhoHV_5x5_10th, na.rm = T),
    rhoHV5 = mean(RhoHV_5x5_50th, na.rm = T),
    rhoHV9 = mean(RhoHV_5x5_90th, na.rm = T),

    zdr   = mean(Zdr, na.rm = T),
    zdr1   = mean(Zdr_5x5_10th, na.rm = T),
    zdr5   = mean(Zdr_5x5_50th, na.rm = T),
    zdr9   = mean(Zdr_5x5_90th, na.rm = T),

    #kdp   = rate_kdp(Kdp, minutes_past),
    maxKdp = max(Kdp, na.rm = T),
    minKdp = min(Kdp, na.rm = T),

    target = log1p(mean(Expected)),
    meanRef = mean(Ref,na.rm = T),
    varRef = var(Ref, na.rm = T),
    sumRef = sum(Ref,na.rm = T),
    maxRef = max(Ref, na.rm = T),
    minRef = min(Ref, na.rm = T),

    yy1 = mean(Ref,na.rm = T) / mean(radardist_km, na.rm = T),
    #yy2 = mean(Ref_5x5_90th, na.rm = T) - mean(Ref_5x5_50th, na.rm = T),
    yy3 = mean(Ref_5x5_90th,na.rm = T) / mean(radardist_km, na.rm = T),
    yy4 = mean(Ref_5x5_50th,na.rm = T) / mean(radardist_km, na.rm = T),
    yy5 = mean(RefComposite,na.rm = T) / mean(radardist_km, na.rm = T),
    yy6 = mean(RefComposite_5x5_90th,na.rm = T) / mean(radardist_km, na.rm = T),
    yy7 = mean(RefComposite_5x5_50th,na.rm = T) / mean(radardist_km, na.rm = T),
    yy8 = ref_diff(Ref),
    yy9 = ref_diff(Ref_5x5_50th),
    yy10 = ref_diff(Ref_5x5_90th),
    yy11 = ref_diff(RefComposite),
    yy12 = ref_diff(RefComposite_5x5_50th),
    yy13 = ref_diff(RefComposite_5x5_90th),
    records = .N,
    naCounts = sum(is.na(Ref))
    #mp50 =  mpalmer(RefComposite_5x5_50th, minutes_past),
    #mp90 =  mpalmer(RefComposite_5x5_90th, minutes_past),
    #mp = log1p(mpalmer(Ref, minutes_past))
    ),Id][records>naCounts,],destination_frame="train.hex")
    
    summary(trainHex)

rfHex<-h2o.randomForest(x=c("dist", "refArea1", "varRefArea1", "refArea5", "varRefArea5", "refArea9", "varRefArea9",
                            "meanRefcomp", "varRefcomp","meanRefcomp1","meanRefcomp5", "varRefcomp5","meanRefcomp9", "varRefcomp9",
                            "maxRefcomp", "minRefcomp",
                            "rhoHV","rhoHV1","rhoHV5","rhoHV9",
                            "zdr", "zdr1", "zdr5", "zdr9",
                            "maxKdp", "minKdp",
                            #"mp50","mp90","mp",
                            "meanRef", "varRef", "sumRef", "records","naCounts", 
                            "maxRef", "minRef",
                            "yy1", "yy3", "yy4", "yy5", "yy6", "yy7","yy8","yy9",
                            "yy10","yy11","yy12","yy13"
                        ),
    y="target",training_frame=trainHex,model_id="rfStarter.hex", ntrees=500, sample_rate = 0.7)
rfHex
h2o.varimp(rfHex)
rm(train)

test<-fread("../test.csv",select=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))

# test$dt <- time_difference(te_raw$minutes_past)
# test$mp5 <- marshall_palmer(te_raw$Ref_5x5_50th)
# test$mp9 <- marshall_palmer(te_raw$Ref_5x5_90th)

#Cut off Ref values < 0
test$Ref_5x5_50th[which(test$Ref_5x5_50th < 0)] <- NA
test$Ref_5x5_90th[which(test$Ref_5x5_90th < 0)] <- NA
test$RefComposite[which(test$RefComposite < 0)] <- NA
test$RefComposite_5x5_10th[which(test$RefComposite_5x5_10th < 0)] <- NA
test$RefComposite_5x5_50th[which(test$RefComposite_5x5_50th < 0)] <- NA
test$RefComposite_5x5_90th[which(test$RefComposite_5x5_90th < 0)] <- NA
test$Ref[which(test$Ref < 0)] <- NA


testHex<-as.h2o(test[,.(
    dist   = mean(radardist_km, na.rm = T),

    refArea1   = mean(Ref_5x5_10th, na.rm = T),
    varRefArea1   = var(Ref_5x5_10th, na.rm = T),
    refArea5   = mean(Ref_5x5_50th, na.rm = T),
    varRefArea5   = var(Ref_5x5_50th, na.rm = T),
    refArea9  = mean(Ref_5x5_90th, na.rm = T),
    varRefArea9  = var(Ref_5x5_90th, na.rm = T),

    meanRefcomp = mean(RefComposite,na.rm=T),
    varRefcomp = var(RefComposite,na.rm = T),
    meanRefcomp1 = mean(RefComposite_5x5_10th,na.rm = T),
    #varRefcomp1 = var(RefComposite_5x5_10th,na.rm = T),
    meanRefcomp5 = mean(RefComposite_5x5_50th,na.rm = T),
    varRefcomp5 = var(RefComposite_5x5_50th,na.rm = T),
    meanRefcomp9 = mean(RefComposite_5x5_90th,na.rm = T),
    varRefcomp9 = var(RefComposite_5x5_90th,na.rm = T),
    maxRefcomp = max(RefComposite, na.rm = T),
    minRefcomp = min(RefComposite, na.rm = T),

    rhoHV = mean(RhoHV, na.rm = T),
    rhoHV1 = mean(RhoHV_5x5_10th, na.rm = T),
    rhoHV5 = mean(RhoHV_5x5_50th, na.rm = T),
    rhoHV9 = mean(RhoHV_5x5_90th, na.rm = T),

    zdr   = mean(Zdr, na.rm = T),
    zdr1   = mean(Zdr_5x5_10th, na.rm = T),
    zdr5   = mean(Zdr_5x5_50th, na.rm = T),
    zdr9   = mean(Zdr_5x5_90th, na.rm = T),

    #kdp   = rate_kdp(Kdp, minutes_past),
    maxKdp = max(Kdp, na.rm = T),
    minKdp = min(Kdp, na.rm = T),
    
    meanRef = mean(Ref,na.rm=T),
    varRef = var(Ref, na.rm=T),
    sumRef = sum(Ref,na.rm=T),
    maxRef = max(Ref, na.rm = T),
    minRef = min(Ref, na.rm = T),

    yy1 = mean(Ref,na.rm = T) / mean(radardist_km, na.rm = T),
    #yy2 = mean(Ref_5x5_90th, na.rm = T) - mean(Ref_5x5_50th, na.rm = T),
    yy3 = mean(Ref_5x5_90th,na.rm = T) / mean(radardist_km, na.rm = T),
    yy4 = mean(Ref_5x5_50th,na.rm = T) / mean(radardist_km, na.rm = T),
    yy5 = mean(RefComposite,na.rm = T) / mean(radardist_km, na.rm = T),
    yy6 = mean(RefComposite_5x5_90th,na.rm = T) / mean(radardist_km, na.rm = T),
    yy7 = mean(RefComposite_5x5_50th,na.rm = T) / mean(radardist_km, na.rm = T),
    yy8 = ref_diff(Ref),
    yy9 = ref_diff(Ref_5x5_50th),
    yy10 = ref_diff(Ref_5x5_90th),
    yy11 = ref_diff(RefComposite),
    yy12 = ref_diff(RefComposite_5x5_50th),
    yy13 = ref_diff(RefComposite_5x5_90th),

    records = .N,
    naCounts = sum(is.na(Ref))
    #mp50 =  mpalmer(RefComposite_5x5_50th, minutes_past),
    #mp90 =  mpalmer(RefComposite_5x5_90th, minutes_past),
    #mp = log1p(mpalmer(Ref, minutes_past))
    ),Id],destination_frame="test.hex")
    
    summary(testHex)

submission<-fread("../sample_solution.csv")
predictions<-as.data.frame(h2o.predict(rfHex,testHex))
submission$Expected<- 0.9 * expm1(predictions$predict) + 0.1 * submission$Expected

#convert expected values to 0.01in values
submission$Expected <- round(submission$Expected / 0.254) * 0.254

summary(submission)
write.csv(submission,"rfv3cn3.csv",row.names=F)


####################################################################################
## Appendix: h2o.randomForest API
####################################################################################
##h2o.randomForest(x, y, training_frame, model_id, validation_frame, checkpoint,
##  mtries = -1, sample_rate = 0.632, build_tree_one_node = FALSE,
##  ntrees = 50, max_depth = 20, min_rows = 1, nbins = 20,
##  nbins_cats = 1024, binomial_double_trees = FALSE,
##  balance_classes = FALSE, max_after_balance_size = 5, seed,
##  offset_column = NULL, weights_column = NULL, nfolds = 0,
##  fold_column = NULL, fold_assignment = c("AUTO", "Random", "Modulo"),
##  keep_cross_validation_predictions = FALSE, ...)
## Arguments

## x	
## A vector containing the names or indices of the predictor variables to use in building the GBM model.

## y	
## The name or index of the response variable. If the data does not contain a header, this is the column index number starting at 1, and increasing from left to right. (The response must be either an integer or a categorical variable).

## training_frame	
## An H2OFrame object containing the variables in the model.

## model_id	
## (Optional) The unique id assigned to the resulting model. If none is given, an id will automatically be generated.

## validation_frame	
## An H2OFrame object containing the variables in the model.

## checkpoint	
## "Model checkpoint (either key or H2ODeepLearningModel) to resume training with."

## mtries	
## Number of variables randomly sampled as candidates at each split. If set to -1, defaults to sqrtp for classification, and p/3 for regression, where p is the number of predictors.

## sample_rate	
## Sample rate, from 0 to 1.0.   (edit: row sampling, per tree)

## build_tree_one_node	
## Run on one node only; no network overhead but fewer cpus used. Suitable for small datasets.

## ntrees	
## A nonnegative integer that determines the number of trees to grow.

## max_depth	
## Maximum depth to grow the tree.

## min_rows	
## Minimum number of rows to assign to teminal nodes.

## nbins	
## For numerical columns (real/int), build a histogram of this many bins, then split at the best point.

## nbins_cats	
## For categorical columns (enum), build a histogram of this many bins, then split at the best point. Higher values can lead to more overfitting.

## binomial_double_trees	
## For binary classification: Build 2x as many trees (one per class) - can lead to higher accuracy.

## balance_classes	
## logical, indicates whether or not to balance training data class counts via over/under-sampling (for imbalanced data)

## max_after_balance_size	
## Maximum relative size of the training data after balancing class counts (can be less than 1.0)

## seed	
## Seed for random numbers (affects sampling) - Note: only reproducible when running single threaded

## offset_column	
## Specify the offset column.

## weights_column	
## Specify the weights column.

## nfolds	
## (Optional) Number of folds for cross-validation. If nfolds >= 2, then validation must remain empty.

## fold_column	
## (Optional) Column with cross-validation fold index assignment per observation

## fold_assignment	
## Cross-validation fold assignment scheme, if fold_column is not specified Must be "AUTO", "Random" or "Modulo"

## keep_cross_validation_predictions	
## Whether to keep the predictions of the cross-validation models
                
