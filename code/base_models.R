
# packages should already be loaded in the main script
# library(caret)
# options(java.parameters = "-Xmx1024m")
# library(bartMachine)
# library(kernlab)
# library(randomForest)
# library(fastAdaboost)
# library(xgboost)
# library(plyr)
# library(lightgbm)
# library(naivebayes)
# library(haven)
# library(readxl)
# library(pROC)

#  ---------------
# function to get summary statistics easily for each model
conf.mat.func <- function(model, data){
  if (model$method == 'bartMachine' ){ # for BART the levels are swapped so we take the second column
    preds <- factor(ifelse(predict(model, newdata=data, type='prob')[,2] >= model$ths, "yes", "no"), levels=c("yes", "no"))
  }else{
    preds <- factor(ifelse(predict(model, newdata=data, type='prob')[,1]  >= model$ths, "yes", "no"), levels=c("yes", "no")) 
  }
  conf.mat <- confusionMatrix(preds, reference=data[['cin3plus']])
  return(round(c(conf.mat$overall[1], conf.mat$byClass[1], conf.mat$byClass[2], conf.mat$byClass[3], conf.mat$byClass[4])*100, 2))
}


# function to find threshold that gives the closest specificity to that of the cytology only program
threshold.finder <- function(model, target = 0.733){
  indx.lower <- 0 
  indx.upper <- 1
  print(c(indx.lower, indx.upper))
  suppressWarnings({
    for (i in c(0.1, 0.01, 0.001)){
      ths.temp <- thresholder(model, seq(indx.lower, indx.upper, i))
      indx.lower <- ths.temp$prob_threshold[max(which(ths.temp[["Specificity"]] <= target))]
      indx.upper <- ths.temp$prob_threshold[min(which(ths.temp[["Specificity"]] >= target))]
      print(c(indx.lower, indx.upper))
    }
    ths.temp <- thresholder(model, seq(indx.lower, indx.upper, 0.0001))
  })
  return(ths.temp$prob_threshold[which.min(abs(ths.temp[["Specificity"]] - target))])
}

# keep the specificity the same as the model with only cytology --> this will result in different sens/referral than current but with same specificity
# cytology only: means refer after abnormal baseline cytology (>=Pap2) ( for Table 2)
cyt.referred <- factor(ifelse(trimws(training$primary.cyt) >= "PAP2", "yes", "no"), levels=c("yes", "no"))
confusionMatrix(cyt.referred, training$cin3plus) # PPV = 0.27458, NPV = 0.97061

cyt.referred.test <- factor(ifelse(trimws(testing$primary.cyt) >= "PAP2", "yes", "no"), levels=c("yes", "no"))
confusionMatrix(cyt.referred.test, testing$cin3plus) # PPV = 0.26928, NPV = 0.96850       

# triage by cytology and hpv16/18 (for Table 2)
cyt.1618.referred <- factor(ifelse(trimws(training$primary.cyt) >= "PAP2" & training$genotype %in% c("HPV16", "HPV18"), "yes", "no"), levels=c("yes", "no"))
confusionMatrix(cyt.1618.referred, training$cin3plus) # PPV = 0.41946, NPV = 0.94088         

cyt.1618.referred.test <- factor(ifelse(trimws(testing$primary.cyt) >= "PAP2" & testing$genotype %in% c("HPV16", "HPV18"), "yes", "no"), levels=c("yes", "no"))
confusionMatrix(cyt.1618.referred.test, testing$cin3plus) # PPV = 0.41066, NPV = 0.93916                 


# (a) logistic regression  ---------------
lm.model <- train(cin3plus~., data=training[,c(keep.cols, "cin3plus")],trControl = train_control, method = "glm", family='binomial')
summary(lm.model)
lm.model$ths <- threshold.finder(lm.model, 0.73310)

conf.mat.func(lm.model, training) # in-sample performance
conf.mat.func(lm.model, testing) # out-of-sample performance

# (b) logistic regression with all interactions   ---------------
lm.int.model <- train(cin3plus~ (age.group + genotype + test.type + cin3.history + recent.screen + ct.value)^2,
                      data=training[,c(keep.cols, "cin3plus")], trControl = train_control, method = "glm", family='binomial')
lm.int.model$ths <- threshold.finder(lm.int.model, 0.73310)

conf.mat.func(lm.int.model, training) # in-sample performance
conf.mat.func(lm.int.model, testing) # out-of-sample performance

# (c) support vector machine classifier  ---------------
svm.grid <- expand.grid(.C = c(0.5, 1, 2), .sigma = 0.1024532)
svm.model <- train(cin3plus~., data=training[,c(keep.cols, "cin3plus")], trControl = train_control,
                   method = "svmLinear2", tuneLength=3, preProc = c("center","scale"), prob.model=T)
svm.model$ths <- threshold.finder(svm.model, 0.73310)

conf.mat.func(svm.model, training) # in-sample performance
conf.mat.func(svm.model, testing) # out-of-sample performance

svm.roc <- roc(testing$cin3plus, predict(svm.model, newdata=testing[,keep.cols], type='prob')[,1])
svm.coords <- coords(svm.roc); svm.coords$Model <- paste0('SVM; AUC: ', round(svm.roc$auc, 3))

# (d) random forest classifier   ---------------
rf.model <- train(x =training[,keep.cols], y =training$cin3plus, trControl = train_control, method = "rf", tuneLength=3)
rf.model$ths <-  threshold.finder(rf.model, 0.73310)

conf.mat.func(rf.model, training) # in-sample performance
conf.mat.func(rf.model, testing) # out-of-sample performance

rf.roc <- roc(testing$cin3plus, predict(rf.model, newdata=testing[,keep.cols], type='prob')[,1])
rf.coords <- coords(rf.roc); rf.coords$Model <- paste0('Random Forest; AUC: ', round(rf.roc$auc, 3))

# (e) AdaBoost classifier  ---------------
ada.grid <- expand.grid(mfinal = (1:3)*3, maxdepth = c(1, 3),coeflearn = c("Breiman", "Freund", "Zhu"))
ada.model <- train(x=training[,keep.cols], y=training$cin3plus,trControl = train_control, method = "AdaBoost.M1", tuneGrid=ada.grid) ; ada.model
ada.model$ths <-  threshold.finder(ada.model, 0.73310) ; ada.model$ths

conf.mat.func(ada.model, training) # in-sample performance
conf.mat.func(ada.model, testing) # out-of-sample performance

# (f) XGBoost classifier  ---------------
xgb.model <- train(cin3plus~., data = training[,c(keep.cols, "cin3plus")], trControl = train_control, method = "xgbTree")
xgb.model$ths <- threshold.finder(xgb.model, 0.73310)

conf.mat.func(xgb.model, training) # in-sample performance
conf.mat.func(xgb.model, testing) # out-of-sample performance

# (g) Naive bayes classifier  ---------------
nb.model <- train(cin3plus~., data = training[, c(keep.cols, "cin3plus")], trControl = train_control, method = "naive_bayes") ; summary(nb.model)
nb.model$ths <- threshold.finder(nb.model, 0.73310) # probability threshold to match specificity of 0.77 

conf.mat.func(nb.model, training) # in-sample performance
conf.mat.func(nb.model, testing) # out-of-sample performance

# (h) BART (bayesian adaptive regression tree) classifier ---------------
# --> make sure to run the 'options(java.parameters = ...)' line before loading the package
set_bart_machine_num_cores(4)

bart.model <- train(cin3plus~., data = training[,c(keep.cols,'cin3plus')], prob_rule_class = 0.11, mem_cache_for_speed = FALSE,
                    trControl = train_control, method = "bartMachine", tuneGrid=data.frame(num_trees=10,k=2,alpha=0.1,beta=0.1, nu=1), serialize = TRUE)
summary(bart.model)
bart.model$ths <- 0.1134 #threshold.finder(bart.model, 0.73310)
conf.mat.func(bart.model, training) # in-sample performance
conf.mat.func(bart.model, testing) # out-of-sample performance

# save models 
save(ada.model, bart.model, lm.model, lm.int.model, nb.model, rf.model, svm.model, xgb.model, file="base_models.RData")

# ROC curve calculations:
# save roc in a dataframe for plotting 
load("base_models.RData")
lm.roc <- roc(testing$cin3plus, predict(lm.model, newdata=testing[,c(keep.cols, "ct.value")], type='prob')[,1])
lm.int.roc <- roc(testing$cin3plus, predict(lm.int.model, newdata=testing[,c(keep.cols, "ct.value")], type='prob')[,1])
ada.roc <- roc(testing$cin3plus, predict(ada.model, newdata=testing[,c(keep.cols, "ct.value")], type='prob')[,1])
bart.roc <- roc(testing$cin3plus, predict(bart.model, newdata=testing[,c(keep.cols, "ct.value")], type='prob')[,1])
nb.roc <- roc(testing$cin3plus, predict(nb.model, newdata=testing[,c(keep.cols, "ct.value")], type='prob')[,1])
xgb.roc <- roc(testing$cin3plus, predict(xgb.model, newdata=testing[,c(keep.cols, "ct.value")], type='prob')[,1])
rf.roc <- roc(testing$cin3plus, predict(rf.model, newdata=testing[,c(keep.cols, "ct.value")], type='prob')[,1])
svm.roc <- roc(testing$cin3plus, predict(svm.model, newdata=testing[,c(keep.cols, "ct.value")], type='prob')[,1])

lm.coords <- coords(lm.roc) ; lm.coords$Model <- paste0('Logistic regression; AUC: ', round(lm.roc$auc, 3))
lm.int.coords <- coords(lm.int.roc) ; lm.int.coords$Model <- paste0('Logistic regression (interactions); AUC: ', round(lm.int.roc$auc, 3))
ada.coords <- coords(ada.roc) ; ada.coords$Model <- paste0('adaBoost; AUC: ', round(ada.roc$auc, 3))
bart.coords <- coords(bart.roc) ; bart.coords$Model <- paste0('BART; AUC: ', round(bart.roc$auc, 3))
nb.coords <- coords(nb.roc) ; nb.coords$Model <- paste0('Naive Bayes; AUC: ', round(nb.roc$auc, 3))
xgb.coords <- coords(xgb.roc) ; xgb.coords$Model <- paste0("XGBoost; AUC: ", round(xgb.roc$auc, 3))
rf.coords <- coords(rf.roc) ; rf.coords$Model <- paste0('Random Forest; AUC: ', round(rf.roc$auc, 3))
svm.coords <- coords(svm.roc) ; svm.coords$Model <- paste0('SVM; AUC: ', round(svm.roc$auc, 3))

roc.coords <- rbind(lm.coords, lm.int.coords, rf.coords, ada.coords, xgb.coords, svm.coords, nb.coords, bart.coords)
save(roc.coords, file="roc.coords.RData")

