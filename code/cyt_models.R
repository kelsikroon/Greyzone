
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

keep.cols.cyt <- c(keep.cols, "primary.cyt")

# (a) logistic regression  ---------------
lm.model.cyt <- train(cin3plus~., data=training[,c(keep.cols.cyt, "cin3plus")],trControl = train_control, method = "glm", family='binomial')
summary(lm.model.cyt)
lm.model.cyt$ths <- threshold.finder(lm.model.cyt, 0.73310)

conf.mat.func(lm.model.cyt, training) # in-sample performance
conf.mat.func(lm.model.cyt, testing) # out-of-sample performance

# (b) logistic regression with all interactions   ---------------
lm.int.model.cyt <- train(cin3plus~ (age.group + genotype + test.type + cin3.history + recent.screen + ct.value + primary.cyt)^2,
                      data=training[,c(keep.cols.cyt, "cin3plus")], trControl = train_control, method = "glm", family='binomial')
lm.int.model.cyt$ths <- threshold.finder(lm.int.model.cyt, 0.73310)

conf.mat.func(lm.int.model.cyt, training) # in-sample performance
conf.mat.func(lm.int.model.cyt, testing) # out-of-sample performance

# (c) support vector machine classifier  ---------------
svm.grid <- expand.grid(.C = c(0.5, 1, 2), .sigma = 0.1024532)
svm.model.cyt <- train(cin3plus~., data=training[,c(keep.cols.cyt, "cin3plus")], trControl = train_control,
                   method = "svmLinear2", tuneLength=3, preProc = c("center","scale"), prob.model.cyt=T)
svm.model.cyt$ths <- threshold.finder(svm.model.cyt, 0.73310)

conf.mat.func(svm.model.cyt, training) # in-sample performance
conf.mat.func(svm.model.cyt, testing) # out-of-sample performance

# (d) random forest classifier   ---------------
rf.model.cyt <- train(x =training[,keep.cols.cyt], y =training$cin3plus, trControl = train_control, method = "rf", tuneLength=3)
rf.model.cyt$ths <-  threshold.finder(rf.model.cyt, 0.73310)

conf.mat.func(rf.model.cyt, training) # in-sample performance
conf.mat.func(rf.model.cyt, testing) # out-of-sample performance

# (e) AdaBoost classifier  ---------------
ada.grid <- expand.grid(mfinal = (1:3)*3, maxdepth = c(1, 3),coeflearn = c("Breiman", "Freund", "Zhu"))
ada.model.cyt <- train(x=training[,keep.cols.cyt], y=training$cin3plus,trControl = train_control, method = "AdaBoost.M1", tuneGrid=ada.grid) ; ada.model.cyt
ada.model.cyt$ths <-  threshold.finder(ada.model.cyt, 0.73310) ; ada.model.cyt$ths

conf.mat.func(ada.model.cyt, training) # in-sample performance
conf.mat.func(ada.model.cyt, testing) # out-of-sample performance

# (f) XGBoost classifier  ---------------
xgb.model.cyt <- train(cin3plus~., data = training[,c(keep.cols.cyt, "cin3plus")], trControl = train_control, method = "xgbTree")
xgb.model.cyt$ths <- threshold.finder(xgb.model.cyt, 0.73310)

conf.mat.func(xgb.model.cyt, training) # in-sample performance
conf.mat.func(xgb.model.cyt, testing) # out-of-sample performance

# (g) Naive bayes classifier  ---------------
nb.model.cyt <- train(cin3plus~., data = training[, c(keep.cols.cyt, "cin3plus")], trControl = train_control, method = "naive_bayes") ; summary(nb.model.cyt)
nb.model.cyt$ths <- threshold.finder(nb.model.cyt, 0.73310) # probability threshold to match specificity of 0.77 

conf.mat.func(nb.model.cyt, training) # in-sample performance
conf.mat.func(nb.model.cyt, testing) # out-of-sample performance

# (h) BART (bayesian adaptive regression tree) classifier ---------------
# --> make sure to run the 'options(java.parameters = ...)' line before loading the package
set_bart_machine_num_cores(4)

bart.model.cyt <- train(cin3plus~., data = training[,c(keep.cols.cyt,'cin3plus')], prob_rule_class = 0.11, mem_cache_for_speed = FALSE,
                    trControl = train_control, method = "bartMachine", tuneGrid=data.frame(num_trees=10,k=2,alpha=0.1,beta=0.1, nu=1), serialize = TRUE)
summary(bart.model.cyt)
bart.model.cyt$ths <- 0.1134 #threshold.finder(bart.model.cyt, 0.73310)
conf.mat.func(bart.model.cyt, training) # in-sample performance
conf.mat.func(bart.model.cyt, testing) # out-of-sample performance

# save models 
save(ada.model.cyt, bart.model.cyt, lm.model.cyt, lm.int.model.cyt, nb.model.cyt, rf.model.cyt, svm.model.cyt, xgb.model.cyt, file="cyt_models.RData")

# ROC curve calculations:
# save roc in a dataframe for plotting 
#load("cyt_models.RData")
lm.roc.cyt <- roc(testing$cin3plus, predict(lm.model.cyt, newdata=testing[,c(keep.cols.cyt, "ct.value")], type='prob')[,1])
lm.int.roc.cyt <- roc(testing$cin3plus, predict(lm.int.model.cyt, newdata=testing[,c(keep.cols.cyt, "ct.value")], type='prob')[,1])
ada.roc.cyt <- roc(testing$cin3plus, predict(ada.model.cyt, newdata=testing[,c(keep.cols.cyt, "ct.value")], type='prob')[,1])
bart.roc.cyt <- roc(testing$cin3plus, predict(bart.model.cyt, newdata=testing[,c(keep.cols.cyt, "ct.value")], type='prob')[,1])
nb.roc.cyt <- roc(testing$cin3plus, predict(nb.model.cyt, newdata=testing[,c(keep.cols.cyt, "ct.value")], type='prob')[,1])
xgb.roc.cyt <- roc(testing$cin3plus, predict(xgb.model.cyt, newdata=testing[,c(keep.cols.cyt, "ct.value")], type='prob')[,1])
rf.roc.cyt <- roc(testing$cin3plus, predict(rf.model.cyt, newdata=testing[,c(keep.cols.cyt, "ct.value")], type='prob')[,1])
svm.roc.cyt <- roc(testing$cin3plus, predict(svm.model.cyt, newdata=testing[,c(keep.cols.cyt, "ct.value")], type='prob')[,1])

lm.coords <- coords(lm.roc.cyt) ; lm.coords$Model <- paste0('Logistic regression; AUC: ', round(lm.roc.cyt$auc, 3))
lm.int.coords <- coords(lm.int.roc.cyt) ; lm.int.coords$Model <- paste0('Logistic regression (interactions); AUC: ', round(lm.int.roc.cyt$auc, 3))
ada.coords <- coords(ada.roc.cyt) ; ada.coords$Model <- paste0('adaBoost; AUC: ', round(ada.roc.cyt$auc, 3))
bart.coords <- coords(bart.roc.cyt) ; bart.coords$Model <- paste0('BART; AUC: ', round(bart.roc.cyt$auc, 3))
nb.coords <- coords(nb.roc.cyt) ; nb.coords$Model <- paste0('Naive Bayes; AUC: ', round(nb.roc.cyt$auc, 3))
xgb.coords <- coords(xgb.roc.cyt) ; xgb.coords$Model <- paste0("XGBoost; AUC: ", round(xgb.roc.cyt$auc, 3))
rf.coords <- coords(rf.roc.cyt) ; rf.coords$Model <- paste0('Random Forest; AUC: ', round(rf.roc.cyt$auc, 3))
svm.coords <- coords(svm.roc.cyt) ; svm.coords$Model <- paste0('SVM; AUC: ', round(svm.roc.cyt$auc, 3))

roc.coords.cyt <- rbind(lm.coords, lm.int.coords, ada.coords, xgb.coords, rf.coords, svm.coords, nb.coords, bart.coords)
save(roc.coords.cyt, file="roc.coords.cyt.RData")
