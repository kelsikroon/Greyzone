# greyzone models
setwd("~/Projects/Grey zone/Analysis")
source("greyzone_functions.R")
#load("cyt_models.RData")
#load("base_models.RData")
keep.cols.cyt <- c(keep.cols, "primary.cyt")

# (a1) Linear model 
lm.results <-  make.greyzone(lm.model, lm.model.cyt, training, plot=F)
lm.results.test <- test.greyzone(lm.model, lm.model.cyt, testing, lm.results$greyzone[1], lm.results$greyzone[2], lm.results$cyt.model.ths)

# (a2) linear model with all interactions 
lm.interactions.results <- make.greyzone(lm.int.model, lm.int.model.cyt, training, plot=F)
lm.interactions.results.test <- test.greyzone(lm.int.model, lm.int.model.cyt, testing, lm.interactions.results$greyzone[1], lm.interactions.results$greyzone[2], lm.interactions.results$cyt.model.ths)

# (b) AdaBoost classifier
ada.results <- make.greyzone(ada.model, ada.model.cyt, training, plot=F)
ada.results.test <- test.greyzone(ada.model, ada.model.cyt, testing, ada.results$greyzone[1], ada.results$greyzone[2], ada.results$cyt.model.ths)

# (c) BART 
bart.results <- make.greyzone(bart.model, bart.model.cyt, training, plot=F)
bart.results.test <- test.greyzone(bart.model, bart.model.cyt, testing, bart.results$greyzone[1], bart.results$greyzone[2], bart.results$cyt.model.ths)

# (d) Random Forest
rf.results <- make.greyzone(rf.model, rf.model.cyt, training, plot=F)
rf.results.test <- test.greyzone(rf.model, rf.model.cyt, testing, rf.results$greyzone[1], rf.results$greyzone[2], rf.results$cyt.model.ths)

# (e) Naive Bayes 
# --> added so code works before running NB models 
testing$nb.pred <- rep(0, 9935)
testing$nb.pred.cyt <- rep(0, 9935)

# --> NB produces very small predictions, so we use a logit model on the predictions to correct them
training$nb.pred <- predict(nb.model, training[,keep.cols], type='prob')[,1]
training$nb.pred.cyt <- predict(nb.model.cyt, training[,keep.cols.cyt], type='prob')[,1]
nb.logit <-  train(cin3plus~nb.pred, data = training, trControl = train_control, method = "glm", family='binomial')
nb.cyt.logit <-  train(cin3plus~nb.pred.cyt, data = training, trControl = train_control, method = "glm", family='binomial')
nb.logit$ths  <- threshold.finder(nb.logit, 0.733)
nb.cyt.logit$ths <- threshold.finder(nb.cyt.logit, 0.733)
nb.results <- make.greyzone(nb.logit, nb.cyt.logit, training, plot=F)

testing$nb.pred <- predict(nb.model, testing[,keep.cols], type='prob')[,1]
testing$nb.pred.cyt <- predict(nb.model.cyt, testing[,keep.cols.cyt], type='prob')[,1]
nb.results.test <- test.greyzone(nb.logit, nb.cyt.logit, testing, nb.results$greyzone[1], nb.results$greyzone[2], nb.results$cyt.model.ths)


# (f) XGBoost
xgb.results <- make.greyzone(xgb.model, xgb.model.cyt, training, plot=F)
xgb.results.test <- test.greyzone(xgb.model, xgb.model.cyt, testing, xgb.results$greyzone[1], xgb.results$greyzone[2], xgb.results$cyt.model.ths)

save(ada.results, bart.results, lm.results, lm.interactions.results, nb.results, rf.results, xgb.results, file="greyzone_models.RData")
save(ada.results.test, bart.results.test, lm.results.test, lm.interactions.results.test, nb.results.test, rf.results.test, xgb.results.test, file="greyzone_models_test.RData")

# 
# rm(lm.model, lm.model.cyt, nb.model, nb.model.cyt, nb.logit, nb.cyt.logit, bart.model, bart.model.cyt,
#    rf.model, rf.model.cyt, ada.model, ada.model.cyt, lm.int.model, lm.int.model.cyt, xgb.model, xgb.model.cyt)