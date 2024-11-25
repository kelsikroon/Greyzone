# greyzone manuscript analysis 

setwd("~/Projects/Grey zone/Analysis")
library(readxl)
library(stringr)
library(lubridate)
library(plyr)
library(dplyr)
library(haven)
library(purrr)
library(caret)
library(pROC)
library(ggplot2)
library(cowplot)
options(java.parameters = "-Xmx1024m")
library(bartMachine) # for bart classifier
#library(kernlab)
library(randomForest) # for random forest classifier
#install.packages("adabag")
# devtools::install_github("souravc83/fastAdaboost") 
library(adabag)
library(xgboost) # for xgboost classifier
library(naivebayes) # for naive bayes classifier
library(e1071) # for support vector machine classifier
library(rpart)
library(DescTools)

set.seed(13062001)


# 1. data cleaning and preparation (make test/train sets) ---------------
source("data_prep.R")

# 2. patient characteristics  ---------------
# check if they're same as what Robine reported in her manuscript
rbind(dat %>% group_by(test.type) %>% summarise(n=n()) %>% mutate(freq = n / sum(n)*100, feature = 'test.type') %>% rename("value" = "test.type") %>% relocate(feature),
      dat %>% group_by(age.group) %>% summarise(n=n()) %>% mutate(freq = n / sum(n)*100, feature = 'age.group')  %>% rename("value" = "age.group") %>% relocate(feature),
      dat %>% group_by(genotype) %>% summarise(n=n()) %>% mutate(freq = n / sum(n)*100, feature = 'genotype') %>% rename("value" = "genotype") %>% relocate(feature),
      dat %>% group_by(ct.group) %>% summarise(n=n()) %>% mutate(freq = n / sum(n)*100, feature = 'ct.group') %>% rename("value" = "ct.group") %>% relocate(feature),
      dat %>% group_by(recent.screen) %>% summarise(n=n()) %>% mutate(freq = n / sum(n)*100, feature = 'recent.screen') %>% rename("value" = "recent.screen") %>% relocate(feature),
      dat %>% group_by(cin3.history) %>% summarise(n=n()) %>% mutate(freq = n / sum(n)*100, feature = 'cin3.history') %>% rename("value" = "cin3.history") %>% relocate(feature))


summary.func <- function(data, columns, outcome){
  summary.df <- data %>% subset(cin3plus == outcome) %>% group_by(cin3plus) %>%  summarise(n=n()) %>% mutate(freq = round(n / dim(data)[1]*100, 2), feature='total') %>% rename("value" = "cin3plus") %>% relocate(feature)
  for (c in columns){
    summary.df <- rbind(summary.df, data %>% subset(cin3plus == outcome) %>% group_by_at(c) %>% summarise(n=n()) %>% mutate(freq = round(n / sum(n)*100, 2), feature=c) %>% rename("value" = c) %>% relocate(feature))
  }
  return(summary.df)
}

# Table 1. patient characteristics 
rbind(c("", NA, rep(c("training", "testing"), each=4)), 
      c("", NA, "CIN3+", "CIN3+",  "<CIN3", "<CIN3","CIN3+", "CIN3+","<CIN3", "<CIN3"), 
      cbind(summary.func(training, c("test.type", "age.group", "genotype", "ct.group", "recent.screen", "cin3.history"), "yes"),
            summary.func(training, c("test.type", "age.group", "genotype", "ct.group", "recent.screen", "cin3.history"), "no")[,c(3,4)],
            summary.func(testing, c("test.type", "age.group", "genotype", "ct.group", "recent.screen", "cin3.history"), "yes")[,c(3,4)],
            summary.func(testing, c("test.type", "age.group", "genotype", "ct.group", "recent.screen", "cin3.history"), "no")[,c(3,4)])) %>% write.csv(., file='table1_sep2024.csv')


# 3. models without baseline cytology  ---------------

# (3.1) fit models
# source("base_models.R") # takes a while to run --> so easier to just load the saved models
load("base_models.RData")

# (3.2) ROC curve
load("roc.coords.RData")
roc.plot <- ggplot(roc.coords, aes(x=specificity, y=sensitivity, color=Model, group=Model)) + geom_vline(xintercept = 0.733) +
  geom_line(lwd=1) + theme_classic()  + scale_x_continuous(trans='reverse') + geom_abline(intercept=1, slope=1, lty=2) +
  geom_point(aes(x=0.733, y=.819), size=3, col='black') +
  theme(legend.position = c(0.75, 0.2),text = element_text(size=10)) 
roc.plot

# 4. models with baseline cytology  ---------------

# (4.1) fit models
# source("cyt_models.R") # takes a while to run --> so easier to just load the saved models
load("cyt_models.RData")

# (4.2) ROC curve
load("roc.coords.cyt.RData")
roc.plot.cyt <- ggplot(roc.coords.cyt, aes(x=specificity, y=sensitivity, color=Model, group=Model)) + 
  geom_line(lwd=1) + theme_classic()  + scale_x_continuous(trans='reverse') + geom_abline(intercept=1, slope=1, lty=2) +
  theme(legend.position = c(0.75, 0.2), text=element_text(size=10)) + geom_point(aes(x=0.733, y=0.819), col='black', size=3) 

# Figure 3. ROC plots 
plot_grid(roc.plot, roc.plot.cyt, labels = c("A", "B")) 
ggsave2("Figure3_sep2024.pdf", width=12, height=7)
rm(roc.plot, roc.coords.cyt)

# 5. greyzone method on all models  ---------------
#source("greyzone_models.R")
load("greyzone_models.RData")
load("greyzone_models_test.Rdata")

{ # Figure 1. Two-curve plot of PPV and NPV against threshold for positive classification used for determining the grey zone with a logistic regression model with interactions. 
  res.model <- lm.interactions.results
  # base model and cyt model --> this makes more sense for deciding threshold
ggplot(res.model$coords.base, aes(x=threshold, y=ppv)) + theme_classic()+
  ylab("PPV/ NPV") + xlab('Threshold for positive classification (probability of CIN3+)') + 
  geom_line(aes(col='PPV: base model'), lwd=1) +
  geom_line(data=res.model$coords.base, aes(x=threshold, y=npv, col='NPV: base model'), lwd=1) +
  geom_vline(xintercept = res.model$greyzone,lty=2) + 
  geom_hline(yintercept = c(0.27458, 0.97061),lty=1) +
  scale_y_continuous(breaks=round(c(0, 0.25, 0.27458, 0.5, 0.75, 0.97061, 1.0), 3)) + 
  scale_x_continuous(breaks=round(sort(c(0, res.model$greyzone, 0.25, 0.50, 0.75, 1.00)), 3)) + 
  theme(legend.position= c(0.9, 0.1), legend.title = element_blank(), text=element_text(size=17))
}
ggsave2("Figure1_sep2024.pdf", width=12, height=7)

# values used in results section paragraph:
# logistic base model: 
# 9290 (31.1%) immediately referred for colposcopy and of these 2188 (23.6%) had CIN3+, but noone got cytology
length(which(predict(lm.model, type='prob', newdata=training)[,1] >= lm.model$ths)) # 9290/29808 ==> 31.2%
length(which(predict(lm.model, type='prob', newdata=training)[,1] >= lm.model$ths & training$cin3plus == 'yes')) # 2188/9290 ==> 23.6%
length(which(predict(lm.model, type='prob', newdata=training)[,1] <= lm.model$ths & training$cin3plus == 'no')) # 19436/20518 ==> 94.7%
BinomCI(9255, 29808, method='wilson')
BinomCI(2181, 9255, method='wilson')
BinomCI(19464, 29808 - 9255, method='wilson')

# logistic cytology model: 
# 12505 (41.9%) immediately referred for colposcopy and of these 2508 (20.1%) had CIN3+, but 100% got cytology
length(which(predict(lm.model.cyt, type='prob', newdata=training)[,1] >= lm.model.cyt$ths)) # 7156/29808 ==> 24.0%
length(which(predict(lm.model.cyt, type='prob', newdata=training)[,1] >= lm.model.cyt$ths & training$cin3plus == 'yes')) # 890/7156 ==> 12.4%
length(which(predict(lm.model.cyt, type='prob', newdata=training)[,1] <= lm.model.cyt$ths & training$cin3plus == 'no')) # 16541/17303 ==> 89.4%
BinomCI(9900, 29808, method='wilson') # prop referred
BinomCI(2791, 9900, method='wilson') # ppv?
BinomCI(19429, 29808 - 9900, method='wilson')

# logistic model grey-zone: 
# Of the 36.4% (10837/29808 patients) inside the base model grey zone, 26.3% (2851/10837) had a predicted probability of CIN3+ 
# equal to or above this cut-off, and would therefore be referred based on their grey-zone prediction.
replaced <- ifelse(predict(lm.model, type='prob', newdata=training)[,1] >= lm.results$greyzone[1] & predict(lm.model, type='prob', newdata=training)[,1] <= lm.results$greyzone[2], 1, 0 ) # 10837/29808
probs <- predict(lm.model, type='prob', newdata=training)[,1]
probs[replaced==1] <- predict(lm.model.cyt, type='prob', newdata=training)[,1][replaced==1]

combined.probs.outcome <-  case_when(
  replaced == 1 & probs >= lm.results$cyt.model.ths ~ "yes",
  replaced == 1 & probs <= lm.results$cyt.model.ths ~ "no",
  replaced == 0 & probs <= lm.results$greyzone[1] ~ "no",
  replaced == 0 & probs >= lm.results$greyzone[2] ~ "yes"
)
sum(replaced==1 & probs >= lm.results$cyt.model.ths)
sum(combined.probs.outcome == 'yes' & training$cin3plus == 'yes')

confusionMatrix(factor(combined.probs.outcome, levels=c("yes", "no")), training$cin3plus)
sum(replaced==0 & probs >= lm.results$greyzone[2]) # 7164 -- immediately referred without cytology


# combine all results into one data frame ---> Training data 
{
results <- list(lm.results, lm.interactions.results, nb.results, rf.results, xgb.results, ada.results, bart.results)
overall <- do.call(rbind, lapply(results, function(x) x$results))
greyzone.lower <- sapply(results, function(x) c(0, 0, x$greyzone[1])) %>% as.numeric()
greyzone.upper <- sapply(results, function(x) c(0, 0, x$greyzone[2])) %>% as.numeric()
cyt.cutoff <- sapply(results, function(x) c(0, 0, x$cyt.model.ths)) %>% as.numeric()
greyzone <- ifelse(greyzone.upper==0, "--", paste0("[", greyzone.lower, " - ", greyzone.upper, "]" )) %>% as.character()
overall <- cbind(model = rep(c("lm", "lm.int", "nb", "rf", "xgb", "ada", "bart"), each=3), greyzone, overall) ; overall
#saveRDS(overall, file="overall_res.RDS")

# combine all results into one data frame ---> Testing data 
results.test <- list(lm.results.test, lm.interactions.results.test, nb.results.test, rf.results.test, xgb.results.test, ada.results.test, bart.results.test)
overall.test <- do.call(rbind, lapply(results.test, function(x) x$results))
overall.test <- cbind(model = rep(c("lm", "lm.interactions", "nb", "rf", "xgb", "ada", "bart"), each=3), overall.test); overall.test
#saveRDS(overall.test, file="overall_test.rds")

table2 <- cbind(overall, overall.test[,3:8])
table2[, 4:15] <- table2[, 4:15]*100
table2 <- table2[, c(1, 3, 2, 4:15)]
table2
}

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


write.csv(table2, file='table2_nov192024.csv')

# 6. sensitivity analysis  ---------------
set.seed(130621)
keep.cols.sens <- c("age.group", "genotype", "cin3.history", "recent.screen", "ct.value")
# (a) clinician-sampling
CS.model <- train(cin3plus~. , data=training[training$test.type =="Uitstrijkje" ,c(keep.cols.sens, "cin3plus")], trControl = train_control, method = "glm", family='binomial') 
CS.model.cyt <- train(cin3plus~., data=training[training$test.type =="Uitstrijkje" ,c(keep.cols.sens, "primary.cyt", "cin3plus")],trControl = train_control, method = "glm", family='binomial')
CS.model$ths <- threshold.finder(CS.model, 0.733)
CS.model.cyt$ths <- threshold.finder(CS.model.cyt, 0.733)

CS.results <-  make.greyzone(CS.model, CS.model.cyt, training[training$test.type =="Uitstrijkje",], plot=F)
CS.results.test <- test.greyzone(CS.model, CS.model.cyt, testing[testing$test.type =="Uitstrijkje",],  CS.results$greyzone[1], CS.results$greyzone[2], CS.results$cyt.model.ths)
CS.results; CS.results.test

# (b) self-sampling
SS.model <-  train(cin3plus~., data=training[training$test.type =="ZAS" ,c(keep.cols.sens, "cin3plus")],trControl = train_control, method = "glm", family='binomial')
SS.model.cyt <-  train(cin3plus~., data=training[training$test.type =="ZAS" ,c(keep.cols.sens,"primary.cyt" ,"cin3plus")],trControl = train_control, method = "glm", family='binomial')
SS.model$ths <- threshold.finder(SS.model, 0.733)
SS.model.cyt$ths <- threshold.finder(SS.model.cyt, 0.733)

SS.results <-  make.greyzone(SS.model, SS.model.cyt, training[training$test.type =="ZAS",], plot=F)
SS.results.test <- test.greyzone(SS.model, SS.model.cyt, testing[testing$test.type =="ZAS",], SS.results$greyzone[1], SS.results$greyzone[2], SS.results$cyt.model.ths)
SS.results; SS.results.test

table3 <- cbind(model = rep(c("CS", "SS"), each=3), 
      greyzone = c("--", "--",  paste0("[", CS.results$greyzone[1], " - ", CS.results$greyzone[2], "]" ), "--", "--", paste0("[", SS.results$greyzone[1], " - ", SS.results$greyzone[2], "]" ) ),
      do.call(rbind, list(CS.results$results, SS.results$results)), do.call(rbind, list(CS.results.test$results[,2:7], SS.results.test$results[,2:7])))
table3[, 4:15] <- table3[, 4:15]*100
table3
write.csv(table3, file='table3_nov192024.csv')
