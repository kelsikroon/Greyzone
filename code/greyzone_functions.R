# greyzone functions
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

train_control <- trainControl(method = "cv",  number = 5, savePredictions=T, classProbs = T,
                              summaryFunction = twoClassSummary, verbose=T)

make.truth.groups <- function(prediction, truth){
  return(case_when(
    prediction =='yes' & truth=='no' ~ 'false.pos',
    prediction =='no' & truth=='yes' ~ 'false.neg',
    prediction =='no' & truth=='no' ~ 'true.neg',
    prediction =='yes' & truth=='yes' ~ 'true.pos',
  ))
}

#--------------------------------------------------------------------------------------------
plot.preds.greyzone <- function(preds, truth, outcome, cutoff.base, cutoff.cyt, replaced,  main = "Density plot of predictions", lower=NULL, upper =NULL, color=T, plot=F){
  
  Grouping <- factor(make.truth.groups(outcome, truth), levels=c("true.pos", "true.neg", "false.pos", "false.neg"))
  plot.data <- data.frame(preds, replaced, truth, Grouping)
  plot.data$grp.alpha <- ifelse(plot.data$Grouping == 'false.pos', 'false.pos', " ")
  
  
  p <- ggplot()+
    geom_histogram(data=plot.data[plot.data$truth=='no',], aes(x=preds, fill=Grouping, y= ..count../sum(..count..)),
                   breaks=seq(0, 1,length.out= 100), col='white', alpha=1, position= 'identity') +
    geom_histogram(data=plot.data[plot.data$truth=='yes',], aes(x=preds, y=after_stat(count)/sum(after_stat(count)), fill=Grouping),
                   breaks=seq(0, 1,length.out= 100), col='white', position='identity' )  +
    geom_histogram(data=plot.data[plot.data$truth=='no',], aes(x=preds, alpha=grp.alpha, fill=Grouping, y= ..count../sum(..count..)),
                   breaks=seq(0, 1,length.out= 100), col='white', position= 'identity') +
    theme_classic() + theme(legend.position=c(0.8, 0.5)) +
    labs(title = main) + ylim(c(0, 0.2)) +
    xlab("Proability of CIN3+ (base model)") + ylab("density (pos & neg each sum to 1)") +
    scale_fill_manual(values = rev(c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"))) +
    scale_alpha_manual(values=c('false.pos'=1, " "= 0)) + guides(alpha='none')
  
  
  
  if (!is.null(lower)){
    p <- p + labs(title = main, subtitle = paste0("Big model cut off: ", cutoff.cyt, ". Grey zone: [", lower, ", ", upper, "]")) +
      geom_vline(xintercept = c(lower, upper), lty=2) + 
      scale_x_continuous(breaks = c(0, lower, upper, 0.25, 0.5, 0.75, 1), labels = c("0.00", as.character(lower), as.character(upper), "0.25", "0.50", "0.75","1.00")) 
  }
  if (!color){
    p <- p + scale_fill_manual(values=c('#F8766D', '#7CAE00', '#7CAE00', '#F8766D'), name='Grouping', labels=c("Truly positive", "Truly negative", "Truly negative", "Truly positive")) 
  }
  if (plot) print(p)
  return(p)
}


#--------------------------------------------------------------------------------------------
# combine predictions from 2 models using given upper and lower boundaries
combo.model.predict <- function(base.model, cyt.model, base.pred, cyt.pred, test.data, lower, upper, cyt.upper){
  #baseline.pred <- predict(base.model, test.data[,c(keep.cols, "nb.pred")], type='prob')[,1]
  #if (base.model$method == "bartMachine") baseline.pred <- 1- baseline.pred
  if (lower ==upper){
    replaced <- rep(0, dim(test.data)[1])
    combined.probs <- base.pred
  }else{
    replaced <- ifelse((base.pred > lower  & base.pred <= upper), 1, 0)
    updated.pred <- cyt.pred[which(replaced==1)] #predict(cyt.model, test.data[which(replaced==1), c(keep.cols.cyt, "nb.pred.cyt")], type='prob')[,1]
    #if (cyt.model$method == 'bartMachine') updated.pred <- 1- updated.pred
    combined.probs <- base.pred
    combined.probs[which(replaced==1)] <- updated.pred
  }
  
  combined.probs.outcome <-  case_when(
    replaced == 1 & combined.probs >= cyt.upper ~ "yes",
    replaced == 1 & combined.probs <= cyt.upper ~ "no",
    replaced == 0 & combined.probs <= lower ~ "no",
    replaced == 0 & combined.probs >= upper ~ "yes"
  )#ifelse(combined.probs >= ifelse(replaced, cyt.model$ths, base.model$ths), "yes", "no")
  return(data.frame(replaced, base.pred, combined.probs, combined.probs.outcome))
}




#--------------------------------------------------------------------------------------------
# function to find the grey-zone using pre-specified PPV & NPV values and return results and plots
make.greyzone <- function(base.model, cyt.model, data, plots = F){
  suppressWarnings({
    base.pred <- predict(base.model, newdata= data, type = 'prob')[,1]
    cyt.pred <- predict(cyt.model, newdata= data, type = 'prob')[,1]
  })
  
  if (base.model$method == 'bartMachine'){
    base.pred = 1 - base.pred
    cyt.pred = 1 - cyt.pred
    base.model$ths <- 0.1134
    cyt.model$ths <- 0.1134
  }
  # find greyzone of base model
  base.model.roc <- roc(ifelse(data$cin3plus=='yes', 1, 0), base.pred, ret='all_coords', direction ='<')
  coords.base <- coords(base.model.roc, x=seq(0, 1, 0.001), input = "threshold", ret='all')[1:1000,]
  base.upper <- min(coords.base$threshold[coords.base$ppv > 0.27458], na.rm=T) 
  base.lower <- max(coords.base$threshold[coords.base$npv > 0.97061], na.rm=T)
  
  # find cutoff of bigger model using ppv and npv of cytology model for those who were replaced
  cyt.model.roc <- roc(ifelse(data$cin3plus=='yes', 1, 0), cyt.pred, ret='all_coords', direction ='<')
  coords.cyt <- coords(cyt.model.roc, x=seq(0, 1, 0.001), input = "threshold", ret='all')[1:1000,]
  cyt.lower <- min(coords.cyt$threshold[coords.cyt$ppv > 0.27458], na.rm=T) 
  cyt.upper <- max(coords.cyt$threshold[coords.cyt$npv > 0.97061], na.rm=T) 
  
  # final results of each model
 # base.model$ths <- coords.base$threshold[which.min(abs(coords.base$specificity - 0.733))]
  base.res <- confusionMatrix(factor(ifelse(base.pred >= base.model$ths, "yes", "no"), levels=c("yes", "no")), data$cin3plus)
  cyt.res <- confusionMatrix(factor(ifelse(cyt.pred >= cyt.model$ths, "yes", "no"), levels=c("yes", "no")), data$cin3plus)
  
  print(c(base.lower, base.upper))
  if (base.lower > base.upper) {
    temp <- base.upper
    base.upper <- base.lower
    base.lower <- temp
  }
  combo.predict <- combo.model.predict(base.model, cyt.model, base.pred, cyt.pred, data, base.lower, base.upper, cyt.upper)
  combo.res <- confusionMatrix(factor(combo.predict$combined.probs.outcome, levels=c("yes", "no")), data$cin3plus)
  percent.cyt <- sum(combo.predict$replaced ==1)/dim(data)[1]
  
  referred <- sum(combo.predict$combined.probs.outcome =="yes")
  # combining results
  results <- data.frame(type= c("base", "cytology", "combined"),
                        prop.cyt = c(0, 1, percent.cyt),
                        acc =  c(base.res$overall[1], cyt.res$overall[1], combo.res$overall[1]),
                        sens =  c(base.res$byClass[1], cyt.res$byClass[1], combo.res$byClass[1]), 
                        spec = c(base.res$byClass[2], cyt.res$byClass[2], combo.res$byClass[2]), 
                        ppv = c(base.res$byClass[3], cyt.res$byClass[3], combo.res$byClass[3]),
                        npv = c(base.res$byClass[4], cyt.res$byClass[4], combo.res$byClass[4]))

  results <- results %>% mutate_if(is.numeric, round, digits=4)
  if (! plots) return(list(greyzone = c(base.lower, base.upper),base.model.ths = base.model$ths,  cyt.model.ths = cyt.upper, results= results, coords.base=coords.base))

  # histogram of predictions coloured by group (e.g., true positive, true negative...) after greyzone model applied
  hist.plot <- plot.preds.greyzone(combo.predict$combined.probs,  data$cin3plus, combo.predict$combined.probs.outcome,
                                   base.model$ths, cyt.upper, combo.predict$replaced, "Grey zone method using PPV/NPV",
                                   base.lower, base.upper)
  return(list(greyzone = c(base.lower, base.upper), base.model.ths = base.model$ths, cyt.model.ths = cyt.upper, results= results, hist.plot = hist.plot))
}

#--------------------------------------------------------------------------------------------

test.greyzone <- function(base.model, cyt.model, data, lower, upper, cyt.cutoff){
  suppressWarnings({
    base.pred <- predict(base.model, newdata= data, type = 'prob')[,1]
    cyt.pred <- predict(cyt.model, newdata= data, type = 'prob')[,1]
  })
  
  if (base.model$method == 'bartMachine'){
    base.pred = 1 - base.pred
    cyt.pred = 1 - cyt.pred
  }
  # else{
  #   base.model$ths <-  threshold.finder(base.model, 0.73310)
  # }
  
  #base.model.roc <- roc(ifelse(data$cin3plus=='yes', 1, 0), base.pred, ret='all_coords', direction ='<')
  #coords.base <- coords(base.model.roc, x=seq(0, 1, 0.001), input = "threshold", ret='all')[1:1000,]
  
  # final results of each model
  #base.model$ths <- coords.base$threshold[which.min(abs(coords.base$specificity - 0.733))]
  base.res <- confusionMatrix(factor(ifelse(base.pred >= base.model$ths, "yes", "no"), levels=c("yes", "no")), data$cin3plus)
  cyt.res <- confusionMatrix(factor(ifelse(cyt.pred >= cyt.model$ths, "yes", "no"), levels=c("yes", "no")), data$cin3plus)
 
  combo.predict <- combo.model.predict(base.model, cyt.model, base.pred, cyt.pred, data, lower, upper, cyt.cutoff)
  combo.res <- confusionMatrix(factor(combo.predict$combined.probs.outcome, levels=c("yes", "no")), data$cin3plus)
  percent.cyt <- sum(combo.predict$replaced ==1)/dim(data)[1]
  
  # combining results
  results <- data.frame(type= c("base", "cytology", "combined"),
                        prop.cyt = c(0, 1, percent.cyt),
                        acc =  c(base.res$overall[1], cyt.res$overall[1], combo.res$overall[1]),
                        sens =  c(base.res$byClass[1], cyt.res$byClass[1], combo.res$byClass[1]), 
                        spec = c(base.res$byClass[2], cyt.res$byClass[2], combo.res$byClass[2]), 
                        ppv = c(base.res$byClass[3], cyt.res$byClass[3], combo.res$byClass[3]),
                        npv = c(base.res$byClass[4], cyt.res$byClass[4], combo.res$byClass[4]))
  
  results <- results %>% mutate_if(is.numeric, round, digits=4)
  return(list(results= results))
}

