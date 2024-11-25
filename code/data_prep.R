# library(readxl)
# library(stringr)
# library(lubridate)
# library(dplyr)
# library(haven)
# library(purrr)
# library(caret)
# set.seed(13062001)
rivm.dat <- read_xlsx("M:/Sectie_DMHE/Projects/RISCC/2. Workpackage 2/2b. Data collection/RIVM/Dataset VUmc.xlsx")
#rivm.dat.robine <- read_sav("M:/Sectie_DMHE/Projects/RISCC/2. Workpackage 2/2b. Data collection/RIVM/HPVpositives_RIVM.sav")

# Women were included for analysis if they did not have a screen recorded in the six months 
# prior to the current visit, had adequate cytology, genotyping results and ct-value available, 
# and did not have a history of cervical cancer in the past 10 years. 
cxca.past <- apply(rivm.dat[,str_detect(colnames(rivm.dat), "-diagnose")], 1, 
                   function(x) any(x %in% c("A", "B", "C", "D", "E", "S", "T")))
prev.6.months <- ifelse(!is.na(rivm.dat$`H1-datum`) & difftime(rivm.dat$`H1-datum`, rivm.dat$`Datum HPV` ) > -6*30, 1, 0)
valid.hpv <- ifelse(rivm.dat$hrHPV !=99 & !is.na(rivm.dat$hrHPV), 1, 0)
inclusion <- ifelse(valid.hpv & !prev.6.months & !cxca.past, T, F)
sum(inclusion)

rivm.dat <- rivm.dat[inclusion, ]
# rivm.dat.robine <- rivm.dat.robine[rivm.dat.robine$HPV==1 & !is.na(rivm.dat.robine$HPV),]
# rivm.dat.robine <- rivm.dat.robine[inclusion==1,]
#-----------------
# feature making 
# - worst history result + date
# - most recent history result + date
# - HPV genotypes indicator + type-specific ct-value
# - screened in the past 7 years 
# - history of CIN3
# - age group (make 60+ age category)
# - sample method

data.cleaning <- function(){
  rivm.code.to.cin <- function(x){
    return(case_when(
      x %in% c("A", "B", "C", "D", "E", "S", "T") ~ "CxCA",
      x %in% c("F", "I", "J", "K", "L") ~ "CIN3",
      x %in% c("M") ~ "CIN2",
      x %in% c("N", "V", "W", "O", "P")~ "<CIN2",
      T ~ ""))
  }
  
  find.worst.past <- function(x){
    diags_raw <- x[str_detect(names(x), "-diagnose")]
    dates <- x[str_detect(names(x), "-datum")]
    if (all(is.na(diags_raw))) return(c(worst = NA, date =NA))
    
    diags <- rivm.code.to.cin(diags_raw[!is.na(diags_raw)])
    dates <- dates[!is.na(dates)]
    
    indx <- which(diags == max(diags))
    worst <- ifelse(diags[indx][1] =="", NA, diags[indx][1])
    date <- ifelse(diags[indx][1] =="", NA, dates[indx][1])
    return(list(worst =worst, date = date))
  }
  
  # - worst previous screening result + date
  worst.past <- matrix(unlist(apply(rivm.dat, 1, function(x) find.worst.past(x))), ncol=2, byrow=T)
  worst.result <- worst.past[,1]
  worst.result.date <- worst.past[,2]
  
  # - history of CIN3 yes/no?
  cin3.history <- ifelse(worst.result != 'CIN3' | is.na(worst.result), 0, 1)
  
  # - most recent screening result + date
  last.result <- rivm.dat$`H1-diagnose`
  last.result.date <- rivm.dat$`H1-datum`
  
  # We combined the results of the different channels into one variable for HPV genotype and divided women
  # into one of three categories: women who are HPV16+, women who are HPV16- and HPV18+, and women who 
  # are HPV16/18- and positive for one of the other hrHPV types, neglecting single or multiple type infections. 
  
  # - HPV16/18 as separate categories: (1) HPV16 or (2) HPV18 or (3) HPVother
  hpv16 <- ifelse(!is.na(rivm.dat$HPV16), 1, 0)
  hpv18 <- ifelse(!is.na(rivm.dat$HPV18) & hpv16!=1, 1, 0)
  hpvother <- ifelse(hpv16==0 & hpv18==0, 1, 0)
  
  genotype <- case_when(
    hpv16 ==1~ "HPV16",
    hpv18==1 ~"HPV18",
    hpvother == 1 ~ "HPVother"
  )
  
  hpv1618 <- ifelse(!is.na(rivm.dat$HPV16) | !is.na(rivm.dat$HPV18), 1, 0)
  
  rivm.dat[,c("HPV16", "HPV18", "HPVother")][is.na(rivm.dat[,c("HPV16", "HPV18", "HPVother")])] <- 41
  
  # - ct value: >= 41 means negative and a lower value means stronger positive result
  # (I took ct-value corresonding to HPV genotype that they were positive for in the hierarchical order 16-->18-->other)
  ct_value <- case_when(
    hpv16==1~ rivm.dat$HPV16,
    hpv18==1~ rivm.dat$HPV18,
    T ~ rivm.dat$hrHPV
  )
  
  # - screened in past 7 years yes/no?
  # We defined women as being recently screened when the time to previous cytology and/or 
  # histology was not more than 7 years
  recent.diff.time <- difftime(rivm.dat$`Datum HPV`, rivm.dat$`H1-datum`)
  recent.screen <- ifelse(recent.diff.time > 7*365.25 | is.na(recent.diff.time), 0, 1)
  
  
  # outcome variable: CIN2+, CIN3+ and cervical cancer
  CIN2 <- ifelse(rivm.dat$`Histologie PO-max` == 'CIN2', 1, 0)
  CIN3 <- ifelse(rivm.dat$`Histologie PO-max` == 'CIN3', 1, 0)
  cxca <- ifelse(rivm.dat$`Histologie PO-max` == 'cxca', 1, 0)
  
  
  # combine all variables into a data set:
  dat <- data.frame(id = rivm.dat$Adminnummer, region=rivm.dat$Screeningsregio, 
                    age.group = rivm.dat$Leeftijdscategorie, test.type = rivm.dat$Monstertype,
                    hpv.test.date = rivm.dat$`Datum HPV`, genotype = genotype, ct.value = ct_value,
                    primary.cyt = rivm.dat$`PO-cytologie`, primary.cyt.date = rivm.dat$`Datum PO-cytologie`,
                    control.cyt = rivm.dat$`CO-cytologie`, control.cyt.date = rivm.dat$`Datum CO-cytologie`,
                    worst.history = worst.result, worst.history.date = worst.result.date,
                    last.result = last.result, last.result.date = last.result.date, 
                    cin3.history = cin3.history, recent.screen= recent.screen,
                    cin2plus = ifelse(CIN2 | CIN3 | cxca, 1, 0), 
                    cin3plus = ifelse(CIN3|cxca, 1, 0), 
                    cxca = cxca)
  dat$age.group[dat$age.group =='65+'] <- '60+'
  dat$age.group[dat$age.group =='60-64'] <- '60+'
  dat$age.group <- factor(dat$age.group, levels =  rev(names(table(dat$age.group))))
  
  dat$ct.group <- case_when(ct_value < 30 ~ "<30", 
                            ct_value >=30 & ct_value <=35~"30-35", 
                            ct_value>35~">35")
  
  dat$recent.screen <- factor(ifelse(dat$recent.screen==1, "yes", "no"), levels=c("yes", "no"))
  dat$cin3.history <-  factor(ifelse(dat$cin3.history==1, "yes", "no"), levels=c("yes", "no"))
  # make outcome variables factors:
  dat$cin2plus[is.na(dat$cin2plus)] <- 0 # they didn't have histology because cytology wasn't bad enough so can say cin2=0 
  
  dat$cin3plus[is.na(dat$cin3plus)] <- 0 # they didn't have histology because cytology wasn't bad enough so can say cin3=0 
  dat$cin3plus <- ifelse(dat$cin3plus == 1, "yes", "no")
  dat$cin3plus <- factor(dat$cin3plus, levels=c("yes", "no"))

  dat$cxca <- ifelse(dat$cxca == 1, "yes", "no")
  dat$cxca <- factor(dat$cxca, levels=c("yes", "no"))
  
  dat$primary.cyt <- factor(dat$primary.cyt)
  
  return(dat)
}

dat <- data.cleaning()

#excludes <- which(ifelse(dat$cin3plus=="no" & dat$cin2plus==1, 1, 0)==1)

#dat <- dat[-excludes,]

# testing and training set
inTraining <- caret::createDataPartition(dat$cin3plus, p = 0.75, list = FALSE)
training <- dat[ inTraining,]
testing  <- dat[-inTraining,]

factor.cols <- c("age.group", "genotype", "test.type", "cin3.history", "recent.screen" )
keep.cols <- c(factor.cols, "ct.value")

training[,factor.cols] <- lapply(training[,factor.cols], factor)
testing[,factor.cols] <- lapply(testing[,factor.cols], factor)


#write.csv(dat, "clean_dat.csv")
