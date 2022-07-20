# Correlation test for each country
# single group CFA -> factor score calculation -> correlation

library(lavaan)
library(psych)
library(tidyverse)
library(dplyr)

# load data
data <- read.csv(("../_cleandata/Final_COVIDiSTRESS_Vol2_cleaned.csv")) 
dat<- data %>% dplyr::select(UserLanguage, #subset df
                                contains("antiexpert"),
                                contains("conspirational_think"),
                                contains("trust"),
                                -contains("NAppl"))
dat<-na.omit(dat)

#extract languages with n>=100
n.langs <- table(dat$UserLanguage)
list.langs <-labels(n.langs)[[1]]
langs.include <-list.langs[n.langs>=100]
n.include <- n.langs[n.langs>=100]

#extract data
for (i in 1:length(langs.include)){
  if (i == 1){
    data.mi <-dat[dat$UserLanguage == langs.include[i],]
  }else{
    current <- dat[dat$UserLanguage == langs.include[i],]
    data.mi <- rbind(data.mi,current)
  }
}


# rename education levels
#data.mi$education[data.mi$education == "University degree (e.g., MA, MSc, BA, BSc)"] <- "University degree"
#data.mi$education[data.mi$education == "Some University or equivalent \r\n(still ongoing, or completed a module or more, but did not graduate)"] <- "Some University"
dat_cor<- data.mi

variables.consp <- c('conspirational_think_1','conspirational_think_2',
                     'conspirational_think_3','conspirational_think_4')
variables.anti <- c('antiexpert_1','antiexpert_2','antiexpert_3')

# first, start with the whole group
dat_whole <- dat_cor


cfa.model.aess <- 'AESS =~ antiexpert_1 + antiexpert_2+antiexpert_3'
cfa.model.cts <- 'CTS =~ conspirational_think_1 + conspirational_think_2+ 
 conspirational_think_3 + conspirational_think_4'

# do cfa
cfa.whole.aess <- cfa(cfa.model.aess,dat_cor,estimator='DWLS')
dat_whole_aess <- lavPredict(cfa.whole.aess, dat_whole)

cfa.whole.cts <- cfa(cfa.model.cts, dat_cor, estimator='DWLS')
dat_whole_cts <- lavPredict(cfa.whole.cts,dat_whole)

dat_whole <- cbind(dat_whole,dat_whole_aess,dat_whole_cts)

# correlation
corr.test(data.matrix(dat_whole[,c(17,16,9:15)]))
