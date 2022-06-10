# do some additional analyses for revisions

# factor score correlation
# whole CFA

if (!require("pacman")){
  install.packages("pacman")
}
pacman::p_load(tidyverse, 
               rio, 
               psych, 
               car,
               lavaan,
               sirt,
               here,
               MASS,
               BayesFactor,
               dplyr)
# set the current directory
here::i_am('README.md')
here::here()

# load files
load(here('1_conspiracy/1_conspiracy.RData'))
load(here('2_antiexpert/2_antiexpert.Rdata'))

# data prep
variables <- c('conspirational_think_1','conspirational_think_2',
               'conspirational_think_3','conspirational_think_4',
               'antiexpert_1','antiexpert_2','antiexpert_3')
dat<- dat %>% dplyr::select(UserLanguage, #subset df
                                 variables,
                                 -contains("NAppl"))


# filter out languages where N < 100
# let's move on to MI test across languages
table(dat$UserLanguage)
# shall be >= 100
# extract  languages where N >= 100
n.langs <- table(dat$UserLanguage)
list.langs <- labels(n.langs)[[1]]
langs.include <- list.langs[n.langs>=100]
n.include <- n.langs[n.langs>=100]


# extract data
for (i in 1:length(langs.include)){
  if (i == 1){
    dat_mi <- dat[dat$UserLanguage == langs.include[i],]
  }else{
    current <- dat[dat$UserLanguage == langs.include[i],]
    dat_mi <- rbind(dat_mi,current)
  }
}

# do whole CFA with to factors
model <- '
CONSP =~ conspirational_think_1 + conspirational_think_2+
  conspirational_think_3 + conspirational_think_4
ANTI =~ antiexpert_1 + antiexpert_2+antiexpert_3

CONSP~~ANTI
'

result <- cfa(model,data=dat_mi,estimator='DWLS')
summary(result,standardized=T)





#### cross-validation
# whether measurement alignment works well across two different groups
# in terms of R2s

# filter out languages where N < 300 for crossvalidation
# let's move on to MI test across languages
table(dat$UserLanguage)
# shall be >= 300
# extract  languages where N >= 300
n.langs <- table(dat$UserLanguage)
list.langs <- labels(n.langs)[[1]]
langs.include <- list.langs[n.langs>=300]
n.include <- n.langs[n.langs>=300]


# extract data
for (i in 1:length(langs.include)){
  if (i == 1){
    dat_mi <- dat[dat$UserLanguage == langs.include[i],]
  }else{
    current <- dat[dat$UserLanguage == langs.include[i],]
    dat_mi <- rbind(dat_mi,current)
  }
}


# randomly pickup 1/2
set.seed(12345)
picked <- sample(seq_len(nrow(dat_mi)),size = floor(.5*nrow(dat_mi)))
dat.train <- dat_mi[picked,]
dat.val <- dat_mi[-picked,]

# consp
# do configural CFA
model.consp <- 'CONSP =~ conspirational_think_1 + conspirational_think_2+
  conspirational_think_3 + conspirational_think_4'
cfa.consp.train <- cfa(model.consp,dat.train,estimator='DWLS',group = 'UserLanguage')
# do alignment
var.consp <- c('conspirational_think_1','conspirational_think_2',
               'conspirational_think_3','conspirational_think_4')
par.consp.train <- invariance_alignment_cfa_config(dat = dat.train[,var.consp], 
                                       group = dat.train$UserLanguage)
mod1.consp.train <- invariance.alignment(lambda = par.consp.train$lambda, nu =
                                           par.consp.train$nu, align.scale = c(0.2, 0.4), 
                                         align.pow = c(0.25, 0.25))
# test performance
mod1.consp.train$es.invariance['R2',]
#  loadings intercepts 
#0.9859013  0.9898674  

# with balidation set
cfa.consp.valid <- cfa(model.consp,dat.val,estimator='DWLS',group = 'UserLanguage')
# do alignment
par.consp.val<- invariance_alignment_cfa_config(dat = dat.val[,var.consp], 
                                                   group = dat.val$UserLanguage)
mod1.consp.val <- invariance.alignment(lambda = par.consp.val$lambda, nu =
                                         par.consp.val$nu, align.scale = c(0.2, 0.4), 
                                         align.pow = c(0.25, 0.25))
# test performance
mod1.consp.val$es.invariance['R2',]
#  loadings intercepts 
#0.980530   0.988533 


# Passed cross-validation!

# now anti expert
var.anti <- c('antiexpert_1','antiexpert_2','antiexpert_3')
model.anti <- 'ANTI =~ antiexpert_1 + antiexpert_2+antiexpert_3 '

# train set
cfa.anti.train <- cfa(model.anti,data=dat.train,estimator='DWLS',group = 'UserLanguage')
par.anti.train <- invariance_alignment_cfa_config(dat = dat.train[,var.anti], 
                                                   group = dat.train$UserLanguage)
mod1.anti.train <- invariance.alignment(lambda = par.anti.train$lambda, nu =
                                          par.anti.train$nu, align.scale = c(0.2, 0.4), 
                                         align.pow = c(0.25, 0.25))
# test performance
mod1.anti.train$es.invariance['R2',]
#  loadings intercepts 
#0.9585151  0.9884734  

# with balidation set
cfa.anti.valid <- cfa(model.anti,dat.val,estimator='DWLS',group = 'UserLanguage')
# do alignment
par.anti.val<- invariance_alignment_cfa_config(dat = dat.val[,var.anti], 
                                                group = dat.val$UserLanguage)
mod1.anti.val <- invariance.alignment(lambda = par.anti.val$lambda, nu =
                                         par.anti.val$nu, align.scale = c(0.2, 0.4), 
                                       align.pow = c(0.25, 0.25))
# test performance
mod1.anti.val$es.invariance['R2',]
#  loadings intercepts 
# 0.9405767  0.9880928 


# print out factor loadings
# consp
load(here('1_conspiracy/1_conspiracy.RData'))


mod1$lambda.aligned

# anti
load(here('2_antiexpert/2_antiexpert.Rdata'))

mod2$lambda.aligned

