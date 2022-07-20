library(psych)
library(lavaan)
library(sirt)
library(MASS)



aligned.factor.scores <- function(lambda,nu,y){
  #calculate inverse matrix
  lambda1 <- ginv((lambda))
  #create matrix for nu
  ns <- nrow(y)
  nus <- matrix(nu,nrow=ns, ncol=length(nu), byrow=T)
  # y - nu
  y_nu <- y - nu
  # f = inv(lambda)*(y-nu)
  F <- lambda1 %*% t(as.matrix(y_nu))
}

#extract compliance_related variable names
# new model: without compliance_3 and 7

compliance <- c('antiexpert_1','antiexpert_2',
                'antiexpert_3')


#Load the cleaned cvs file

#load data
load('../2_antiexpert/2_antiexpert.RData')
data<-dat


DATA <- data


#extract languages with n>=100
n.langs <- table(data$UserLanguage)
list.langs <-labels(n.langs)[[1]]
langs.include <-list.langs[n.langs>=100]
n.include <- n.langs[n.langs>=100]

#extract data
for (i in 1:length(langs.include)){
  if (i == 1){
    data.mi <-data[data$UserLanguage == langs.include[i],]
 }else{
   current <- data[data$UserLanguage == langs.include[i],]
   data.mi <- rbind(data.mi,current)
 }
}

#set and examine fitmeasures
fits <- c('rmsea','srmr','cfi','tli')


#####
# CFA

#general CFA
cfa.model.cle <- 'CLE =~ antiexpert_1 + antiexpert_2+antiexpert_3'
cfa.whole.cle <-cfa(model=cfa.model.cle, data=data.mi,estimator='DWLS', group=
                      'UserLanguage')
fitMeasures(cfa.whole.cle)[fits]
#msea.sclaed       srmr   cfi.scaled    tli.scaled
#   0.07030909 0.03697213 0.99303018 0.97909053  
# acceptable configural invariance


# then metric
cfa.metric.cle <-cfa(model=cfa.model.cle, data=data.mi,estimator='DWLS', group=
                      'UserLanguage',group.equal ='loadings')
fitMeasures(cfa.metric.cle)[fits]
#rmsea.scaled         srmr   cfi.scaled   tli.scaled 
#0.08146663 0.05849418 0.97710762 0.97192760 
# unacceptable. Alignment.

# extract parameters
par.cle <- invariance_alignment_cfa_config(dat = data.mi[,compliance], 
                                            group = data.mi$UserLanguage)
# do alignment
mod1.cle <- invariance.alignment(lambda = par.cle$lambda, nu =
                                   par.cle$nu, align.scale = c(0.2, 0.4), align.pow = c(0.25, 0.25))
# test performance
mod1.cle$es.invariance['R2',]
#  loadings intercepts 
# 0.9540601  0.9800864  
# all ≥ 75%. well addressed.

# monte carlo. IRT test (optional)
#--- find parameter constraints for prespecified tolerance
cmod1.cle <- sirt::invariance_alignment_constraints(model=mod1.cle, nu_parm_tol=.7,
                                                lambda_parm_tol=.7 )

# simulation function
simulation_CLE <- function(times,n,data,n.include,seed=1){
  # get data
  data.mi <- data
  compliance <- c('antiexpert_1','antiexpert_2',
                  'antiexpert_3')
  
  # list for return
  cor.mean <- rep(0,times)
  cor.var <- rep(0,times)
  R2.loading <- rep(0,times)
  R2.intercept <- rep(0,times)
  
  # simulation replication
  # repeat
  for (j in 1:times){
    set.seed(seed)
    G <- n.include # groups
    I <- 3 # items
    
    # lambda, nu, and error_var
    err_var.cle <- matrix(1, nrow=G,ncol=I)
    # simulate data
    # enter group mu and sigma
    data.mi$COMP <- rowMeans(data.mi[,compliance])
    mu<-scale(aggregate(x=data.mi$COMP,
              by = list(data.mi$UserLanguage),
              FUN=mean, na.rm=T)[,2])[,1]
    sigma <- (aggregate(x=data.mi$COMP,
                       by = list(data.mi$UserLanguage),
                       FUN=sd, na.rm=T)[,2])
    N <- rep(n,G)
    dat <- invariance_alignment_simulate(
      par.cle$nu,par.cle$lambda,err_var.cle,mu,sigma,N
    )
    par.simul <- invariance_alignment_cfa_config(dat = dat[,compliance], 
                                               group = dat$group)
    mod1.simul <- invariance.alignment(lambda = par.simul$lambda, nu =
                                         par.simul$nu, align.scale = c(0.2, 0.4), align.pow = c(0.25, 0.25))
    
    # true vs aligned scores
    cfa.model.simul <- 'CLE =~ antiexpert_1 + antiexpert_2+antiexpert_3'
    
     
    
    cfa.simul <- cfa(cfa.model.simul,dat,estimator='DWLS',group='group',
                     group.equal=c('loadings','intercepts'),meanstructure=T)
    
    # get group mean
    params.simul <- parameterEstimates(cfa.simul)
    alpha.simul <- params.simul[(params.simul$op=='~1')&(params.simul$lhs=='CLE'),'est']
    
    # group mean correlation (Muthen 2018)
    correlation <- corr.test(alpha.simul,mod1.simul$pars$alpha0,method='spearman')$r
    
    
    
    psi.simul <- params.simul[(params.simul$op=='~~')&(params.simul$lhs=='CLE'),'est']
    correlation.psi <- corr.test(psi.simul,mod1.simul$pars$psi0,method='spearman')$r
    
    cor.mean[j] <- correlation
    cor.var[j] <- correlation.psi
    
    R2.loading[j] <- mod1.simul$es.invariance['R2',1]
    R2.intercept[j] <- mod1.simul$es.invariance['R2',2]
  }
  # make matrix
  to.return <-cbind(cor.mean,cor.var,R2.loading,R2.intercept)
  to.return <- data.matrix(to.return)
  
  message  (sprintf('%d/%d Done',j,times))
  
  return(to.return)
  
}

# use five cores

library(foreach)
library(parallel)
library(doParallel)

cores <- 5
times <- 500

cl <- parallel::makeCluster(cores,type='FORK')
doParallel::registerDoParallel(cl)

# n = 100

start_100 <-Sys.time()
now <- foreach (i = seq(1,times)) %dopar%{
  # get result
  simulation_CLE(1,100,data.mi,28,i)
  #  message(sprintf('%d',i))
}
end_100<-Sys.time()
elapsed_100 <- end_100 - start_100
# merge result
for (i in 1:times){
  if (i == 1){
    simulate_100 <- now[[1]]
  }else{
    simulate_100 <- rbind(simulate_100,now[[i]])
  }
}

# save n = 100
write.csv(data.frame(simulate_100),file='simulate_aess_100.csv',row.names = F)

# n=200


# multicore processing for n=200
start_200 <-Sys.time()
now <- foreach (i = seq(1,times)) %dopar%{
  # get result
  simulation_CLE(1,200,data.mi,28,i)
#  message(sprintf('%d',i))
}
end_200<-Sys.time()
elapsed_200 <- end_200 - start_200
# merge result
for (i in 1:times){
  if (i == 1){
    simulate_200 <- now[[1]]
  }else{
    simulate_200 <- rbind(simulate_200,now[[i]])
  }
}
# save n = 200
write.csv(data.frame(simulate_200),file='simulate_aess_200.csv',row.names = F)

# n= 500
start_500 <-Sys.time()
now <- foreach (i = seq(1,times)) %dopar%{
  # get result
  simulation_CLE(1,500,data.mi,28,i)
  #  message(sprintf('%d',i))
}
end_500<-Sys.time()
elapsed_500 <- end_500 - start_500
# merge result
for (i in 1:times){
  if (i == 1){
    simulate_500 <- now[[1]]
  }else{
    simulate_500 <- rbind(simulate_500,now[[i]])
  }
}
# save n = 500
write.csv(data.frame(simulate_500),file='simulate_aess_500.csv',row.names = F)
# stop cluster
parallel::stopCluster(cl)
