rm(list=ls())
library(h2o)
library(nse)
library(coda)
library(TSA)
library(ggplot2)
library(reshape2)
library(glmnetRcpp)

source("InfluenceFunctions.R")
source("SampleMomentConditionFunctions.R")
source("LongRunVarianceEstimators.R")
source("Measures.R")

CI.length = function(p, mean = 0, StdDev = 1){
  return(qnorm(1-p/2, mean, StdDev) - qnorm(p/2, mean, StdDev))
}


#################################### start of testing code ################
#### This file is testing estimators on Sharpe Ratio

# h2o.init(nthreads = -1)
h2o.init(nthreads = -1, max_mem_size = "12g")
mydate=Sys.time()
myseed=format(mydate,"%Y%m%d%H")
myseed=as.integer(myseed)
set.seed(myseed)

## Max order for the polynomial in GLM.LASSO estimator
d.GLM.LASSO=7
isTrunc=FALSE
trunc=0.1
trunc1=0.25

## Max order for the spectrum0 estimator
d.spectrum0=2

## K and d as in H&W1981
K.HW1981=25
d.HW1981=2

## Parameters for generating sample AR1 series
ndata=100
n.mc=100 ### Number of MC runs to compute one SE
n.start=100
nsim=2
mu=0
sigma=1

# Names of the measures to be studied
measure.names = c("Mean", "StdDev", "VaR", "ES", "SR", "SoR")
measure.values = c(0, 1, 0, 0, 0, 0)
measure.funs = list(mean, sd, VaR.hist, ES.hist, SR, SoR.const)
measure.IFs = list(mu.IF, SD.IF, VaR.IF, ES.IF, SR.IF, SoR.const.IF)

# measure.names = c("ES", "SR", "SoR")
# measure.values = c(0, 0, 0)
# measure.funs = list(ES.hist, SR, SoR.const)
# measure.IFs = list(ES.IF, SR.IF, SoR.const.IF)

for(measure.iter in 1:length(measure.names)){
  ### Name of the measure being studied
  measure.name = measure.names[measure.iter]
  measure.value = measure.values[measure.iter]
  Current.Fun=function(data){
    return(measure.funs[[measure.iter]](data))
  }
  ### The influence function of the measure being studied
  IF.fun=function(data){
    return(measure.IFs[[measure.iter]](data))
  }
  
  ### The methods being used
  if(isTrunc){
    method.names=c("Monte Carlo","Newey-West",
                   "HW1981","GLM.Elastic",paste("IF.SE.GLM.ELASTIC.trunc",trunc),
                   paste("IF.SE.GLM.ELASTIC.trunc",trunc1))
  } else {
    method.names=c("Monte Carlo","Newey-West", "nse.Spec0.AR.0",
                   "HW1981","GLM.Elastic",
                   "GLM.BFGS")
  }
  
  phis = seq(0,0.9,by=0.1)
  theta= 0
  param.name = "phi"
  
  ####################
  
  res.list = list()
#  res.mat=matrix(0,nrow=length(phis),ncol=length(method.names))
  
  ptm=proc.time()
  
  for(phi.iter in 1:length(phis)){
    phi=phis[phi.iter]
    res=matrix(0,nrow=nsim,ncol=length(method.names)-1)
    for(sim.iter in 1:nsim){
      
      data.raw=as.numeric(arima.sim(model=list(ar=phi,ma=theta),n=ndata,n.start=n.start,sd=sqrt(1-phi^2)*sigma))+mu
      data=IF.fun(data.raw)
      if(isTrunc){
        tmp=c(nse.nw(data),
              SE.HW1981(data,K=K.HW1981,d=d.HW1981),
              SE.GLM.LASSO(data,d=d.GLM.LASSO,alpha=0.5),
              SE.GLM.LASSO(data,d=d.GLM.LASSO,keep=trunc,alpha=0.5),
              SE.GLM.LASSO(data,d=d.GLM.LASSO,keep=trunc1,alpha=0.5))
      } else {
        tmp=c(nse.nw(data),
              nse.spec0(x = data, type = "ar", lag.prewhite = 0),
#              coda::spectrum0.ar(x = data)$spec,
#              nse.spec0(x = data, type = "glm", lag.prewhite = NULL),

              SE.HW1981(data,K=K.HW1981,d=d.HW1981),
              SE.GLM.LASSO(data,d=d.GLM.LASSO,alpha=0.5),
              SE.GLM.BFGS(data, d=d.GLM.LASSO,alpha = 0.5)
              )
      }
      res[sim.iter,]=sqrt(tmp*length(data))
    }
    data.raw=as.numeric(arima.sim(model=list(ar=phi,ma=theta),n=ndata*n.mc*nsim,n.start=n.start,sd=sqrt(1-phi^2)*sigma))+mu
    data.raw=array(data.raw,c(ndata,n.mc,nsim))
    se.mc=mean(apply(data.raw,3,function(x) sd(apply(x,2,Current.Fun))))*sqrt(ndata)
    # variance.analytical=1/(1-phi)+1/2/(1-phi^2)*(mu/sigma)^2
#    res=apply(res,2,mean)
#    res=c(se.mc, res)
    res = cbind(se.mc, res)
    colnames(res)=method.names
    res.list[[paste(param.name,phi,sep = "=")]] = res
    saveRDS(res.list, file=paste("./data/",measure.name,"_SE_",myseed,"_",nsim,"_",ndata,"_",measure.value,"_d=",d.GLM.LASSO,"_AR1",".rds",sep=""))
    # h2o.shutdown()
#    res.mat[phi.iter,]=res
#    write.csv(res.mat,file=paste("./data/",measure.name,"_SE_",myseed,"_",nsim,"_",ndata,"_",measure.value,"_d=",d.GLM.LASSO,"_AR1",".csv",sep=""),row.names = FALSE)
    
  }
  
  ### Add PureMC results to the first column of res.mat
  
  
  running.time=proc.time()-ptm
  
  ### organize results
  saveRDS(res.list, file=paste("./data/",measure.name,"_SE_",myseed,"_",nsim,"_",ndata,"_",measure.value,"_d=",d.GLM.LASSO,"_AR1",".rds",sep=""))
}

# shutdown h2o server
h2o.shutdown(prompt = FALSE)
