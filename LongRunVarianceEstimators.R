#### Collection of long run variance estimation functions

library(coda)
library(TSA)

sampleVariance.cor=function(data){
  mu.hat=mean(data)
  sum=0
  N=length(data)
  tmp=matrix(0,nrow=N,ncol=N)
  for(i in 1:N){
    for(j in 1:N){
     tmp[i,j]=(data[i]-mu.hat)*(data[j]-mu.hat)
    }
  }
  return(tmp)
}

#' Compute the periodogram as defined in H&W 1981
#'
#' @param data Vector of data
#' @param max.freq Maximum frequency to be computed
#'
#' @return list of frequencies and corresponding periodograms
#' @export
#'
#' @examples
#' myperiodogram(rnorm(10))
myperiodogram=function(data,max.freq=0.5){
  ## data.fft=myfft(data) This is very slow
  data.fft=fft(data)
  N=length(data)
  #   tmp=1:N
  #   inset=tmp[(1:N)<floor(N/2)]
  tmp = Mod(data.fft[1:floor(N/2)])^2/N
  tmp = sapply(tmp, function(x) max(0.01,x))
  return(list(spec=tmp,freq=((0:(floor(N/2)-1))/N)))
}


#' Exact implementation of H&W 1981, compute the variance of the sample mean of the data
#'
#' @param data Vector of data
#' @param K Number of periodograms to keep
#' @param d Order of the polynomial fit
#'
#' @return Variance of the sample mean. To get asymptotic variance, multiply by length(data)
#' @export
#'
#' @examples
#' SE.HW1981(rnorm(10))
SE.HW1981=function(data,K=25,d=2){
  N=length(data)
  
  # Step 1: compute I(n/N) for n=1..2K and J(fn) for n=1..K 
  my.periodogram=myperiodogram(data)
  my.freq=my.periodogram$freq
  my.periodogram=my.periodogram$spec
  Is=my.periodogram[1:(2*K)]
  ns=1:K
  fn=(ns*4-1)/2/N
  Js=rep(0,K)
  for(n in 1:K){
    Js[n]=log((Is[2*n-1]+Is[2*n])/2)
  }
  
  # Step 2: use OLS to fit polynomial
  
  # create 1, x, x^2, ..., x^d
  x.mat=rep(1,K)
  for(col.iter in 1:d){
    x.mat=cbind(x.mat,fn^col.iter)
  }
  x.df=as.data.frame(x.mat)
  colnames(x.df)=paste0("V",0:d)
  Js=Js+0.27
  x.df=cbind(Js,x.df)
  lm.spec=lm(Js~.-1,data=x.df)
  
  a0.hat=as.numeric(coef(lm.spec)[1])
  # Step 3: adjust for bias
  
  s=solve(t(x.mat)%*%x.mat)
  c1=exp(-0.645*s[1,1]/2)
  
  # Step 4: return the estimated variance
  return(c1*exp(a0.hat)/N)
}



###### function to compute the standard error of sample mean of the data, to get asymptotic SE, multiply the SE by sqrt(length(data))
###### Using GLM with L1 regularization to fit a polynomial from the periodogram
###### then the intercept of the fitted polynomial is the spectral density p(0) at frequency 0
###### the variance of of the sample mean is p(0)/length(data)
###### and the asymptotic variance of the sample mean is p(0)
###### K is the proportion of frequencies to keep, 1 means keep all

#' Compute the variance of sample mean of the data
#'
#' Using GLM with L1 regularization to fit a polynomial from the periodogram then the intercept of the fitted polynomial is the spectral density p(0) at frequency 0
#' 
#' @param data Vector of data
#' @param d Maximum order of the polynomial
#' @param alpha Weight of regularization. alpha=1 LASSO, alpha=0 Ridge
#' @param keep Percentage of periodograms to use for regression
#'
#' @return variance of the sample mean. To get the asymptotic variance, multiply by length(data)
#' @export
#'
#' @examples
#' SE.GLM.LASSO(rnorm(10))
SE.GLM.LASSO=function(data,d=5,alpha=1,keep=1){
  
  N=length(data)
  # Step 1: compute the periodograms 
  my.periodogram=myperiodogram(data)
  my.freq=my.periodogram$freq
  my.periodogram=my.periodogram$spec
  
  # remove values of frequency 0 as it does not contain information about the variance
  my.freq=my.freq[-1]
  my.periodogram=my.periodogram[-1]
  
  # implement cut-off
  nfreq=length(my.freq)
  my.freq=my.freq[1:floor(nfreq*keep)]
  my.periodogram=my.periodogram[1:floor(nfreq*keep)]
  
  # Step 2: use GLM with L1 regularization to fit polynomial
  
  # create 1, x, x^2, ..., x^d
  x.mat=rep(1,length(my.freq))
  for(col.iter in 1:d){
    x.mat=cbind(x.mat,my.freq^col.iter)
  }
  x.df=as.data.frame(x.mat)
  colnames(x.df)=paste0("V",0:d)
  x.df=cbind(my.periodogram,x.df)
  x.h2o.df=as.h2o(x.df)
  
  # fit GLM with L1 regularization
  
  my.glm.lasso=h2o.glm(x=colnames(x.h2o.df)[-1],y=colnames(x.h2o.df)[1],
                       training_frame = x.h2o.df, family = "gamma",
                       link="log",lambda_search = TRUE,alpha=alpha,
                       ignore_const_cols = FALSE)
  my.glm.lasso@model$coefficients
  
  # predict with new data V0=1 and all other variables=0
  newx.mat=matrix(c(1,rep(0,d)),nrow=1)
  newx.df=as.data.frame(newx.mat)
  newx.h2o.df=as.h2o(newx.df)
  h2o::colnames(newx.h2o.df)=paste0("V",0:d)
  p0.hat=h2o.predict(my.glm.lasso,newx.h2o.df)[1,1]
  
  
  # Step 4: return the estimated variance
  return(p0.hat/N)
}

#' Title log-likelihood of sample data from exponential distribution
#' Given the link function g(z)=log(z), the independent variables x 
#' and the coefficients beta
#' \eqn{g(\mu_i)=X\beta}.
#'
#' @param X n by p matrix of the independent variables
#' @param b p by 1 matrix of the coefficients
#' @param y vector of sample data, length is n
#'
#' @return the log-likelihood
#' @export

LL.exp = function(b, y, X, lambda.en = 0, alpha.lasso = 0.5, link = "identity", normalize = FALSE){
  #  if(link == "identity") theta = 
  theta = exp(-X%*%b)
  tmp = sum(-theta*y+log(theta)) -  # notice the minus sign because we are maximizing LL 
    lambda.en * (alpha.lasso * sum(abs(b)) + (1-alpha.lasso)/2 * sum(b^2))
  return(tmp)
}

#' Title Compute the gradient of the log-likelihood of the sample response for exponential distribution
#' vectorized version
#'
#' @param b the coefficients
#' @param y the sample responses
#' @param X the matrix of independent variables
#'
#' @return the gradient
#' @export
#'

LL.exp.grad = function(b, y, X, lambda.en = 0, alpha.lasso = 0.5, normalize = FALSE){
  res = t(X) %*% (exp(-X%*%b) * y - 1) - 
    lambda.en *(alpha.lasso * sign(b) + (1 - alpha.lasso) * b)
  if (normalize){
    return(normalize_vector(res))
  } else{
    return(res)
  }
}

#' Title Compute the gradient of the log-likelihood of the sample response for exponential distribution
#' slow with loops
#'
#' @param b the coefficients
#' @param y the sample responses
#' @param X the matrix of independent variables
#'
#' @return the gradient
#' @export
#'

LL.exp.grad_loop = function(b, y, X, lambda.en = 0, alpha.lasso = 0.5){
  res = matrix(0, nrow = length(b), ncol = 1)
  for(j in 1:length(b)){
    tmp = 0
    for(i in 1:length(y)){
      tmp = tmp + (y[i]*exp(-sum(X[i,]*b))-1)*X[i,j]
    }
    tmp = tmp - lambda.en * (alpha.lasso * sign(b[j]) + (1 - alpha.lasso) * 2 * b[j])
    res[j,1] = tmp
  }
  return(res)
}

SE.GLM.BFGS=function(data, d=7, 
                     lambda.en = 0, alpha=0.5, 
                     normalize = TRUE, keep=1){
  
  N=length(data)
  # Step 1: compute the periodograms 
  my.periodogram=myperiodogram(data)
  my.freq=my.periodogram$freq
  my.periodogram=my.periodogram$spec
  
  # remove values of frequency 0 as it does not contain information about the variance
  my.freq=my.freq[-1]
  my.periodogram=my.periodogram[-1]
  
  # implement cut-off
  nfreq=length(my.freq)
  my.freq=my.freq[1:floor(nfreq*keep)]
  my.periodogram=my.periodogram[1:floor(nfreq*keep)]
  
  # Step 2: use GLM with BFGS optimization
  
  # create 1, x, x^2, ..., x^d
  x.mat=rep(1,length(my.freq))
  for(col.iter in 1:d){
    x.mat=cbind(x.mat,my.freq^col.iter)
  }
  
  b0 = rnorm(d + 1)
  
  res.grad = optim(b0, LL.exp, gr = LL.exp.grad, 
                   y = my.periodogram, X = x.mat, 
                   lambda.en = lambda.en, normalize = normalize,
                   method = "BFGS", control = list(fnscale=-1))
  
  # Step 3: return the estimated variance
  return(exp(res.grad$par[1])/N)
}

normalize_vector = function(x){
  return(x/sqrt(sum(x^2)))
}

# use crossvalidation to choose lambda.en

cv.SE.GLM.BFGS = function(data, d=7, 
                          lambda.en = 0, alpha=0.5, 
                          normalize = TRUE,
                          nfolds = 5,
                          cv.results = FALSE,
                          keep=1){
  N=length(data)
  # Step 1: compute the periodograms 
  my.periodogram=myperiodogram(data)
  my.freq=my.periodogram$freq
  my.periodogram=my.periodogram$spec
  
  # remove values of frequency 0 as it does not contain information about the variance
  my.freq=my.freq[-1]
  my.periodogram=my.periodogram[-1]
  
  # implement cut-off
  nfreq=length(my.freq)
  my.freq=my.freq[1:floor(nfreq*keep)]
  my.periodogram=my.periodogram[1:floor(nfreq*keep)]
  
  # Step 2: use GLM with BFGS optimization
  
  # create 1, x, x^2, ..., x^d
  x.mat=rep(1,length(my.freq))
  for(col.iter in 1:d){
    x.mat=cbind(x.mat,my.freq^col.iter)
  }
  
  res = list()
  for (cv.iter in 1:nfolds){
    b0 = rnorm(d + 1)
    # generate the indices for test set
    test_idx = sample(1:length(my.periodogram), 
                      length(my.periodogram)/nfolds)
    
    # train model excluding test set
    res.grad = optim(b0, LL.exp, gr = LL.exp.grad, 
                     y = my.periodogram[-test_idx], 
                     X = x.mat[-test_idx,], 
                     lambda.en = lambda.en, 
                     normalize = normalize,
                     method = "BFGS", 
                     control = list(fnscale=-1))
    # compute prediction using test set
    y.pred = exp(x.mat[test_idx,]%*%res.grad$par)
    
    # compute the RMSE
    rmse = sqrt(sum((my.periodogram[test_idx] - y.pred)^2))
    
    # store the results
    res[[cv.iter]] = list(par = res.grad$par, rmse = rmse)
  }
  
  # find the solution with the smallest rmse and return it
  min_i = 1
  for(i in 2:length(res)){
    if(res[[i]]$rmse < res[[min_i]]$rmse){
      min_i = i
    }
  }
  
  # if need detailed cv results,
  # return only the results for the coefficients
  # if not, return the estimated variance using the best coefs
  if (cv.results){
    return(res)
  } else {
    return(exp(res[[min_i]]$par[1])/N)
  }
}



