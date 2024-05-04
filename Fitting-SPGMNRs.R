if(!require("mixtools")){install.packages("mixtools")}else{library(mixtools)}
if(!require("locpol")){install.packages("locpol")}else{library(locpol)}
library(splines)
library(MASS)
###Kernel function
Kern<-function(x,x0,h){
  z=(x-x0)/h
  #f=ifelse(abs(z)<=1, 0.75*(1 - z^2),0)/h
  f=dnorm(z)
  out=f
  if(sum(f)>0){out=f/sum(f)};
  return(out)
}

Dist=function(x){
  k=ncol(x)
  f=apply(x,2,sum)
  idx=order(f,decreasing = T)
  y=x[,idx]
  res=apply(y,1,which.max)
  return(list(switch=any(res!=1),id=idx))
}

mix.linear.fit=function(x,y,k,init.model0=NULL){
  library(mixtools)
  mix.prop0=init.model0$mix.prop
  mix.sigma0=sqrt(init.model0$mix.sigma2)
  model=regmixEM(y,x,k=k,lambda=mix.prop,sigma=mix.sigma0)
  r=model$posterior
  mix.beta=model$beta
  mix.pi=model$lambda
  mix.sigma2=(model$sigma)^2
  LogLik=model$loglik
  mix.mu=(model$x)%*%mix.beta
  df=2*k-1+2*k
  BIC0=-2*LogLik+log(n)*df
  R2=Rsquared(y,mix.mu,r)
  c.dist=cond_dist(y,mix.mu,mix.pi,mix.sigma2)
  est.cdf=CDFmixture(x,y,mix.pi,mix.sigma2,mix.mu)$Fy
  res=Dist(mix.mu);mix.mu=mix.mu[,res$id];mix.pi=mix.pi[res$id];mix.sigma2=mix.sigma2[res$id]
  return(list(resp=r,mix.mu=mix.mu,mix.prop=mix.pi,mix.sigma2=mix.sigma2,df=df,R2=R2,BIC=BIC0,est.dist=c.dist,est.cdf=est.cdf,LL=LogLik))
}

ColStat=function(x){
  apply(x,2,function(x) c(mean(x),sd(x)))
}

ColStat1=function(x){
  c(mean(x),sd(x))
}
##Normalizing functions

Normalize=function(x){
  x/sum(x)
}
minmaxscaler=function(x){
  (x-min(x))/(max(x)-min(x))
}

standardizer=function(x){
  z=(x-mean(x))/sd(x)
  return(z)
}

##Conditional distribution of y|x
cond_dist=function(y,mix.mu,mix.prop,mix.sigma2){
  k=ncol(mix.mu)
  rowSums(sapply(1:k,function(j) mix.prop[j]*dnorm(y-mix.mu[,j],0,sqrt(mix.sigma2[j]))))
}

mse=function(est.model,true.model){
  library(combinat)
  mu.est=est.model$mix.mu;mu.true=true.model$mu
  mix.prop=est.model$mix.prop;true.prop=true.model$prop
  mix.sigma2=est.model$mix.sigma2;true.sigma2=true.model$sigma2
  n=nrow(mu.est);K=length(mix.prop)
  out=permn(K)
  res=sapply(1:length(out),function(i){
    idx=out[[i]];
    sqrt(sum(rowSums(mu.est[,idx]-mu.true)^2)/n)
  })
  RASE=min(res)
  MSE_pi=(mix.prop[out[[which.min(res)]]]-true.prop)^2
  MSE_sigma2=(mix.sigma2[out[[which.min(res)]]]-true.sigma2)^2
  return(c(MSE_pi[1],MSE_sigma2,RASE))
}

perf=function(est,true){
  n=length(est)
  rase=sqrt(sum((est-true)^2)/n)
  mae=max(abs(est-true))
  kl=sum(true*log(true/est))
  chi2=sum(((true-est)^2)/est)
  ks=max(abs(est-true))
  return(c(rase,kl,chi2,ks))
}
##
BIC=function(x,bw,K,LogLik,dfn){
  n=length(x)
  rk=2.5375;ck=0.7737;
  dfm=(rk*abs(diff(range(x)))*ck)/bw
  dfp=2*K-1
  df=K*dfm+2*K-1
  BIC1=-2*LogLik+log(n)*df
  BIC2=-2*LogLik+log(n)*dfn
  return(c(BIC1,BIC2))
}

##Calculating the R-Squared
Rsquared=function(y,mh,gn){
  K=ncol(gn);n=length(y)
  ybar_k=colSums(y*gn)/colSums(gn)
  EWSS=sum(sapply(1:K,function(j) sum(gn[,j]*(mh[,j]-ybar_k[j])^2)))
  RWSS0=sapply(1:K,function(j) sum(gn[,j]*(y-mh[,j])^2))
  RWSS=sum(RWSS0)
  WSS=EWSS+RWSS
  RMSE=sqrt(sum(sapply(1:K,function(j){nk=sum(gn[,j]);RWSS0[j]/nk})))
  return(list(RMSE=RMSE,R2=round(EWSS/WSS,2)*100))
}

###Kernel polynomial smoother
polynomial.smoother=function(x,y,xgrid,d,W){
  ##d: denotes the degree of the polynomial
  n=length(x);ngrid=length(xgrid)
  fit=sapply(1:ngrid,function(i){
    x0=xgrid[i];W=diag(W[,i])+1e-100*diag(n)
    z=(x-x0);Z=t(t(matrix(z,n,d+1))^(0:d))
    (solve(t(Z)%*%W%*%Z)%*%t(Z)%*%W%*%y)[1]
  })
  return(fit)
}

###Kernel polynomial smoother matrix
polynomial.smoother.matrix=function(x,x.grid,d,W){
  n=length(x)
  S=sapply(1:length(x.grid),function(i){
    w=diag(W[,i]);
    X=t(t(matrix((x-x.grid[i]),n,d+1))^(0:d))
    (solve(t(X)%*%w%*%X)%*%t(X)%*%w)[1,]
  })
  return(S)
}

##Backfitting function
backfit=function(y,x,xgrid,d,k,mh,pi_init,sigma2_init,bw,backfit=TRUE){
  n=length(y)
  ##Estimating the global parameters given the non-parametric estimates
  pi0=pi_init
  sigma20=sigma2_init
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-mh[,j],0,sqrt(sigma20)[j])))))
  diff=1e10
  count=0
  while(diff>1e-10){
    #E-Step
    g=sapply(1:k,function(j) pi0[j]*dnorm(y-mh[,j],0,sqrt(sigma20[j]))+1e-300)
    gn=g/rowSums(g)
    #M-Step
    pi1=colSums(gn)/n
    sigma21=NULL
    for(j in 1:k){
      sigma21=c(sigma21,sum(gn[,j]*(y-mh[,j])^2)/sum(gn[,j]))
    }
    #Evaluate for convergence
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mh[,j],0,sqrt(sigma21[j]))))))
    diff=abs(LogLik0-LogLik1)
    LogLik0=LogLik1
    sigma20=sigma21
    pi0=pi1
    count=count+1
    if(count==5e2) diff=1e-100
  }
  
  #Re-estimating the non-parametric functions given the global parameter estimates
  ##Initialiaze the algorithm using the responsibilities from the previous stage
  ##Since the parameter estimates are well-labelled
  mu=mh
  mu0=mh
  #mu0=sapply(1:k,function(j) approx(x,mu0[,j],xgrid,rule=2)$y)
  if(backfit){
    LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mh[,j],0,sqrt(sigma21)[j])))))
    Kh=sapply(xgrid,function(x0) Kern(x,x0,bw))
    ngrid=length(xgrid)
    diff=1e10
    count=0
    while(diff>1e-10){
      ##local E-step
      #gn=lapply(1:ngrid,function(t){g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu0[t,j],0,sqrt(sigma21[j]))+1e-100);gn=g/rowSums(g)})
      g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu0[,j],0,sqrt(sigma21[j]))+1e-100);gn=g/rowSums(g)
      ##local M-step
      mu1=sapply(1:k,function(j){
        W=gn[,j]
        local.polynomial.smoother(x,y,xgrid,bw,d,W)[,2]
      })
      mu=sapply(1:k,function(j) approx(xgrid,mu1[,j],x,rule=2)$y)
      ##Evaluating convergence
      LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mu[,j],0,sqrt(sigma21)[j])))))
      diff=abs(LogLik0-LogLik1)
      LogLik0=LogLik1
      mu0=mu
      count=count+1
      if(count==5e2) diff=1e-100
    }
  }
  out=list(mu=mu,pi1=pi1,sigma21=sigma21)
}

##Local polynomial smoother
local.polynomial.smoother=function(x,y,xgrid,bw,d,W){
  library(locpol)
  n=length(y)
  g=locPolSmootherC(x,y,xgrid,bw,d,gaussK,weig=W)
  return(g)
}

###Computing the roughness of a curve
Rough_curve<-function(x,f=NULL){
  if(!is.null(f)) f=f[order(x)]
  x=sort(x)
  n=length(x);
  N=diag(n)
  ##Checking whether the x's are distinct
  if(length(unique(x))!=n){z=unique(x);n=length(z);N=sapply(z,function(u) as.numeric(x==u));  x=z}
  h<-NULL
  for(i in (1:(n-1))){h=c(h,x[i+1]-x[i])}
  Q<-matrix(0,nrow=n,ncol=n-2)
  for(j in (2:(n-1))){Q[j-1,j-1]=1/h[j-1];Q[j,j-1]=-(1/h[j-1]+1/h[j]);Q[j+1,j-1]=1/h[j]}
  R<-matrix(0,n-2,n-2)
  for(i in (2:(n-1))){R[i-1,i-1]=(h[i-1]+h[i])/3;if(i<=(n-2)){R[i,i-1]=R[i-1,i]=h[i]/6}}
  K=Q%*%solve(R)%*%t(Q)
  Rh=NULL
  if(!is.null(f)){
    f=apply(N,2,function(x) mean(f[which(x==1)])); 
    Rh=t(f)%*%K%*%f} ##Roughness value
  return(list(Rh=Rh,K=K,N=N))
}

GCV.kernel=function(x,y,xgrid=NULL,k,h.range,model=1,d=0){
  ##model=1, local EM algorithm
  ##model=2, global EM algorithm
  n=length(y)
  if(missing(xgrid)) xgrid=x
  m0=initialize.model(x,y,k=k,method=1)
  #m=Kernel_Mix_EM4(x,y,k,0.4,0,xgrid,init.model = m0)
  #m0$mu0=m$mix.mu
  plot(x,y);for(j in 1:k) lines(x,m0$mu0[,j])
  out1=out2=NULL
  for(h in h.range){
    if(model==1) fit=Kernel_Mix_EM_loc(x,y,k=k,bw=h,d=d,xgrid,init.model=m0)
    if(model==2) fit=Kernel_Mix_EM_g(x,y,k=k,bw=h,d=d,xgrid,init.model=m0)
    ##Calculating the GCV error
    r=fit$resp;
    Kh=sapply(x,function(u) Kern(x,u,h))
    GCV1=NULL;df=0
    for(j in 1:k){
      nk=sum(r[,j])
      Wk=r[,j]*Kh
      mk=(fit$mix.mu)[,j]
      Sk=polynomial.smoother.matrix(x,x,d,Wk)
      dfk=sum(diag(Sk))
      yh=mk;df=df+dfk
      SSEk=sum(r[,j]*(y-yh)^2)
      GCV1=c(GCV1,(SSEk/nk)/(1-(dfk/nk))^2) ###The one in use
    }
    yh2=rowSums(r*fit$mix.mu);
    GCV2=(sum((y-yh2)^2)/n)/(1-df/n)^2
    out1=rbind(out1,c(h,df,round(sum(GCV1),4),GCV1))
    out2=rbind(out2,c(h,df,round(GCV2,4)))
  }
  return(list(GCV1=out1,GCV2=out2))
}

###Initialization function
initialize.model=function(x,y,k,method=NULL,true.init=NULL,p=1){
  n=length(y)
  BIC=1e6
  if(method=="1"){##Mixtures of Regression splines
    for(j in 1:1e3){
      m=list(BIC=1e6,sw=FALSE,mu=0)
      try({m=mix.reg.splines(x,y,k)},silent=T)
      if(m$BIC<BIC & !m$sw){init.model=m$init.model0;BIC=m$BIC}
    }
  }
  if(method=="2"){##Mixtures of polynomial regressions
    for(j in 1:1e3){
      m=list(BIC=1e6)
      try({m=mix.poly(x,y,k,p)})
      if(m$BIC<BIC){init.model=m$init.model;BIC=m$BIC}
    }
  }
  if(method=="3"){##True values
    m0=true.init
    init.model=list(mix.mu=m0$mu,mix.sigma2=m0$sigma2,mix.prop=m0$rho)
  }
  return(init.model)
}

##mixture of polynomial regressions (to initialize)
mix.poly=function(x,y,k,p){
  n=length(y)
  Xd=t(t(matrix(x,n,p+1))^(0:p))
  ##Initial state
  pi0=rep(1/k,k)
  sigma20=rgamma(k,1,1)^2
  Beta0=matrix(rnorm(k*(p+1)),p+1,k)
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-Xd%*%Beta0[,j],0,sqrt(sigma20[j]))))))
  diff=1e6
  count=0
  while(diff>1e-10){
    ##E-Step
    g=sapply(1:k,function(j) pi0[j]*dnorm(y-Xd%*%Beta0[,j],0,sqrt(sigma20[j]))+1e-100)
    gn=g/rowSums(g)
    ##M-Step
    pi1=colSums(gn)/n
    Beta1=sigma21=NULL
    for(j in 1:k){
      W=diag(gn[,j])
      Beta1=cbind(Beta1,solve(t(Xd)%*%W%*%Xd)%*%t(Xd)%*%W%*%y)
      sigma21=c(sigma21,sum(gn[,j]*(y-Xd%*%Beta1)^2)/sum(gn[,j]))
    }
    #Evaluate for convergence
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-Xd%*%Beta1[,j],0,sqrt(sigma21[j]))))))
    diff=abs(LogLik0-LogLik1)
    LogLik0=LogLik1
    Beta0=Beta1
    sigma20=sigma21
    pi0=pi1
    count=count+1
    if(count==1e3) diff=1e-100
  }
  df=(p+1)*k+2*k-1
  BIC=-2*LogLik0+log(n)*df
  mu=Xd%*%Beta1
  r0=pi1*dnorm(y-Xd%*%Beta1,0,sqrt(sigma21));r=r0/rowSums(r0)
  if(which.max(pi1)==2) {pi1=pi1[2:1];sigma21=sigma21[2:1];Beta1=Beta1[,2:1]} 
  
  model0=list(r=r,Beta0=t(as.matrix(Beta1[-1,])),pi0=pi1,sigma20=sigma21,mu0=mu)
  return(list(init.model0=model0,BIC=BIC))
}

###Mixture of regression splines 
mix.reg.splines=function(x,y,k){
  library(splines)
  n=length(y)
  Xd=bs(x,knots=quantile(unique(x),probs=seq(0.2,0.8,0.2)),intercept=T)##B-spline basis matrix
  d=ncol(Xd);
  ##Initial state
  #pi0=rep(1/k,k)
  pi0=c(0.65,0.35)
  sigma20=c(0.09,0.16)
  #sigma20=rgamma(k,1,1)^2
  Beta0=matrix(rnorm(k*d),d,k)
  LogLik0=sum(log(rowSums(pi0*dnorm(y-Xd%*%Beta0,0,sqrt(sigma20)))))
  diff=1e6
  ll=NULL
  count=0
  while(diff>1e-10){
    ##E-Step
    g=sapply(1:k,function(j) pi0[j]*dnorm(y-Xd%*%Beta0[,j],0,sqrt(sigma20[j]))+1e-300)
    gn=g/rowSums(g)
    ##M-Step
    pi1=colSums(gn)/n
    Beta1=NULL;
    sigma21=NULL
    for(j in 1:k){
      W=diag(gn[,j])
      Beta1=cbind(Beta1,solve(t(Xd)%*%W%*%Xd)%*%t(Xd)%*%W%*%y)
      sigma21=c(sigma21,sum(gn[,j]*(y-Xd%*%Beta1[,j])^2)/sum(gn[,j]))
    }
    #Evaluate for convergence
    LogLik1=sum(log(rowSums(pi1*dnorm(y-Xd%*%Beta1,0,sqrt(sigma21)))))
    diff=abs(LogLik0-LogLik1)
    LogLik0=LogLik1
    Beta0=Beta1
    sigma20=sigma21
    pi0=pi1
    count=count+1
    if(count==1e3) diff=1e-100
  }
  df=k*d+2*k-1
  df_reg=k*d
  BIC=-2*LogLik0+df*log(n)
  mu1=Xd%*%Beta1
  res=Dist(mu1);mu1=mu1[,res$id];pi1=pi1[res$id];sigma21=sigma21[res$id]
  r0=pi1*dnorm(y-Xd%*%Beta1,0,sqrt(sigma21));r=r0/rowSums(r0)
  R2=Rsquared(y,mu1,r)
  c.dist=cond_dist(y,mu1,pi1,sigma21)
  est.cdf=CDFmixture(x,y,pi1,sigma21,mu1)$Fy
  model0=list(resp=r,mix.prop=pi1,mix.mu=mu1,mix.sigma2=sigma21,df=df,LL=LogLik1,R2=R2,BIC=BIC,est.dist=c.dist,est.cdf=est.cdf,sw=res$switch)
  return(list(r=r,BIC=BIC,mu=mu1,Beta=Beta1,pi=pi1,sigma2=sigma21,init.model0=model0,df_reg=df_reg,sw=res$switch))
}

GMM=function(y,mix.mu,mix.prop,mix.sigma){
  k=length(mix.mu)
  out=rowSums(sapply(1:k,function(j) mix.prop[j]*dnorm(y-mix.mu[j],0,mix.sigma[j])))
  return(out)
}
                     
##Semi-parametric mixtures of non-parametric regressions(Local EM)
Kernel_Mix_EM_loc=function(x,y,k,bw,d,xgrid,init.model){
  n=length(y)
  ngrid=length(xgrid)
  ##Initial state
  pi0=init.model$mix.prop
  mu0=init.model$mix.mu;
  sigma20=init.model$mix.sigma2
  ##
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[j])))))
  mu0=sapply(1:k,function(j) approx(x,mu0[,j],xgrid,rule=2)$y)
  pi0=matrix(pi0,ngrid,k,byrow=T)
  sigma20=matrix(sigma20,ngrid,k,byrow=T)
  Kh=sapply(xgrid,function(x0) Kern(x,x0,bw))
  diff=1e6
  count=0
  tol=NULL
  while(diff>1e-10){
    ##E-step
    g=sapply(1:k,function(j) pi0[,j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[,j])+1e-100);gn=g/rowSums(g)
    ##M-step
    mu1=NULL;
    pi1=sigma21=NULL
    for(j in 1:k){
      W=gn[,j]
      pi1=cbind(pi1,local.polynomial.smoother(x,rep(1,n),xgrid,bw,0,W)[,2])
      mh=local.polynomial.smoother(x,y,xgrid,bw,d,W)[,2]
      mu1=cbind(mu1,approx(xgrid,mh,x,rule=2)$y)
      sigma21=cbind(sigma21,local.polynomial.smoother(x,(y-mu1[,j])^2,xgrid,bw,0,W)[,2])
    }
    pi1=sapply(1:k,function(j) approx(xgrid,pi1[,j],x,rule=2)$y)
    sigma21=sapply(1:k,function(j) approx(xgrid,sigma21[,j],x,rule=2)$y)
    ##Evaluating convergence
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[,j]*dnorm(y-mu1[,j],0,sqrt(sigma21)[,j])))))
    diff=abs(LogLik0-LogLik1)
    LogLik0=LogLik1
    mu0=mu1
    pi0=pi1
    sigma20=sigma21
    count=count+1
    tol=c(tol,LogLik1)
    if(count==1e3) diff=1e-100
  }
  #plot.ts(tol)
 out=backfit(y,x,xgrid,d,k,mu1,colMeans(pi1),colMeans(sigma21),bw);
  mu1=out$mu;pi1=out$pi1;sigma21=out$sigma21;
  res=Dist(mu1);mu1=mu1[,res$id];pi1=init.model$mix.prop[res$id];sigma21=init.model$mix.sigma2[res$id]
  g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu1[,j],0,sqrt(sigma21[j])));gn=g/rowSums(g)
  R2=Rsquared(y,mu1,gn)
  LL1=sum(log(rowSums(g)))
  Kh=sapply(x,function(x0) Kern(x,x0,bw))
  df=sum(sapply(1:k,function(j){w=gn[,j]*Kh;S=polynomial.smoother.matrix(x,x,d,w);sum(diag(S))}))
  BIC=BIC(x,bw,k,LL1,df)
  out=list(resp=gn,mix.prop=pi1,mix.mu=mu1,mix.sigma2=sigma21,df=df,LL=LL1,BIC=BIC[2])
  return(out)
}

##NaiveEM algorithm
Kernel_Mix_EM.Naive=function(x,y,k,bw,d,xgrid,init.model){
  n=length(y)
  ngrid=length(xgrid)
  ##Initial state
  pi0=init.model$mix.prop
  mu0=init.model$mix.mu;
  sigma20=init.model$mix.sigma2
  ##
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[j])))))
  mu0=sapply(1:k,function(j) approx(x,mu0[,j],xgrid,rule=2)$y)
  pi0=matrix(pi0,ngrid,k,byrow=T)
  sigma20=matrix(sigma20,ngrid,k,byrow=T)
  Kh=sapply(xgrid,function(x0) Kern(x,x0,bw))
  diff=1e6
  count=0
  tol=NULL
  while(diff>1e-10){
    ##local E-step
    g=lapply(1:ngrid,function(i){g=sapply(1:k,function(j) pi0[i,j]*dnorm(y-mu0[i,j],0,sqrt(sigma20[i,j])));gn=g/rowSums(g)})
    ##local M-step
    mu1=t(sapply(1:ngrid,function(t){
      sapply(1:k,function(j){
        w=g[[t]][,j]
        mh=local.polynomial.smoother(x,y,xgrid[t],bw,0,w)[,2]
      })
    }))
    mu=sapply(1:k,function(j) approx(xgrid,mu1[,j],x,rule=2)$y)
    pi1=t(sapply(1:ngrid,function(t){
      sapply(1:k,function(j){
        w=g[[t]][,j]
        mh=local.polynomial.smoother(x,rep(1,n),xgrid[t],bw,0,w)[,2]
      })
    }))
    prop=sapply(1:k,function(j) approx(xgrid,pi1[,j],x,rule=2)$y)
    sigma21=t(sapply(1:ngrid,function(t){
      sapply(1:k,function(j){
        w=g[[t]][,j]
        mh=local.polynomial.smoother(x,(y-mu1[t,j])^2,xgrid[t],bw,0,w)[,2]
      })
    }))
    sigma2=sapply(1:k,function(j) approx(xgrid,sigma21[,j],x,rule=2)$y)
    ##Evaluating convergence
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) prop[,j]*dnorm(y-mu[,j],0,sqrt(sigma2)[,j])))))
    diff=abs(LogLik0-LogLik1)
    LogLik0=LogLik1
    mu0=mu1
    pi0=pi1
    sigma20=sigma21
    count=count+1
    tol=c(tol,LogLik1)
    if(count==1e3) diff=1e-100
  }
  mu1=mu
  out=backfit(y,x,xgrid,d,k,mu1,colMeans(pi1),colMeans(sigma21),bw,backfit = F);
  mu2=out$mu;pi1=out$pi1;sigma21=out$sigma21;
  res=Dist(mu2);mu2=mu2[,res$id];pi1=pi1[res$id];sigma21=sigma21[res$id]
  g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu2[,j],0,sqrt(sigma21)[j]));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  Kh=sapply(x,function(x0) Kern(x,x0,bw))
  df=sum(sapply(1:k,function(j){w=gn[,j]*Kh;S=polynomial.smoother.matrix(x,x,d,w);sum(diag(S))}))
  BIC=BIC(x,bw,k,LL1,df)
  out=list(resp=gn,mix.prop=pi1,mix.mu=mu2,mix.sigma2=sigma21,df=df,LL=LL1,BIC=BIC[2])
  return(out)
}

###Model-based ECM (MB-ECM) algorithm
Kernel_Mix_MB_ECM=function(x,y,k,bw,d,xgrid,init.model,lmd_0=1e-5){
  n=length(y)
  ngrid=length(xgrid)
  grid=1:ngrid
  xgrid0=xgrid
  ##Initial state
  lmd0=rep(1/ngrid,ngrid)
  pi0=init.model$mix.prop
  mu0=init.model$mix.mu;
  sigma20=init.model$mix.sigma2
  ##
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[j])))))
  mu0=sapply(1:k,function(j) approx(x,mu0[,j],xgrid,rule=2)$y)
  pi0=matrix(pi0,ngrid,k,byrow=T)
  sigma20=matrix(sigma20,ngrid,k,byrow=T)
  Kh=sapply(xgrid,function(x0) Kern(x,x0,bw))
  diff=dif=1e6
  tol=NULL
  count=0
  while(diff>1e-10){
    ##E-Step
    v=sapply(1:length(grid),function(t) lmd0[t]*GMM(y,mu0[t,],pi0[t,],sqrt(sigma20)[t,]))+1e-300;vh=v/rowSums(v)
    zh=lapply(1:length(grid),function(t){
      g=sapply(1:k,function(j) pi0[t,j]*dnorm(y,mu0[t,j],sqrt(sigma20[t,j])))+1e-300
      ;gh=g/rowSums(g)})
    ##M-Step
    lmd1=colSums(vh)/n
    grid=which(lmd1>=lmd_0)
    lmd1=lmd1[grid]
    xgrid=xgrid[grid]
    Kh=Kh[,grid]
    ##mu
    mu=sapply(1:k,function(j){
      sapply(1:length(grid),function(t){
        id=grid[t]
        W=lmd1[t]*zh[[id]][,j]
        mh=local.polynomial.smoother(x,y,xgrid[t],bw,d,W)[,2]
        #sum(W*y)/sum(W)
      })
    })
    ##sigma2
    sigma2=sapply(1:k,function(j){
      sapply(1:length(grid),function(t){
        id=grid[t]
        W=lmd1[t]*zh[[id]][,j]*Kh[,t]
        yh=mu[t,j]
        res2=c((y-yh)^2)
        sig2=sum(W*res2)/sum(W)
      })
    })
    
    prop=sapply(1:k,function(j){
      sapply(1:length(grid),function(t){
        id=grid[t]
        W=zh[[id]][,j]*Kh[,t]
        rho=sum(W)/sum(lmd1[t]*Kh[,t])
      })
    })
    prop=prop/rowSums(prop)
    ##Evaluating convergence
    mix.mu=mix.prop=mix.sigma2=NULL
    for(j in 1:k){
      mix.prop=cbind(mix.prop,approx(xgrid,prop[,j],xout=x,rule=2)$y)
      mix.sigma2=cbind(mix.sigma2,approx(xgrid,sigma2[,j],xout=x,rule=2)$y)
      mix.mu=cbind(mix.mu,approx(xgrid,mu[,j],xout=x,rule=2)$y)
    }
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) mix.prop[,j]*dnorm(y-mix.mu[,j],0,sqrt(mix.sigma2)[,j])))))
    diff=abs(LogLik1-LogLik0)
    dif=max(abs(lmd0[grid]-lmd1))
    lmd0=lmd1
    sigma20=sigma2
    pi0=prop
    mu0=mu
    LogLik0=LogLik1
    count=count+1
    tol=c(tol,LogLik1)
    if(count==1e3) diff=1e-100
  }
  mu1=mix.mu
  out=backfit(y,x,xgrid,d,k,mu1,colMeans(pi0),colMeans(sigma20),bw);
  mu2=out$mu;pi1=out$pi1;sigma21=out$sigma21;
  res=Dist(mu2);mu2=mu2[,res$id];pi1=pi1[res$id];sigma21=sigma21[res$id]
  g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu2[,j],0,sqrt(sigma21)[j]));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  Kh=sapply(x,function(x0) Kern(x,x0,bw))
  df=sum(sapply(1:k,function(j){w=gn[,j]*Kh;S=polynomial.smoother.matrix(x,x,d,w);sum(diag(S))}))
  BIC=BIC(x,bw,k,LL1,df)
  out=list(resp=gn,mix.prop=pi1,mix.mu=mu2,mix.sigma2=sigma21,df=df,LL=LL1,BIC=BIC[2])
  return(out)
}

###Model-based EM  (MB-EM) algorithm
Kernel_Mix_MB_EM=function(x,y,k,bw,d,xgrid,init.model){
  n=length(y)
  ngrid=length(xgrid)
  ##Initial state
  lmd0=rep(1/ngrid,ngrid)
  pi0=init.model$mix.prop
  mu0=init.model$mix.mu;
  sigma20=init.model$mix.sigma2
  ##
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[j])))))
  mu0=sapply(1:k,function(j) approx(x,mu0[,j],xgrid,rule=2)$y)
  pi0=matrix(pi0,ngrid,k,byrow=T)
  sigma20=matrix(sigma20,ngrid,k,byrow=T)
  Kh=sapply(xgrid,function(x0) Kern(x,x0,bw))
  diff=dif=1e6
  count=0
  tol=NULL
  while(diff>1e-10){
    ##E-Step
    v=sapply(1:ngrid,function(t) lmd0[t]*GMM(y,mu0[t,],pi0[t,],sqrt(sigma20)[t,]))+1e-300;vh=v/rowSums(v)
    zh=lapply(1:ngrid,function(t){
      g=sapply(1:k,function(j) pi0[t,j]*dnorm(y,mu0[t,j],sqrt(sigma20[t,j])))+1e-300
      ;gh=g/rowSums(g)})
    ##M-Step
    lmd1=colSums(vh)/n
    ##mu
    mu=sapply(1:k,function(j){
      sapply(1:ngrid,function(t){
        W=lmd1[t]*zh[[t]][,j]*Kh[,t]
        #mh=local.polynomial.smoother(x,y,xgrid[t],bw,d,W)[,2]
        sum(W*y)/sum(W)
      })
    })
    ##sigma2
    sigma2=sapply(1:k,function(j){
      sapply(1:ngrid,function(t){
        W=lmd1[t]*zh[[t]][,j]*Kh[,t]
        yh=mu[t,j]
        res2=c((y-yh)^2)
        sig2=sum(W*res2)/sum(W)
      })
    })
    
    prop=sapply(1:k,function(j){
      sapply(1:ngrid,function(t){
        W=lmd1[t]*zh[[t]][,j]*Kh[,t]
        rho=sum(W)/sum(lmd1[t]*Kh[,t])
      })
    })
    prop=prop/rowSums(prop)
    ##Evaluating convergence
    mix.mu=mix.prop=mix.sigma2=NULL
    for(j in 1:k){
      mix.prop=cbind(mix.prop,approx(xgrid,prop[,j],xout=x,rule=2)$y)
      mix.sigma2=cbind(mix.sigma2,approx(xgrid,sigma2[,j],xout=x,rule=2)$y)
      mix.mu=cbind(mix.mu,approx(xgrid,mu[,j],xout=x,rule=2)$y)
    }
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) mix.prop[,j]*dnorm(y-mix.mu[,j],0,sqrt(mix.sigma2)[,j])))))
    diff=abs(LogLik1-LogLik0)
    dif=max(abs(lmd0-lmd1))
    lmd0=lmd1
    sigma20=sigma2
    pi0=prop
    mu0=mu
    LogLik0=LogLik1
    count=count+1
    tol=c(tol,LogLik1)
    if(count==1e3) diff=1e-100
  }
  mu1=mix.mu
  out=backfit(y,x,xgrid,d,k,mu1,colMeans(pi0),colMeans(sigma20),bw);
  mu2=out$mu;pi1=out$pi1;sigma21=out$sigma21;
  res=Dist(mu2);mu2=mu2[,res$id];pi1=pi1[res$id];sigma21=sigma21[res$id]
  g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu2[,j],0,sqrt(sigma21)[j]));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  Kh=sapply(x,function(x0) Kern(x,x0,bw))
  df=sum(sapply(1:k,function(j){w=gn[,j]*Kh;S=polynomial.smoother.matrix(x,x,d,w);sum(diag(S))}))
  BIC=BIC(x,bw,k,LL1,df)
  out=list(resp=gn,mix.prop=pi1,mix.mu=mu2,mix.sigma2=sigma21,df=df,LL=LL1,BIC=BIC[2])
  return(out)
}

