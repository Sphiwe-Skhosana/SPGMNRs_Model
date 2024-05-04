simulator=function(r,t){
  library(mclust)
  samp=c(2.5e2,5e2,1e3,2e3)
  band0=c(0.05,0.045,0.04,0.035) ##Local constant #tdf= approx. 20 and 25, respectively.
  band1=c(0.055,0.05,0.045,0.03) ##Local linear
  LL1=LL0=out.prop=out.sigma2=out20=out23=out3=out4=out5=out50=out6=NULL
  outN=count=0
  out24=rep(0,2)
  n=samp[t];
  band=rbind(band0,band1)[,t]
  h0=band[1];h1=band[2];
  lmd0=1e-5
  while(count<1){
    fit1=fit1.1=fit2=fit2.1=fit3=fit3.1=fit4=NULL
    x=sort(runif(n))
    k=2
    mix.prop=c(0.65,0.35)
    #m=cbind(1-cos(2*pi*x),exp(2*x)) ##Scenario 1
    #m=cbind(1+cos(2*pi*x),exp(2*x)) ###Scenario 2
    #m=cbind(-cos(2*pi*x),exp(2*x^2)) ##Scenario 3
    z=sample(1:2,n,replace=T,prob=mix.prop)
    y=ifelse(z==1,m[,1]+rnorm(n,0,sqrt(0.09)),m[,2]+rnorm(n,0,sqrt(0.16)))
    true.dist=cond_dist(y,m,mix.prop,c(0.09,0.16))
    true.cdf=CDFmixture(x,y,mix.prop,c(0.09,0.16),m)$Fy
    true.model=list(prop=mix.prop,sigma2=c(0.09,0.16),mu=m)
    truth=list(mu=m,rho=mix.prop,sigma2=c(0.09,0.16))
    xgrid=seq(min(x),max(x),length.out=100)
    try({
      init.model0=initialize.model(x,y,2,method=1,p=1,true.init = truth)
      par(mfrow=c(1,2))
      s0=Sys.time()
      fit1=Kernel_Mix_EM_MB3(x,y,2,h0,0,xgrid,init.model0,lmd_0 = lmd0)
      init.model01=init.model0
      init.model01$mix.mu=fit1$mix.mu;init.model01$mix.prop=fit1$mix.prop;init.model01$mix.sigma2=fit1$mix.sigma2
      fit1.1=Kernel_Mix_EM_MB3(x,y,2,h1,1,xgrid,init.model01,lmd_0 = lmd0)
      fit2=Kernel_Mix_EM.Naive(x,y,2,h0,0,xgrid,init.model0)
      init.model02=init.model0
      init.model02$mix.mu=fit2$mix.mu;init.model02$mix.prop=fit2$mix.prop;init.model02$mix.sigma2=fit2$mix.sigma2
      fit2.1=Kernel_Mix_EM.Naive(x,y,2,h1,1,xgrid,init.model02)
      fit3=Kernel_Mix_EM_loc(x,y,2,h0,0,xgrid,init.model0)
      init.model05=init.model0
      init.model05$mix.mu=fit3$mix.mu;init.model05$mix.prop=fit3$mix.prop;init.model05$mix.sigma2=fit3$mix.sigma2
      fit3.1=Kernel_Mix_EM_loc(x,y,2,h1,1,xgrid,init.model05)
      fit4=init.model0
      print(Sys.time()-s0)
    })
    if(!is.null(fit1)&!is.null(fit1.1)&!is.null(fit2)&!is.null(fit2.1)&!is.null(fit3)&!is.null(fit3.1)&!is.null(fit4)){
      count=count+1
      for(d in 0:1){
        if(d==0){fit1=fit1;fit2=fit2;fit3=fit3}else{fit1=fit1.1;fit2=fit2.1;fit3=fit3.1}
        ##
        out20=rbind(out20,c(perf(fit1$est.dist,true.dist)[1],perf(fit2$est.dist,true.dist)[1],perf(fit3$est.dist,true.dist)[1]))
        out23=rbind(out23,c(perf(fit1$est.cdf,true.cdf)[4],perf(fit2$est.cdf,true.cdf)[4],perf(fit3$est.cdf,true.cdf)[4]))
        out3=rbind(out3,c(mse(fit1,true.model)[1],mse(fit2,true.model)[1],mse(fit3,true.model)[1]))
        out4=rbind(out4,c(mse(fit1,true.model)[2],mse(fit2,true.model)[2],mse(fit3,true.model)[2]))
        out5=rbind(out5,c(mse(fit1,true.model)[3],mse(fit2,true.model)[3],mse(fit3,true.model)[3]))
        out50=rbind(out50,c(mse(fit1,true.model)[4],mse(fit2,true.model)[4],mse(fit3,true.model)[4]))
        out6=rbind(out6,c(adjustedRandIndex(z,ifelse(fit1$resp[,1]>=0.5,1,2)),adjustedRandIndex(z,ifelse(fit2$resp[,1]>=0.5,1,2)),adjustedRandIndex(z,ifelse(fit3$resp[,1]>=0.5,1,2))))
        LL0=rbind(LL0,fit1$LL);LL1=rbind(LL1,fit2$LL)
        pdf(paste0(n,"/d=",d,"/plot-",r,".pdf"),width=10,height=10)
        par(mfrow=c(2,2))
        plot(x,y,pch=19,type="n",ylab=bquote(m[k](x)),main="MBEM")
        for(j in 1:2) lines(x,m[,j],lwd=3,col="black")
        for(j in 1:2) lines(x,fit1$mix.mu[,j],col="red",lwd=3)
        legend("topleft",legend=c("True-function","Estimate"),col=c("black","red"),lwd=rep(3,2),cex=0.9,lty=rep(1,2))
        plot(x,y,pch=19,type="n",ylab=bquote(m[k](x)),main="NaiveEM")
        for(j in 1:2) lines(x,m[,j],lwd=3,col="black")
        for(j in 1:2) lines(x,fit2$mix.mu[,j],col="red",lwd=3)
        legend("topleft",legend=c("True-function","Estimate"),col=c("black","red"),lwd=rep(3,2),cex=0.9,lty=rep(1,2))
        plot(x,y,pch=19,type="n",ylab=bquote(m[k](x)),main="LEM")
        for(j in 1:2) lines(x,m[,j],lwd=3,col="black")
        for(j in 1:2) lines(x,fit3$mix.mu[,j],col="red",lwd=3)
        legend("topleft",legend=c("True-function","Estimate"),col=c("black","red"),lwd=rep(3,2),cex=0.9,lty=rep(1,2))
        plot(x,y,pch=19,type="n",ylab=bquote(m[k](x)),main="Initial")
        for(j in 1:2) lines(x,m[,j],lwd=3,col="black")
        for(j in 1:2) lines(x,fit4$mix.mu[,j],col="red",lwd=3)
        legend("topleft",legend=c("True-function","Estimate"),col=c("black","red"),lwd=rep(2,2),cex=0.9,lty=rep(1,2))
        dev.off()
      }
      
    }
  }
  return(list(out.prop=out.prop,out.sigma2=out.sigma2,out20=out20,out23=out23,out3=out3,
              out4=out4,out5=out5,out50=out50,out6=out6,LL_Naive=LL1))
}

##Computations
for(t in 1:4){
  R=500
  n=c(2.5e2,5e2,1e3,2e3)[t]
  res=lapply(1:R,function(r) simulator(r=r,t=t))
  for(d in 0:1){
  ##d is the degree of the polynomial
  # d=0 : local constant or kernel estimator
  # d=1 : local linear estimator
  out20=out23=out3=out4=out5=out50=out6=matrix(0,R,3)
  LL0=LL1=cbind(1:R,rep(0,R))
  if(d==0){j=1}else{j=2}
  for(i in 1:R){
    out20[i,]=(res[[i]]$out20)[j,];out23[i,]=(res[[i]]$out23)[j,];out3[i,]=(res[[i]]$out3)[j,];out4[i,]=(res[[i]]$out4)[j,];out5[i,]=(res[[i]]$out5)[j,];
    out50[i,]=(res[[i]]$out50)[j,];out6[i,]=(res[[i]]$out6)[j,]
    LL1[i,2]=(res[[i]]$LL_MB)[j,]
  }
  ##Performance measures
  RES=round(rbind(ColStat(out20),ColStat(out23),ColStat(out3),ColStat(out4),ColStat(out5),ColStat(out50),ColStat(out6)),3);
  rownames(RES)=c("RASE_f","","KS","","MSE_prop","","MSE_sig1","","MSE_sig2","","RASE_m","","ARIndx","")
  colnames(RES)=c("MB-EM","NaiveEM","LEM")
  write.csv(RES,paste0(n,"/d=",d,"/R.csv"))
}
}
