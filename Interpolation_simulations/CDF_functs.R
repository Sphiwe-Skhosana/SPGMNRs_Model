##Calculating the CDF of a conditional probability distribution
CDFmixture=function(x,y,prop,sigma2,mu){
  n=length(y)
  k=length(prop)
  mu=mu[order(y),]
  y=sort(y);
  
  PDF=rep(0,n)
  for(j in 1:k){
    PDF=PDF+prop[j]*pnorm(y-mu[,j],0,sqrt(sigma2[j]))
  }
  PDF=PDF/sum(PDF)
  return(list(x=y,Fy=cumsum(PDF),fy=PDF))
}

##Binning a variable based on conditional probabilities

Binmixture=function(y,prop,sigma2,mu,nbins=5){
  n=length(y)
  mu=mu[order(y),]
  y=sort(y);
  k=length(prop)
  PDF=rep(0,n)
  for(j in 1:k){
    PDF=PDF+prop[j]*pnorm(y,mu[,j],sqrt(sigma2[j]))
  }
  PDF=PDF/sum(PDF)
  f=cut(y,nbins); bins=unique(f);
  freq=sapply(1:nbins,function(i) n*sum(PDF[which(f==bins[i])]))
  return(list(bins=bins,y=y,freq=freq))
}
