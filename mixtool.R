########## comparing mixtool with our model#####
###### simulation #######
for(i in 1:nsim){
  
  mix_data <- rnormmix(n=n, lambda=lamb, mu=med, sigma=sig)###generate the data
  
  ours<-GMM(y=mix_data)
  
  out<-normalmixEM(mix_data, arbvar = FALSE, epsilon = 1e-04,
                   ECM = FALSE)
  
  meds[i,]=c(ours[[2]],out$mu)
  mix[i,]=c(ours[[1]], out$lambda)
  sig2[i,]=c(ours$sigma, out$sigma[1])
}

mmed=matrix(c(apply(meds,2, mean),med), ncol=2, byrow = TRUE)
mixm=matrix(c(apply(mix,2, mean)[c(1,3)],lamb[1]))
sig2m=matrix(c(apply(sig2,2, mean), sig[1]), nrow = 3)
resul=cbind.data.frame(mmed, mixm, sig2m)


colnames(resul)=c("M1","M2","pi","sig2")
rownames(resul)=c("ours","Mixtools", "simu")


write.csv(resul,"mixtool.csv")

####### mixtool data ######
rm(list=ls())
data(faithful)
attach(faithful)

GMM<-function(y, mu, ve, conv.crit){
  
  y=y
  y=as.matrix(y)
  n=length(y)
  
  P=matrix(0.5,n,2)
  
  if(missing(ve)){ve=var(y)/4}
  
  if(missing("mu")){
    mu=c(0,0)  
    mu[1]=quantile(y,0.4)
    mu[2]=quantile(y,0.6)
  }
  
  if(missing("conv.crit")){conv.crit=1e-04}
  
  P1=dnorm(y,mu[1],sqrt(ve))
  P2=dnorm(y,mu[2],sqrt(ve))
  P[,1]=P1/(P1+P2)
  P[,2]=P2/(P1+P2)
  
  
  iter=0
  maxiter=300
  
  repeat{
    W1=P
    JJ=crossprod(P,P)
    Jy=crossprod(P,y)
    C1=JJ
    
    
    C2=solve(C1)
    sol=C2%*%Jy
    
    m1=sol[1]
    m2=sol[2]
    
    
    ma=m1
    mb=m2
    
    pred=P%*%c(m1,m2)
    e=y-pred
    
    dd=W1%*%C2%*%t(W1)
    
    ve1= (t(e)%*%e+sum(diag(dd*c(ve))))/n ####duvida aqui
    
    
    dif=max(abs(ve-ve1))
    
    ve=ve1
    iter=iter+1
    
    pi=sum(P[,1])/n
    P1=pi*dnorm(y,ma,sqrt(ve))
    P2=(1-pi)*dnorm(y,mb,sqrt(ve))
    
    P[,1]=P1/(P1+P2)
    P[,2]=P2/(P1+P2)  
    
    ###print(c(m1,m2,iter))
    
    ##if (iter>maxiter) dif=0
    if((dif<var(y)*conv.crit)|(iter==maxiter)) break
  }
  
  colnames(P)<-c("prob_subP1","prob_subP2")
  prop<-c(pi, (1-pi))
  med<-c(m1,m2)
  
  saida<-list(lmabda=unlist(prop), mu=unlist(med), sigma=unlist(sqrt(ve)))
  return(saida)
}


out1<-normalmixEM(waiting, arbvar = FALSE, epsilon = 1e-03)
out2<-GMM(waiting)

summary(out1)
out2
