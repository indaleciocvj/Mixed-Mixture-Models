

########## Mixture Mixed Models: Biennial Growth as a Latent Variable in Coffee Bean Progenies #######

## Code to run simulation analysis as showed on table 4
rm(list = ls())

library(mixtools)
library(lme4)

###### Function for mixture model

GMM_Uni <- function(gen = gen, rep = rep, response = yield, maxiter = 2000, conv.crit = 1e-5){
  
  maxiter <- maxiter   
  gen <- gen
  rep <- rep
  y   <- response
  n   <- length(y)  
  
  x1    <- factor(rep)
  z     <- factor(gen)
  nefix <- as.numeric(max(levels(x1)))
  negen <- as.numeric(length(levels(z)))
  x1    <- model.matrix(~x1 -1)
  z     <- model.matrix(~z -1)
  
  TT    <- t(c(0,0,rep(1,ncol(x1)),rep(0,ncol(z))))
  TT    <- crossprod(TT,TT)
  
  XX    <- crossprod(x1)
  XZ    <- crossprod(x1,z)
  ZZ    <- crossprod(z)
  Xy    <- crossprod(x1,y)
  Zy    <- crossprod(z,y)
  P     <- matrix(0.5,n,2)
  ve    <- var(y)
  va    <- var(y)/2
  
  m1    <- quantile(y,0.4)
  m2    <- quantile(y,0.6)
  P1    <- dnorm(y,m1,sqrt(ve))
  P2    <- dnorm(y,m2,sqrt(ve))
  P[,1] <- P1/(P1+P2)
  P[,2] <- P2/(P1+P2)
  ####W1    <- cbind(P,x1,z)
  
  iter  <- 0
  
  
  repeat{
    
    W1    <- cbind(P,x1,z)
    
    G <- diag(ncol(z))*c(ve/va)
    
    JJ <- crossprod(P)
    XJ <- crossprod(x1,P)
    ZJ <- crossprod(z,P)
    
    Jy <- crossprod(P,y)
    
    C1 <- cbind(JJ,t(XJ),t(ZJ))
    C2 <- cbind(XJ,XX,XZ)
    C3 <- cbind(ZJ,t(XZ),ZZ+G)
    
    C  <- rbind(C1,C2,C3)+TT
    so  <- as.matrix(c(Jy,Xy,Zy))
    C2  <- solve(C)
    sol <- C2%*%so
    
    m1 <- sol[1]
    m2 <- sol[2]
    beta <- sol[3:(2+ncol(x1))]
    a <- sol[(2+ncol(x1)+1):(2+ncol(x1)+ncol(z))]
    ma <- m1+x1%*%beta+z%*%a
    mb <- m2+x1%*%beta+z%*%a
    
    pred <- P%*%c(m1,m2)+x1%*%beta+z%*%a
    e <- y-pred
    dd <- W1%*%C2%*%t(W1)
    
    ve1 <- (t(e)%*%e+sum(diag(dd*c(ve))))/n
    va1 <- (t(a)%*%a+sum(diag(C2[(2+ncol(x1)+1):(2+ncol(x1)+ncol(z)),(2+ncol(x1)+1):(2+ncol(x1)+ncol(z))]*c(ve))))/ncol(z)
    
    dif <- max(abs(ve-ve1),abs(va-va1))
    
    iter <- iter+1
    
    va <- va1
    ve <- ve1
    pi <- sum(P[,1])/n
    P1 <- pi*dnorm(y,ma,sqrt(ve))
    
    P2    <- (1-pi)*dnorm(y,mb,sqrt(ve))
    P[,1] <- P1/(P1+P2)
    P[,2] <- P2/(P1+P2)
    
    cat("Iteration:", iter,", Mean_1:", m1,", Mean_2:", m2, "\n")
    
    if(iter>maxiter)dif <- 0
    if((dif<var(y)*conv.crit)|(iter==maxiter)) break
    
  }
  
  #Va <- unlist(va); names(Va) <- "Va"
  colnames(P) <- c("prob group1","prob group2")
  Prob <- P
  blup<-a
  SE<-sqrt(diag(C2[(2+ncol(x1)+1):(2+ncol(x1)+ncol(z)),(2+ncol(x1)+1):(2+ncol(x1)+ncol(z))]*c(ve)))
  meds<-c(m1,m2)
  MixPar <- pi
  Out <- list(Va = unlist(va), Ve = unlist(ve), Prob = Prob, blup=blup,
                SE.blup=SE, means=meds, MixPar = MixPar)
  return(Out)
  
}



######### Create simulation scenarios

mean1 <- c(5)
mean2 <- c(7, 20)
piMix <- c(0.2, 0.5)
ngen <- 100

scen <- as.data.frame(expand.grid(mean1, piMix, mean2))
colnames(scen) <- c("Means1", "Pmix", "Means2")

Rep <- c(2, -2)##replication effect


SigG <- 2 ###genetic variance
Sige <- ((1-0.5)/0.5)*SigG ##error variance


########### simulating the phenotype
nsim <-1000

Out_Mixture <- Out_Unimod <- Out_Simu <- c()
for(i in 1:nsim){
  for(s in 1:nrow(scen)){
    
  
  ##### Simulating phenotype
  blups.Sim <- rnorm(ngen, 0, sd=sqrt(SigG))
  
  pop1 <- scen[s,"Means1"]+ blups.Sim[1:round(ngen*(scen[s,2]))] ## population 1
  pop2 <- scen[s,"Means2"]+ blups.Sim[round(ngen*(scen[s,2])+1):ngen] ## population 2
  
  pop <- c(pop1, pop2)
  
  y1 <- data.frame(geno= c(1:ngen), Rep=1, Pheno= pop+Rep[1])
  y2 <- data.frame(geno= c(1:ngen), Rep=2, Pheno= pop+Rep[2])
  
  SimPheno <- rbind(y1,y2)
  err <- rnorm(2*ngen, 0, sd=sqrt(Sige))
  SimPheno$Pheno <- SimPheno$Pheno+ err
  
  SimPheno$Rep <- as.factor(SimPheno$Rep)
  SimPheno$geno <- as.factor(SimPheno$geno)
  
  #### run models
  
  M.Mixture <- GMM_Uni(gen=SimPheno$geno, rep=SimPheno$Rep, response=SimPheno$Pheno, maxiter = 2000, conv.crit = 1e-6)
  Unimod <- lmer(Pheno ~ Rep+ (1|geno), data=SimPheno)
  
  rm(SimPheno)
  
  ######Heritabilaires
  
  Vcomp_Unimod <- as.data.frame(VarCorr(Unimod))
  
  h2_Mixture <- M.Mixture$Va/c(M.Mixture$Va+ M.Mixture$Ve)
  h2_Unimod <- Vcomp_Unimod[1,4]/c(Vcomp_Unimod[1,4]+ Vcomp_Unimod[2,4])
  h2_Sim <- var(blups.Sim)/(var(blups.Sim)+ var(err))
  #######summarize results
  
  Out_Mixture <- rbind(Out_Mixture, data.frame(means= t(M.Mixture$means), GenVar = M.Mixture$Va, MixPi = c(M.Mixture$MixPar),
                                       crr.blup= cor(M.Mixture$blup, blups.Sim) ,Scen = paste0("Scen",s)))
  
  Out_Unimod <- rbind(Out_Unimod, data.frame(GenVar = Vcomp_Unimod[1,4], h2= h2_Unimod, crr.blup = cor(ranef(Unimod)$geno,blups.Sim),
                                   Scen = paste0("Scen",s)))
  
  Out_Simu <- rbind(Out_Simu, data.frame(GenVar= var(blups.Sim), h2=h2_Sim, Scen = paste0("Scen",s)))
  
  }
}


####################### Mixture Models results
### Mean of paramer estimates
Resul_Mixture <- do.call(rbind, lapply(split(Out_Mixture, Out_Mixture$Scen), FUN = function(i) apply(i[,1:5], 2, mean)))
Resul_Mixture <-  as.data.frame(Resul_Mixture); Resul_Mixture$Scen <- rownames(Resul_Mixture)


### Standard deviations for parameter estimates on the simalation
Sd_Mixture <- do.call(rbind, lapply(split(Out_Mixture, Out_Mixture$Scen), FUN = function(i) apply(i[,1:5], 2, sd)))
Sd_Mixture <-  as.data.frame(Sd_Mixture); Sd_Mixture$Scen <- rownames(Sd_Mixture)  
  
##save results
#write.csv(Resul_Mixture, "Resul_Mixture.csv", row.names = FALSE)
#write.csv(Sd_Mixture, "Sd_Mixture.csv", row.names = FALSE)



##################### Mixed model
### Mean of parameter estiamtes
Resul_Unimod <- do.call(rbind, lapply(split(Out_Unimod, Out_Unimod$Scen), FUN = function(i) apply(i[,1:3], 2, mean)))
Resul_Unimod <- as.data.frame(Resul_Unimod); Resul_Unimod$Scen <- rownames(Resul_Unimod)

### Standard deviations for parameter estimates on the simalation
Sd_Unimod <- do.call(rbind, lapply(split(Out_Unimod, Out_Unimod$Scen), FUN = function(i) apply(i[,1:3], 2, sd)))
Sd_Unimod <- as.data.frame(Sd_Unimod); Sd_Unimod$Scen <- rownames(Sd_Unimod)

##save results
#write.csv(Resul_Unimod, "Resul_Unimod.csv", row.names = FALSE)
#write.csv(Sd_Unimod, "Sd_Unimod.csv", row.names = FALSE)




######## Realized simulated values
Real_Simu <- do.call(rbind, lapply(split(Out_Simu, Out_Simu$Scen), FUN = function(i) apply(i[,1:2], 2, mean)))
Real_Simu <- as.data.frame(Real_Simu); Real_Simu$Scen <- rownames(Real_Simu)

##save results
#write.csv(Resul_Mixture, "Real_Simu.csv", row.names = FALSE)



