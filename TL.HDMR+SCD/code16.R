################### varying vNKdt_rho ##################################

library(MASS)
library(ncvreg)
library(mrcovreg)
library(glmnet)

library(foreach)
library(doParallel)
library(doMC)

source('Compared_Methods.R')
source('DataGeneration.R')
source('cisMRncvreg.R')
source('Main.R')
source('SelectAK.R')

# par in source dataset
Nx <- 500
Ny <- 5000
P <- 100
hx2 <- 0.05
Pcausal <- P*0.1
cv.lambd <- exp(c(seq(-20,5,0.1)))
mc <- 10

# causal effect
bl <- 0.05
bu <- 0.05
hxy2 <- c(runif(Pcausal,bl,bu)^2,rep(0,(P-Pcausal)))
sqrt(hxy2)

g <- 120
rho <- 0.8
#vXcor <- c(0,0.2,0.5,0.8)
Xcor <- 0

# par in aux dataset
# delta is a sparse high-dimensional vector whose L1-norm is no larger than h

# h<s*sqrt(log(p/n0))
#Pcausal*sqrt(log(P/g))
NKdt_distorb_prop <- c(rep(0.1,5),rep(0.5,5))
distorb_value <- c(rep(0.05,5),rep(0.2,5))^2
#sqrt(vdistorb_value)*P*NKdt_distorb_prop #delta

NKdt <- 10
Nx_kb <- floor(runif(NKdt,500,1000))
Ny_kb <- floor(runif(NKdt,5000,10000))

# correlation between multiple datasets
NKdt_rho <- 0

r.A0 <- 3
vC0 <- seq(0.1,2,0.1)

#i1=1


for(i1 in 1:length(vC0)){
  
  C0 <- vC0[i1]
  
  cat('vC0=',vC0[i1],'\n')
  
  # par in aux dataset
  g_kb <- floor(runif(NKdt,50,120))
  rho_kb <- rnorm(NKdt,rho,0.1)
  rho_kb[rho_kb<0] <- 0
  rho_kb[rho_kb>0.9] <- 0.9
  Xcor_kb <- rnorm(NKdt,Xcor,0.1)
  Xcor_kb[Xcor_kb<0] <- 0
  Xcor_kb[Xcor_kb>0.9] <- 0.9
  
  
  result <- NULL
  cl <- makeCluster(mc)
  registerDoParallel(cl)
  #registerDoMC(mc)
  result <- foreach(ii=1:200,
                    .packages = c("ncvreg",
                                  "MASS",
                                  "mrcovreg",'glmnet'),
                    .errorhandling =  "remove"
  ) %dopar% {
    Main(ii,Nx,Ny,g,rho,P,hx2,hxy2,cv.lambd,Xcor,
         NKdt,NKdt_distorb_prop,distorb_value,
         Nx_kb,Ny_kb,g_kb,rho_kb,Xcor_kb,NKdt_rho,r.A0,C0)
  }
  stopCluster(cl)
  save(result,file=paste0('code16_g_',g,'_P_',P,'_rho_',rho,'_Xcor_',Xcor,
                          '_bl_',bl,'_bu_',bu,'_NKdt_',NKdt, '_NKdt_rho_',
                          NKdt_rho,'_C0_',C0,'_0115.Rdata'))
  
  rm(result)
  gc()

}
