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

# par in source dataset
Nx <- 500
Ny <- 5000
P <- 100
hx2 <- 0.05
Pcausal <- P*0.1
cv.lambd <- exp(c(seq(-15,5,0.1)))
mc <- 10

# causal effect
bl <- 0.05
bu <- 0.05
hxy2 <- c(runif(Pcausal,bl,bu)^2,rep(0,(P-Pcausal)))
sqrt(hxy2)

vg <- c(50,80,120)
vrho <- c(0.2,0.5,0.8)
#vXcor <- c(0,0.2,0.5,0.8)
vXcor <- 0

# par in aux dataset
# delta is a sparse high-dimensional vector whose L1-norm is no larger than h
# h<s*sqrt(log(p/n0))
Pcausal*sqrt(log(P/vg))
NKdt_distorb_prop <- 0.1
distorb_value <- 0.05^2
sqrt(distorb_value)*P*NKdt_distorb_prop #delta
NKdt <- 3
Nx_kb <- floor(runif(NKdt,500,1000))
Ny_kb <- floor(runif(NKdt,5000,10000))

# correlation between multiple datasets
vNKdt_rho <- c(0.2,0.5,0.8)

# i1=3
# i2=2
# i3=i4=1

for(i1 in 1:length(vNKdt_rho)){
  
  NKdt_rho <- vNKdt_rho[i1]
  
  for(i2 in 1:length(vg)){
    
    g <- vg[i2]
    
    for(i3 in 1:length(vrho)){
      
      rho <- vrho[i3]
      
      for(i4 in 1:length(vXcor)){
        
        cat('vNKdt_rho=',vNKdt_rho[i1],' vg=',vg[i2],' vrho=',vrho[i3],' vXcor=',vXcor[i4],'\n')
        
        Xcor <- vXcor[i4]
        
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
                                        "mrcovreg",'glmnet')
                          #.errorhandling =  "remove"
        ) %dopar% {
          Main(ii,Nx,Ny,g,rho,P,hx2,hxy2,cv.lambd,Xcor,
               NKdt,NKdt_distorb_prop,distorb_value,
               Nx_kb,Ny_kb,g_kb,rho_kb,Xcor_kb,NKdt_rho)
        }
        stopCluster(cl)
        save(result,file=paste0('code11_g_',g,'_P_',P,'_rho_',rho,'_Xcor_',Xcor,
                                '_bl_',bl,'_bu_',bu,'_distorb_value_',sqrt(distorb_value),
                                '_NKdt_',NKdt, '_NKdt_rho_',NKdt_rho,'_0114.Rdata'))
        
        rm(result)
        gc()
        
        
      }
      
    }

  }

}


################### varying vXcor ##################################

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

# par in source dataset
Nx <- 500
Ny <- 5000
P <- 100
hx2 <- 0.05
Pcausal <- P*0.1
cv.lambd <- exp(c(seq(-15,5,0.1)))
mc <- 20

# causal effect
bl <- 0.05
bu <- 0.05
hxy2 <- c(runif(Pcausal,bl,bu)^2,rep(0,(P-Pcausal)))
sqrt(hxy2)

vg <- c(50,80,120)
vrho <- c(0.2,0.5,0.8)
vXcor <- c(0,0.2,0.5,0.8)
#vXcor <- 0

# par in aux dataset
# delta is a sparse high-dimensional vector whose L1-norm is no larger than h
# h<s*sqrt(log(p/n0))
Pcausal*sqrt(log(P/vg))
NKdt_distorb_prop <- 0.1
distorb_value <- 0.05^2
sqrt(distorb_value)*P*NKdt_distorb_prop #delta
NKdt <- 3
Nx_kb <- floor(runif(NKdt,500,1000))
Ny_kb <- floor(runif(NKdt,5000,10000))

# correlation between multiple datasets
# vNKdt_rho <- c(0,0.2,0.5,0.8)
vNKdt_rho <- 0

# i1=3
# i2=2
# i3=i4=1

for(i1 in 1){
  
  NKdt_rho <- vNKdt_rho[i1]
  
  for(i2 in 1:length(vg)){
    
    g <- vg[i2]
    
    for(i3 in 1:length(vrho)){
      
      rho <- vrho[i3]
      
      for(i4 in 1:length(vXcor)){
        
        cat('vNKdt_rho=',vNKdt_rho[i1],' vg=',vg[i2],' vrho=',vrho[i3],' vXcor=',vXcor[i4],'\n')
        
        Xcor <- vXcor[i4]
        
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
                                        "mrcovreg",'glmnet')
                          #.errorhandling =  "remove"
        ) %dopar% {
          Main(ii,Nx,Ny,g,rho,P,hx2,hxy2,cv.lambd,Xcor,
               NKdt,NKdt_distorb_prop,distorb_value,
               Nx_kb,Ny_kb,g_kb,rho_kb,Xcor_kb,NKdt_rho)
        }
        stopCluster(cl)
        save(result,file=paste0('code11_g_',g,'_P_',P,'_rho_',rho,'_Xcor_',Xcor,
                                '_bl_',bl,'_bu_',bu,'_distorb_value_',sqrt(distorb_value),
                                '_NKdt_',NKdt, '_NKdt_rho_',NKdt_rho,'_0114.Rdata'))
        
        rm(result)
        gc()
        
        
      }
      
    }
    
  }
  
}