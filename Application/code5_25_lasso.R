############### transfer learning #####################
library(glmnet)
library(ncvreg)
library(mrcovreg)

library(foreach)
library(doParallel)
library(doMC)

source('TransCondi.R')
source('SelectDt.R')
source('cisMRncvreg.R')
source('Compared_Methods.R')
source('RemoveIV.R')
load('code41.fdata_NK1all_top.Rdata')

Ancesname <- names(fdata_NK.Listall)
outname <- names(fdata_NK.Listall[[1]])

resultall <- list()
IV.presso <- list()

oua=2
ouy=5

######

fdata1 <- fdata_NK1.Listall[[oua]][[ouy]]

fdata_NK1 <- list()
ok=1
for(oj in 1:length(Ancesname)){
  if(oj==oua) next
  fdata_NK1[[ok]] <- fdata_NK1.Listall[[oj]][[ouy]]
  ok <- ok+1
}  

cv.lambd <- exp(c(seq(-10,-2,1)))



##################### Lasso remove invalid IV #############
cat('Lasso remove invalid IV...','\n')
fdata1j <- RemoveIV(fdata1,
                    bx_corMatrxi=cor(fdata1$betaGX_new),
                    methth='Lasso',
                    cv.lambd)
fdata1j <- fdata1j$fdata2
fdata_NK1o <- list()
for(kpi in 1:length(fdata_NK1)){
  onon <- RemoveIV(fdata_NK1[[kpi]],
                   bx_corMatrxi=cor(fdata_NK1[[kpi]]$betaGX_new),
                   methth='Lasso',cv.lambd)
  fdata_NK1o[[kpi]] <- onon$fdata2
}
## Transfer learning+removeIV
cat('Transfer learning...','\n')
fdata_condi <- DATA_condi(fdata1j,fdata_NK1o)
cisMVMR.Trans.lasso.presso <- Oracle.Trans.lasso(fdata_condi$betaGX_A0,fdata_condi$betaGY_A0,
                                                 fdata_condi$betaGX_AK,fdata_condi$betaGY_AK,
                                                 cv.lambd,Methodname='Lasso')
IV.Lasso <- list(fdata1j=fdata1j,
                 fdata_NK1o=fdata_NK1o)

save(cisMVMR.Trans.lasso.presso,IV.Lasso,file=paste0('code5_resultall_oua_',oua,'_ouy_',ouy,'_transLASSO.Rdata'))


