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

oua=1
ouy=1

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

cat('selecting datasets...','\n')
res.select.all <- list()
res.select <- SelectDt.AK(fdata1=fdata1,
                          fdata_NK1=fdata_NK1,
                          cv.lambd=cv.lambd,
                          Methodname='MCP',
                          r.A0=3,
                          C0=1.6)
res.select.all[[1]] <- res.select$Ind.loc

res.select <- SelectDt.AK(fdata1=fdata1,
                          fdata_NK1=fdata_NK1,
                          cv.lambd,
                          Methodname='Lasso',
                          r.A0=3,
                          C0=1.6)
res.select.all[[2]] <- res.select$Ind.loc

res.select <- SelectDt.AK(fdata1=fdata1,
                          fdata_NK1=fdata_NK1,
                          cv.lambd,
                          Methodname='ElasticNet',
                          r.A0=3,
                          C0=1.6)
res.select.all[[3]] <- res.select$Ind.loc


save(res.select.all,file=paste0('code5_resultall_oua_',oua,'_ouy_',ouy,'_selectIV.Rdata'))
