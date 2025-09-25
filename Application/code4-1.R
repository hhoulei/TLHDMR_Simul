
################### 减少IV个数 ###############################

library(R.utils)
library(data.table)
dt1 <- fread('PMID36635386_studies_export.tsv')

ancestry <- strsplit(dt1$discoverySampleAncestry,' ')
N <- lapply(ancestry, function(x) x[1])
N <- as.numeric(unlist(N))

ancestry1 <- lapply(ancestry, function(x) paste0(unlist(x[-1]),collapse = ' '))
ancestry1 <- unlist(ancestry1)

dt1$ancestry <- ancestry1
dt1$N <- N

table(dt1$ancestry)

l1 <- table(dt1$reportedTrait)
sum(l1==4)
choname <- names(l1)[l1==4]
dt1 <- dt1[dt1$reportedTrait %in% choname, ]
table(dt1$ancestry)

#################### South Asian ###########################

SA.loc <- dt1$ancestry=='South Asian'
dt_ori <- dt1[SA.loc,]
load('code21.SNP_top10_SA.Rdata')

load('code4.SA.552.Rdata')
dt_SA <- dt2

expname <- dt_SA$reportedTrait

SNP_dt_uni <- list()
for(i in 1:length(expname)){
  
  cat('i=',i,'\n')
  
  once.name <- expname[i]
  once.loc <- which(dt_ori$reportedTrait==once.name)
  
  once.SNP <- SNP_finaldt_SA_out[[i]]
  once.SNP <- once.SNP[once.SNP$rsid %in% SNP_top10_SA[[once.loc]]$rsid,]
  once.SNP$p_value <- as.numeric(once.SNP$p_value)
  once.SNP <- once.SNP[which.min(once.SNP$p_value),]
  
  SNP_dt_uni[[i]] <- once.SNP
}
names(SNP_dt_uni) <- expname
all(unlist(lapply(SNP_dt_uni, function(x) nrow(x)))==1)

dt_SA$finalSNP <- unlist(lapply(SNP_dt_uni, function(x) x$rsid))

finalSNP.all <- unique(dt_SA$finalSNP)

SNP_finaldt_SA_out <- lapply(SNP_finaldt_SA_out, function(x) x[x$rsid %in% finalSNP.all,])
Stroke_finaldt_SA_out <- lapply(Stroke_finaldt_SA_out,function(x) x[x$rsid %in% finalSNP.all,])

save(dt_SA,SNP_finaldt_SA_out,Stroke_finaldt_SA_out,file='code41top.SA.552.Rdata')



#################### African unspecified ###########################

SA.loc <- dt1$ancestry=='African unspecified'
dt_ori <- dt1[SA.loc,]
load('code21.SNP_top10_AFR.Rdata')

load('code4.AFR.552.Rdata')
dt_AFR <- dt2

expname <- dt_AFR$reportedTrait

SNP_dt_uni <- list()
for(i in 1:length(expname)){
  
  cat('i=',i,'\n')
  
  once.name <- expname[i]
  once.loc <- which(dt_ori$reportedTrait==once.name)
  
  once.SNP <- SNP_finaldt_AFR_out[[i]]
  once.SNP <- once.SNP[once.SNP$rsid %in% SNP_top10_AFR[[once.loc]]$rsid,]
  once.SNP$p_value <- as.numeric(once.SNP$p_value)
  once.SNP <- once.SNP[which.min(once.SNP$p_value),]
  
  SNP_dt_uni[[i]] <- once.SNP
}
names(SNP_dt_uni) <- expname
all(unlist(lapply(SNP_dt_uni, function(x) nrow(x)))==1)

dt_AFR$finalSNP <- unlist(lapply(SNP_dt_uni, function(x) x$rsid))

finalSNP.all <- unique(dt_AFR$finalSNP)

SNP_finaldt_AFR_out <- lapply(SNP_finaldt_AFR_out, function(x) x[x$rsid %in% finalSNP.all,])
Stroke_finaldt_AFR_out <- lapply(Stroke_finaldt_AFR_out,function(x) x[x$rsid %in% finalSNP.all,])

save(dt_AFR,SNP_finaldt_AFR_out,Stroke_finaldt_AFR_out,file='code41top.AFR.552.Rdata')




################### East Asian ###########################

SA.loc <- dt1$ancestry=='East Asian'
dt_ori <- dt1[SA.loc,]
load('code21.SNP_top10_EA.Rdata')

load('code4.EA.552.Rdata')
dt_EA <- dt2

expname <- dt_EA$reportedTrait

SNP_dt_uni <- list()
for(i in 1:length(expname)){
  
  cat('i=',i,'\n')
  
  once.name <- expname[i]
  once.loc <- which(dt_ori$reportedTrait==once.name)
  
  once.SNP <- SNP_finaldt_EA_out[[i]]
  once.SNP <- once.SNP[once.SNP$rsid %in% SNP_top10_EA[[once.loc]]$rsid,]
  once.SNP$p_value <- as.numeric(once.SNP$p_value)
  once.SNP <- once.SNP[which.min(once.SNP$p_value),]
  
  SNP_dt_uni[[i]] <- once.SNP
}
names(SNP_dt_uni) <- expname
all(unlist(lapply(SNP_dt_uni, function(x) nrow(x)))==1)

dt_EA$finalSNP <- unlist(lapply(SNP_dt_uni, function(x) x$rsid))

finalSNP.all <- unique(dt_EA$finalSNP)

SNP_finaldt_EA_out <- lapply(SNP_finaldt_EA_out, function(x) x[x$rsid %in% finalSNP.all,])
Stroke_finaldt_EA_out <- lapply(Stroke_finaldt_EA_out,function(x) x[x$rsid %in% finalSNP.all,])

save(dt_EA,SNP_finaldt_EA_out,Stroke_finaldt_EA_out,file='code41top.EA.552.Rdata')



################### European ###########################

SA.loc <- dt1$ancestry=='European'
dt_ori <- dt1[SA.loc,]
load('code21.SNP_top10_EUR.Rdata')

load('code4.EUR.552.Rdata')
dt_EUR <- dt2

expname <- dt_EUR$reportedTrait

SNP_dt_uni <- list()
for(i in 1:length(expname)){
  
  cat('i=',i,'\n')
  
  once.name <- expname[i]
  once.loc <- which(dt_ori$reportedTrait==once.name)
  
  once.SNP <- SNP_finaldt_EUR_out[[i]]
  once.SNP <- once.SNP[once.SNP$rsid %in% SNP_top10_EUR[[once.loc]]$rsid,]
  once.SNP$p_value <- as.numeric(once.SNP$p_value)
  once.SNP <- once.SNP[which.min(once.SNP$p_value),]
  
  SNP_dt_uni[[i]] <- once.SNP
}
names(SNP_dt_uni) <- expname
all(unlist(lapply(SNP_dt_uni, function(x) nrow(x)))==1)

dt_EUR$finalSNP <- unlist(lapply(SNP_dt_uni, function(x) x$rsid))

finalSNP.all <- unique(dt_EUR$finalSNP)

SNP_finaldt_EUR_out <- lapply(SNP_finaldt_EUR_out, function(x) x[x$rsid %in% finalSNP.all,])
Stroke_finaldt_EUR_out <- lapply(Stroke_finaldt_EUR_out,function(x) x[x$rsid %in% finalSNP.all,])

save(dt_EUR,SNP_finaldt_EUR_out,Stroke_finaldt_EUR_out,file='code41top.EUR.552.Rdata')


########################################################

load('code41top.SA.552.Rdata')
load('code41top.AFR.552.Rdata')
load('code41top.EA.552.Rdata')
load('code41top.EUR.552.Rdata')

loc <- order(dt_SA$reportedTrait)
dt_SA <- dt_SA[loc,]
SNP_finaldt_SA_out <- SNP_finaldt_SA_out[loc]

loc <- order(dt_AFR$reportedTrait)
dt_AFR <- dt_AFR[loc,]
SNP_finaldt_AFR_out <- SNP_finaldt_AFR_out[loc]

loc <- order(dt_EA$reportedTrait)
dt_EA <- dt_EA[loc,]
SNP_finaldt_EA_out <- SNP_finaldt_EA_out[loc]

loc <- order(dt_EUR$reportedTrait)
dt_EUR <- dt_EUR[loc,]
SNP_finaldt_EUR_out <- SNP_finaldt_EUR_out[loc]

sum(dt_SA$reportedTrait==dt_AFR$reportedTrait)
sum(dt_SA$reportedTrait==dt_EA$reportedTrait)
sum(dt_SA$reportedTrait==dt_EUR$reportedTrait)

sum(names(SNP_finaldt_SA_out)==names(SNP_finaldt_AFR_out))
sum(names(SNP_finaldt_SA_out)==names(SNP_finaldt_EA_out))
sum(names(SNP_finaldt_SA_out)==names(SNP_finaldt_EUR_out))


aa <- SNP_finaldt_SA_out[[1]]$rsid
all(unlist(lapply(SNP_finaldt_SA_out,function(x) sum(x$rsid==aa)))==523)
all(unlist(lapply(Stroke_finaldt_SA_out,function(x) sum(x$rsid==aa)))==523)

aa <- SNP_finaldt_AFR_out[[1]]$rsid
all(unlist(lapply(SNP_finaldt_AFR_out,function(x) sum(x$rsid==aa)))==538)
all(unlist(lapply(Stroke_finaldt_AFR_out,function(x) sum(x$rsid==aa)))==538)

aa <- SNP_finaldt_EA_out[[1]]$rsid
all(unlist(lapply(SNP_finaldt_EA_out,function(x) sum(x$rsid==aa)))==536)
all(unlist(lapply(Stroke_finaldt_EA_out,function(x) sum(x$rsid==aa)))==536)

aa <- SNP_finaldt_EUR_out[[1]]$rsid
all(unlist(lapply(SNP_finaldt_EUR_out,function(x) sum(x$rsid==aa)))==411)
all(unlist(lapply(Stroke_finaldt_EUR_out,function(x) sum(x$rsid==aa)))==411)


GXbetaList <- list()
seGXbetaList <- list()
GYbetaList <- list()
seGYbetaList <- list()

GXbetaList[[1]] <- matrix(unlist(lapply(SNP_finaldt_SA_out,function(x) x$beta)),ncol = length(SNP_finaldt_SA_out))
seGXbetaList[[1]] <- matrix(unlist(lapply(SNP_finaldt_SA_out,function(x) x$standard_error)),ncol = length(SNP_finaldt_SA_out))
GYbetaList[[1]] <- matrix(unlist(lapply(Stroke_finaldt_SA_out,function(x) x$beta)),ncol = length(Stroke_finaldt_SA_out))
seGYbetaList[[1]] <- matrix(unlist(lapply(Stroke_finaldt_SA_out,function(x) x$standard_error)),ncol = length(Stroke_finaldt_SA_out))

GXbetaList[[2]] <- matrix(unlist(lapply(SNP_finaldt_AFR_out,function(x) x$beta)),ncol = length(SNP_finaldt_AFR_out))
seGXbetaList[[2]] <- matrix(unlist(lapply(SNP_finaldt_AFR_out,function(x) x$standard_error)),ncol = length(SNP_finaldt_AFR_out))
GYbetaList[[2]] <- matrix(unlist(lapply(Stroke_finaldt_AFR_out,function(x) x$beta)),ncol = length(Stroke_finaldt_AFR_out))
seGYbetaList[[2]] <- matrix(unlist(lapply(Stroke_finaldt_AFR_out,function(x) x$standard_error)),ncol = length(Stroke_finaldt_AFR_out))

GXbetaList[[3]] <- matrix(unlist(lapply(SNP_finaldt_EA_out,function(x) x$beta)),ncol = length(SNP_finaldt_EA_out))
seGXbetaList[[3]] <- matrix(unlist(lapply(SNP_finaldt_EA_out,function(x) x$standard_error)),ncol = length(SNP_finaldt_EA_out))
GYbetaList[[3]] <- matrix(unlist(lapply(Stroke_finaldt_EA_out,function(x) x$beta)),ncol = length(Stroke_finaldt_EA_out))
seGYbetaList[[3]] <- matrix(unlist(lapply(Stroke_finaldt_EA_out,function(x) x$standard_error)),ncol = length(Stroke_finaldt_EA_out))

GXbetaList[[4]] <- matrix(unlist(lapply(SNP_finaldt_EUR_out,function(x) x$beta)),ncol = length(SNP_finaldt_EUR_out))
seGXbetaList[[4]] <- matrix(unlist(lapply(SNP_finaldt_EUR_out,function(x) x$standard_error)),ncol = length(SNP_finaldt_EUR_out))
GYbetaList[[4]] <- matrix(unlist(lapply(Stroke_finaldt_EUR_out,function(x) x$beta)),ncol = length(Stroke_finaldt_EUR_out))
seGYbetaList[[4]] <- matrix(unlist(lapply(Stroke_finaldt_EUR_out,function(x) x$standard_error)),ncol = length(Stroke_finaldt_EUR_out))

names(GXbetaList) <- names(seGXbetaList) <- 
  names(GYbetaList) <- names(seGYbetaList) <- c('SA','AFR','EA','EUR')

colnames(GXbetaList[[1]]) <- colnames(GXbetaList[[2]]) <- 
  colnames(GXbetaList[[3]]) <- colnames(GXbetaList[[4]]) <- dt_SA$reportedTrait

colnames(seGXbetaList[[1]]) <- colnames(seGXbetaList[[2]]) <- 
  colnames(seGXbetaList[[3]]) <- colnames(seGXbetaList[[4]]) <- dt_SA$reportedTrait

colnames(GYbetaList[[1]]) <- colnames(GYbetaList[[2]]) <- 
  colnames(GYbetaList[[3]]) <- colnames(GYbetaList[[4]]) <- names(Stroke_finaldt_SA_out)

colnames(seGYbetaList[[1]]) <- colnames(seGYbetaList[[2]]) <- 
  colnames(seGYbetaList[[3]]) <- colnames(seGYbetaList[[4]]) <- names(Stroke_finaldt_SA_out)

rownames(GXbetaList[[1]]) <- rownames(seGXbetaList[[1]]) <- SNP_finaldt_SA_out[[1]]$rsid
rownames(GXbetaList[[2]]) <- rownames(seGXbetaList[[2]]) <- SNP_finaldt_AFR_out[[1]]$rsid
rownames(GXbetaList[[3]]) <- rownames(seGXbetaList[[3]]) <- SNP_finaldt_EA_out[[1]]$rsid
rownames(GXbetaList[[4]]) <- rownames(seGXbetaList[[4]]) <- SNP_finaldt_EUR_out[[1]]$rsid

save(GXbetaList,seGXbetaList,GYbetaList,seGYbetaList,file='code41top.summMatrix.Rdata')

#################### LD matrix #########################################
library(data.table)
load('code41top.summMatrix.Rdata')

for(i in 1:4){
  
  cat('i=',i,'\n')
  
  snp_once <- rownames(GXbetaList[[i]])
  
  cat('length(snp_once)=',  length(snp_once),'\n')
  
  if(i==1){
    Ances <- 'SAS'
  }else if(i==2){
    Ances <- 'AFR'
  }else if(i==3){
    Ances <- 'EAS'
  }else if(i==4){
    Ances <- 'EUR'
  }
  bim <- fread(paste0("H:/LD_REFERENCE_PANEL/1kg.v3/",Ances,".bim"), header = FALSE)
  ref_snps <- bim$V2
  missing_snps <- setdiff(snp_once, ref_snps)
  snp_once1 <- snp_once[!(snp_once %in% missing_snps)]
  write.table(snp_once1,file='LDmatrix/snp_once.txt',quote =F,row.names =F,col.names =F)
  
  COMMD <- paste0('H:/PLINK/plink_win64_20241022/plink ',
                  '--bfile H:/LD_REFERENCE_PANEL/1kg.v3/',Ances,
                  ' --extract LDmatrix/snp_once.txt',
                  ' --r square --out LDmatrix/',Ances,'_ld_matrix_top')
  system(COMMD)
  ld_matrix <- fread(paste0('LDmatrix/',Ances,'_ld_matrix_top.ld'))
  ld_matrix <- as.matrix(ld_matrix)
  rownames(ld_matrix) <- colnames(ld_matrix) <- snp_once1
  ld_matrix1 <- matrix(0,nrow=length(snp_once),ncol=length(snp_once))
  rownames(ld_matrix1) <- colnames(ld_matrix1) <- snp_once
  for(o1 in 1:length(snp_once)){
    for(o2 in 1:length(snp_once)){
      if((snp_once[o1] %in% snp_once1) & (snp_once[o2] %in% snp_once1)){
        ld_matrix1[snp_once[o1],snp_once[o2]] <- ld_matrix[snp_once[o1],snp_once[o2]]
      }
      
    }
  }
  
  cat('length(snp_once1)=',  length(snp_once1),'\n')
  cat('nrow(ld_matrix)=',nrow(ld_matrix),'\n')
  cat('nrow(ld_matrix1)=',nrow(ld_matrix1),'\n')
  
  save(ld_matrix1,file=paste0('LDmatrix/',Ances,'_ld_matrix_top.Rdata'))
  
}

#######
library(MASS)

load('LDmatrix/EUR_ld_matrix_top.Rdata')
EUR_ld_matrix <- ld_matrix1
load('LDmatrix/SAS_ld_matrix_top.Rdata')
SAS_ld_matrix <- ld_matrix1
load('LDmatrix/AFR_ld_matrix_top.Rdata')
AFR_ld_matrix <- ld_matrix1
load('LDmatrix/EAS_ld_matrix_top.Rdata')
EAS_ld_matrix <- ld_matrix1

SAS_ld_matrix <- as.matrix(SAS_ld_matrix)
AFR_ld_matrix <- as.matrix(AFR_ld_matrix)
EUR_ld_matrix <- as.matrix(EUR_ld_matrix)
EAS_ld_matrix <- as.matrix(EAS_ld_matrix)

SAS_ld_matrix <- SAS_ld_matrix^2
AFR_ld_matrix <- AFR_ld_matrix^2
EUR_ld_matrix <- EUR_ld_matrix^2
EAS_ld_matrix <- EAS_ld_matrix^2

SAS_ld_matrix_inv <- ginv(SAS_ld_matrix)
AFR_ld_matrix_inv <- ginv(AFR_ld_matrix)
EUR_ld_matrix_inv <- ginv(EUR_ld_matrix)
EAS_ld_matrix_inv <- ginv(EAS_ld_matrix)

save(SAS_ld_matrix,AFR_ld_matrix,EUR_ld_matrix,EAS_ld_matrix,
     SAS_ld_matrix_inv,AFR_ld_matrix_inv,EUR_ld_matrix_inv,EAS_ld_matrix_inv,
     file='LDmatrix/ld_matrix_all_top.Rdata')

rm(list=ls())
gc()

################### marginal association to conditional association ########
library(MASS)
load('code41top.summMatrix.Rdata')
load('LDmatrix/ld_matrix_all_top.Rdata')
source('TransCondi.R')

Ancesname <- names(GXbetaList)
outname <- colnames(GYbetaList[[1]])

fdata_NK.Listall <- list()
fdata_NK1.Listall <- list()

for(i in 1:length(Ancesname)){
  
  if(Ancesname[i] == 'SA'){
    LD_mat <- SAS_ld_matrix
    LD_inv <- SAS_ld_matrix_inv
  }else if(Ancesname[i] == 'AFR'){
    LD_mat <- AFR_ld_matrix
    LD_inv <- AFR_ld_matrix_inv
  }else if(Ancesname[i] == 'EUR'){
    LD_mat <- EUR_ld_matrix
    LD_inv <- EUR_ld_matrix_inv
  }else if(Ancesname[i] == 'EA'){
    LD_mat <- EAS_ld_matrix
    LD_inv <- EAS_ld_matrix_inv
  }
  
  fdata_NK.List <- list()
  fdata_NK1.List <- list()
  
  for(j in 1:length(outname)){
    
    cat('i=',i,' j=',j,'\n')
    
    once <- list(Bx_obs=GXbetaList[[i]],
                 By_obs=GYbetaList[[i]][,j],
                 Sigy_condi=seGYbetaList[[i]][,j],
                 LD_mat=LD_mat,
                 LD_inv=LD_inv)
    
    fdata_NK.List[[j]] <- once
    fdata_NK1.List[[j]] <- TransCondi(once)
    
  }
  
  names(fdata_NK.List) <- names(fdata_NK1.List) <- outname
  
  fdata_NK.Listall[[i]] <- fdata_NK.List
  fdata_NK1.Listall[[i]] <- fdata_NK1.List
  
}

names(fdata_NK.Listall) <- names(fdata_NK1.Listall) <- Ancesname

save(fdata_NK.Listall,fdata_NK1.Listall,file='code41.fdata_NK1all_top.Rdata')



