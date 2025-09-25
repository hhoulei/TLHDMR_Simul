
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
table(dt1$ancestry) #899 metabolites

# SNP_sig <- list()
# 
# for(i in 1:nrow(dt1)){
#   cat('i=',i,'\n')
# 
#   deskfile <- paste0('H:/metabolite899/',dt1$accessionId[i],'.h.tsv.gz')
#   dtonce <- fread(deskfile)
#   dtonce$p_value <- as.numeric(dtonce$p_value)
#   dtonce1 <- dtonce[dtonce$p_value<5e-8,]
#   head(dtonce1)
#   if(nrow(dtonce1)==0){
#     dtonce1 <- dtonce[which.min(dtonce$p_value),]
#   }
#   SNP_sig[[i]] <- dtonce1
#   rm(dtonce1,dtonce)
# }
# save(SNP_sig,file='code2.SNP_sig.Rdata')

load('code2.SNP_sig.Rdata')

snp_num <- lapply(SNP_sig, function(x) nrow(x))
snp_num <- unlist(snp_num)
dt1$snp_num <- snp_num

# #### 每个种族分开找SNP，每个性状只选top snp ######
################# South Asian ##########################
SA.loc <- dt1$ancestry=='South Asian'
dt2 <- dt1[SA.loc,]
summary(dt2$snp_num)

SNP_TOP_SA <- lapply(SNP_sig[SA.loc], function(x) x$rsid[which.min(x$p_value)])
SNP_TOP_SA <- unlist(SNP_TOP_SA)
dt2$SNP_TOP_SA <- SNP_TOP_SA
SNP_TOP_SA <- unique(SNP_TOP_SA) #853个SNP

# SNP_topdt_SA <- list()
# 
# for(i in 1:nrow(dt2)){
#   cat('i=',i,'\n')
# 
#   deskfile <- paste0('H:/metabolite899/',dt2$accessionId[i],'.h.tsv.gz')
#   dtonce <- fread(deskfile)
# 
#   dtonce1 <- dtonce[dtonce$rsid %in% SNP_TOP_SA,]
#   dtonce1 <- dtonce1[order(dtonce1$rsid),]
#   dtonce1 <- dtonce1[!duplicated(dtonce1$rsid),]
# 
#   SNP_topdt_SA[[i]] <- dtonce1
#   rm(dtonce1,dtonce)
# }
# save(SNP_topdt_SA,file='code2.SNP_topdt_SA.Rdata')

##
load('code2.SNP_topdt_SA.Rdata')
snp_uni <- lapply(SNP_topdt_SA, function(x) unique(x$rsid))
snp_uni <- unlist(snp_uni)
snp_uni <- table(snp_uni)
snp_uni <- names(snp_uni)[snp_uni==899] ##614
### 只有614个SNP在所有的数据库中存在

dt2$SNP_TOP_allexist <- ifelse(dt2$SNP_TOP %in% snp_uni,1,0)
sum(dt2$SNP_TOP_allexist) #650

## 剩余 899-650=249个性状找top10个的SNP
# SNP_top10_SA <- list()
# 
# for(i in 1:nrow(dt2)){
#   cat('i=',i,'\n')
# 
#   if(dt2$SNP_TOP_allexist[i]==1) next
# 
#   deskfile <- paste0('H:/metabolite899/',dt2$accessionId[i],'.h.tsv.gz')
#   dtonce <- fread(deskfile)
# 
#   dtonce$p_value <- as.numeric(dtonce$p_value)
#   dtonce <- dtonce[order(dtonce$p_value),]
#   dtonce <- dtonce[1:10,]
# 
#   SNP_top10_SA[[i]] <- dtonce
#   rm(dtonce)
# }
# save(SNP_top10_SA,file='code2.SNP_top10_SA.Rdata')

## 看看这top10的SNP哪几个在所有的数据库中都存在
load('code2.SNP_top10_SA.Rdata')
SNP_top10_uni <- lapply(SNP_top10_SA,function(x) x$rsid)
SNP_top10_uni <- unique(unlist(SNP_top10_uni)) #2418

# SNP_top10dt_SA <- list()
# 
# for(i in 1:nrow(dt2)){
#   cat('i=',i,'\n')
# 
#   #if(dt2$SNP_TOP_allexist[i]==1) next
# 
#   deskfile <- paste0('H:/metabolite899/',dt2$accessionId[i],'.h.tsv.gz')
#   dtonce <- fread(deskfile)
# 
#   dtonce1 <- dtonce[dtonce$rsid %in% SNP_top10_uni,]
#   dtonce1 <- dtonce1[order(dtonce1$rsid),]
# 
#   SNP_top10dt_SA[[i]] <- dtonce1
#   rm(dtonce,dtonce1)
# }
# save(SNP_top10dt_SA,file='code2.SNP_top10dt_SA.Rdata')

## 加一步，看看这top10的SNP哪几个在结局的数据库中都存在
load('code2.SNP_top10dt_SA.Rdata')
outdt1 <- read.csv('dt_out_stroke.csv')
unique(outdt1$ancestry)
outdt1 <- outdt1[outdt1$ancestry=='South Asian',]

SNP_top10dt_SA_out <- SNP_top10dt_SA
for(i in 1:nrow(outdt1)){
  cat('i=',i,'\n')

  #if(dt2$SNP_TOP_allexist[i]==1) next

  deskfile <- paste0('H:/metabolite899/',outdt1$accessionId[i],'.h.tsv.gz')
  dtout_once <- fread(deskfile)
  
  for(j in 1:length(SNP_top10dt_SA_out)){
    once.exp <- SNP_top10dt_SA_out[[j]]
    once.exp <- once.exp[once.exp$rsid %in% dtout_once$rsid,]
    SNP_top10dt_SA_out[[j]] <- once.exp
  }
  
  rm(dtout_once,once.exp)
}

outdt2 <- read.csv('dt_out_lung.csv')
unique(outdt2$ancestry)
outdt2 <- outdt2[outdt2$ancestry=='South Asian',]

for(i in 1:nrow(outdt2)){
  cat('i=',i,'\n')
  
  #if(dt2$SNP_TOP_allexist[i]==1) next
  
  deskfile <- paste0('H:/metabolite899/',outdt2$accessionId[i],'.tsv.gz')
  dtout_once <- fread(deskfile)
  
  for(j in 1:length(SNP_top10dt_SA_out)){
    once.exp <- SNP_top10dt_SA_out[[j]]
    once.exp <- once.exp[once.exp$rsid %in% dtout_once$rsid,]
    SNP_top10dt_SA_out[[j]] <- once.exp
  }
  
  rm(dtout_once,once.exp)
}

save(SNP_top10dt_SA_out,file='code2.SNP_top10dt_SA_out.Rdata')


## 统计频率，看看哪些SNP在所有数据库中都有
load('code2.SNP_top10dt_SA_out.Rdata')
SNP_top10dt_SA <- SNP_top10dt_SA_out

snpal1 <- lapply(SNP_top10dt_SA, function(x) unique(x$rsid))
snpal1 <- unlist(snpal1)
snpal1 <- table(snpal1)

sum(dt2$SNP_TOP_allexist==0)
choSNP <- names(snpal1)[snpal1==899] #1459个

choSNP_loc <- lapply(SNP_top10_SA,function(x) which(x$rsid %in% choSNP))
choSNP_loc <- lapply(choSNP_loc, function(x) if(length(x)>0){x[1]}else{NA})

for(i in 1:length(choSNP_loc)){
  if(is.na(choSNP_loc[[i]])) next

  loc.once <-  choSNP_loc[[i]]
  choSNP_loc[[i]] <- SNP_top10_SA[[i]]$rsid[loc.once]
}

dt2$SNP_TOP10 <- unlist(choSNP_loc)
dt2$SNP_TOP10_allexist <- ifelse((dt2$SNP_TOP_allexist==1) | (!is.na(dt2$SNP_TOP10)),1,0)
sum(dt2$SNP_TOP10_allexist) #887
dt2$final_SNP <- ifelse(dt2$SNP_TOP10_allexist==0,NA,ifelse(dt2$SNP_TOP_allexist==1,dt2$SNP_TOP_SA,dt2$SNP_TOP10))

trait_loc <- which(!is.na(dt2$final_SNP))
dt2 <- dt2[trait_loc,]

SNP_final_uni <- unique(dt2$final_SNP) #836

SNP_dtfinal_SA <- list()
for(i in 1:nrow(dt2)){
  cat('i=',i,'\n')

  deskfile <- paste0('H:/metabolite899/',dt2$accessionId[i],'.h.tsv.gz')
  dtonce <- fread(deskfile)

  dtonce1 <- dtonce[dtonce$rsid %in% SNP_final_uni,]
  dtonce1 <- dtonce1[order(dtonce1$rsid),]
  dtonce1 <- dtonce1[!duplicated(dtonce1$rsid),]

  SNP_dtfinal_SA[[i]] <- dtonce1

  rm(dtonce,dtonce1)
}

save(dt2,SNP_dtfinal_SA,file='code2.SNP_dtfinal_SA.Rdata')

all(unlist(lapply(SNP_dtfinal_SA,function(x) nrow(x)))==836)

###############################################################
################# African unspecified ##########################

SA.loc <- dt1$ancestry=='African unspecified'
dt3 <- dt1[SA.loc,]
summary(dt3$snp_num)

SNP_TOP_AFR <- lapply(SNP_sig[SA.loc], function(x) x$rsid[which.min(x$p_value)])
SNP_TOP_AFR <- unlist(SNP_TOP_AFR)
dt3$SNP_TOP_AFR <- SNP_TOP_AFR
SNP_TOP_AFR <- unique(SNP_TOP_AFR) #875个SNP

# SNP_topdt_AFR <- list()
# 
# for(i in 1:nrow(dt3)){
#   cat('i=',i,'\n')
# 
#   deskfile <- paste0('H:/metabolite899/',dt3$accessionId[i],'.h.tsv.gz')
#   dtonce <- fread(deskfile)
# 
#   dtonce1 <- dtonce[dtonce$rsid %in% SNP_TOP_AFR,]
#   dtonce1 <- dtonce1[order(dtonce1$rsid),]
#   dtonce1 <- dtonce1[!duplicated(dtonce1$rsid),]
# 
#   SNP_topdt_AFR[[i]] <- dtonce1
#   rm(dtonce1,dtonce)
# }
# save(SNP_topdt_AFR,file='code2.SNP_topdt_AFR.Rdata')

##
load('code2.SNP_topdt_AFR.Rdata')
snp_uni <- lapply(SNP_topdt_AFR, function(x) unique(x$rsid))
snp_uni <- unlist(snp_uni)
snp_uni <- table(snp_uni)
snp_uni <- names(snp_uni)[snp_uni==899] ##634
### 只有634个SNP在所有的数据库中存在

dt3$SNP_TOP_allexist <- ifelse(dt3$SNP_TOP %in% snp_uni,1,0)
sum(dt3$SNP_TOP_allexist) #651

## 剩余 899-651=248个性状找top10个的SNP
# SNP_top10_AFR <- list()
# 
# for(i in 1:nrow(dt3)){
#   cat('i=',i,'\n')
# 
#   if(dt3$SNP_TOP_allexist[i]==1) next
# 
#   deskfile <- paste0('H:/metabolite899/',dt3$accessionId[i],'.h.tsv.gz')
#   dtonce <- fread(deskfile)
# 
#   dtonce$p_value <- as.numeric(dtonce$p_value)
#   dtonce <- dtonce[order(dtonce$p_value),]
#   dtonce <- dtonce[1:10,]
# 
#   SNP_top10_AFR[[i]] <- dtonce
#   rm(dtonce)
# }
# save(SNP_top10_AFR,file='code2.SNP_top10_AFR.Rdata')

## 看看这top10的SNP哪几个在所有的数据库中都存在
load('code2.SNP_top10_AFR.Rdata')
SNP_top10_uni <- lapply(SNP_top10_AFR,function(x) x$rsid)
SNP_top10_uni <- unique(unlist(SNP_top10_uni)) #2430

# SNP_top10dt_AFR <- list()
# 
# for(i in 1:nrow(dt3)){
#   cat('i=',i,'\n')
# 
#   deskfile <- paste0('H:/metabolite899/',dt3$accessionId[i],'.h.tsv.gz')
#   dtonce <- fread(deskfile)
# 
#   dtonce1 <- dtonce[dtonce$rsid %in% SNP_top10_uni,]
#   dtonce1 <- dtonce1[order(dtonce1$rsid),]
# 
#   SNP_top10dt_AFR[[i]] <- dtonce1
#   rm(dtonce,dtonce1)
# }
# save(SNP_top10dt_AFR,file='code2.SNP_top10dt_AFR.Rdata')

## 统计频率，看看哪些SNP在所有数据库中都有
load('code2.SNP_top10dt_AFR.Rdata')

snpal1 <- lapply(SNP_top10dt_AFR, function(x) unique(x$rsid))
snpal1 <- unlist(snpal1)
snpal1 <- table(snpal1)

sum(dt3$SNP_TOP_allexist==0)
choSNP <- names(snpal1)[snpal1==899] #1559个

choSNP_loc <- lapply(SNP_top10_AFR,function(x) which(x$rsid %in% choSNP))
choSNP_loc <- lapply(choSNP_loc, function(x) if(length(x)>0){x[1]}else{NA})
choSNP_loc[[899]] <- NA

for(i in 1:length(choSNP_loc)){
  if(is.na(choSNP_loc[[i]])) next
  
  loc.once <-  choSNP_loc[[i]]
  choSNP_loc[[i]] <- SNP_top10_AFR[[i]]$rsid[loc.once]
}

dt3$SNP_TOP10 <- unlist(choSNP_loc)
dt3$SNP_TOP10_allexist <- ifelse((dt3$SNP_TOP_allexist==1) | (!is.na(dt3$SNP_TOP10)),1,0)
sum(dt3$SNP_TOP10_allexist) #895
dt3$final_SNP <- ifelse(dt3$SNP_TOP10_allexist==0,NA,ifelse(dt3$SNP_TOP_allexist==1,dt3$SNP_TOP_AFR,dt3$SNP_TOP10))

trait_loc <- which(!is.na(dt3$final_SNP))
dt3 <- dt3[trait_loc,]

SNP_final_uni <- unique(dt3$final_SNP) #871

SNP_dtfinal_AFR <- list()
for(i in 1:nrow(dt3)){
  cat('i=',i,'\n')
  
  deskfile <- paste0('H:/metabolite899/',dt3$accessionId[i],'.h.tsv.gz')
  dtonce <- fread(deskfile)
  
  dtonce1 <- dtonce[dtonce$rsid %in% SNP_final_uni,]
  dtonce1 <- dtonce1[order(dtonce1$rsid),]
  dtonce1 <- dtonce1[!duplicated(dtonce1$rsid),]
  
  SNP_dtfinal_AFR[[i]] <- dtonce1
  
  rm(dtonce,dtonce1)
}

save(dt3,SNP_dtfinal_AFR,file='code2.SNP_dtfinal_AFR.Rdata')


all(unlist(lapply(SNP_dtfinal_AFR,function(x) nrow(x)))==871)

###############################################################
################# East Asian ##########################

EA.loc <- dt1$ancestry=='East Asian'
dt4 <- dt1[EA.loc,]
summary(dt4$snp_num)

SNP_TOP_EA <- lapply(SNP_sig[EA.loc], function(x) x$rsid[which.min(x$p_value)])
SNP_TOP_EA <- unlist(SNP_TOP_EA)
dt4$SNP_TOP_EA <- SNP_TOP_EA
SNP_TOP_EA <- unique(SNP_TOP_EA) #858个SNP

# SNP_topdt_EA <- list()
# 
# for(i in 1:nrow(dt4)){
#   cat('i=',i,'\n')
# 
#   deskfile <- paste0('H:/metabolite899/',dt4$accessionId[i],'.h.tsv.gz')
#   dtonce <- fread(deskfile)
# 
#   dtonce1 <- dtonce[dtonce$rsid %in% SNP_TOP_EA,]
#   dtonce1 <- dtonce1[order(dtonce1$rsid),]
#   dtonce1 <- dtonce1[!duplicated(dtonce1$rsid),]
# 
#   SNP_topdt_EA[[i]] <- dtonce1
#   rm(dtonce1,dtonce)
# }
# save(SNP_topdt_EA,file='code2.SNP_topdt_EA.Rdata')

##
load('code2.SNP_topdt_EA.Rdata')
snp_uni <- lapply(SNP_topdt_EA, function(x) unique(x$rsid))
snp_uni <- unlist(snp_uni)
snp_uni <- table(snp_uni)
snp_uni <- names(snp_uni)[snp_uni==899] ##607
### 只有607个SNP在所有的数据库中存在

dt4$SNP_TOP_allexist <- ifelse(dt4$SNP_TOP %in% snp_uni,1,0)
sum(dt4$SNP_TOP_allexist) #633

## 剩余 899-633=266个性状找top10个的SNP
# SNP_top10_EA <- list()
# 
# for(i in 1:nrow(dt4)){
#   cat('i=',i,'\n')
# 
#   if(dt4$SNP_TOP_allexist[i]==1) next
# 
#   deskfile <- paste0('H:/metabolite899/',dt4$accessionId[i],'.h.tsv.gz')
#   dtonce <- fread(deskfile)
# 
#   dtonce$p_value <- as.numeric(dtonce$p_value)
#   dtonce <- dtonce[order(dtonce$p_value),]
#   dtonce <- dtonce[1:10,]
# 
#   SNP_top10_EA[[i]] <- dtonce
#   rm(dtonce)
# }
# save(SNP_top10_EA,file='code2.SNP_top10_EA.Rdata')

## 看看这top10的SNP哪几个在所有的数据库中都存在
load('code2.SNP_top10_EA.Rdata')
SNP_top10_uni <- lapply(SNP_top10_EA,function(x) x$rsid)
SNP_top10_uni <- unique(unlist(SNP_top10_uni)) #2522

# SNP_top10dt_EA <- list()
# 
# for(i in 1:nrow(dt4)){
#   cat('i=',i,'\n')
#   
#   deskfile <- paste0('H:/metabolite899/',dt4$accessionId[i],'.h.tsv.gz')
#   dtonce <- fread(deskfile)
#   
#   dtonce1 <- dtonce[dtonce$rsid %in% SNP_top10_uni,]
#   dtonce1 <- dtonce1[order(dtonce1$rsid),]
#   
#   SNP_top10dt_EA[[i]] <- dtonce1
#   rm(dtonce,dtonce1)
# }
# save(SNP_top10dt_EA,file='code2.SNP_top10dt_EA.Rdata')

## 统计频率，看看哪些SNP在所有数据库中都有
load('code2.SNP_top10dt_EA.Rdata')

snpal1 <- lapply(SNP_top10dt_EA, function(x) unique(x$rsid))
snpal1 <- unlist(snpal1)
snpal1 <- table(snpal1)

sum(dt4$SNP_TOP_allexist==0)
choSNP <- names(snpal1)[snpal1==899] #1491个

choSNP_loc <- lapply(SNP_top10_EA,function(x) which(x$rsid %in% choSNP))
choSNP_loc <- lapply(choSNP_loc, function(x) if(length(x)>0){x[1]}else{NA})

for(i in 1:length(choSNP_loc)){
  if(is.na(choSNP_loc[[i]])) next
  
  loc.once <-  choSNP_loc[[i]]
  choSNP_loc[[i]] <- SNP_top10_EA[[i]]$rsid[loc.once]
}

dt4$SNP_TOP10 <- unlist(choSNP_loc)
dt4$SNP_TOP10_allexist <- ifelse((dt4$SNP_TOP_allexist==1) | (!is.na(dt4$SNP_TOP10)),1,0)
sum(dt4$SNP_TOP10_allexist) #890
dt4$final_SNP <- ifelse(dt4$SNP_TOP10_allexist==0,NA,ifelse(dt4$SNP_TOP_allexist==1,dt4$SNP_TOP_EA,dt4$SNP_TOP10))

trait_loc <- which(!is.na(dt4$final_SNP))
dt4 <- dt4[trait_loc,]

SNP_final_uni <- unique(dt4$final_SNP) #853

SNP_dtfinal_EA <- list()
for(i in 1:nrow(dt4)){
  cat('i=',i,'\n')
  
  deskfile <- paste0('H:/metabolite899/',dt4$accessionId[i],'.h.tsv.gz')
  dtonce <- fread(deskfile)
  
  dtonce1 <- dtonce[dtonce$rsid %in% SNP_final_uni,]
  dtonce1 <- dtonce1[order(dtonce1$rsid),]
  dtonce1 <- dtonce1[!duplicated(dtonce1$rsid),]
  
  SNP_dtfinal_EA[[i]] <- dtonce1
  
  rm(dtonce,dtonce1)
}

save(dt4,SNP_dtfinal_EA,file='code2.SNP_dtfinal_EA.Rdata')

all(unlist(lapply(SNP_dtfinal_EA,function(x) nrow(x)))==853)


###############################################################
################# European ##########################

EUR.loc <- dt1$ancestry=='European'
dt5 <- dt1[EUR.loc,]
summary(dt5$snp_num)

SNP_TOP_EUR <- lapply(SNP_sig[EUR.loc], function(x) x$rsid[which.min(x$p_value)])
SNP_TOP_EUR <- unlist(SNP_TOP_EUR)
dt5$SNP_TOP_EUR <- SNP_TOP_EUR
SNP_TOP_EUR <- unique(SNP_TOP_EUR) #658个SNP

# SNP_topdt_EUR <- list()
# 
# for(i in 1:nrow(dt5)){
#   cat('i=',i,'\n')
#   
#   deskfile <- paste0('H:/metabolite899/',dt5$accessionId[i],'.h.tsv.gz')
#   dtonce <- fread(deskfile)
#   
#   dtonce1 <- dtonce[dtonce$rsid %in% SNP_TOP_EUR,]
#   dtonce1 <- dtonce1[order(dtonce1$rsid),]
#   dtonce1 <- dtonce1[!duplicated(dtonce1$rsid),]
#   
#   SNP_topdt_EUR[[i]] <- dtonce1
#   rm(dtonce1,dtonce)
# }
# save(SNP_topdt_EUR,file='code2.SNP_topdt_EUR.Rdata')

##
load('code2.SNP_topdt_EUR.Rdata')
snp_uni <- lapply(SNP_topdt_EUR, function(x) unique(x$rsid))
snp_uni <- unlist(snp_uni)
snp_uni <- table(snp_uni)
snp_uni <- names(snp_uni)[snp_uni==899] ##646
### 只有646个SNP在所有的数据库中存在

dt5$SNP_TOP_allexist <- ifelse(dt5$SNP_TOP %in% snp_uni,1,0)
sum(dt5$SNP_TOP_allexist) #887

## 剩余 899-887=12个性状找top10个的SNP
# SNP_top10_EUR <- list()
# 
# for(i in 1:nrow(dt5)){
#   cat('i=',i,'\n')
#   
#   if(dt5$SNP_TOP_allexist[i]==1) next
#   
#   deskfile <- paste0('H:/metabolite899/',dt5$accessionId[i],'.h.tsv.gz')
#   dtonce <- fread(deskfile)
#   
#   dtonce$p_value <- as.numeric(dtonce$p_value)
#   dtonce <- dtonce[order(dtonce$p_value),]
#   dtonce <- dtonce[1:10,]
#   
#   SNP_top10_EUR[[i]] <- dtonce
#   rm(dtonce)
# }
# save(SNP_top10_EUR,file='code2.SNP_top10_EUR.Rdata')

## 看看这top10的SNP哪几个在所有的数据库中都存在
load('code2.SNP_top10_EUR.Rdata')
SNP_top10_uni <- lapply(SNP_top10_EUR,function(x) x$rsid)
SNP_top10_uni <- unique(unlist(SNP_top10_uni)) #119

# SNP_top10dt_EUR <- list()
# 
# for(i in 1:nrow(dt5)){
#   cat('i=',i,'\n')
#   
#   deskfile <- paste0('H:/metabolite899/',dt5$accessionId[i],'.h.tsv.gz')
#   dtonce <- fread(deskfile)
#   
#   dtonce1 <- dtonce[dtonce$rsid %in% SNP_top10_uni,]
#   dtonce1 <- dtonce1[order(dtonce1$rsid),]
#   
#   SNP_top10dt_EUR[[i]] <- dtonce1
#   rm(dtonce,dtonce1)
# }
# save(SNP_top10dt_EUR,file='code2.SNP_top10dt_EUR.Rdata')

## 统计频率，看看哪些SNP在所有数据库中都有
load('code2.SNP_top10dt_EUR.Rdata')

snpal1 <- lapply(SNP_top10dt_EUR, function(x) unique(x$rsid))
snpal1 <- unlist(snpal1)
snpal1 <- table(snpal1)

sum(dt5$SNP_TOP_allexist==0)
choSNP <- names(snpal1)[snpal1==899] #53个

choSNP_loc <- lapply(SNP_top10_EUR,function(x) which(x$rsid %in% choSNP))
choSNP_loc <- lapply(choSNP_loc, function(x) if(length(x)>0){x[1]}else{NA})

for(i in 1:length(choSNP_loc)){
  if(is.na(choSNP_loc[[i]])) next
  
  loc.once <-  choSNP_loc[[i]]
  choSNP_loc[[i]] <- SNP_top10_EUR[[i]]$rsid[loc.once]
}
choSNP_loc[817:899] <- NA

dt5$SNP_TOP10 <- unlist(choSNP_loc)
dt5$SNP_TOP10_allexist <- ifelse((dt5$SNP_TOP_allexist==1) | (!is.na(dt5$SNP_TOP10)),1,0)
sum(dt5$SNP_TOP10_allexist) #895
dt5$final_SNP <- ifelse(dt5$SNP_TOP10_allexist==0,NA,ifelse(dt5$SNP_TOP_allexist==1,dt5$SNP_TOP_EUR,dt5$SNP_TOP10))

trait_loc <- which(!is.na(dt5$final_SNP))
dt5 <- dt5[trait_loc,]

SNP_final_uni <- unique(dt5$final_SNP) #654

SNP_dtfinal_EUR <- list()
for(i in 1:nrow(dt5)){
  cat('i=',i,'\n')
  
  deskfile <- paste0('H:/metabolite899/',dt5$accessionId[i],'.h.tsv.gz')
  dtonce <- fread(deskfile)
  
  dtonce1 <- dtonce[dtonce$rsid %in% SNP_final_uni,]
  dtonce1 <- dtonce1[order(dtonce1$rsid),]
  dtonce1 <- dtonce1[!duplicated(dtonce1$rsid),]
  
  SNP_dtfinal_EUR[[i]] <- dtonce1
  
  rm(dtonce,dtonce1)
}

save(dt5,SNP_dtfinal_EUR,file='code2.SNP_dtfinal_EUR.Rdata')

all(unlist(lapply(SNP_dtfinal_EUR,function(x) nrow(x)))==654)

################################## combine ####################################

load('code2.SNP_dtfinal_SA.Rdata')
load('code2.SNP_dtfinal_AFR.Rdata')
load('code2.SNP_dtfinal_EA.Rdata')
load('code2.SNP_dtfinal_EUR.Rdata')

uni_trait <- dt2$reportedTrait
uni_trait <- intersect(uni_trait,dt3$reportedTrait)
uni_trait <- intersect(uni_trait,dt4$reportedTrait)
uni_trait <- intersect(uni_trait,dt5$reportedTrait)

loc2 <- dt2$reportedTrait %in% uni_trait
dt2 <- dt2[loc2,]
SNP_dtfinal_SA <- SNP_dtfinal_SA[loc2]

loc3 <- dt3$reportedTrait %in% uni_trait
dt3 <- dt3[loc3,]
SNP_dtfinal_AFR <- SNP_dtfinal_AFR[loc3]

loc4 <- dt4$reportedTrait %in% uni_trait
dt4 <- dt4[loc4,]
SNP_dtfinal_EA <- SNP_dtfinal_EA[loc4]

loc5 <- dt5$reportedTrait %in% uni_trait
dt5 <- dt5[loc5,]
SNP_dtfinal_EUR <- SNP_dtfinal_EUR[loc5]

##

loc2 <- order(dt2$reportedTrait)
dt2 <- dt2[loc2,]
SNP_dtfinal_SA <- SNP_dtfinal_SA[loc2]

loc3 <- order(dt3$reportedTrait)
dt3 <- dt3[loc3,]
SNP_dtfinal_AFR <- SNP_dtfinal_AFR[loc3]

loc4 <- order(dt4$reportedTrait)
dt4 <- dt4[loc4,]
SNP_dtfinal_EA <- SNP_dtfinal_EA[loc4]

loc5 <- order(dt5$reportedTrait)
dt5 <- dt5[loc5,]
SNP_dtfinal_EUR <- SNP_dtfinal_EUR[loc5]

##

sum(dt2$reportedTrait==dt3$reportedTrait)
sum(dt3$reportedTrait==dt4$reportedTrait)
sum(dt4$reportedTrait==dt5$reportedTrait)

###

trait_name <- dt2$reportedTrait

save(dt2,dt3,dt4,dt5,trait_name,
     SNP_dtfinal_SA,SNP_dtfinal_AFR,SNP_dtfinal_EA,SNP_dtfinal_EUR,
     file='code2.final870.Rdata')

