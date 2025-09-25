
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

# #### 每个种族分开找SNP，每个性状选top 10 snp ######
################# South Asian ##########################
SA.loc <- dt1$ancestry=='South Asian'
dt2 <- dt1[SA.loc,]
summary(dt2$snp_num)

## 899个性状找top10个的SNP
# SNP_top10_SA <- list()
# 
# for(i in 1:nrow(dt2)){
#   cat('i=',i,'\n')
#   top <- 10
#   
#   if(dt2$snp_num[i]>10){top <- dt2$snp_num[i]}
# 
#   deskfile <- paste0('H:/metabolite899/',dt2$accessionId[i],'.h.tsv.gz')
#   dtonce <- fread(deskfile)
# 
#   dtonce$p_value <- as.numeric(dtonce$p_value)
#   dtonce <- dtonce[order(dtonce$p_value),]
#   dtonce <- dtonce[1:top,]
# 
#   SNP_top10_SA[[i]] <- dtonce
#   rm(dtonce)
# }
# save(SNP_top10_SA,file='code21.SNP_top10_SA.Rdata')

## 看看这top10的SNP哪几个在所有的数据库中都存在
load('code21.SNP_top10_SA.Rdata')
SNP_top10_uni <- lapply(SNP_top10_SA,function(x) x$rsid)
SNP_top10_uni <- unique(unlist(SNP_top10_uni)) #9780

# SNP_top10dt_SA <- list()
# 
# for(i in 1:nrow(dt2)){
#   cat('i=',i,'\n')
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
# save(SNP_top10dt_SA,file='code21.SNP_top10dt_SA.Rdata')

## 加一步，看看这top10的SNP哪几个在结局的数据库中都存在
# load('code21.SNP_top10dt_SA.Rdata')
# outdt1 <- read.csv('dt_out_stroke.csv')
# unique(outdt1$ancestry)
# outdt1 <- outdt1[outdt1$ancestry=='South Asian',]
# SNP_top10dt_SA_out <- SNP_top10dt_SA
# for(i in 1:nrow(outdt1)){
#   cat('i=',i,'\n')
#   
#   deskfile <- paste0('H:/metabolite899/',outdt1$accessionId[i],'.h.tsv.gz')
#   dtout_once <- fread(deskfile)
#   
#   for(j in 1:length(SNP_top10dt_SA_out)){
#     once.exp <- SNP_top10dt_SA_out[[j]]
#     once.exp <- once.exp[once.exp$rsid %in% dtout_once$rsid,]
#     SNP_top10dt_SA_out[[j]] <- once.exp
#   }
#   
#   rm(dtout_once,once.exp)
# }
# 
# save(SNP_top10dt_SA_out,file='code21.SNP_top10dt_SA_stroke.Rdata')


## 统计频率，看看哪些SNP在所有数据库中都有
load('code21.SNP_top10dt_SA_stroke.Rdata')

snpal1 <- lapply(SNP_top10dt_SA_out, function(x) unique(x$rsid))
snpal1 <- unlist(snpal1)
snpal1 <- table(snpal1)

choSNP <- names(snpal1)[snpal1==899] #2564个

choSNP_loc <- lapply(SNP_top10dt_SA_out,function(x) x[which(x$rsid %in% choSNP),])
choSNP_loc <- lapply(choSNP_loc,function(x) x[order(x$rsid),])

all(unlist(lapply(choSNP_loc,function(x) x$rsid==choSNP)))

dt2$keep_trait <- 0
for(i in 1:899){
  once <- sum(choSNP %in% SNP_top10_SA[[i]]$rsid)
  if(once>0){dt2$keep_trait[i] <- 1}
}
sum(dt2$keep_trait) #664

loc.trait <- which(dt2$keep_trait==1)
SNP_finaldt_SA_out <- choSNP_loc[loc.trait]
dt2 <- dt2[dt2$keep_trait==1,]
names(SNP_finaldt_SA_out) <- dt2$reportedTrait

outdt1 <- read.csv('dt_out_stroke.csv')
unique(outdt1$ancestry)
outdt1 <- outdt1[outdt1$ancestry=='South Asian',]
Stroke_finaldt_SA_out <- list()
for(i in 1:nrow(outdt1)){
  cat('i=',i,'\n')

  deskfile <- paste0('H:/metabolite899/',outdt1$accessionId[i],'.h.tsv.gz')
  dtout_once <- fread(deskfile)
  dtout_once <- dtout_once[dtout_once$rsid %in% choSNP,]
  dtout_once <- dtout_once[order(dtout_once$rsid),]
  
  Stroke_finaldt_SA_out[[i]] <- dtout_once
}
names(Stroke_finaldt_SA_out) <- outdt1$reportedTrait

all(unlist(lapply(Stroke_finaldt_SA_out,function(x) x$rsid==choSNP)))

save(dt2,SNP_finaldt_SA_out,Stroke_finaldt_SA_out,file='code21.SNP_dtfinal_SA_stroke.Rdata')

###############################################################
################# African unspecified ##########################


SA.loc <- dt1$ancestry=='African unspecified'
dt2 <- dt1[SA.loc,]
summary(dt2$snp_num)

## 899个性状找top10个的SNP
SNP_top10_AFR <- list()

for(i in 1:nrow(dt2)){
  cat('i=',i,'\n')
  top <- 10

  if(dt2$snp_num[i]>10){top <- dt2$snp_num[i]}

  deskfile <- paste0('H:/metabolite899/',dt2$accessionId[i],'.h.tsv.gz')
  dtonce <- fread(deskfile)

  dtonce$p_value <- as.numeric(dtonce$p_value)
  dtonce <- dtonce[order(dtonce$p_value),]
  dtonce <- dtonce[1:top,]

  SNP_top10_AFR[[i]] <- dtonce
  rm(dtonce)
}
save(SNP_top10_AFR,file='code21.SNP_top10_AFR.Rdata')

## 看看这top10的SNP哪几个在所有的数据库中都存在
#load('code21.SNP_top10_AFR.Rdata')
SNP_top10_uni <- lapply(SNP_top10_AFR,function(x) x$rsid)
SNP_top10_uni <- unique(unlist(SNP_top10_uni)) #15173

SNP_top10dt_AFR <- list()

for(i in 1:nrow(dt2)){
  cat('i=',i,'\n')

  deskfile <- paste0('H:/metabolite899/',dt2$accessionId[i],'.h.tsv.gz')
  dtonce <- fread(deskfile)

  dtonce1 <- dtonce[dtonce$rsid %in% SNP_top10_uni,]
  dtonce1 <- dtonce1[order(dtonce1$rsid),]

  SNP_top10dt_AFR[[i]] <- dtonce1
  rm(dtonce,dtonce1)
}
save(SNP_top10dt_AFR,file='code21.SNP_top10dt_AFR.Rdata')

## 加一步，看看这top10的SNP哪几个在结局的数据库中都存在
#load('code21.SNP_top10dt_AFR.Rdata')
outdt1 <- read.csv('dt_out_stroke.csv')
unique(outdt1$ancestry)
outdt1 <- outdt1[outdt1$ancestry=='African American or Afro-Caribbean',]
SNP_top10dt_AFR_out <- SNP_top10dt_AFR
for(i in 1:nrow(outdt1)){
  cat('i=',i,'\n')

  deskfile <- paste0('H:/metabolite899/',outdt1$accessionId[i],'.h.tsv.gz')
  dtout_once <- fread(deskfile)

  for(j in 1:length(SNP_top10dt_AFR_out)){
    once.exp <- SNP_top10dt_AFR_out[[j]]
    once.exp <- once.exp[once.exp$rsid %in% dtout_once$rsid,]
    SNP_top10dt_AFR_out[[j]] <- once.exp
  }

  rm(dtout_once,once.exp)
}

save(SNP_top10dt_AFR_out,file='code21.SNP_top10dt_AFR_stroke.Rdata')


## 统计频率，看看哪些SNP在所有数据库中都有
#load('code21.SNP_top10dt_AFR_stroke.Rdata')

snpal1 <- lapply(SNP_top10dt_AFR_out, function(x) unique(x$rsid))
snpal1 <- unlist(snpal1)
snpal1 <- table(snpal1)

choSNP <- names(snpal1)[snpal1==899] #6647个

choSNP_loc <- lapply(SNP_top10dt_AFR_out,function(x) x[which(x$rsid %in% choSNP),])
choSNP_loc <- lapply(choSNP_loc,function(x) x[order(x$rsid),])

all(unlist(lapply(choSNP_loc,function(x) x$rsid==choSNP)))

dt2$keep_trait <- 0
for(i in 1:899){
  once <- sum(choSNP %in% SNP_top10_AFR[[i]]$rsid)
  if(once>0){dt2$keep_trait[i] <- 1}
}
sum(dt2$keep_trait) #667

loc.trait <- which(dt2$keep_trait==1)
SNP_finaldt_AFR_out <- choSNP_loc[loc.trait]
dt2 <- dt2[dt2$keep_trait==1,]
names(SNP_finaldt_AFR_out) <- dt2$reportedTrait

outdt1 <- read.csv('dt_out_stroke.csv')
unique(outdt1$ancestry)
outdt1 <- outdt1[outdt1$ancestry=='African American or Afro-Caribbean',]
Stroke_finaldt_AFR_out <- list()
for(i in 1:nrow(outdt1)){
  cat('i=',i,'\n')
  
  deskfile <- paste0('H:/metabolite899/',outdt1$accessionId[i],'.h.tsv.gz')
  dtout_once <- fread(deskfile)
  dtout_once <- dtout_once[dtout_once$rsid %in% choSNP,]
  dtout_once <- dtout_once[order(dtout_once$rsid),]
  
  Stroke_finaldt_AFR_out[[i]] <- dtout_once
}
names(Stroke_finaldt_AFR_out) <- outdt1$reportedTrait

all(unlist(lapply(Stroke_finaldt_AFR_out,function(x) x$rsid==choSNP)))

save(dt2,SNP_finaldt_AFR_out,Stroke_finaldt_AFR_out,file='code21.SNP_dtfinal_AFR_stroke.Rdata')


###############################################################
################# East Asian ##########################

SA.loc <- dt1$ancestry=='East Asian'
dt2 <- dt1[SA.loc,]
summary(dt2$snp_num)

## 899个性状找top10个的SNP
SNP_top10_EA <- list()

for(i in 1:nrow(dt2)){
  cat('i=',i,'\n')
  top <- 10
  
  if(dt2$snp_num[i]>10){top <- dt2$snp_num[i]}
  
  deskfile <- paste0('H:/metabolite899/',dt2$accessionId[i],'.h.tsv.gz')
  dtonce <- fread(deskfile)
  
  dtonce$p_value <- as.numeric(dtonce$p_value)
  dtonce <- dtonce[order(dtonce$p_value),]
  dtonce <- dtonce[1:top,]
  
  SNP_top10_EA[[i]] <- dtonce
  rm(dtonce)
}
save(SNP_top10_EA,file='code21.SNP_top10_EA.Rdata')

## 看看这top10的SNP哪几个在所有的数据库中都存在
#load('code21.SNP_top10_EA.Rdata')
SNP_top10_uni <- lapply(SNP_top10_EA,function(x) x$rsid)
SNP_top10_uni <- unique(unlist(SNP_top10_uni)) #9780

SNP_top10dt_EA <- list()

for(i in 1:nrow(dt2)){
  cat('i=',i,'\n')
  
  deskfile <- paste0('H:/metabolite899/',dt2$accessionId[i],'.h.tsv.gz')
  dtonce <- fread(deskfile)
  
  dtonce1 <- dtonce[dtonce$rsid %in% SNP_top10_uni,]
  dtonce1 <- dtonce1[order(dtonce1$rsid),]
  
  SNP_top10dt_EA[[i]] <- dtonce1
  rm(dtonce,dtonce1)
}
save(SNP_top10dt_EA,file='code21.SNP_top10dt_EA.Rdata')

## 加一步，看看这top10的SNP哪几个在结局的数据库中都存在
#load('code21.SNP_top10dt_EA.Rdata')
outdt1 <- read.csv('dt_out_stroke.csv')
unique(outdt1$ancestry)
outdt1 <- outdt1[outdt1$ancestry=='East Asian',]
SNP_top10dt_EA_out <- SNP_top10dt_EA
for(i in 1:nrow(outdt1)){
  cat('i=',i,'\n')
  
  deskfile <- paste0('H:/metabolite899/',outdt1$accessionId[i],'.h.tsv.gz')
  dtout_once <- fread(deskfile)
  
  for(j in 1:length(SNP_top10dt_EA_out)){
    once.exp <- SNP_top10dt_EA_out[[j]]
    once.exp <- once.exp[once.exp$rsid %in% dtout_once$rsid,]
    SNP_top10dt_EA_out[[j]] <- once.exp
  }
  
  rm(dtout_once,once.exp)
}

save(SNP_top10dt_EA_out,file='code21.SNP_top10dt_EA_stroke.Rdata')


## 统计频率，看看哪些SNP在所有数据库中都有
#load('code21.SNP_top10dt_EA_stroke.Rdata')

snpal1 <- lapply(SNP_top10dt_EA_out, function(x) unique(x$rsid))
snpal1 <- unlist(snpal1)
snpal1 <- table(snpal1)

choSNP <- names(snpal1)[snpal1==899] #2564个

choSNP_loc <- lapply(SNP_top10dt_EA_out,function(x) x[which(x$rsid %in% choSNP),])
choSNP_loc <- lapply(choSNP_loc,function(x) x[order(x$rsid),])

all(unlist(lapply(choSNP_loc,function(x) x$rsid==choSNP)))

dt2$keep_trait <- 0
for(i in 1:899){
  once <- sum(choSNP %in% SNP_top10_EA[[i]]$rsid)
  if(once>0){dt2$keep_trait[i] <- 1}
}
sum(dt2$keep_trait) #664

loc.trait <- which(dt2$keep_trait==1)
SNP_finaldt_EA_out <- choSNP_loc[loc.trait]
dt2 <- dt2[dt2$keep_trait==1,]
names(SNP_finaldt_EA_out) <- dt2$reportedTrait

outdt1 <- read.csv('dt_out_stroke.csv')
unique(outdt1$ancestry)
outdt1 <- outdt1[outdt1$ancestry=='East Asian',]
Stroke_finaldt_EA_out <- list()
for(i in 1:nrow(outdt1)){
  cat('i=',i,'\n')
  
  deskfile <- paste0('H:/metabolite899/',outdt1$accessionId[i],'.h.tsv.gz')
  dtout_once <- fread(deskfile)
  dtout_once <- dtout_once[dtout_once$rsid %in% choSNP,]
  dtout_once <- dtout_once[order(dtout_once$rsid),]
  
  Stroke_finaldt_EA_out[[i]] <- dtout_once
}
names(Stroke_finaldt_EA_out) <- outdt1$reportedTrait

all(unlist(lapply(Stroke_finaldt_EA_out,function(x) x$rsid==choSNP)))

save(dt2,SNP_finaldt_EA_out,Stroke_finaldt_EA_out,file='code21.SNP_dtfinal_EA_stroke.Rdata')

###############################################################
################# European ##########################

SA.loc <- dt1$ancestry=='European'
dt2 <- dt1[SA.loc,]
summary(dt2$snp_num)

## 899个性状找top10个的SNP
SNP_top10_EUR <- list()

for(i in 1:nrow(dt2)){
  cat('i=',i,'\n')
  top <- 10
  
  deskfile <- paste0('H:/metabolite899/',dt2$accessionId[i],'.h.tsv.gz')
  dtonce <- fread(deskfile)
  
  dtonce$p_value <- as.numeric(dtonce$p_value)
  
  num_cut <- sum(dtonce$p_value<5e-10)
  if(num_cut>10){top <- num_cut}
  
  dtonce <- dtonce[order(dtonce$p_value),]
  dtonce <- dtonce[1:top,]
  
  SNP_top10_EUR[[i]] <- dtonce
  rm(dtonce)
}
save(SNP_top10_EUR,file='code21.SNP_top10_EUR.Rdata')

## 看看这top10的SNP哪几个在所有的数据库中都存在
#load('code21.SNP_top10_EUR.Rdata')
SNP_top10_EUR <- lapply(SNP_top10_EUR,function(x) if(nrow(x)>50){x <- x[1:50,]}else{x})
SNP_top10_uni <- lapply(SNP_top10_EUR,function(x) x$rsid)
SNP_top10_uni <- unique(unlist(SNP_top10_uni)) #12554

SNP_top10dt_EUR <- list()

for(i in 1:nrow(dt2)){
  cat('i=',i,'\n')
  
  deskfile <- paste0('H:/metabolite899/',dt2$accessionId[i],'.h.tsv.gz')
  dtonce <- fread(deskfile)
  
  dtonce1 <- dtonce[dtonce$rsid %in% SNP_top10_uni,]
  dtonce1 <- dtonce1[order(dtonce1$rsid),]
  
  SNP_top10dt_EUR[[i]] <- dtonce1
  rm(dtonce,dtonce1)
}
save(SNP_top10dt_EUR,file='code21.SNP_top10dt_EUR.Rdata')

## 加一步，看看这top10的SNP哪几个在结局的数据库中都存在
#load('code21.SNP_top10dt_EA.Rdata')
outdt1 <- read.csv('dt_out_stroke.csv')
unique(outdt1$ancestry)
outdt1 <- outdt1[outdt1$ancestry=='European',]
SNP_top10dt_EUR_out <- SNP_top10dt_EUR
for(i in 1:nrow(outdt1)){
  cat('i=',i,'\n')
  
  deskfile <- paste0('H:/metabolite899/',outdt1$accessionId[i],'.h.tsv.gz')
  dtout_once <- fread(deskfile)
  
  for(j in 1:length(SNP_top10dt_EUR_out)){
    once.exp <- SNP_top10dt_EUR_out[[j]]
    once.exp <- once.exp[once.exp$rsid %in% dtout_once$rsid,]
    SNP_top10dt_EUR_out[[j]] <- once.exp
  }
  
  rm(dtout_once,once.exp)
}

save(SNP_top10dt_EUR_out,file='code21.SNP_top10dt_EUR_stroke.Rdata')


## 统计频率，看看哪些SNP在所有数据库中都有
#load('code21.SNP_top10dt_EUR_stroke.Rdata')

snpal1 <- lapply(SNP_top10dt_EUR_out, function(x) unique(x$rsid))
snpal1 <- unlist(snpal1)
snpal1 <- table(snpal1)

choSNP <- names(snpal1)[snpal1==899] #8466个

choSNP_loc <- lapply(SNP_top10dt_EUR_out,function(x) x[which(x$rsid %in% choSNP),])
choSNP_loc <- lapply(choSNP_loc,function(x) x[order(x$rsid),])

all(unlist(lapply(choSNP_loc,function(x) x$rsid==choSNP)))

dt2$keep_trait <- 0
for(i in 1:899){
  once <- sum(choSNP %in% SNP_top10_EUR[[i]]$rsid)
  if(once>0){dt2$keep_trait[i] <- 1}
}
sum(dt2$keep_trait) #794

loc.trait <- which(dt2$keep_trait==1)
SNP_finaldt_EUR_out <- choSNP_loc[loc.trait]
dt2 <- dt2[dt2$keep_trait==1,]
names(SNP_finaldt_EUR_out) <- dt2$reportedTrait

outdt1 <- read.csv('dt_out_stroke.csv')
unique(outdt1$ancestry)
outdt1 <- outdt1[outdt1$ancestry=='European',]
Stroke_finaldt_EUR_out <- list()
for(i in 1:nrow(outdt1)){
  cat('i=',i,'\n')
  
  deskfile <- paste0('H:/metabolite899/',outdt1$accessionId[i],'.h.tsv.gz')
  dtout_once <- fread(deskfile)
  dtout_once <- dtout_once[dtout_once$rsid %in% choSNP,]
  dtout_once <- dtout_once[order(dtout_once$rsid),]
  
  Stroke_finaldt_EUR_out[[i]] <- dtout_once
}
names(Stroke_finaldt_EUR_out) <- outdt1$reportedTrait

all(unlist(lapply(Stroke_finaldt_EUR_out,function(x) x$rsid==choSNP)))

save(dt2,SNP_finaldt_EUR_out,Stroke_finaldt_EUR_out,file='code21.SNP_dtfinal_EUR_stroke.Rdata')


################################## combine ####################################

load('code21.SNP_dtfinal_SA_stroke.Rdata')
dt_SA <- dt2
load('code21.SNP_dtfinal_AFR_stroke.Rdata')
dt_AFR <- dt2
load('code21.SNP_dtfinal_EA_stroke.Rdata')
dt_EA <- dt2
load('code21.SNP_dtfinal_EUR_stroke.Rdata')
dt_EUR <- dt2

uni_trait <- dt_SA$reportedTrait
uni_trait <- intersect(uni_trait,dt_AFR$reportedTrait)
uni_trait <- intersect(uni_trait,dt_EA$reportedTrait)
uni_trait <- intersect(uni_trait,dt_EUR$reportedTrait) #552个性状

save(uni_trait,file='code21.uni_trait.Rdata')


