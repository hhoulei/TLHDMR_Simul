
library(utils)
library(data.table)
dt1 <- fread('PMID36180795_studies_export.tsv')

ancestry <- strsplit(dt1$discoverySampleAncestry,' ')
N <- lapply(ancestry, function(x) x[1])
N <- as.numeric(unlist(N))

ancestry1 <- lapply(ancestry, function(x) paste0(unlist(x[-1]),collapse = ' '))
ancestry1 <- unlist(ancestry1)

dt1$ancestry <- ancestry1
dt1$N <- N

table(dt1$ancestry)

dt1 <- dt1[dt1$ancestry %in% c('African American or Afro-Caribbean',
                               'East Asian',
                               'European',
                               'South Asian'),]
table(dt1$ancestry)

write.csv(dt1,'dt_out_stroke.csv')

for(i in 5:nrow(dt1)){
  cat('i=',i,'\n')
  
  options(timeout = 60000)
  link <- paste0(dt1$summaryStatistics[i],'/harmonised/',dt1$accessionId[i],'.h.tsv.gz')
  deskfile <- paste0('H:/metabolite899/',dt1$accessionId[i],'.h.tsv.gz')
  download.file(link,deskfile)
}

##############################################################################
rm(list = ls())
gc()


library(utils)
library(data.table)
dt1 <- fread('PMID36914875_studies_export.tsv')

ancestry <- strsplit(dt1$discoverySampleAncestry,' ')
N <- lapply(ancestry, function(x) x[1])
N <- as.numeric(unlist(N))

ancestry1 <- lapply(ancestry, function(x) paste0(unlist(x[-1]),collapse = ' '))
ancestry1 <- unlist(ancestry1)

dt1$ancestry <- ancestry1
dt1$N <- N

table(dt1$ancestry)

dt1 <- dt1[dt1$ancestry %in% c('African unspecified',
                               'East Asian',
                               'European',
                               'South Asian'),]
table(dt1$ancestry)

dt1 <- dt1[order(dt1$reportedTrait),]
dt1 <- dt1[!((dt1$ancestry=='European') & (!(dt1$N %in% c(475645,393161)))),]

write.csv(dt1,'dt_out_lung.csv')

for(i in 1:nrow(dt1)){
  cat('i=',i,'\n')
  
  options(timeout = 60000)
  link <- paste0(dt1$summaryStatistics[i],'/',dt1$accessionId[i],'.tsv.gz')
  deskfile <- paste0('H:/metabolite899/',dt1$accessionId[i],'.tsv.gz')
  download.file(link,deskfile)
}


