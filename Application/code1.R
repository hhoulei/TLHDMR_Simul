
library(utils)
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
 
for(i in 319:nrow(dt1)){
  cat('i=',i,'\n')
  
  options(timeout = 60000)
  link <- paste0(dt1$summaryStatistics[i],'/harmonised/',dt1$accessionId[i],'.h.tsv.gz')
  deskfile <- paste0('H:/metabolite899/',dt1$accessionId[i],'.h.tsv.gz')
  download.file(link,deskfile)
}
