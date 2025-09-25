
RemoveIV <- function(fdata1,bx_corMatrxi,methth,trueX,cv.lambd){
  
  giv <- nrow(fdata1$betaGX_new)
  Piv <- ncol(fdata1$betaGX_new)
  Pvalue <- NULL
  
  
  mc=50
  # cl <- makeCluster(mc)
  # registerDoParallel(cl)
  registerDoMC(mc)
  result <- foreach(gi=1:giv,
                    .packages = c("ncvreg",
                                  "MASS",
                                  "mrcovreg",'glmnet')
                    #.errorhandling =  "remove"
  ) %dopar% {
    
    fdata_once <- fdata1
    
    fdata_once$betaGX_new <- fdata_once$betaGX_new[-gi,]
    fdata_once$betaGY_new <- fdata_once$betaGY_new[-gi]
    fdata_once$sebetaGX_new <- fdata_once$sebetaGX_new[-gi,]
    fdata_once$sebetaGY_new <- fdata_once$sebetaGY_new[-gi]
    
    if(methth=='MCP'){
      res_cisMVMR.MCP <- cisMRncvreg(fdata_once,trueX,cv.lambd,penalname='MCP')
      resj <- res_cisMVMR.MCP
    }else if(methth=='SCAD'){
      res_cisMVMR.SCAD <- cisMRncvreg(fdata_once,trueX,cv.lambd,penalname='SCAD')
      resj <- res_cisMVMR.SCAD
    }else if(methth=='Lasso'){
      res.LDAlasso <-MR.LDA.Lasso(fdata_once,cv.lambd,trueX,alpha=1)
      resj <- res.LDAlasso
    }else if(methth=='ElasticNet'){
      res.elas <- MR.LDA.Lasso(fdata_once,cv.lambd,trueX,alpha=0.5)
      resj <- res.elas
    }
    
    
    thtaj <- resj$betall 
    predgammaj <- sum(fdata1$betaGX_new[gi]*thtaj)
    RSSobsj <- fdata1$betaGY_new[gi]-predgammaj
    
    Aj <- fdata1$betaGX_new[gi,]
    seAj <- fdata1$sebetaGX_new[gi,]
    bx_cor_Matrxionce <- bx_corMatrxi
    for(ro in 1:Piv){
      for(co in 1:Piv){
        bx_cor_Matrxionce[ro,co] <- seAj[ro]*seAj[co]*bx_corMatrxi[ro,co]
      }
    }
    
    nboot <- 100
    Pvalj <- 0
    for(kp in 1:nboot){
      Ajrandom <- mvrnorm(1,Aj,bx_cor_Matrxionce)
      gammarandom <- rnorm(1,predgammaj,fdata1$sebetaGY_new[gi])
      RSSexpjk <- gammarandom-sum(Ajrandom*thtaj)
      Pvalj <- Pvalj + as.numeric(RSSexpjk>RSSobsj)
    }
    Pvalj <- Pvalj/nboot
    
    Pvalj

  }
  
  #stopCluster(cl)
  
  Pvalue <- unlist(result)
  
  # for(gi in 1:giv){
  #   
  #   fdata_once <- fdata1
  #   
  #   fdata_once$betaGX_new <- fdata_once$betaGX_new[-gi,]
  #   fdata_once$betaGY_new <- fdata_once$betaGY_new[-gi]
  #   fdata_once$sebetaGX_new <- fdata_once$sebetaGX_new[-gi,]
  #   fdata_once$sebetaGY_new <- fdata_once$sebetaGY_new[-gi]
  #   
  #   if(methth=='MCP'){
  #     res_cisMVMR.MCP <- cisMRncvreg(fdata_once,trueX,cv.lambd,penalname='MCP')
  #     resj <- res_cisMVMR.MCP
  #   }else if(methth=='SCAD'){
  #     res_cisMVMR.SCAD <- cisMRncvreg(fdata_once,trueX,cv.lambd,penalname='SCAD')
  #     resj <- res_cisMVMR.SCAD
  #   }else if(methth=='Lasso'){
  #     res.LDAlasso <-MR.LDA.Lasso(fdata_once,cv.lambd,trueX,alpha=1)
  #     resj <- res.LDAlasso
  #   }else if(methth=='ElasticNet'){
  #     res.elas <- MR.LDA.Lasso(fdata_once,cv.lambd,trueX,alpha=0.5)
  #     resj <- res.elas
  #   }
  # 
  #   
  #   thtaj <- resj$betall 
  #   predgammaj <- sum(fdata1$betaGX_new[gi]*thtaj)
  #   RSSobsj <- fdata1$betaGY_new[gi]-predgammaj
  #   
  #   Aj <- fdata1$betaGX_new[gi,]
  #   seAj <- fdata1$sebetaGX_new[gi,]
  #   bx_cor_Matrxionce <- bx_corMatrxi
  #   for(ro in 1:Piv){
  #     for(co in 1:Piv){
  #       bx_cor_Matrxionce[ro,co] <- seAj[ro]*seAj[co]*bx_corMatrxi[ro,co]
  #     }
  #   }
  #   
  #   Pvalj <- 0
  #   for(kp in 1:1000){
  #     Ajrandom <- mvrnorm(1,Aj,bx_cor_Matrxionce)
  #     gammarandom <- rnorm(1,predgammaj,fdata1$sebetaGY_new[gi])
  #     RSSexpjk <- gammarandom-sum(Ajrandom*thtaj)
  #     Pvalj <- Pvalj + as.numeric(RSSexpjk>RSSobsj)
  #   }
  #   Pvalj <- Pvalj/1000
  #   
  #   Pvalue <- c(Pvalue,Pvalj)
  # }

  which(Pvalue<0.05)
  which(Pvalue<(0.05/giv))
  
  Pcut <- (0.05/giv)
  Pcut <- (0.05)
  
  fdata2 <- fdata1
  fdata2$betaGX <- fdata2$betaGX[Pvalue>=Pcut,]
  fdata2$betaGY <- fdata2$betaGY[Pvalue>=Pcut]
  fdata2$betaGX_new <- fdata2$betaGX_new[Pvalue>=Pcut,]
  fdata2$betaGY_new <- fdata2$betaGY_new[Pvalue>=Pcut]
  fdata2$sebetaGX_new <- fdata2$sebetaGX_new[Pvalue>=Pcut,]
  fdata2$sebetaGY_new <- fdata2$sebetaGY_new[Pvalue>=Pcut]
  
  return(list(fdata2=fdata2,
              removeIV=which(Pvalue<Pcut)))
  
}
