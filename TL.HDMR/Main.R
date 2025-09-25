
Metrics <- function(choX,trueX){
  
  TP <- sum(choX==1 & trueX==1)
  allP <- sum(trueX==1)
  
  FP <- sum(choX==1 & trueX==0)
  allN <- sum(trueX==0)
  
  res <- c(TP/allP,FP/allN)
  names(res) <- c('TPR','FPR')
  
  return(res)
  
}


Main <- function(ii,Nx,Ny,g,rho,P,hx2,hxy2,cv.lambd,Xcor,
                 NKdt,NKdt_distorb_prop,distorb_value,
                 Nx_kb,Ny_kb,g_kb,rho_kb,Xcor_kb,NKdt_rho){
  
  # dimn <- P*(NKdt+1)
  # multi_bx_cov <- ar1_cor(dimn,NKdt_rho)
  # diag(multi_bx_cov) <- NKdt_rho+0.1
  # 
  # if(NKdt_rho==0) diag(multi_bx_cov) <- NKdt_rho
  # 
  # sssig <- ar1_cor(P,rho)
  # multi_bx_cov[1:P,1:P] <- sssig
  # 
  # for(ink in 2:(NKdt+1)){
  #   sssig <- ar1_cor(P,Xcor_kb[ink-1])
  #   locc <- ((ink-1)*P+1):(ink*P)
  #   multi_bx_cov[locc,locc] <- sssig
  # }
  # # multi_bx_cov最左上角为source的matrix，依次是NKdt个matrix
  # 
  # gall <- max(max(g_kb),g)
  # multi_bx_cond <- mvrnorm(gall,rep(0,P*(NKdt+1)),multi_bx_cov)
  # 
  # bx_condall <- multi_bx_cond[,1:P]
  # if(nrow(bx_condall)!=g){
  #   bx_condall <- bx_condall[sample(1:nrow(bx_condall),g),]
  # }
  # 
  # multi_bx_condall <- list()
  # for(ink in 2:(NKdt+1)){
  #   locc <- ((ink-1)*P+1):(ink*P)
  #   multi_bx_condall[[ink-1]] <- multi_bx_cond[,locc]
  # }
  ####
  
  if(NKdt_rho==0){
    
    sssig <- ar1_cor(P,Xcor)
    bx_condall <- mvrnorm(g,rep(0,P),sssig)
    
    multi_bx_condall <- list()
    for(ink in 2:(NKdt+1)){
      sssig0 <- ar1_cor(P,Xcor_kb[ink-1])
      bx_condall0 <- mvrnorm(g_kb[ink-1],rep(0,P),sssig0) 
      multi_bx_condall[[ink-1]] <- bx_condall0
    }
    
  }else{
    
    dimn <- NKdt+1
    multi_bx_cov <- ar1_cor(dimn,NKdt_rho)
    
    bx_condall <- NULL
    multi_bx_condall <- vector('list',NKdt)
    for(kp in 1:P){
      gall <- max(max(g_kb),g)
      once <- mvrnorm(gall,rep(0,dimn),multi_bx_cov) 
      
      bx_condall <- cbind(bx_condall,once[,1])
      
      for(i_dim in 1:NKdt){
        multi_bx_condall[[i_dim]] <- cbind(multi_bx_condall[[i_dim]],once[,i_dim+1])
      }
      
    }
    
    
  }
  
  if(nrow(bx_condall)!=g){
    bx_condall <- bx_condall[sample(1:nrow(bx_condall),g),]
  }

  ###############
  fdata <- DataGeneration(Nx,Ny,g,rho,P,hx2,hxy2,Xcor,bx_condall)
  trueXe <- fdata$theta
  trueX <- as.numeric(fdata$theta!=0)
  
  fdata1 <- TransCondi(fdata)
  
  fdata_NK <- Multi_DataGeneration(Nx,Ny,g,rho,P,hx2,hxy2,Xcor,
                                   NKdt,NKdt_distorb_prop,distorb_value,
                                   Nx_kb,Ny_kb,g_kb,rho_kb,Xcor_kb,multi_bx_condall)
  
  fdata_NK1 <- list()
  for(olk in 1:length(fdata_NK)){
    fdata_NK1[[olk]] <- TransCondi(fdata_NK[[olk]])
  }
  
  ################# compared MR methods ##################################

  
  # covreg
  res.covreg <- MRcovreg(fdata,cv.lambd,trueX)

  # mvivw  lasso
  res.lasso <- MRLasso(fdata,cv.lambd,trueX,alpha=1)
  res.lasso$cv

  # mvivw LDA lasso
  res.LDAlasso <-MR.LDA.Lasso(fdata1,cv.lambd,trueX,alpha=1)
  res.LDAlasso$cv

  # mvivw  elastic-net
  res.elas <- MRLasso(fdata,cv.lambd,trueX,alpha=0.5)

  # mvivw  LDA elastic-net
  res.LDAelas <- MR.LDA.Lasso(fdata1,cv.lambd,trueX,alpha=0.5)

  # MR BMA
  # too slow

  # our method
  res_cisMVMR.MCP <- cisMRncvreg(fdata1,trueX,cv.lambd,penalname='MCP')
  res_cisMVMR.SCAD <- cisMRncvreg(fdata1,trueX,cv.lambd,penalname='SCAD')
  
  ## Transfer learning
  cisMVMR.Trans.lasso <- Oracle.Trans.lasso(fdata1,fdata_NK1,cv.lambd,Methodname='Lasso',trueX)
  cisMVMR.Trans.elas <- Oracle.Trans.lasso(fdata1,fdata_NK1,cv.lambd,Methodname='ElasticNet',trueX)
  cisMVMR.Trans.MCP <- Oracle.Trans.lasso(fdata1,fdata_NK1,cv.lambd,Methodname='MCP',trueX)
  cisMVMR.Trans.SCAD <- Oracle.Trans.lasso(fdata1,fdata_NK1,cv.lambd,Methodname='SCAD',trueX)

  return(list(res.covreg=res.covreg,
              res.lasso=res.lasso,
              res.LDAlasso=res.LDAlasso,
              res.elas=res.elas,
              res.LDAelas=res.LDAelas,
              res_cisMVMR.MCP=res_cisMVMR.MCP,
              res_cisMVMR.SCAD=res_cisMVMR.SCAD,
              cisMVMR.Trans.lasso=cisMVMR.Trans.lasso,
              cisMVMR.Trans.elas=cisMVMR.Trans.elas,
              cisMVMR.Trans.MCP=cisMVMR.Trans.MCP,
              cisMVMR.Trans.SCAD=cisMVMR.Trans.SCAD,
              trueX=trueX,
              trueXe=trueXe
              ))
  
}