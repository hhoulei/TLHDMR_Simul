
Metrics <- function(choX,trueX){
  
  TP <- sum(choX==1 & trueX==1)
  allP <- sum(trueX==1)
  
  FP <- sum(choX==1 & trueX==0)
  allN <- sum(trueX==0)
  
  res <- c(TP/allP,FP/allN)
  names(res) <- c('TPR','FPR')
  
  return(res)
  
}

DATA_condi <- function(fdata1,fdata_NK0){
  
  betaGX_A0 <- fdata1$betaGX_new
  betaGY_A0 <- fdata1$betaGY_new
  
  betaGX_AK <- betaGX_A0
  betaGY_AK <- betaGY_A0
  for(oak in 1:length(fdata_NK0)){
    betaGX_AK <- rbind(betaGX_AK,fdata_NK0[[oak]]$betaGX_new)
    betaGY_AK <- c(betaGY_AK,fdata_NK0[[oak]]$betaGY_new)
  }
  
  return(list(betaGX_A0=betaGX_A0,
              betaGY_A0=betaGY_A0,
              betaGX_AK=betaGX_AK,
              betaGY_AK=betaGY_AK))
}


Main <- function(ii,Nx,Ny,g,rho,P,hx2,hxy2,cv.lambd,Xcor,
                 NKdt,NKdt_distorb_prop,distorb_value,
                 Nx_kb,Ny_kb,g_kb,rho_kb,Xcor_kb,NKdt_rho,r.A0,C0){
  
  ###################
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

  ################# data generation ##################################
  fdata <- DataGeneration(Nx,Ny,g,rho,P,hx2,hxy2,Xcor,bx_condall)
  trueXe <- fdata$theta
  trueX <- as.numeric(fdata$theta!=0)
  
  fdata1 <- TransCondi(fdata)
  
  fdata_NK <- Multi_DataGeneration(Nx,Ny,g,rho,P,hx2,hxy2,Xcor,
                                   NKdt,NKdt_distorb_prop,distorb_value,
                                   Nx_kb,Ny_kb,g_kb,rho_kb,Xcor_kb,multi_bx_condall)
  
  fdata_NK0 <- list()
  for(olk in 1:length(fdata_NK)){
    fdata_NK0[[olk]] <- TransCondi(fdata_NK[[olk]])
  }
  
  ################# select datasets ##################################
  
  #Hh <- sum(hxy2!=0)*sqrt(log(P/g))
  
  # res.select <- SelectDt.AK(fdata1,fdata_NK0,cv.lambd,
  #                           Methodname='MCP',
  #                           r.A0,C0)
  # res.select$diff.loss
  # res.select$loss.A0.sd
  # res.select$Ind.loc
  
  
  ################# compared MR methods ##################################

  # mvivw LDA lasso
  res.LDAlasso <-MR.LDA.Lasso(fdata1,cv.lambd,trueX,alpha=1)
  #res.LDAlasso$cv

  # mvivw  LDA elastic-net
  res.LDAelas <- MR.LDA.Lasso(fdata1,cv.lambd,trueX,alpha=0.5)
  #res.LDAelas$cv

  # our method
  res_cisMVMR.MCP <- cisMRncvreg(fdata1,trueX,cv.lambd,penalname='MCP')
  #res_cisMVMR.MCP$cv
  res_cisMVMR.SCAD <- cisMRncvreg(fdata1,trueX,cv.lambd,penalname='SCAD')
  #res_cisMVMR.SCAD$cv
  
  
  ## Transfer learning
  fdata_NK1 <- DATA_condi(fdata1,fdata_NK0)
  cisMVMR.Trans.lasso.nosel <- Oracle.Trans.lasso(fdata_NK1$betaGX_A0,fdata_NK1$betaGY_A0,
                                            fdata_NK1$betaGX_AK,fdata_NK1$betaGY_AK,
                                            cv.lambd,Methodname='Lasso',trueX)
  cisMVMR.Trans.elas.nosel <- Oracle.Trans.lasso(fdata_NK1$betaGX_A0,fdata_NK1$betaGY_A0,
                                           fdata_NK1$betaGX_AK,fdata_NK1$betaGY_AK,
                                           cv.lambd,Methodname='ElasticNet',trueX)
  cisMVMR.Trans.MCP.nosel <- Oracle.Trans.lasso(fdata_NK1$betaGX_A0,fdata_NK1$betaGY_A0,
                                          fdata_NK1$betaGX_AK,fdata_NK1$betaGY_AK,
                                          cv.lambd,Methodname='MCP',trueX)
  cisMVMR.Trans.SCAD.nosel <- Oracle.Trans.lasso(fdata_NK1$betaGX_A0,fdata_NK1$betaGY_A0,
                                           fdata_NK1$betaGX_AK,fdata_NK1$betaGY_AK,
                                           cv.lambd,Methodname='SCAD',trueX)
  rm(fdata_NK1)
  
  ## Transfer learning Lasso
  res.select <- SelectDt.AK(fdata1,fdata_NK0,cv.lambd,
                            Methodname='Lasso',trueX,
                            r.A0,C0)
  res.select$diff.loss
  C0*max(res.select$loss.A0.sd,0.01)
  sel.NK.lasso <- res.select$Ind.loc
  fdata_NK1 <- fdata_NK0[sel.NK.lasso]
  fdata_NK1 <- DATA_condi(fdata1,fdata_NK1)
  cisMVMR.Trans.lasso <- Oracle.Trans.lasso(fdata_NK1$betaGX_A0,fdata_NK1$betaGY_A0,
                                            fdata_NK1$betaGX_AK,fdata_NK1$betaGY_AK,
                                            cv.lambd,Methodname='Lasso',trueX)
  rm(fdata_NK1)
  
  ## Transfer learning ElasticNet
  res.select <- SelectDt.AK(fdata1,fdata_NK0,cv.lambd,
                            Methodname='ElasticNet',trueX,
                            r.A0,C0)
  res.select$diff.loss
  C0*max(res.select$loss.A0.sd,0.01)
  sel.NK.elas <- res.select$Ind.loc
  fdata_NK1 <- fdata_NK0[sel.NK.elas]
  fdata_NK1 <- DATA_condi(fdata1,fdata_NK1)
  cisMVMR.Trans.elas <- Oracle.Trans.lasso(fdata_NK1$betaGX_A0,fdata_NK1$betaGY_A0,
                                           fdata_NK1$betaGX_AK,fdata_NK1$betaGY_AK,
                                           cv.lambd,Methodname='ElasticNet',trueX)
  rm(fdata_NK1)
  
  ## Transfer learning MCP
  res.select <- SelectDt.AK(fdata1,fdata_NK0,cv.lambd,
                            Methodname='MCP',trueX,
                            r.A0,C0)
  res.select$diff.loss
  C0*max(res.select$loss.A0.sd,0.01)
  sel.NK.MCP <- res.select$Ind.loc
  fdata_NK1 <- fdata_NK0[sel.NK.MCP]
  fdata_NK1 <- DATA_condi(fdata1,fdata_NK1)
  cisMVMR.Trans.MCP <- Oracle.Trans.lasso(fdata_NK1$betaGX_A0,fdata_NK1$betaGY_A0,
                                          fdata_NK1$betaGX_AK,fdata_NK1$betaGY_AK,
                                          cv.lambd,Methodname='MCP',trueX)
  rm(fdata_NK1)
  
  ## Transfer learning SCAD
  res.select <- SelectDt.AK(fdata1,fdata_NK0,cv.lambd,
                            Methodname='SCAD',trueX,
                            r.A0,C0)
  sel.NK.SCAD <- res.select$Ind.loc
  fdata_NK1 <- fdata_NK0[sel.NK.SCAD]
  fdata_NK1 <- DATA_condi(fdata1,fdata_NK1)
  cisMVMR.Trans.SCAD <- Oracle.Trans.lasso(fdata_NK1$betaGX_A0,fdata_NK1$betaGY_A0,
                                           fdata_NK1$betaGX_AK,fdata_NK1$betaGY_AK,
                                           cv.lambd,Methodname='SCAD',trueX)

  return(list(res.LDAlasso=res.LDAlasso,
              res.LDAelas=res.LDAelas,
              res_cisMVMR.MCP=res_cisMVMR.MCP,
              res_cisMVMR.SCAD=res_cisMVMR.SCAD,
              cisMVMR.Trans.lasso.nosel=cisMVMR.Trans.lasso.nosel,
              cisMVMR.Trans.elas.nosel=cisMVMR.Trans.elas.nosel,
              cisMVMR.Trans.MCP.nosel=cisMVMR.Trans.MCP.nosel,
              cisMVMR.Trans.SCAD.nosel=cisMVMR.Trans.SCAD.nosel,
              cisMVMR.Trans.lasso=cisMVMR.Trans.lasso,
              cisMVMR.Trans.elas=cisMVMR.Trans.elas,
              cisMVMR.Trans.MCP=cisMVMR.Trans.MCP,
              cisMVMR.Trans.SCAD=cisMVMR.Trans.SCAD,
              sel.NK.lasso=sel.NK.lasso,
              sel.NK.elas=sel.NK.elas,
              sel.NK.MCP=sel.NK.MCP,
              sel.NK.SCAD=sel.NK.SCAD,
              trueX=trueX,
              trueXe=trueXe
              ))
  
}