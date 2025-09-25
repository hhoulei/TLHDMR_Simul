

cisMRncvreg <- function(fdata1,trueX,cv.lambd,penalname){
  
  #penalname=MCP,SCAD
  
  betaGX_new <- fdata1$betaGX_new
  betaGY_new <- fdata1$betaGY_new
  
  ################# ncvreg #####################
  
  ######
  cvfit1 <- ncvreg(betaGX_new, betaGY_new,
                   family = "gaussian",
                   penalty = penalname,
                   lambda=cv.lambd) 
  loc1 <- cvfit1$convex.min
  Estbeta1 <- cvfit1$beta[,loc1]
  
  cv.est <- cvfit1$beta
  cv.est <- as.matrix(cv.est)
  cv.est <- cv.est[-1,]
  cv1  <- apply(cv.est,2,function(x) Metrics(as.numeric(x!=0),trueX))

  ######
  resall <- Estbeta1[-1]
  
  P <- ncol(betaGX_new)
  names(resall) <- paste0('X',1:P)
  
  return(list(betall=resall,
               cv=cv1))
}



# cisMRncvreg.Permutation <- function(Est_ori,fdata1,nPerm){
#   
#   betaGX_new <- fdata1$betaGX_new
#   betaGY_new <- fdata1$betaGY_new
#   
#   res_all <- NULL
#   for(pm in 1:nPerm){
#     
#     cat('cisMRncvreg Permutation: ',pm,'\n')
#     
#     ns <- sample(1:length(betaGY_new),length(betaGY_new))
#     betaGY_new1 <- betaGY_new[ns]
#     
#     resonce <- cisMRncvreg(fdata1)
#     res_all <- rbind(res_all,resonce)
#     
#   }
#   
#   #calculate the p-value
#   p_val = rep(1,length(Est_ori))
#   for(i in 1:length(Est_ori)){
#     p_val[i]=(sum(res_all[,i]>Est_ori[i])+1)/(length(res_all[,i])+1)
#   }
#   
#   names(p_val) <-  names(Est_ori)
#   
#   return(p_val)
# }


Oracle.Trans.lasso <- function(fdata1,fdata_NK1,cv.lambd,Methodname,trueX){
  
  betaGX_A0 <- fdata1$betaGX_new
  betaGY_A0 <- fdata1$betaGY_new
  
  betaGX_AK <- betaGX_A0
  betaGY_AK <- betaGY_A0
  for(oak in 1:length(fdata_NK1)){
    betaGX_AK <- rbind(betaGX_AK,fdata_NK1[[oak]]$betaGX_new)
    betaGY_AK <- c(betaGY_AK,fdata_NK1[[oak]]$betaGY_new)
  }
  
  p <- ncol(betaGX_AK)
  g <- nrow(betaGX_AK)
  g0 <- nrow(betaGX_A0)
  
  if(Methodname=='Lasso'){
    
    ### all data containing the source data
    cv.init <- cv.glmnet(betaGX_AK, 
                         betaGY_AK, 
                         alpha = 1,
                         #lambda=seq(1,0.1,length.out=10)*sqrt(2*log(p)/g)
                         lambda=cv.lambd*sqrt(2*log(p)/g))
    cv.lam.const <- cv.init$lambda/sqrt(2*log(p)/g)
    cv.w.kA <- as.matrix(cv.init$glmnet.fit$beta)
    cv.w.kA<- t(apply(cv.w.kA,1,function(x) x*(abs(x)>=cv.init$lambda)))
    
    
    lam.const <- cv.init$lambda.min/sqrt(2*log(p)/g)
    loc.min <- which(cv.init$lambda==cv.init$lambda.min)
    w.kA <- as.numeric(glmnet(betaGX_AK,
                              betaGY_AK,
                              alpha = 1,
                              lambda=lam.const*sqrt(2*log(p)/g))$beta)
    w.kA<-w.kA*(abs(w.kA)>=lam.const*sqrt(2*log(p)/g))
    
    ### only source data
    cv.delta.kA <- glmnet(x=betaGX_A0,
                          y=betaGY_A0-betaGX_A0%*%w.kA, 
                          alpha = 1,
                          lambda=cv.lam.const*sqrt(2*log(p)/g0))
    cv.delta.kA <- as.matrix(cv.delta.kA$beta)
    cv.delta.kA<- t(apply(cv.delta.kA,1,function(x) x*(abs(x)>=cv.lam.const*sqrt(2*log(p)/g0))))
    
    delta.kA <- as.numeric(glmnet(x=betaGX_A0,
                                  y=betaGY_A0-betaGX_A0%*%w.kA, 
                                  alpha = 1,
                                  lambda=lam.const*sqrt(2*log(p)/g0))$beta)
    delta.kA<-delta.kA*(abs(delta.kA)>=lam.const*sqrt(2*log(p)/g0))
    
    # combine
    beta.kA <- w.kA + delta.kA
    cv.beta.kA <- cv.w.kA + cv.delta.kA
    
  }else 
    if(Methodname=='ElasticNet'){
    
    ### all data containing the source data
    cv.init <- cv.glmnet(betaGX_AK, 
                         betaGY_AK, 
                         alpha = 0.5,
                         #lambda=seq(1,0.1,length.out=10)*sqrt(2*log(p)/g)
                         lambda=cv.lambd*sqrt(2*log(p)/g))
    cv.lam.const <- cv.init$lambda/sqrt(2*log(p)/g)
    cv.w.kA <- as.matrix(cv.init$glmnet.fit$beta)
    cv.w.kA<- t(apply(cv.w.kA,1,function(x) x*(abs(x)>=cv.init$lambda)))
    
    
    lam.const <- cv.init$lambda.min/sqrt(2*log(p)/g)
    loc.min <- which(cv.init$lambda==cv.init$lambda.min)
    w.kA <- as.numeric(glmnet(betaGX_AK,
                              betaGY_AK,
                              alpha = 0.5,
                              lambda=lam.const*sqrt(2*log(p)/g))$beta)
    w.kA<-w.kA*(abs(w.kA)>=lam.const*sqrt(2*log(p)/g))
    
    ### only source data
    cv.delta.kA <- glmnet(x=betaGX_A0,
                          y=betaGY_A0-betaGX_A0%*%w.kA, 
                          alpha = 0.5,
                          lambda=cv.lam.const*sqrt(2*log(p)/g0))
    cv.delta.kA <- as.matrix(cv.delta.kA$beta)
    cv.delta.kA<- t(apply(cv.delta.kA,1,function(x) x*(abs(x)>=cv.lam.const*sqrt(2*log(p)/g0))))
    
    delta.kA <- as.numeric(glmnet(x=betaGX_A0,
                                  y=betaGY_A0-betaGX_A0%*%w.kA, 
                                  alpha = 0.5,
                                  lambda=lam.const*sqrt(2*log(p)/g0))$beta)
    delta.kA<-delta.kA*(abs(delta.kA)>=lam.const*sqrt(2*log(p)/g0))
    
    # combine
    beta.kA <- w.kA + delta.kA
    cv.beta.kA <- cv.w.kA + cv.delta.kA
    
  }else 
    if(Methodname=='MCP' | Methodname=='SCAD'){
    
      penalname <- Methodname

      ### all data containing the source data
      cv.init <- ncvreg(betaGX_AK,
                        betaGY_AK,
                        family = "gaussian",
                        penalty = penalname,
                        lambda=cv.lambd*sqrt(2*log(p)/g))
      
      cv.lam.const <- cv.init$lambda/sqrt(2*log(p)/g)
      cv.w.kA <- as.matrix(cv.init$beta)
      cv.w.kA <- cv.w.kA[-1,]
      cv.w.kA<- t(apply(cv.w.kA,1,function(x) x*(abs(x)>=cv.init$lambda)))
      
      loc.min <- cv.init$convex.min
      lam.const <- cv.init$lambda[loc.min]
      w.kA <- cv.w.kA[,loc.min]
      
      ### only source data
      cv.delta.kA <- ncvreg(betaGX_A0,
                            betaGY_A0-betaGX_A0%*%w.kA, 
                            family = "gaussian",
                            penalty = penalname,
                            lambda=cv.lam.const*sqrt(2*log(p)/g0))
      cv.delta.kA <- as.matrix(cv.delta.kA$beta)
      cv.delta.kA <- cv.delta.kA[-1,]
      cv.delta.kA<- t(apply(cv.delta.kA,1,function(x) x*(abs(x)>=cv.lam.const*sqrt(2*log(p)/g0))))
      
      delta.kA <- cv.delta.kA[,loc.min]
      
      # combine
      beta.kA <- w.kA + delta.kA
      cv.beta.kA <- cv.w.kA + cv.delta.kA

    }
  
  cv.beta.kA1  <- apply(cv.beta.kA,2,function(x) Metrics(as.numeric(x!=0),trueX))
  
  return(list(betall=beta.kA,
              cv=cv.beta.kA1))
  
}

