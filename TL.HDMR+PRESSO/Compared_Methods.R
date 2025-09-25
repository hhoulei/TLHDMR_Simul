
MRcovreg <- function(fdata,cv.lambd,trueX){
  
  P <- ncol(fdata$Bx_obs)
  # require the number of SNP above 10
  res1 <- mr_covreg(bx=fdata$Bx_obs[,1], 
                    bw=fdata$Bx_obs[,-1], 
                    by=fdata$By_obs, 
                    S=diag(1/fdata$sy_obs))
  Estbeta.covreg <- c(c(res1$thest),res1$a)
  names(Estbeta.covreg) <- paste0('X',1:P)
  # Metrics.covreg <- Metrics(as.numeric(Estbeta.covreg!=0),trueX)
  # MSE.covreg <- mean((Estbeta.covreg-trueXe)^2,na.rm=T)
  
  cvfit <- cv.mr_covreg(bx=fdata$Bx_obs[,1], 
                        bw=fdata$Bx_obs[,-1], 
                        by=fdata$By_obs, 
                        S=diag(1/fdata$sy_obs),
                        lambda=cv.lambd
  )
  
  cv.est <- cvfit$glmnet.fit$beta
  cv.est <- as.matrix(cv.est)
  cv.est <- rbind(1,cv.est)
  cv.covreg <- apply(cv.est,2,function(x) Metrics(as.numeric(x!=0),trueX))
  
  return(list(betall=Estbeta.covreg,
              cv=cv.covreg))
  
}


MRLasso <- function(fdata,cv.lambd,trueX,alpha){
  
  cvfit <- cv.glmnet(fdata$Bx_obs, fdata$By_obs, 
                     alpha = alpha, 
                     weights = fdata$sy_obs^(-2),
                     lambda=cv.lambd)
  cv.est <-  cvfit$glmnet.fit$beta
  cv.est <- as.matrix(cv.est)
  cv.elas  <- apply(cv.est,2,function(x) Metrics(as.numeric(x!=0),trueX))
  
  fit <- glmnet(fdata$Bx_obs, fdata$By_obs, 
                weights = fdata$sy_obs^(-2),
                alpha = alpha, 
                lambda = cvfit$lambda.min)
  Estbeta.elas <- as.vector(fit$beta)
  names(Estbeta.elas) <- paste0('X',1:P)
  # Metrics.elas <- Metrics(as.numeric(Estbeta.elas!=0),trueX)
  # MSE.elas <- mean((Estbeta.elas-trueXe)^2,na.rm=T)

  return(list(betall=Estbeta.elas,
              cv=cv.elas))
  
}

MR.LDA.Lasso <- function(fdata1,cv.lambd,trueX,alpha){
  
  cvfit <- cv.glmnet(fdata1$betaGX_new, fdata1$betaGY_new, 
                     alpha = alpha, 
                     lambda=cv.lambd)
  cv.est <-  cvfit$glmnet.fit$beta
  cv.est <- as.matrix(cv.est)
  cv.elas  <- apply(cv.est,2,function(x) Metrics(as.numeric(x!=0),trueX))
  
  fit <- glmnet(fdata1$betaGX_new, fdata1$betaGY_new, 
                alpha = alpha, 
                lambda = cvfit$lambda.min)
  Estbeta.elas <- as.vector(fit$beta)
  names(Estbeta.elas) <- paste0('X',1:P)
  
  return(list(betall=Estbeta.elas,
              cv=cv.elas))
  
}

