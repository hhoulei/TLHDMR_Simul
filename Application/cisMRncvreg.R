

cisMRncvreg <- function(fdata1,cv.lambd,penalname){
  
  #penalname=MCP,SCAD
  
  betaGX_new <- fdata1$betaGX_new
  betaGY_new <- fdata1$betaGY_new
  
  ################# ncvreg #####################
  
  ######
  cvfit1 <- ncvreg(betaGX_new, betaGY_new,
                   family = "gaussian",
                   penalty = penalname,
                   lambda=cv.lambd,
max.iter = 10000000) 
  loc1 <- cvfit1$convex.min
  Estbeta1 <- cvfit1$beta[,loc1]
  
  # cv.est <- cvfit1$beta
  # cv.est <- as.matrix(cv.est)
  # cv.est <- cv.est[-1,]
  # cv1  <- apply(cv.est,2,function(x) Metrics(as.numeric(x!=0),trueX))

  ######
  resall <- Estbeta1[-1]
  
  P <- ncol(betaGX_new)
  names(resall) <- paste0('X',1:P)
  
  return(resall)
}

