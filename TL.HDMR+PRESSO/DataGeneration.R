
######################## Data Generation ############################
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}


Multi_DataGeneration <- function(Nx,Ny,g,rho,P,hx2,hxy2,Xcor,NKdt,NKdt_distorb_prop,distorb_value,
                                 Nx_kb,Ny_kb,g_kb,rho_kb,Xcor_kb,multi_bx_condall){
  
  
  multi_fdata <- list()
  
  for(k in 1:NKdt){
    
    NKdt_distorb <- runif(P,0,distorb_value)
    loc <- sample(1:P,P*NKdt_distorb_prop)
    hxy2_k <- hxy2
    hxy2_k[loc] <- NKdt_distorb[loc]
    
    Nx_k <- Nx_kb[k]
    Ny_k <- Ny_kb[k]
    g_k <- g_kb[k]
    rho_k <- rho_kb[k]
    Xcor_k <- Xcor_kb[k]
    
    oncebx_condall <- multi_bx_condall[[k]]
    if(nrow(oncebx_condall)!=g_k){
      oncebx_condall <- oncebx_condall[sample(1:nrow(oncebx_condall),g_k),]
    }
    
    multi_fdata[[k]] <- DataGeneration(Nx_k,Ny_k,g_k,rho_k,P,hx2,hxy2_k,Xcor,oncebx_condall)
    
  }
  
 return(multi_fdata)

}

DataGeneration <- function(Nx,Ny,g,rho,P,hx2,hxy2,Xcor,bx_condall){
  
  LD_mat <- ar1_cor(g,rho)
  
  LD_inv <- solve(LD_mat)
  L <- chol(LD_mat)
  #summary(c(LD_mat))
  
  # sssig <- diag(P)
  # sssig[upper.tri(sssig)] <- runif((P-1)*P/2,0.2,0.5)
  # sssig[lower.tri(sssig)] <- t(sssig)[lower.tri(t(sssig))]
  
  # sssig <- ar1_cor(P,Xcor)
  # bx_condall <- mvrnorm(g,rep(0,P),sssig)
  
  Bx_condi <- NULL
  Bx_obs <- NULL
  for(op in 1:P){
    
    #bx_cond <- rnorm(g,0,1)
    bx_cond <- bx_condall[,op]
    sigmax <- sqrt(hx2/c(t(bx_cond) %*% LD_mat %*% bx_cond))
    bx_cond <- sigmax*bx_cond
    #bx_cond[sample(1:g,10)] <- 0
    Bx_condi <- cbind(Bx_condi,bx_cond)
    
    epix <- rnorm(g,0,(1-hx2)/Nx)
    bx_obs <- LD_mat %*% bx_cond + t(L) %*% epix
    Bx_obs <- cbind(Bx_obs,bx_obs)
  }
  theta <- sqrt(hxy2)
  
  pleio <- c(rep(0.05,floor(g*0.2)),rep(0,g-floor(g*0.2)))
  
  By_condi <- Bx_condi %*% theta + pleio + rnorm(g,0,0.005)
  epiy <- rnorm(g,0,(1-hx2*sum(hxy2))/Ny)
  By_obs <- LD_mat %*% By_condi + t(L) %*% epiy
  
  colnames(Bx_condi) <- paste0('X',1:P)
  colnames(Bx_obs) <- paste0('X',1:P)
  By_condi <- c(By_condi)
  By_obs <- c(By_obs)
  
  sx_obs <- matrix(rep(sqrt((1-hx2)/Nx),g*P),ncol=P,nrow=g)
  Sigx_condi <- LD_mat %*% crossprod(t(sx_obs))
  sy_obs <- rep(sqrt((1-hx2*sum(hxy2))/Ny),g)
  Sigy_condi <- LD_mat %*% crossprod(t(sy_obs))

  fdata <- list(By_condi=By_condi,
                By_obs=By_obs,
                Bx_condi=Bx_condi,
                Bx_obs=Bx_obs,
                sx_obs=sx_obs,
                sy_obs=sy_obs,
                Sigx_condi=Sigx_condi,
                Sigy_condi=Sigy_condi,
                Nx=Nx,
                Ny=Ny,
                theta=theta,
                LD_mat=LD_mat,
                LD_inv=LD_inv
                )
  return(fdata)
}

Matrix12 <- function(Sigma){
  #Sigma <- sebetaGX
  resegi <- eigen(Sigma)
  lamda <- solve(resegi$vectors)%*%Sigma%*%(resegi$vectors)
  lamda_sqrt <- diag(sqrt(diag(lamda)))
  Sigma_sqrt <- (resegi$vectors)%*%lamda_sqrt%*%solve(resegi$vectors)
  Sigma_12 <- solve(Sigma_sqrt)
  return(Sigma_12)
}

TransCondi <- function(fdata){

  ################# transfer marginal to conditional #####################
  
  betaGX <- fdata$LD_inv %*% fdata$Bx_obs
  betaGY <- fdata$LD_inv %*% fdata$By_obs
  #sebetaGX <-  LD_inv %*% (LD_mat * crossprod(t(Sigx_condi))) %*% LD_inv
  sebetaGY <- fdata$LD_inv %*% (fdata$LD_mat * crossprod(t(fdata$Sigy_condi))) %*% fdata$LD_inv
  
  sebetaGY_new <- Matrix12(sebetaGY)
  betaGY_new <- sebetaGY_new %*% betaGY
  betaGX_new <- sebetaGY_new %*% betaGX
  
  sebetaGX_new <- matrix(1,nrow=nrow(betaGX_new),ncol=ncol(betaGX_new))
  sebetaGY_new <- matrix(1,nrow=nrow(betaGY_new),ncol=ncol(betaGY_new))
  
  colnames(sebetaGX_new) <- paste0('X',1:P)
  colnames(betaGX_new) <- paste0('X',1:P)
  
  return(list(betaGX=betaGX,
              betaGY=betaGY,
              betaGX_new=betaGX_new,
              betaGY_new=betaGY_new,
              sebetaGX_new=sebetaGX_new,
              sebetaGY_new=sebetaGY_new
              ))
  
}

