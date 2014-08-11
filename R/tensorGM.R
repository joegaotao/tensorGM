###########################################################################
#  Tensor Gaussian Graphical Model with l1 penalty
##########################################################################
#' @x tensor class by asTensor
#' @lambda a list of tuning parameters(scalar), or a list of matrices for adaptive penalty
#' @maxIter maximum number of iteration
#' @tol the tolerance value

#####main optimization function#######
glasso.tensor <- function(x, lambda, maxIter = 10, tol = 0.01){
  if(!inherits(x[[1]], "tensor"))
    stop("please input list of tensor data!")
  N <- length(x)
  dim <- x[[1]]@dim
  iter <- 0
  maxDel <- tol + 1  #maximal change of the estimators
  Dprod <- prod(dim)  #m=m1*m2***mK
  K <- length(dim)
  #initialization
  OmegaHat <- list()
  for(i in 1:K){
    OmegaHat <- c(OmegaHat, list(diag(dim[i])))
  } 
  while(maxDel > tol && iter < maxIter){ #condition to continue optimization
    iter <- iter+1
    maxDel <- 0
    for(kk in 1:K){ #update the kk-th mode
      S <- matrix(0, dim[kk], dim[kk])
      for(n in 1:N){
        for(i in 1:K){ #remove correlation
          root <- chol(OmegaHat[[i]])
          if(i != kk) 
            sample@data <- tensorProd(x[[n]], root, i)#####
        }
        ss <- toMatrix(sample, kk)
        S <- S + ss %*% t(ss) / (Dprod / dim[kk])
      }
      S <- S/N
      res <- glasso(S, rho = lambda[[kk]], penalize.diagonal = FALSE) #optimize on the kk-th mode
      oldOmega <- OmegaHat[[kk]]
      OmegaHat[[kk]] <- res$wi
      maxDel <- max(abs(OmegaHat[[kk]] - oldOmega), maxDel)
    }# update the kk(th) Omega
    print(paste("iter=",iter, ";maxDel=", round(maxDel,4)))
  }# total iteration
  
  #output
  for(i in 1:(K-1)){
    b <- (OmegaHat[[i]])[1,1]
    OmegaHat[[i]] <- OmegaHat[[i]] / b
    OmegaHat[[K]] <- OmegaHat[[K]] * b
  }
  ## log-likelihood
  # loglike = 
  # BIC = 
  return(list(OmegaHat = OmegaHat, BIC = 0))
}

