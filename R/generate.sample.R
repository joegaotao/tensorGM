##########simulation#######
geneOmega <- function(d,p){
  K <- length(d)  #dimension length
  AllOmega <- list()
  AllSigmaRoot <- list()
  for(i in 1:K){
    curE <- d[i] * (d[i] - 1) / 2  #the number of entries for the matrix
    curM <- matrix(0, d[i], d[i]) #new matrix
    pos <- rgeom(curE * p[i] * 2,p[i]) + 1
    for(j in 2:length(pos))
      pos[j]=pos[j] + pos[j-1]
    pos <- pos[pos <= curE]
    sign <- ifelse(runif(length(pos),-1,1)>0,1,-1)
    upper <- rep(0,curE)
    upper[pos] <- sign*runif(length(pos),0.2,0.8)
    curM[upper.tri(curM)] <- upper
    curM <- curM + t(curM)
    curM <- curM %*% diag(1 / colSums(1.2 * abs(curM) + 0.001))
    curM <- (curM + t(curM)) / 2
    diag(curM) <- 1
    AllOmega <- c(AllOmega,list(curM))   # store the concentration matrices
    AllSigmaRoot <- c(AllSigmaRoot, list(solve(chol(curM)))) # store their "squre root", which will be used in generating samples
  }
  list(Omega = AllOmega, Root = AllSigmaRoot)
}


#######################
#########################
####generate n samples for estimation
geneSample <- function(N, d, AllSigmaRoot){
  sample <- new("tensor",dim = d)
  D <- prod(d)
  K <- length(d)
  sample.data <- matrix(0, nr = N, nc = D)
  for(n in 1:N){
    sample@data <- rnorm(D)
    for(kk in 1:K){
      sample@data<-tensorProd(sample,AllSigmaRoot[[kk]],kk)
    }
    sample.data[n, ] <- sample@data
  }
  sample.data <- array(sample.data, dim = c(N, d))
  return(sample.data)
}


