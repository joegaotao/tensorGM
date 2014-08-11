###############################################################
# nonparanormal methods
###############################################################
#' @x tensor data
#' @npn.func three nonparanormal transformation method, "shrinkage", "truncation", "skeptic"
#' @positive only for skeptic method, project transformed matrix into positive definite matrix
#' @npn.thresh default NULL
#' @verbose default print runing message


## Main function
tensor.npn = function(x, npn.func = "shrinkage", positive = TRUE, npn.thresh = NULL, verbose = TRUE){
  gcinfo(FALSE)
  n = nrow(x)
  d = ncol(x)
  x.col = colnames(x)
  x.row = rownames(x)
  
  # Shrinkaage transformation
  if(npn.func == "shrinkage"){
    if(verbose) cat("Conducting the nonparanormal (npn) transformation via shrunkun ECDF....")
    
    x = qnorm(apply(x,2,rank)/(n+1))
    x = x/sd(x[,1])
    
    if(verbose) cat("done.\n")
    rm(n,d,verbose)
    gc()	
    colnames(x) = x.col
    rownames(x) = x.row
  }
  
  # Truncation transformation
  if(npn.func == "truncation"){
    if(verbose) cat("Conducting nonparanormal (npn) transformation via truncated ECDF....")
    if(is.null(npn.thresh)) npn.thresh = 1/(4*(n^0.25)*sqrt(pi*log(n)))
    
    x = qnorm(pmin(pmax(apply(x,2,rank)/n, npn.thresh), 1-npn.thresh))
    x = x/sd(x[,1])
    
    if(verbose) cat("done.\n")
    rm(n,d,npn.thresh,verbose)
    gc()
    colnames(x) = x.col
    rownames(x) = x.row
  }
  
  if(npn.func == "skeptic"){
    if(verbose) cat("Conducting nonparanormal (npn) transformation via skeptic....")
    x = 2 * sin(pi/6 * cor(x, method = "spearman"))
    ### project x into positive definite matrix
    
    if(verbose) cat("done.\n")
    rm(n,d,verbose)
    gc()
    colnames(x) = x.col
    rownames(x) = x.col
  }
  
  return(x)
}