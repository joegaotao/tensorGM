#difine the tensor class and its basic operation
dyn.load("C:\\Users\\gaotao01\\Documents\\tensorGM\\src\\product.dll")
dyn.load("C:\\Users\\gaotao01\\Documents\\tensorGM\\src\\expansion.dll")
library(glasso)


#define the tensor class, 
#@data stores a vector, e.g., for a three mode tensor: a(1,1,1) a(1,1,2)...a(1,1,m_3),a(1,2,1)..
#This a little bit different from tensor vectorization, because on indices on the right runs faster
setClass("tensor", representation(data = "numeric",
                                 dim = "numeric"))


setMethod("initialize","tensor",
          function(.Object,dim){
            .Object@data <- rep(0, prod(dim))
            .Object@dim <- dim
            .Object
          })

#define the product between tensor X and matrix y on the n-th mode
.tensorProd <- function(x, y, n){
  sizeC <- prod(x@dim) * dim(y)[1] / x@dim[n]
  .C("product",
     as.double(x@data),
     as.double(as.vector(t(y))),  #byRow
     as.integer(x@dim),
     as.integer(dim(y)),
     as.integer(length(x@dim)),
     cc = double(sizeC),
     as.integer(n)
     )$cc
}


setGeneric("tensorProd", function(x, y, n) standardGeneric("tensorProd"))

setMethod("tensorProd", c("tensor", "matrix", "numeric"), .tensorProd)

#expand a tensor x to matrix by the n-th mode
.toMatrix<-function(x, n){
  if(n > length(x@dim)) 
    stop("n should not be larger than then dimension of the tensor")
  sizeC <- prod(x@dim)
  res <- .C("expansion",
            as.double(x@data),
            as.integer(x@dim),
            as.integer(length(x@dim)),
            cc = double(sizeC),
            as.integer(n))
  matrix(res$cc, nrow = x@dim[n], ncol = sizeC / x@dim[n], byrow=TRUE)
}

setGeneric("toMatrix",function(x,n) standardGeneric("toMatrix"))
setMethod("toMatrix", c("tensor","numeric"), .toMatrix)


### index claims the sample-way
.asTensor <- function(x, index){
  D <- dim(x)
  D.len <- length(D)
  N <- D[index]
  left.dim <- D[-index]
  nprod <- prod(left.dim)
  
  sample.tensor <- list()
  #sample.id <- sapply(seq_len(N), function(x) seq(x, by = N, len = nprod))
  for(i in seq_len(N)){
    sample.tensor[[i]] <- new("tensor", dim = left.dim)
    id <- paste(c("x[", rep(",", index - 1), i, rep(",", D.len - index), "]"), collapse = "")
    sample.tensor[[i]]@data <- as.vector(eval(parse(text = id)))
  }
  return(sample.tensor)
}

setGeneric("asTensor",function(x, index) standardGeneric("asTensor"))
setMethod("asTensor", c("array", "numeric"), .asTensor)


