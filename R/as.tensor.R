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


