N=100
#d=c(30,30,30,30)  
d=c(10,10)  
# the dimension of the tensor, it can be quite arbitrary.
# e.g., d=(20,40,50) for the tensor of dimension (20,40,50)

#p=c(0.05,0.075,0.1,0.2) # the probability of non-zero entry for the concentration matrices
p = c(0.05, 0.075)
set.seed(200)
Omega = geneOmega(d, p)   #generate concentration matrices
sample.data <- geneSample(N, d, Omega$Root)  ## generate samples for estimation

## transform into tensor class
sample.data <- asTensor(sample.data, 1)

lambda <- list(0.01,0.01,0.01)  ## tuning parameter for l1 norm
res2 <- glasso.tensor(sample.data, lambda, maxIter=10,tol=0.01)  #estimation

########################
###evaluation of accuaracy based on matrix norms, specificity and sensitivity
accu<-matrix(0,nrow=length(d),ncol=7)
print("1Frob Norm 2Inf Norm 3Spectral Norm 4entry-wise Inf Norm 5SPE 6SEN 7MCC")
for(i in 1:length(d)){
  A=Omega$Omega[[i]]
  B=res2$OmegaHat[[i]]
  diff=A-B
  accu[i,3]=max(sqrt(eigen(t(diff)%*%diff)$values))
  accu[i,1]=sqrt(sum((diff)^2))
  diff=abs(diff)
  accu[i,2]=sum(max(rowSums(diff)))
  accu[i,4]=max(diff)
  posA=(abs(A)>0)
  posB=(abs(B)>0)
  tp=sum(posA & posB)-dim(A)[1]
  fp=sum((!posA) & posB)
  tn=sum((!posA) & (!posB))
  fn=sum((posA) & (!posB))
  accu[i,5]=tn/(tn+fp)
  accu[i,6]=tp/(tp+fn)
  dden=sqrt(tp+fp)*sqrt(tp+fn)*sqrt(tn+fp)*sqrt(tn+fn)
  if(dden>0){accu[i,7]=(tp*tn-fp*fn)/dden}
}
print("1Frob Norm 2Inf Norm 3Spectral Norm 4entry-wise Inf Norm 5SPE 6SEN 7MCC")
print(accu)
