# Simulate CDF of doubly singularbeta ensemble
# and compare to the exact formula
library(rootWishart)
library(rWishart)
library(corpcor)
# set parameteres of doubly singular beta ensemble (DSB)
p=500
m=100
q=6
n=5000
#generate singular wishart matrices with Sigma = I_p
A=rSingularWishart(n, m, diag(1, p))
B=rSingularWishart(n, q, diag(1, p))
#compute \lambda_max for DSB samples
lmax<-rep(0,n)
for(i in 1:n){
  svdA <- svd(A[,,i])
  rankVr <- corpcor::rank.condition(A[,,i])$rank
  eigVecA <- svdA$v[, 1:rankVr]
  eigValAInv <- 1/svdA$d[1:rankVr]
  Xp <- eigVecA %*% diag(eigValAInv)
  C <- crossprod(Xp, B[,,i] %*% Xp)
  svdC <- corpcor::fast.svd(C)
  Xpp <- svdC$u
  singWeights <- Xp %*% Xpp
  largestRoot <- max(crossprod(singWeights,  B[,,i] %*% singWeights))
  lmax[i] <- largestRoot
}
lmax <- sort(lmax)
#rescale to doubleWishart package setting:
lval <- sort(lmax/(1+lmax))
#plot empirical cdf:
DW<-doubleWishart(lval, p=q, n=m, m = p-m+q,  type = "multiple")
# plot on the original scale:
plot(ecdf(lmax))
points(lmax, DW, col='red', pch='.')





