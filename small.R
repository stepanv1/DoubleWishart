###################################################################
library(rootWishart)
library(rWishart)
library(corpcor)
# Simulate CDF of lambda_min of doubly singular beta ensemble
# and compare to the exact formula

p=20
m=19
q=2
n=500
#generate singular wishart matrices with Sigma = I_p
A=rSingularWishart(n, m, diag(1, p))
B=rSingularWishart(n, q, diag(1, p))

lmin<-rep(0,n)
rnks<-rep(0,n)
for(i in 1:n){
  svdA <- eigen(A[,,i])
  rankVr <- corpcor::rank.condition(A[,,i])$rank
  rnks[i] <- rankVr
  eigVecA <- svdA$vectors[, 1:rankVr]
  eigValAInv <- sqrt(1/svdA$values[1:rankVr])
  Xp <- eigVecA %*% diag(eigValAInv)
  C <- crossprod(Xp, B[,,i] %*% Xp)
  svdC <- corpcor::fast.svd(C)
  Xpp <- svdC$u
  singWeights <- Xp %*% Xpp
  smallestRoot <- (crossprod(singWeights,  B[,,i] %*% singWeights))[q,q]
  lmin[i] <- smallestRoot
}
lmin <- sort(lmin)
#rescale to doubleWishart package setting:
lvalmin <- sort(lmin/(1+lmin))
#plot empirical cdf:
DW<-doubleWishart(lvalmin, p=q, n=p-m+q, m = m,  type = "multiple")
# plot on the original scale:
plot(ecdf(lmin), pch='.')
points(lmin, DW, col='red', pch='.')





