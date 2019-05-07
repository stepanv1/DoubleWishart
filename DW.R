# Simulations illustrating the main result
# of 'Exact Largest Eigenvalue Distribution for Doubly Singular Beta Ensemble'
# https://arxiv.org/abs/1905.01774
# Simulate CDF of lambda_max of doubly singular beta ensemble
# and compare to the exact formula
library(rootWishart)
library(rWishart)
library(corpcor)
# set parameteres of doubly singular beta ensemble (DSB)
p=120
m=100
q=6
n=500
#generate singular wishart matrices with Sigma = I_p
A=rSingularWishart(n, m, diag(1, p))
B=rSingularWishart(n, q, diag(1, p))
#compute \lambda_max for DSB samples
lmax<-rep(0,n)
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
  largestRoot <- max(crossprod(singWeights,  B[,,i] %*% singWeights))
  lmax[i] <- largestRoot
}
lmax <- sort(lmax)
#rescale to doubleWishart package setting:
lval <- sort(lmax/(1+lmax))
#plot empirical cdf:
DW<-doubleWishart(lval, p=q, n=m, m = p-m+q,  type = "multiple")
# plot on the original scale:
plot(ecdf(lmax), pch='.')
points(lmax, DW, col='red', pch='.')
