# Simulations illustrating the main result
# of 'Exact Largest Eigenvalue Distribution for Doubly Singular Beta Ensemble'
# https://arxiv.org/abs/1905.01774
# Simulate CDF of lambda_max of doubly singular beta ensemble
# and compare to the exact formula
library(rootWishart)
library(rWishart)
library(corpcor)
# set parameteres of doubly singular beta ensemble (DSB)
p=150
m=100
q=35
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

# Tracy-Widom distribution
hist(lmax,500)
hist(lval,500)
sqrt(150)
library(RMTstat)
# create Johstone TW approximation 
pt = q
nt = m
mt = p-m+q 
gamma <- 2*asin(sqrt((min(pt,nt)-0.5)/(mt+nt-1)))
phi <- 2*asin(sqrt((max(pt,nt)-0.5 )/(mt+nt-1)))
mu <- 2 * log(tan( 0.5*(phi+gamma)))
sigma <- (16/(mt+nt-1)^2)^(1/3)*(sin(phi+gamma)^2*sin(phi)*sin(gamma))^((-1)*(1/3))
Wp <- log(lval/(1-lval))
Z1 <- (Wp - mu)/sigma
TW <- ptw(Z1)
plot(ecdf(Z1))
points(Z1, TW, col='red', pch='.')



