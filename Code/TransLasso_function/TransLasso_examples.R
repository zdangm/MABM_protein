###An example
##coefficient generating functions used in the simulation
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
Coef.gen<- function(s, h,q=30, size.A0, M, sig.beta,sig.delta1, sig.delta2, p, exact=T){
  beta0 <- c(rep(sig.beta,s), rep(0, p - s))
  W <- rep.col(beta0,  M)# ten prior estimates
  W[1,]<-W[1,]-2*sig.beta
  for(k in 1:M){
    if(k <= size.A0){
      if(exact){
        samp0<- sample(1:p, h, replace=F)
        W[samp0,k] <-W[samp0,k] + rep(-sig.delta1, h)
      }else{
        W[1:100,k] <-W[1:100,k] + rnorm(100, 0, h/100)
      }
    }else{
      if(exact){
        samp1 <- sample(1:p, q, replace = F)
        W[samp1,k] <- W[samp1,k] + rep(-sig.delta2,q)
      }else{
        W[1:100,k] <-W[1:100,k] + rnorm(100, 0, q/100)
      }
    }
  }
  
  return(list(W=W, beta0=beta0))
}


source("~/TransLasso_function/TransLasso.R")
set.seed(123)
p = 500
s = 16
M = 20
sig.beta = 0.3

sig.z <- 1
n0 <- 150
M = 20
n.vec <- c(n0, rep(100, M))
Sig.X <- diag(1, p)
Niter = 200
l1=T

size.A0 = 12 
h=6
A0 = 1:size.A0
beta0<- 
coef.all <-Coef.gen( s, h = h, q = 2*s, size.A0 = size.A0,  M = M,   sig.beta = sig.beta,
             sig.delta1 = sig.beta, sig.delta2 = sig.beta+0.2, p = p, exact=F)
B <- cbind(coef.all$beta0, coef.all$W)
beta0 <- coef.all$beta0

###generate the data
X <- NULL
y <- NULL
for (k in 1:(M + 1)) {
  X <- rbind(X, rmvnorm(n.vec[k], rep(0, p), Sig.X))
  ind.k <- ind.set(n.vec, k)
  y <- c(y, X[ind.k, ] %*% B[, k] + rnorm (n.vec[k], 0, 1))
}
###compute init beta ####
beta.init <-
  as.numeric(glmnet(X[1:n.vec[1], ], y[1:n.vec[1]], lambda = sqrt(2 * log(p) / n.vec[1]))$beta)
###########Trans-Lasso#############
prop.re1 <- Trans.lasso(X, y, n.vec, I.til = 1:50, l1 = l1)
prop.re2 <- Trans.lasso(X, y, n.vec, I.til = 101:n.vec[1], l1=l1)
if(size.A0 > 0 & size.A0< M){ #Rank.re characterizes the performance of the sparsity index Rk
  Rank.re<- (sum(prop.re1$rank.pi[1:size.A0]<=size.A0) +
                     sum(prop.re2$rank.pi[1:size.A0]<=size.A0))/2/size.A0
}else{ Rank.re <- 1 }
beta.prop <- (prop.re1$beta.hat + prop.re2$beta.hat) / 2





