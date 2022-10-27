library(Rcpp); library(MASS); sourceCpp("~/Library/CloudStorage/Dropbox/ongoing_works/B RRR/BRRR with missing Y/Rcodes/testabcd.cpp")
sig2 = 1
# random data generation Y = XBZ , (nl) = (np)*(pk)*(kl)
n = 100    # samples in Y, X
l = 10    # response in Y,Z
k = 20    # latent feature samples in Z
p = 100    # predictors in X
r = 2     # true rank

X = matrix(rnorm(n*p), ncol =p) 
Z = matrix(rnorm(k*l), ncol =l) 
B = matrix(rnorm(p*r),nr=p) %*% t(matrix(rnorm(k*r),nr=k ))  *2 + rnorm(p*k, sd=0.1)
Y = X%*%B%*%Z  + rnorm(n*l)

# 2 cases of ols fail: when k >l l and when n< p
XX = t(X)%*%X       # p x p 
XYZ = t(X) %*% Y %*%t(Z)  
ZZ = Z%*%t(Z)       # k x k

ols = ginv(XX) %*% XYZ %*% ginv(ZZ) ; c( mean((ols - B )^2) ,mean((ols - B )^2)/mean(B^2) )

SVD = svd(ols); plot(SVD$d,type = 'b');grid(lwd = 2)
rank_selected = 4
SVD$d [ rank_selected : min(p,k) ] = 0
Mhat = SVD$u %*% diag(SVD$d) %*% t(SVD$v) ; mean((Mhat - B )^2)

### MALA
lam = 1  # in the prior
ystar = diag(p)*lam^2
Iters = 10000
burnin = 1000

Bm = matrix(data=0,nr=p,nc=k)
M = ols
h = 1/(p*k)^1.7
a = 0  ; start_time <- Sys.time()
for(s in 1:Iters){
  MtM = tcrossprod(M)
  tam1 = solve(ystar + MtM, M )
  hxtyxm = h*(XYZ - mymatmul(XX,M,ZZ) )
  tam = M + hxtyxm/sig2 + h*(p+k+2)*tam1 + sqrt(2*h)*rnorm(p*k)
  
  tattam = tcrossprod(tam)
  pro.tam = -sum((Y- mymatmul(X,tam,Z) )^2)/(2*sig2) - 0.5*(p+k+2)*determinant(ystar + tattam)[[1]][1]
  pro.M = -sum((Y- mymatmul(X,M,Z) )^2)/(2*sig2) - 0.5*(p+k+2)*determinant(ystar + MtM)[[1]][1]
  
  tam2 = solve(ystar + tattam,tam)
  tran.m = -sum((M-tam-h*(XYZ - mymatmul(XX,tam,ZZ) )/sig2 - h*(p+k+2)*tam2 )^2)/(4*h)
  tran.tam = -sum((tam-M- hxtyxm/sig2 - h*(p+k+2)*tam1 )^2)/(4*h)
  
  pro.trans = pro.tam+tran.m-pro.M-tran.tam
  if(log(runif(1)) <= pro.trans){
    M = tam
    print(s)
    a = a+1
  } 
  if (s>burnin)Bm = Bm + M/(Iters-burnin)
}
end_time <- Sys.time();  end_time - start_time; a/Iters

#Langevin MC for BRRR
M1 = ols  # initialize
h = h/n/n
M_lmc = matrix(data=0, nr=p,nc=k)
start_time <- Sys.time()
for(s in 1:Iters){
  tam = solve(ystar + tcrossprod(M1),M1)
  M1 = M1 + h*(XYZ - mymatmul(XX,M1,ZZ) )/sig2 + h*(p+k+2)*tam + sqrt(2*h)*rnorm(p*k)
  if(s>burnin) M_lmc =  M_lmc + M1/(Iters-burnin)
}
end_time <- Sys.time();  end_time - start_time

c(mean((X%*% (M_lmc - B) %*%Z )^2), mean((M_lmc - B )^2) ,mean((M_lmc - B )^2)/mean(B^2) )
c(mean((X%*% (ols - B) %*%Z )^2), mean((ols - B )^2) ,mean((ols - B )^2)/mean(B^2) )
c(mean((X%*% (Bm - B) %*%Z )^2), mean((Bm - B )^2) ,mean((Bm - B )^2)/mean(B^2) )
c(mean((X%*% (Mhat - B) %*%Z )^2),mean((Mhat - B )^2) ,mean((Mhat - B )^2)/mean(B^2) )

a/Iters


