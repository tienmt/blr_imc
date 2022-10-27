library(Rcpp); library(MASS) ;library(softImpute)
sourceCpp("/Users/thetm/Dropbox/ongoing_works/B RRR/BRRR with missing Y/Rcodes/testabcd.cpp")
P_Omega = function(a,entri){ a[entri] = 0; return(a) }; P_Omega <- compiler::cmpfun(P_Omega)

# random data generation Y = XBZ , (nl) = (np)*(pk)*(kl)
n = 100    # samples in Y, X
l = 10    # response in Y,Z
k = 20    # latent feature samples in Z
p = 10    # predictors in X
r = 2     # true rank
sig2 = 1

X = matrix(rnorm(n*p), ncol =p) ; tX = t(X)
Z = matrix(rnorm(k*l), ncol =l) ; tZ = t(Z)
B = matrix(rnorm(p*r),nr=p) %*% t(matrix(rnorm(k*r),nr=k )) * 2 + 0.1*rnorm(p*k)
Y = X%*%B%*%Z  + rnorm(n*l)

missrate = 0.3
imiss = sample(n*l,n*l*missrate )

Y_Omega = P_Omega(Y,imiss)
Yna = Y
Yna[imiss] <-NA

fit1 = softImpute(Yna,maxit = 1000,type = 'als')
yimp = complete(Yna,fit1); mean(( (yimp - Y)[imiss])^2) 

# 2 cases of ols fail: when k >l l and when n< p
XX = t(X)%*%X       # p x p 
XYZ_im = t(X) %*% yimp %*%t(Z) 
XYZ_Omega = t(X) %*% Y_Omega %*%t(Z) 
ZZ = Z%*%t(Z)       # k x k

ols_imp = ginv(XX) %*% XYZ_im %*% ginv(ZZ) ; mean((ols_imp - B)^2) 
SVD = svd(ols_imp) ;rank_selected = 3 ;SVD$d [ rank_selected : min(p,k) ] = 0 ; 
Mhat = SVD$u %*% diag(SVD$d) %*% t(SVD$v);  mean((Mhat - B)^2)

### MALA
lam = 1  # in the prior
ystar = diag(p)*lam^2
Iters = 10000
burnin = 1000

Bm = matrix(data=0,nr=p,nc=k)
M = ols_imp
h = 1/(p*k)^1.76
a = 0  ; start_time <- Sys.time()
for(s in 1:Iters){
  MtM = tcrossprod(M)
  tam1 = solve(ystar + MtM, M )
  hxtyxm = h*mymatmul(tX ,P_Omega(Y-mymatmul(X,M,Z),imiss), tZ)
  tam = M + hxtyxm/sig2 + h*(p+k+2)*tam1 + sqrt(2*h)*rnorm(p*k)
  
  tattam = tcrossprod(tam)
  pro.tam = -sum( P_Omega(Y-mymatmul(X,tam,Z),imiss)^2)/(2*sig2) - 0.5*(p+k+2)*determinant(ystar + tattam)[[1]][1]
  pro.M = -sum(P_Omega(Y-mymatmul(X,M,Z),imiss)^2) /(2*sig2) - 0.5*(p+k+2)*determinant(ystar + MtM)[[1]][1]
  
  tam2 = solve(ystar + tattam,tam)
  tran.m = -sum((M-tam-h*mymatmul(tX ,P_Omega(Y-mymatmul(X,tam,Z),imiss), tZ)/sig2 - h*(p+k+2)*tam2 )^2)/(4*h)
  tran.tam = -sum((tam-M- hxtyxm/sig2 - h*(p+k+2)*tam1 )^2)/(4*h)
  
  pro.trans = pro.tam+tran.m-pro.M-tran.tam
  if(log(runif(1)) <= pro.trans){
    M = tam;  print(s);  a = a+1
  } 
  if (s>burnin)Bm = Bm + M/(Iters-burnin)
}
end_time <- Sys.time();  end_time - start_time; a/Iters

#Langevin MC for BRRR
M1 = ols_imp  # initialize
h = h/n^4
M_lmc = matrix(data=0, nr=p,nc=k)
for(s in 1:Iters){
  tam = solve(ystar + tcrossprod(M1),M1)
  M1 = M1 + h*mymatmul(tX ,P_Omega(Y-mymatmul(X,tam,Z),imiss), tZ)/sig2 + h*(p+l+2)*tam + sqrt(2*h)*rnorm(p*k)
  if(s>burnin) M_lmc =  M_lmc + M1/(Iters-burnin)
}

c(mean((X%*% (M_lmc - B) %*%Z )^2), mean((M_lmc - B )^2) ,mean((M_lmc - B )^2)/mean(B^2) ,mean(abs( (X%*% M_lmc %*%Z -Y)[imiss] ) ) )
c(mean((X%*% (ols_imp - B) %*%Z )^2), mean((ols_imp - B )^2) ,mean((ols_imp - B )^2)/mean(B^2),mean(abs( (X%*% ols_imp %*%Z -Y)[imiss] ) ) )
c(mean((X%*% (Bm - B) %*%Z )^2), mean((Bm - B )^2) ,mean((Bm - B )^2)/mean(B^2), mean(abs( (X%*% Bm %*%Z -Y)[imiss] )  ))
c(mean((X%*% (Mhat - B) %*%Z )^2),mean((Mhat - B )^2) ,mean((Mhat - B )^2)/mean(B^2), mean(abs( (X%*% Mhat %*%Z -Y)[imiss] )  ))
a/Iters



