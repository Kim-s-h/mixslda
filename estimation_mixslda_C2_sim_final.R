#############################################
#
# MixsLDA Estimation
#
#############################################
library(MASS)           
library(gtools)
library(Rcpp)
library(RcppArmadillo)

rm(list=ls())


nDocs = 500
nW = 50
r = 1
#----------------------------------------------------------------------------------#
# Load data
#----------------------------------------------------------------------------------#

load(file=paste0('data/D', nDocs, '/L', nW,  '/data_C2_', r, '.Rdata'))

set.seed(452)

K=4
burnIn = 10000
nSamples = 40000

nIter = nSamples


#--------Cpp functions ----------------
sourceCpp('Cpp/z_step_ver2.cpp') 
sourceCpp('Cpp/theta_step_ver2_c2.cpp') 

# --- Specify hyper-parameters: Beta ------------
beta_p = array(1, dim=c(K,V))
for (v in 1:100){
  beta_p[1,v] = 2
}
for (v in 101:200){
  beta_p[2,v] = 2
}
for (v in 201:300){
  beta_p[3,v] = 2
}
for (v in 301:400){
  beta_p[4,v] = 2
}


# --- Specify hyper-parameters: Mean intercept and slope ---------
beta_0 = rep(0,K)
Sigma_0 = 20^2*diag(K)

#--- Specify hyper-parameters: Theta -----------------------------
theta_p = array(2, dim=c(nDocs,K))

#--- Specify hyper-parameters: Xi --------------------------------
nC = 2
xi0 = c(3,7)


#---- Initialization: theta ---------------------
theta.hist = array(NA, dim=c(nDocs,K,nSamples))   # storage for theta
theta0 = matrix(NA, nrow=nDocs, ncol=K)
for(d in 1:nDocs) theta0[d,] = rdirichlet(1,alpha=theta_p[d,])


#--- Initialization: phi -----------------------------------
phi.hist = array(NA, dim=c(K,V,nSamples))   # storage for phi
phi0 = tmp_phi0 = matrix(NA, nrow=K, ncol=V)
for(k in 1:K) phi0[k,] = rdirichlet(1,alpha=beta_p[k,])


#--- Initialization: Z  -----------------------

Z0 = sample(c(1:K), sum(N.D), replace=T )

#y_kv
y_kv = array(0, dim=c(K,V))
wz = as.data.frame(cbind(W[,1], Z0))
y.tmp = table(wz$V1, wz$Z0)
for(k in 1:K){
    y_kv[k,] = y.tmp[,k]
}

#x_dk
zw = as.data.frame(cbind(W[,2], Z0))
x.tmp = aggregate(zw, by=list(zw$V1), function(x) table(factor(x, levels=1:K)))
x_dk = x.tmp[,3]



#--- Initialization: Membership  ------------------
mem = c()
mem.hist = array(NA, dim=c(nDocs,nSamples))
mem[1] = 1
for(i in 2:nDocs){
    mem[i] = sample(c(1:nC), 1)
}

#--- Initialization:Beta_Coefficient  -------------
b1.hist = array(NA, dim=c(K,nSamples))
b1 = mvrnorm(1,mu=t.coef1, Sigma=50*diag(4))
b2.hist = array(NA, dim=c(K,nSamples))
b2 = mvrnorm(1,mu=t.coef2, Sigma=50*diag(4))


ve = 3
se = 60

sig1 = sig2 = ceiling(abs(rnorm(1,0,10^2)))
sig1.hist = c()

pi.hist = array(NA, dim=c(nSamples,nC))

dec = array(NA, dim=c(nDocs, nSamples))

# ---- MCMC -----------------------------

ptm = proc.time()
b=0
for(m in 1:nIter){
  #-----phi Step-------------------------------------
  post.beta = y_kv + beta_p
  for(k in 1:K) phi0[k,] = rdirichlet(1,alpha=post.beta[k,])

  

  #-----Theta Step-------------------------------------
  theta_n = matrix(NA, nrow=nDocs, ncol=K)
  for(d in 1:nDocs){theta_n[d,] = rdirichlet(1,alpha=100*theta0[d,])}
  new.theta = theta_step(sig1, sig2, b1, b2, resp, mem, theta_p, theta_n, theta0, x_dk)
  dec0 = -1*((theta0==new.theta)[,1]*1 - 1)
  theta0 = new.theta

  #---Z Step -----------------------------------------
  #Remember to update Beta_t
  zstep = Zstep(phi0, W, Z0, theta0,  y_kv, x_dk )
  Z0 = zstep[[1]]
  y_kv = zstep[[2]]
  x_dk = zstep[[3]]


  #--- Pi step ------------------------------------------
  nc = as.vector(table(factor(mem, levels=1:nC)))
  pi.est = rdirichlet(1, alpha=nc + xi0)


  #--- C(membership) step ------------------------------------------
  for(d in 2:nDocs){
    pb.mem = c()
    pb.mem[1] = pi.est[1]*dnorm(resp[d], mean=t(b1)%*%theta0[d,], sd=sqrt(sig1))
    pb.mem[2] = pi.est[2]*dnorm(resp[d], mean=t(b2)%*%theta0[d,], sd=sqrt(sig2))
    mem[d] = sample(c(1:nC),1,prob = pb.mem, replace=F)
  }

  #--- Beta Step 1 -----------------------------------------
  theta_c1 = theta0[mem==1,]
  Vb1 = solve(solve(Sigma_0) + (1/sig1) * t(theta_c1)%*%theta_c1)
  mu_b1 = Vb1%*%(solve(Sigma_0)%*%beta_0 + (1/sig1)*t(theta_c1)%*%resp[mem==1])
  b1 = mvrnorm(1, mu=mu_b1, Sigma=Vb1)

  #--- Beta Step 2 -----------------------------------------
  theta_c2 = theta0[mem==2,]
  Vb2 = solve(solve(Sigma_0) + (1/sig2) * t(theta_c2)%*%theta_c2)
  mu_b2 = Vb2%*%(solve(Sigma_0)%*%beta_0 + (1/sig2)*t(theta_c2)%*%resp[mem==2])
  b2 = mvrnorm(1, mu=mu_b2, Sigma=Vb2)
  

  #-----sig1 (=sig2) Step----------------------------------------
  ae1 = 0.5*nDocs + ve
  be1 = se + 0.5*t(resp[mem==1] - theta_c1%*%b1)%*%(resp[mem==1] - theta_c1%*%b1)+0.5*t(resp[mem==2] - theta_c2%*%b2)%*%(resp[mem==2] - theta_c2%*%b2)
  sig2 = sig1 = 1/rgamma(1,shape=ae1, rate=be1)



  b = b+1
  phi.hist[,,b] = phi0
  theta.hist[,,b] = theta0
  b1.hist[,b] = b1
  b2.hist[,b] = b2
  sig1.hist[b] = sig1
  dec[,b] = dec0
  mem.hist[,b] = mem
  pi.hist[b,] = pi.est

}


