pdf(file = "./Relation_fitness_QTL.pdf",width = 10,height = 5)

############################################################
#Eyre-Walker function PNAS 2010 (function relationship s,a)
# Relation fitness versus a QTL. Only negative s values. Stabilising selection
############################################################
z_function <- function(delta,Ne,s,eps,tau) {
  z <- delta*(4*Ne*s)^(tau) * (1+eps)
}

#number of negative mutations and Ne
nit_neg <- 10000
Ne <- 1000

#gamma disribution for negative s mutations
shape_s <- 0.2
mean_s  <- 3000
s <- rgamma(nit_neg,shape=shape_s,rate=shape_s/mean_s)

#distribution of positive and negative values in phenotype
skew.a <- 1#0.05 #to make assymetrical (directional selection). change mean etc.
delta <- sample(x=c(-1,skew.a),size=nit_neg,replace=T,prob=c(1,skew.a))

#parameters for the Eyre-Waalker function:
tau_values <- c(0,0.25,0.5,0.75,1)
tau_values <- c(0,0.25,0.5,0.75,1)
sigma_values <- c(0.25,0.75)

par(mfrow=c(2,3))
for(sigma in sigma_values) {
  eps <- rnorm(nit_neg,mean=0,sd=sigma)
  for(tau in tau_values[c(1,3,5)]) {
    z <- NULL
    for(i in 1:nit_neg) {
      z <- c(z,z_function(delta=delta[i],Ne=Ne,s=s[i],eps=eps[i],tau=tau))
    }
    plot(-s,z,main=sprintf("tau=%.3f\nsd=%.3f\ncorr(s,z)=%.3f",tau,sigma,cor(s,abs(z))))
  }
}

############################################################
#Caballero et al. Genetics 2015: Bivariate gamma distribution 
# Relation fitness versus a QTL. Only negative s values. Stabilising selection
############################################################
if(!require(simstudy)){
  install.packages("simstudy")
}
library(simstudy)

#number of negative mutations and Ne
nit_neg <- 10000
Ne <- 1000

#gamma disribution for negative s mutations 
shape_s <- 0.2
mean_s  <- 4000
#gamma disribution for phenotype values
shape_a <- 0.2
mean_a  <- 1

#distribution of positive and negative values in phenotype
skew.a <- 1#0.05 #to make assymetrical (directional selection). change mean etc.
delta <- sample(x=c(-1,skew.a),size=nit_neg,replace=T,prob=c(1,skew.a))
cor_values <- c(0,0.50,0.75,0.90,0.95,0.995)

par(mfrow=c(2,3))
for(corsz in cor_values) {
  res <- genCorGen(n=nit_neg, nvars=2, params1=c(mean_s,mean_a), params2=c(1/shape_s,1/shape_a), 
                   dist = "gamma", rho = corsz,  corstr = "cs", wide = T)
  s <- res$V1
  a <- delta * res$V2 # mean(abs(a))
  plot(-s,a,main=sprintf("rho=%.3f\ncor=%.3f\n",corsz,cor(s,abs(a))))
}

############################################################
#Eyre-Walker function PNAS 2010 (function relationship s,a)
# Relation fitness versus a QTL.
# include positive s values correlated with a QTL.
# Also include skew par to mimmic directional selection
############################################################
z_function <- function(delta,Ne,s,eps,tau) {
  z <- delta*(4*Ne*s)^(tau) * (1+eps)
}

#number and proportion of negative and positive mutations and Ne
nit_neg <- 10000
nit_pos <- nit_neg * 0.01
Ne <- 1000

#the fitness has two distributions:  positive (exponential) and  negative (gamma).
shape_s_neg <- 0.2
mean_s_neg  <- 4000
mean_s_pos  <- 1 #exponential

#the QTL has gamma distribution
shape_a_neg <- 0.2
mean_a  <- 1

#obtain s values from distributions:
s_neg <- -rgamma(nit_neg,shape=shape_s_neg,rate=shape_s_neg/mean_s_neg)
s_pos <- rexp(nit_pos,1/mean_s_pos)
s <- c(s_neg,s_pos)

#distribution of positive and negative values in phenotype
skew.a <- 1#0.05 #to make assymetrical (directional selection). change mean etc.
delta <- sample(x=c(-1,skew.a),size=nit_neg+nit_pos,replace=T,prob=c(1,skew.a))

#parameters for the Eyre-Waalker function:
tau_values <- c(0,0.25,0.5,0.75,1)
sigma_values <- c(0.25,0.75)

par(mfrow=c(2,3))
for(sigma in sigma_values) {
  eps <- rnorm(nit_neg+nit_pos,mean=0,sd=sigma)
  for(tau in tau_values[c(1,3,5)]) {
    z <- NULL
    #for negative values, do the same than before.
    for(i in 1:nit_neg) {
      z <- c(z,z_function(delta=delta[i],Ne=Ne,s=-s[i],eps=eps[i],tau=tau))
    }
    #for positive, invert the selective effect (1/s) to have higher z close to 0
    for(i in (nit_neg+1):(nit_neg+nit_pos)) {
      z <- c(z,z_function(delta=delta[i],Ne=Ne,s=1/s[i],eps=eps[i],tau=tau))
    }
    plot(s,z,main=sprintf("tau=%.3f\nsd=%.3f",tau,sigma),xlim=c(0,10)) #only positive
    #plot(s,z,main=sprintf("tau=%.3f\nsd=%.3f\ncorr(s,z)=%.3f",tau,sigma,cor(s,abs(z))))
  }
}

############################################################
#Caballero et al. Genetics 2015: Bivariate gamma distribution
# Relation fitness versus a QTL.
# include positive s values (using exponential) correlated with a QTL.
# Also include skew par to mimmic directional selection
# BIVARIATE GAMMA CORRELATED distributions:
############################################################
if(!require(simstudy)){
  install.packages("simstudy")
}
library(simstudy)

#number and proportion of negative and positive mutations and Ne
nit_neg <- 10000
nit_pos <- nit_neg * 0.01
Ne <- 1000

#the fitness has two distributions:  positive (exponential) and  negative (gamma).
shape_s_neg <- 0.2
mean_s_neg  <- 4000
mean_s_pos  <- 1 #exponential distribution

#the QTL has gamma distribution
shape_a_neg <- 0.2
mean_a  <- 1

#distribution of positive and negative values in phenotype
skew.a <- 1#0.05 #to make assymetrical (directional selection). change mean etc.
delta <- sample(x=c(-1,skew.a),size=nit_neg+nit_pos,replace=T,prob=c(1,skew.a))

#parameter to correlate s versus a:
cor_values <- c(0,0.50,0.75,0.90,0.95,0.995)

par(mfrow=c(2,3))
for(corsz in cor_values) {
  #for negative values, do the same than before.
  res.neg <- genCorGen(n=nit_neg, nvars=2, params1=c(mean_s_neg,mean_a), params2=c(1/shape_s_neg,1/shape_a_neg), 
                       dist = "gamma", rho = corsz,  corstr = "cs", wide = T)
  s <- -res.neg$V1
  a <- delta[1:nit_neg] * (res.neg$V2) # mean(abs(a))
  
  #for positive, make a negative correlation s vs a
  res.pos <- genCorGen(n=nit_pos, nvars=2, params1=c(mean_s_pos,mean_a), params2=c(mean_s_pos,mean_a), 
                       dist = "gamma", rho = -corsz,  corstr = "cs", wide = T)
  s <- c(s,res.pos$V1)
  a <- c(a,delta[(nit_neg+1):(nit_neg+nit_pos)] * (res.pos$V2)) # mean(abs(a))
  
  plot(s,a,main=sprintf("rho=%.3f",corsz),xlim=c(0,10)) #only positive
  #plot(s,a,main=sprintf("rho=%.3f\ncor=%.3f\n",corsz,cor(s,abs(a))))
}

############################################################
# include positive s values correlated with a QTL.
# Also include negative s values, correlated with QTL.
#TRy using Eyre-Walker function:
############################################################
z_function <- function(delta,Ne,s,eps,tau) {
  z <- delta*(4*Ne*s)^(tau) * (1+eps)
}

#number and proportion of negative and positive mutations and Ne
nit_neg <- 10000
nit_pos <- nit_neg * c(0.01,0.0001)
Ne <- 1000

#the fitness has two distributions:  positive (exponential) and  negative (gamma).
shape_s_neg <- 0.2
mean_s_neg  <- 4000
mean_s_pos  <- c(1,100) #exponential

#the QTL has gamma distribution
shape_a_neg <- 0.2
mean_a  <- 1

#obtain s values from distributions:
lims <- c(10,100)
par(mfrow=c(2,3))
for(ss in c(1:2)) {
  s_neg <- -rgamma(nit_neg,shape=shape_s_neg,rate=shape_s_neg/mean_s_neg)
  s_pos <- rexp(nit_pos[ss],1/mean_s_pos[ss])
  s <- c(s_neg,s_pos)
  
  #distribution of positive and negative values in phenotype
  skew.a <- 1#0.05 #to make assymetrical (directional selection). change mean etc.
  delta <- sample(x=c(-1,skew.a),size=nit_neg+nit_pos[ss],replace=T,prob=c(1,skew.a))
  
  #parameters for the Eyre-Waalker function:
  tau_values <- c(0,0.10,0.5,1)
  sigma.values <- c(0.1,0.25,0.5,0.75) 
  comb_tau_sig <- matrix(c(tau_values[1],sigma.values[4],
                           #tau_values[2],sigma.values[3],
                           tau_values[2],sigma.values[2],
                           #tau_values[3],sigma.values[3],
                           #tau_values[3],sigma.values[2],
                           tau_values[4],sigma.values[1]),
                           #nrow=6,ncol=2,byrow=T)
                           nrow=3,ncol=2,byrow=T)
  
  for(comb in 1:length(comb_tau_sig[,1])) {
    eps <- rnorm(nit_neg+nit_pos[ss],mean=0,sd=comb_tau_sig[comb,2])
    z <- NULL
    #for negative values, do the same than before.
    for(i in 1:nit_neg) {
      z <- c(z,z_function(delta=delta[i],Ne=Ne,s=-s[i],eps=eps[i],tau=comb_tau_sig[comb,1]))
    }
    #for positive, invert the selective effect (1/s) to have higher z close to 0
    for(i in (nit_neg+1):(nit_neg+nit_pos[ss])) {
      z <- c(z,z_function(delta=delta[i],Ne=Ne,s=1/s[i],eps=eps[i],tau=comb_tau_sig[comb,1]))
    }
    #plot(s,z,main=sprintf("tau=%.3f\nsd=%.3f",tau,sigma),xlim=c(0,10))
    plot(s,z,main=sprintf("s_beneficial=%.0f\ntau=%.3f; sd=%.3f\ncorr(s,z)=%.3f",mean_s_pos[ss],comb_tau_sig[comb,1],comb_tau_sig[comb,2],cor(s,abs(z))))
  }
}
############################################################
#number and proportion of negative and positive mutations and Ne
nit_neg <- 10000
nit_pos <- nit_neg * c(0.01,0.0001)
Ne <- 1000

#the fitness has two distributions:  positive (exponential) and  negative (gamma).
shape_s_neg <- 0.2
mean_s_neg  <- 4000
mean_s_pos  <- c(1,100) #exponential

#the QTL has gamma distribution
shape_a_neg <- 0.2
mean_a  <- 1

#obtain s values from distributions:
lims <- c(10,100)
par(mfrow=c(2,3))
for(ss in c(1:2)) {
  s_neg <- -rgamma(nit_neg,shape=shape_s_neg,rate=shape_s_neg/mean_s_neg)
  s_pos <- rexp(nit_pos[ss],1/mean_s_pos[ss])
  s <- c(s_neg,s_pos)
  
  #distribution of positive and negative values in phenotype
  skew.a <- 1#0.05 #to make assymetrical (directional selection). change mean etc.
  delta <- sample(x=c(-1,skew.a),size=nit_neg+nit_pos[ss],replace=T,prob=c(1,skew.a))
  
  #parameters for the Eyre-Waalker function:
  tau_values <- c(0,0.10,0.5,1)
  sigma.values <- c(0.1,0.25,0.5,0.75) 
  comb_tau_sig <- matrix(c(tau_values[1],sigma.values[4],
                           #tau_values[2],sigma.values[3],
                           tau_values[2],sigma.values[2],
                           #tau_values[3],sigma.values[3],
                           #tau_values[3],sigma.values[2],
                           tau_values[4],sigma.values[1]),
                         #nrow=6,ncol=2,byrow=T)
                         nrow=3,ncol=2,byrow=T)
  
  for(comb in 1:length(comb_tau_sig[,1])) {
    eps <- rnorm(nit_neg+nit_pos[ss],mean=0,sd=comb_tau_sig[comb,2])
    z <- NULL
    #for negative values, do the same than before.
    for(i in 1:nit_neg) {
      z <- c(z,z_function(delta=delta[i],Ne=Ne,s=-s[i],eps=eps[i],tau=comb_tau_sig[comb,1]))
    }
    #for positive, invert the selective effect (1/s) to have higher z close to 0
    for(i in (nit_neg+1):(nit_neg+nit_pos[ss])) {
      z <- c(z,z_function(delta=delta[i],Ne=Ne,s=1/s[i],eps=eps[i],tau=comb_tau_sig[comb,1]))
    }
    plot(s,z,main=sprintf("s_beneficial=%.0f\ntau=%.3f; sd=%.3f",mean_s_pos[ss],tau,sigma),xlim=c(0,lims[ss]))
    #plot(s,z,main=sprintf("tau=%.3f\nsd=%.3f\ncorr(s,z)=%.3f",comb_tau_sig[comb,1],comb_tau_sig[comb,2],cor(s,abs(z))))
  }
}
############################################################
dev.off()