############################################################
# Eyre-Walker function phenotype-genotype
# Relation fitness versus a QTL.
# include positive s values (using exponential) correlated with a QTL.
# Also include skew par to mimmic directional selection
############################################################
if(!require(simstudy)){
  install.packages("simstudy")
}
library(simstudy)

z_function <- function(delta,Ne,s,eps,tau,tau2,max_s) {
  z <- delta * (4*Ne*s)^(tau) * (1 + eps)  / (4*Ne*max_s)^(tau2)
}

n.mutations <- 1e4
Ne <- 1000 #
#the fitness has two distributions:  positive (exponential) and  negative (gamma).
shape_s_neg <- 0.2 #
mean_s_neg  <- 0.03*4*Ne #
shape_s_pos <- 1 #
#the QTL has gamma distribution
shape_a <- 1 #
mean_a  <- 1 #
#the correlation value
#h_mean
h.mean <- 0.36
#output pathname file
output.file <- "GammBiv_EW_plot_Fig1"

skew_values <- 1#c(1,0.25)
################
cor_values <- c(0,0.80,0.95)
################
tau_values <- c(0.00,0.25,0.50)
sd_values <- c(0.50,0.40,0.30)
################
prop.positive_values <- c(0,0.05,0.005)
mean_s_pos_values <- c(0,1,50)
################

pdf(sprintf("%s.pdf",output.file))
par(mfrow=c(2,3))#par(mfrow=c(2,3))

for(ps in c(3:3)) {
  #number and proportion of negative and positive mutations and Ne
  nit_neg <- round(n.mutations * (1 - prop.positive_values[ps]),0)
  nit_pos <- round(n.mutations * (    prop.positive_values[ps]),0)
  
  #gamma distribution for s
  sn <- rgamma(nit_neg,shape=shape_s_neg,rate=shape_s_neg/mean_s_neg) / (4*Ne)
  sp <- rgamma(nit_pos,shape=shape_s_pos,rate=shape_s_pos/mean_s_pos_values[ps]) / (4*Ne)
  s <- c(-sn,sp)
  
  for(skew.a in skew_values) {
    #distribution of positive and negative values in phenotype: SCENARIO 1
    delta.a <- sample(x=c(-1,skew.a),size=nit_neg+nit_pos,replace=T,prob=c(1,skew.a))
    
    for(ts in 1:3) {
    #for(tau in tau_values) {				
	    #for(sigma in sd_values) {
	      #epsilon value from a normal
	      eps <- rnorm(nit_neg+nit_pos,mean=0,sd=sd_values[ts])#sigma)
	      
	      #for negative values, do a positive correlation s vs a
	      zn <- NULL
	      for(i in 1:nit_neg) {
	        zn <- c(zn,z_function(delta=delta.a[i],Ne=Ne,s=sn[i],eps=eps[i],tau=tau_values[ts],tau2=tau_values[ts],max_s=max(abs(s))))
	      }
	      #for positive, make a negative correlation s vs a
	      if(prop.positive_values[ps]>0) {
	        zp <- NULL
	        for(i in 1:nit_pos) {
	          zp <- c(zp,z_function(delta=delta.a[nit_neg+i],Ne=Ne,s=(sp[i]),eps=eps[nit_neg+i],tau=-(tau_values[ts]),tau2=(tau_values[ts]),max_s=max(abs(s))))
	        }
          z <- c(zn,zp)
	        
		      #dominance
	        alpha_n <- shape_s_neg / -mean_s_neg
	        alpha_p <- shape_s_pos / mean_s_pos_values[ps]			  
	        Kn = alpha_n * ((2.0*h.mean)^((-1.0)/shape_s_neg)-1.0);
	        Kp = alpha_p * ((2.0*h.mean)^((-1.0)/shape_s_pos)-1.0);			  
	        hn <- runif(n=length(sn),0,exp(-sn*Kn))
	        hp <- runif(n=length(sp),0,exp(-sp*Kp))			  
	        h <- c(hn,hp)
	        #s <- c(-sn,sp)
	      } else {
	        z <- zn
	        #dominance
	        alpha_n <- shape_s_neg / -mean_s_neg
	        Kn = alpha_n * ((2.0*h.mean)^((-1.0)/shape_s_neg)-1.0);
	        hn <- runif(n=length(sn),0,exp(-sn*Kn))
	        h <- hn
	        #s <- -sn
	      }
	      
	      #plot
	      s[s < -1] <- -1
	      if(skew.a==1) {xlup <- 10} 
	      else {xlup <- 2}
	      plot(s,z,pch=20,main=sprintf("tau=%.3f sd=%.3f\n cor(sel,phe): %.2f",tau_values[ts],sd_values[ts],cor(abs(s),abs(z),method="spearman")),xlab="s",ylab="Phen.Eff",xlim=c(-1,0.1))#,ylim=c(-10,xlup))
	     	#plot(s/(4*Ne),a,pch=20,main=sprintf("rho: %.2f skew: %.2f",corsz,skew.a),xlab="Sel.Coeff",ylab="Phen.Eff",ylim=c(-10,xlup),xlim=c(-1,0.1))#,xlim=c(-1.5e5,1e3)) #4Ns(+): %.0f perc(+): %.0e \n ,mean_s_pos_values[ps],prop.positive_values[ps]
	      abline(v = 0)
	      #plot(s,a,pch=20,main=sprintf("4Ns(+): %.0f perc(+): %.0e \n rho: %.2f skew: %.2f",mean_s_pos_values[ps],prop.positive_values[ps],corsz,skew.a),xlab="Sel.Coeff",ylab="Phen.Eff",xlim=c(-10,max(s)))
	      #abline(v = 0)
	      #plot(s,h,pch=20,main=sprintf("Relation s versus dominance.\n mean_sn: %.3e mean_sp: %.3e \n kn: %.3e kp: %.3e",mean_s_neg,mean_s_pos,Kn,Kp),xlab="Sel.Coeff",ylab="dominance",xlim=c(quantile(sn,probs=0.1),c(quantile(sp,probs=0.9)))) 
	    #}
	  }
	}
}
mean_s_neg  <- -mean_s_neg #
for(skew.a in skew_values) {
  for(ps in c(3:3)) {
    for(corsz in cor_values) {				
      #number and proportion of negative and positive mutations and Ne
      nit_neg <- round(n.mutations * (1 - prop.positive_values[ps]),0)
      nit_pos <- round(n.mutations * (    prop.positive_values[ps]),0)
      
      #distribution of positive and negative values in phenotype: SCENARIO 1
      delta.a <- sample(x=c(-1,skew.a),size=nit_neg+nit_pos,replace=T,prob=c(1,skew.a))
      
      #for negative values, do a positive correlation s vs a
      res.neg <- genCorGen(n=nit_neg, nvars=2, params1=c(-mean_s_neg,mean_a), params2=c(1/shape_s_neg,1/shape_a), 
                           dist = "gamma", rho = corsz,  corstr = "cs", wide = T)
      sn <- -res.neg$V1
      a <- delta.a[1:nit_neg] * (res.neg$V2) # mean(abs(a))  
      thresh.a0 <- c(min(a[sn/(4*Ne)>-0.01]),max(a[sn/(4*Ne)>-0.01]))  #max values for a in positive s
      show(thresh.a0)
      
      #for positive, make a negative correlation s vs a
      if(prop.positive_values[ps]>0) {
        sp <- NULL
        ap <- NULL
        while(length(sp)<nit_pos) {
          res.pos <- genCorGen(n=nit_pos, nvars=2, params1=c(mean_s_pos_values[ps],mean_a), params2=c(1/shape_s_pos,1/shape_a), 
                               dist = "gamma", rho = -corsz,  corstr = "cs", wide = T)
          sp <- c(sp,res.pos$V1[res.pos$V1-mean_s_pos_values[ps]>0 & res.pos$V2>thresh.a0[1] & res.pos$V2<thresh.a0[2]])
          ap <- c(ap,res.pos$V2[res.pos$V1-mean_s_pos_values[ps]>0 & res.pos$V2>thresh.a0[1] & res.pos$V2<thresh.a0[2]])
        }
        min.sp <- min(sp)
        sp <- sp[1:nit_pos]#-min.sp
        ap <- ap[1:nit_pos]
        # mean(abs(a))
        #sp <- res.pos$V1
        #ap <- res.pos$V2
        s <- c(sn,sp)
        a <- c(a,delta.a[(nit_neg+1):(nit_neg+nit_pos)] * ap)
        #dominance
        alpha_n <- shape_s_neg / -mean_s_neg
        alpha_p <- shape_s_pos / mean_s_pos_values[ps]			  
        Kn = alpha_n * ((2.0*h.mean)^((-1.0)/shape_s_neg)-1.0);
        Kp = alpha_p * ((2.0*h.mean)^((-1.0)/shape_s_pos)-1.0);			  
        hn <- runif(n=length(sn),0,exp( sn*Kn))
        hp <- runif(n=length(sp),0,exp(-sp*Kp))			  
        h <- c(hn,hp)
      } else {
        s <- sn
        alpha_n <- shape_s_neg / -mean_s_neg
        Kn = alpha_n * ((2.0*h.mean)^((-1.0)/shape_s_neg)-1.0);
        hn <- runif(n=length(sn),0,exp( sn*Kn))
        h <- hn
      }
      
      #plot
      s[s/(4*Ne) < -1] <- -4*Ne
      if(skew.a==1) {xlup <- 10} 
      else {xlup <- 2}
      plot(s/(4*Ne),a,pch=20,main=sprintf("rho: %.2f \n cor(sel,phe): %.2f",corsz,cor(abs(s),abs(a),method="spearman")),xlab="Sel.Coeff",ylab="Phen.Eff",xlim=c(-1,0.1))#,ylim=c(-10,xlup)#,xlim=c(-1.5e5,1e3)) #4Ns(+): %.0f perc(+): %.0e \n ,mean_s_pos_values[ps],prop.positive_values[ps]
      abline(v = 0)
      #plot(s,a,pch=20,main=sprintf("4Ns(+): %.0f perc(+): %.0e \n rho: %.2f skew: %.2f",mean_s_pos_values[ps],prop.positive_values[ps],corsz,skew.a),xlab="Sel.Coeff",ylab="Phen.Eff",xlim=c(-10,max(s)))
      #abline(v = 0)
      #plot(s,h,pch=20,main=sprintf("Relation s versus dominance.\n mean_sn: %.3e mean_sp: %.3e \n kn: %.3e kp: %.3e",mean_s_neg,mean_s_pos,Kn,Kp),xlab="Sel.Coeff",ylab="dominance",xlim=c(quantile(sn,probs=0.1),c(quantile(sp,probs=0.9)))) 
    }
  }
}
dev.off()

