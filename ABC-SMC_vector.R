##Parameter estimation
transmission_model <- function(N, Gm, Gf, V, rho, r, alpha_go, beta_go, alpha_os, beta_os){
  record <- matrix(0, N, 7)
  M <- rbinom(N,Gm,V/1000)
  F <- rbinom(N,Gf,V/1000)
  n <- vector()
  t <- vector()
  for (i in 1:N) {
    n[i] <- rbinom(1,round(min(rho*M[i],F[i]),0),r)
    p <- runif(1,0,1)
    t[i] <- uniroot.all(function(x){1-p-exp(-0.0043140*x*exp(0.0874910*x))},c(0,100))
  }
  record[,c(1,4,5,6)] <- c(n,t,M,F)
  for (i in 1:N) {
    if (n[i]!=0){
      oocy <- rgamma(n[i],alpha_go,beta_go)
      sporo <- rgamma(n[i],alpha_os,beta_os)
      p <- runif(1,0,1)
      t <- uniroot.all(function(x){1-p-exp(-0.0043140*x*exp(0.0874910*x+0.3604109))},c(0,100))
      record[i,2:4] <- c(min(oocy),min(oocy+sporo),t)
      record[i,7] <- max(oocy+sporo)
    }
  }
  return(record)
}


oocyst_prev <- function(N, Gm, Gf, V, rho, r, alpha_go, beta_go, alpha_os, beta_os){
  a <- transmission_model(N, Gm, Gf, V, rho, r, alpha_go, beta_go, alpha_os, beta_os)
  survival <- matrix(0,N,240)
  for (i in 1:N) {
    if (ceiling(a[i,4]*24)<240){
      survival[i,(ceiling(a[i,4])*24):240] <- 1
    }
  }
  oocyst <- matrix(0,N,240)
  for (i in 1:N) {
    if (ceiling(a[i,2]*24)<240 & a[i,2]!=0 ){
      oocyst[i, (ceiling(a[i,2]*24)):(ceiling(a[i,7]*24))] <- 1
    }
  }
  prev_oocyst <- vector()
  for (i in 1:240) {
    prev_oocyst[i] <- length(which(oocyst[,i]==1 & survival[,i]==0))/length(which(survival[,i]==0))
  }
  return(prev_oocyst)
}


library(foreach)
library(doParallel)
library(mvtnorm)
cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)
parallel::stopCluster(cl)


time.step=20
n=1000
d_total <- matrix(0,time.step,n)
weights <- matrix(1/n,time.step,n) 
para_rec <- matrix(0,n,3)
times <- matrix(0,time.step,n)
para <- list()
epsilon <- rep(0,time.step)
epsilon[1] <- 0.3699087

for (t in 1:14){
  if(t==1){
    rec <- matrix(0,n,4)
    result <- foreach(i=1:n,.packages=c("tmvtnorm","rootSolve")) %dopar% {
      d <- epsilon[t]+1
      time <- 0
      while(d > epsilon[t]){
        set.seed(n*time+i)
        proposed.mu <- c(runif(1,0,8),runif(1,0,1))
        #simulation <- try(sto_model(proposed.mu[1]*4000,proposed.mu[2],proposed.mu[3],proposed.mu[4],kd(pk101,proposed.mu[5],proposed.mu[6],proposed.mu[7]),proposed.mu[8],proposed.mu[9],proposed.mu[10],proposed.mu[11],proposed.mu[12],42,25,672),silent=T)
        oocyst_record <- vector()
        for (j in 1:33){
          simulation <- oocyst_prev(30,round(as.numeric(data[j,2]),0),round(as.numeric(data[j,3]),0),0.87,proposed.mu[1],proposed.mu[2],13.6, 2.7, 17.8, 2.2)
          oocyst_record[j] <- simulation[85]
        }
        time <- time+1
        d <- sqrt(sum((oocyst_record-as.numeric(data[,4])/as.numeric(data[,5]))^2))
        }
      c(proposed.mu,d,time)
    }
    for (i in 1:n){
      rec[i,] <- result[[i]]
    }
    para[[t]] <- rec[,1:2]
    times[t,] <- rec[,4]
    d_total[t,] <- rec[,3]
    
  }else{
    rec <- matrix(0,n,5)
    epsilon[t] <- as.numeric(quantile(d_total[t-1,],0.5))
    
    mean.prev <- colSums(para[[t-1]]*weights[t-1,])
    
    index <- which(d_total[t-1,]<as.numeric(quantile(d_total[t-1,],0.5)))
    notin_para <- para[[t-1]][-index,]
    notin_weights <- weights[t-1,-index]
    hat_para <- para[[t-1]][index,]
    hat_weights <- weights[t-1,index]/sum(weights[t-1,index])
    
    cov_matr <- matrix(0,2,2)
    for (k in 1:dim(hat_para)[1]){
      for (l in 1:length(notin_weights)){
        A <- as.matrix(hat_para[k,],1,2)
        B <- as.matrix(notin_para[l,],1,2)
        cov_matr <- cov_matr + hat_weights[k]*notin_weights[l]*(A-B)%*%t(A-B)
      }
    }
    
    para_rec <- matrix(0,n,2)
    result <- foreach (i=1:n,.packages=c("truncnorm","tmvtnorm","rootSolve")) %dopar% {
      d <- epsilon[t]+1
      time <- 0
      while(d > epsilon[t]){
        sample.particle <- sample(n,1,prob=weights[t-1,])
        proposed.mu0 <- para[[t-1]][sample.particle,]
        set.seed(n*time+i)
        proposed.mu <- as.numeric(rtmvnorm(1,proposed.mu0,cov_matr,c(0,0) ,c(8,1)))
        oocyst_record <- vector()
        for (j in 1:33){
          simulation <- oocyst_prev(30,round(as.numeric(data[j,2]),0),round(as.numeric(data[j,3]),0),0.87,proposed.mu[1],proposed.mu[2],13.7, 2.7, 17.8, 2.2)
          oocyst_record[j] <- simulation[85]
        }
        time <- time+1
        d <- sqrt(sum((oocyst_record-as.numeric(data[,4])/as.numeric(data[,5]))^2))
        #d <- sqrt(sum((oocyst_record-as.numeric(data[,4])/as.numeric(data[,5]))^2))
      }
      mu.weights.denominator=0
      for (k in 1:n){
        mu.weights.denominator <- mu.weights.denominator+weights[t-1,k]*dtmvnorm(proposed.mu,as.numeric(para[[t-1]][k,]),cov_matr,c(0,0) ,c(8,1))
      }
      mu.weights.numerator = prod(dunif(proposed.mu,c(0,0) ,c(8,1)))
      weights_current = mu.weights.numerator/mu.weights.denominator
      c(proposed.mu,d,time,weights_current)
    }
    for (i in 1:n){
      rec[i,] <- result[[i]]
    }
    para[[t]] <- rec[,1:2]
    times[t,] <- rec[,4]
    d_total[t,] <- rec[,3]
    weights[t,] <- rec[,5]
    weights[t,] <- weights[t,]/sum(weights[t,])
  }
}
result <- as.data.frame(para[[14]])
names(result) <- c("rho","p_f")

###Figure 4
hist(para[[14]][,1],breaks=50,main="",xlab=expression(rho),cex.axis=2,cex.lab=2,yaxt="n")
axis(2,at=c(0,25,50),cex.axis=2)
abline(v=median(para[[14]][,1]),col="red")
legend("topright",lty=1,col="red",legend="median",cex=2,bty="n")
hist(para[[14]][,2],breaks=55,main="",xlab=expression(p[r]),ylim=c(0,50),cex.axis=2,cex.lab=2,yaxt="n")
axis(2,at=c(0,25,50),cex.axis=2)
abline(v=median(para[[14]][,2]),col="red")
legend("topright",lty=1,col="red",legend="median",cex=2,bty="n")

#correlation
plot(para[[14]][,1],para[[14]][,2],xlab=expression(rho),ylab=expression(p[r]),main="",,cex.axis=1.5,cex.lab=1.5)


