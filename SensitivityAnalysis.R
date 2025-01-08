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

sporozoite_prev <- function(N, Gm, Gf, V, rho, r, alpha_go, beta_go, alpha_os, beta_os){
  a <- transmission_model(N, Gm, Gf, V, rho, r, alpha_go, beta_go, alpha_os, beta_os)
  survival <- matrix(0,N,720)
  for (i in 1:N) {
    if (ceiling(a[i,4]*24)<720){
      survival[i,(ceiling(a[i,4])*24):720] <- 1
    }
  }
  sporozoite <- matrix(0,N,720)
  for (i in 1:N) {
    if (ceiling(a[i,3]*24)<720 & a[i,3]!=0){
      sporozoite[i, ceiling(a[i,3]*24):720] <- 1
    }
  }
  prev_sporozoite <- vector()
  for (i in 1:720) {
    prev_sporozoite[i] <- length(which(sporozoite[,i]==1 & survival[,i]==0))/length(which(survival[,i]==0))
  }
  return(prev_sporozoite)
}
library(truncnorm)
library(foreach)
library(doParallel)
cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)
parallel::stopCluster(cl)

####Model overview (Figure 3)
benchmark <- c(0.08,7,5.5,71,0.006,0.03,0.0004,0.007,0.02,0.2,0.87,2.4,0.085,13.6, 2.7, 17.8, 2.2)
game_male_record <- matrix(0,100,480)
game_female_record <- matrix(0,100,480)
total_record <- matrix(0,100,480)
for (i in 1:100) {
  simulation <- try(sto_model(benchmark[1]*4000,benchmark[2],benchmark[3],rep(benchmark[4],2880),matrix(0,2880,3),benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],al=42,as=25,480),silent = TRUE)
  total_record[i,] <- ifelse(rowSums(simulation[,c(1:25,60:83,88,89)])==0,1,rowSums(simulation[,c(1:25,60:83,88,89)]))/4000
  game_male_record[i,] <- as.numeric(simulation[,88])/4000
  game_female_record[i,] <- as.numeric(simulation[,89])/4000
}
sporozoite_all <- list()
sporozoite_record <- matrix(0,480,720)

for (j in 1:100) {
  result <- foreach(i=1:480,.packages=c("rootSolve")) %dopar% {
    set.seed(i)
    sporozoite_record <- sporozoite_prev(1000,round(game_male_record[j,i],0), round(game_female_record[j,i],0),0.87,2.42,0.2,13.6, 2.7, 17.8, 2.2)
    sporozoite_record
  }
  for (i in 1:480) {
    sporozoite_record[i,] <- result[[i]]
  }
  sporozoite_all[[j]] <- sporozoite_record
  print(j)
}

prob_sporozoite <- matrix(0,100,480)
for (i in 1:100) {
  prob_sporozoite[i,] <- sporozoite_all[[i]][,15*24]
}


par(mar = c(5.1, 5.1, 4.1, 1.1))
total_median <- vector()
total_upper <- vector()
total_lower <- vector()
male_median <- vector()
male_upper <- vector()
male_lower <- vector()
female_median <- vector()
female_upper <- vector()
female_lower <- vector()
for (i in 1:480) {
  total_median[i] <- median(total_record[,i])
  total_upper[i] <- quantile(total_record[,i],0.975)
  total_lower[i] <- quantile(total_record[,i],0.025)
  male_median[i] <- median(game_male_record[,i])
  male_upper[i] <- quantile(game_male_record[,i],0.975)
  male_lower[i] <- quantile(game_male_record[,i],0.025)
  female_median[i] <- median(game_female_record[,i])
  female_upper[i] <- quantile(game_female_record[,i],0.975)
  female_lower[i] <- quantile(game_female_record[,i],0.025)
}
top95 <- vector()
bot95 <- vector()
sporo_median <- vector()
for (i in 1:480) {
  top95[i] <- quantile(prob_sporozoite[,i],0.975)
  bot95[i] <- quantile(prob_sporozoite[,i],0.025)
  sporo_median[i] <- quantile(prob_sporozoite[,i],0.5)
}
par(mar = c(5.1, 5.1, 4.1, 1.1))
plot(1:480/24,total_median,type="n",log = 'y',xlab = "Time post infection (Days)",ylab="Parasitemia/ml",ylim=c(1,1e+13),xlim = c(10,20),cex.axis=2,cex.lab=1.8,cex.main=2,yaxt="n")
axis(2,at=c(1,10^4,10^8,10^12) ,lab=c(expression(1),expression(10^4),expression(10^8),expression(10^12)),cex.axis=2)
#title(ylab= "", line = 4,cex.lab=2)
polygon(c(1:480/24, rev(1:480/24)), c(total_lower, rev(total_upper)),
        col = rgb(192, 192, 192, maxColorValue=255, alpha=150), border = NA)  

polygon(c(1:480/24, rev(1:480/24)), c(female_lower, rev(female_upper)),
        col = rgb(255, 0, 0, maxColorValue=255, alpha=100), border = NA)

polygon(c(1:480/24, rev(1:480/24)), c(male_lower, rev(male_upper)),
        col = rgb(0, 0, 255, maxColorValue=255, alpha=100), border = NA)
lines(1:480/24,total_median)
lines(1:480/24,male_median,col="blue")
lines(1:480/24,female_median,col="red")
legend("topleft",col=c("black","red","blue"),legend=c("Total parasitemia /ml","Female gametocytemia /ml","Male gametocytemia /ml"),lty=1,bty="n",cex=1.5)


plot(1:480/24,sporo_median,type="l",xlim=c(10,20),ylim=c(0,1),xlab = "Time post infection (Days)",ylab = expression(p[i]),cex.axis=2,cex.lab=2,cex.main=2,yaxt="n")
axis(2,at=c(0,0.25,0.5,0.75,1),cex.axis=2)

polygon(c(1:480/24, rev(1:480/24)), c(top95, rev(bot95)),
        col = rgb(128, 128, 128, maxColorValue=255, alpha=100), border = NA)  






#PRCC
library(lhs)
library(rootSolve)
library(doParallel)
set.seed(1)
A <- randomLHS(1000, 17)
lhs_para <- matrix(0,1000,17)
### +-50% of benchmark values
### +-75%: change 0.5 and 1.5 to 1.25 and 1.75
for (i in 1:17) {
  lhs_para[,i] <- qunif(A[,i], min = 0.5*benchmark[i], max = 1.5*benchmark[i])
}

game_male_record <- matrix(0,1000,720)
game_female_record <- matrix(0,1000,720)
sporozoite_record1 <- matrix(0,1000,720)
for (i in 1:1000) {
    simulation <- try(sto_model((lhs_para[i,1])*4000,lhs_para[i,2],lhs_para[i,3],rep(lhs_para[i,4],2880),matrix(0,2880,3),lhs_para[i,5],lhs_para[i,6],lhs_para[i,7],lhs_para[i,8],lhs_para[i,8],lhs_para[i,9],lhs_para[i,10],al=42,as=25,720),silent = TRUE)
    game_male_record[i,] <- as.numeric(simulation[,88])/4000
    game_female_record[i,] <- as.numeric(simulation[,89])/4000
    for (j in 1:720){
      if (game_male_record[i,j]>10000&game_female_record[i,j]>50000){
        sporozoite_inter <- 1
      }else{
        if(game_male_record[i,j]<500&game_female_record[i,j]<1000){
          sporozoite_inter <- 0
        }else{
          sporozoite <- sporozoite_prev(1000,round(game_male_record[i,j],0), round(game_female_record[i,j],0),lhs_para[i,11],lhs_para[i,12],lhs_para[i,13],lhs_para[i,14],lhs_para[i,15],lhs_para[i,16],lhs_para[i,17])
          sporozoite_inter <- sporozoite[15*24]
        }
      }
      sporozoite_record[i,j] <- sporozoite_inter
    }
  print(i)
}

days_Ip <- vector() ##calculate the time when the transmission probability reaches 0.5
for (i in 1:1000) {
  days_Ip[i] <- min(which(sporozoite_record[i,]>0.5))
}

dataset <- cbind(lhs_para,days_Ip)

rank_trans_data <- matrix(0,1000,18)
for (i in 1:18){
  rank_trans_data[,i] <- rank(dataset[,i])
}
rank_trans_data <- as.data.frame(rank_trans_data)
names(rank_trans_data) <- c(expression(P[init]),expression(mu),expression(sigma),expression(r[p]),expression(f),expression(delta[p]),expression(delta[g]),expression(delta[gc]),expression(m),expression(f[g]),expression(V),expression(rho),expression(p[f]),
                            expression(alpha[go]),expression(beta[go]),expression(alpha[os]),expression(beta[os]),"Ip")
                     

prcc_result <- vector()
for (i in 1:17){
  lm_p <- lm(Ip~. , data = rank_trans_data[,-c(i)])
  lm_para <- lm(rank_trans_data[,i]~. ,data = rank_trans_data[,-c(i,18)])
  prcc_result[i] <- cor(lm_p$residuals,lm_para$residuals)
}
###Figure5 and S1
col_names <- c(expression(P[init]),expression(mu),expression(sigma),expression(r[p]),expression(f),expression(delta[p]),expression(delta[g]),expression(delta[gc]),expression(m),expression(f[g]),expression(V),expression(rho),expression(p[f]),
               expression(alpha[zo]),expression(beta[zo]),expression(alpha[os]),expression(beta[os]),"Ip")

prcc_result_order <- prcc_result[c(10+order(prcc_result[11:17]),order(prcc_result[1:10])), drop=TRUE]
col_names_order <- col_names[c(10+order(prcc_result[11:17]),order(prcc_result[1:10])), drop=TRUE]
par(mar = c(4, 3.5, 2, 1.1))
dotchart(prcc_result_order, labels=col_names_order,pch=19, xlim=c(-1,1),
         xlab = "PRCC values",cex=1.5)
abline(v=0,lty=2)


###sensitivity analysis
cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)
parallel::stopCluster(cl)
###Here we give an example for testing r_p
###Change parameter value e.g., change benchmark[4] to test how r_p affects the result
game_male_record <- matrix(0,100,480)
game_female_record <- matrix(0,100,480)
for (i in 1:100) {
  simulation <- try(sto_model((benchmark[1])*4000,benchmark[2],benchmark[3],rep(benchmark[4],2880),matrix(0,2880,3),benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],al=42,as=25,480),silent = TRUE)
  game_male_record[i,] <- as.numeric(simulation[,88])/4000
  game_female_record[i,] <- as.numeric(simulation[,89])/4000
}

sporozoite_all <- list()
sporozoite_record <- matrix(0,480,720)

for (j in 1:100) {
  result <- foreach(i=1:480,.packages=c("rootSolve")) %dopar% {
    set.seed(i)
    sporozoite_record <- sporozoite_prev(1000,round(game_male_record[j,i],0), round(game_female_record[j,i],0),benchmark[11],benchmark[12],benchmark[13],benchmark[14], benchmark[15], benchmark[16], benchmark[17])
    sporozoite_record
  }
  for (i in 1:480) {
    sporozoite_record[i,] <- result[[i]]
  }
  sporozoite_all[[j]] <- sporozoite_record
  print(j)
}

prob_sporozoite <- matrix(0,100,480)
for (i in 1:100) {
  prob_sporozoite[i,] <- sporozoite_all[[i]][,15*24]
}

time_0.5 <- vector()
for (i in 1:100) {
  time_0.5[i] <- min(which(prob_sporozoite[i,]>=0.5))
}

prob_median <- vector()
prob_upper <- vector()
prob_lower <- vector()

for (i in 1:480) {
  prob_median[i] <- median(prob_sporozoite[,i])
  prob_upper[i] <- quantile(prob_sporozoite[,i],0.975)
  prob_lower[i] <- quantile(prob_sporozoite[,i],0.025)
}

prob_f <- matrix(0,3,480)
prob_f[1,] <- prob_median
prob_f[2,] <- prob_upper
prob_f[3,] <- prob_lower

write.csv(time_0.5,"time_rp30.csv",row.names = FALSE) ##save the results
write.csv(prob_f,"prob_rp30.csv",row.names = FALSE)

###Figure 6 and 7 (one example r_p)
par(mar = c(5.1, 5.1, 4.1, 1.1))
prob_rp30 <-read.csv("prob_rp30.csv")
prob_rp50 <-read.csv("prob_rp50.csv")
prob_rp70 <-read.csv("prob_rp70.csv")
prob_rp90 <-read.csv("prob_rp90.csv")
time_rp30 <-read.csv("time_rp30.csv")$x
time_rp50 <-read.csv("time_rp50.csv")$x
time_rp70 <-read.csv("time_rp70.csv")$x
time_rp90 <-read.csv("time_rp90.csv")$x
plot(1:576/24,prob_rp30[1,],type="l",ylim=c(0,1),xlab = "Time post infection (Days)",ylab = expression(p[i]),col = rgb(255, 153, 204, maxColorValue=255),xlim=c(10,23),cex.axis=2,cex.lab=2,cex.main=2,yaxt="n")
axis(2,at=c(0,0.25,0.5,0.75,1),cex.axis=2)
polygon(c(1:576/24, rev(1:576/24)), c(prob_rp30[2,], rev(prob_rp30[3,])),
        col = rgb(255, 153, 204, maxColorValue=255, alpha=150), border = NA)  

lines(1:576/24,prob_rp50[1,],col = rgb(255,51,153, maxColorValue=255))
polygon(c(1:576/24, rev(1:576/24)), c(prob_rp50[2,], rev(prob_rp50[3,])),
        col = rgb(255,51,153, maxColorValue=255, alpha=150), border = NA)

lines(1:576/24,prob_rp70[1,],col = rgb(204, 0, 102, maxColorValue=255))
polygon(c(1:576/24, rev(1:576/24)), c(prob_rp70[2,], rev(prob_rp70[3,])),
        col = rgb(204, 0, 102, maxColorValue=255, alpha=150), border = NA)
lines(1:576/24,prob_rp90[1,],col = rgb(102,0,51, maxColorValue=255))
polygon(c(1:576/24, rev(1:576/24)), c(prob_rp90[2,], rev(prob_rp90[3,])),
        col = rgb(102,0,51, maxColorValue=255, alpha=150), border = NA)
legend("topleft",col=c(rgb(102,0,51, maxColorValue=255,alpha=150),rgb(204, 0, 102, maxColorValue=255,alpha=150),rgb(255,51,153, maxColorValue=255,alpha=150),rgb(255, 153, 204, maxColorValue=255,alpha=150)),legend=c(expression(r[p] == 90 (PMF==25.4)),expression(r[p] == 70 (PMF==19.7)),expression(r[p] == 50 (PMF==14.1)),expression(r[p] == 30 (PMF==8.5))),lty=1,lwd=5,bty="n")


plot(c(30,50,70,90),c(median(time_rp30)/24,median(time_rp50)/24,median(time_rp70)/24,median(time_rp90)/24),pch=16,xlab=expression(r[p]),ylab=expression(t[p]),cex.axis=2,cex.lab=1.8,cex.main=2,ylim=c(15,21))
lines(c(30,50,70,90),c(quantile(time_rp30,0.975)/24,quantile(time_rp50,0.975)/24,quantile(time_rp70,0.975)/24,quantile(time_rp90,0.975)/24),lty=2)
lines(c(30,50,70,90),c(quantile(time_rp30,0.025)/24,quantile(time_rp50,0.025)/24,quantile(time_rp70,0.025)/24,quantile(time_rp90,0.025)/24),lty=2)



###subpatent patient
game_male_record <- matrix(0,1000,2880)
game_female_record <- matrix(0,1000,2880)
ring_record <- matrix(0,1000,2880)
total_record <- matrix(0,1000,2880)

asymp_time <- c(195, 234, 248, 286)
#10^4 195  10^5 234  10^6 248  10^7 286
#benchmark value of r_p-67.38 such that the parasitemia could be stabilized at a level of interest

simulation <- try(sto_model((benchmark[1])*4000,benchmark[2],benchmark[3],c(rep(benchmark[4],asymp_time[1]),rep(benchmark[4]-67.38,2880)),matrix(0,2880,3),benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],al=42,as=25,2880),silent = TRUE)
ring <- ifelse(rowSums(simulation[,1:25])==0,1,rowSums(simulation[,1:25]))/4000

total <- ifelse(rowSums(simulation[,c(1:25,60:83,88,89)])==0,1,rowSums(simulation[,c(1:25,60:83,88,89)]))/4000
plot(1:2880/24,ring,type="l",log = 'y',xlab = "Time post infection (Days)",ylab = "Parasitemia/ml",ylim=c(1,1e+8),xlim = c(1,120),cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(1:2880/24,simulation[,88]/4000)
abline(h=1e+05)


j=2
for (i in 1:1000) {
  simulation <- try(sto_model((benchmark[1])*4000,benchmark[2],benchmark[3],c(rep(benchmark[4],asymp_time[j]),rep(benchmark[4]-67.38,2880)),matrix(0,2880,3),benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],al=42,as=25,2880),silent = TRUE)
  ring_record[i,] <- ifelse(rowSums(simulation[,1:25])==0,1,rowSums(simulation[,1:25]))/4000
  total_record[i,] <- ifelse(rowSums(simulation[,c(1:25,60:83,88,89)])==0,1,rowSums(simulation[,c(1:25,60:83,88,89)]))/4000
  game_male_record[i,] <- as.numeric(simulation[,88])/4000
  game_female_record[i,] <- as.numeric(simulation[,89])/4000
  
}

library(foreach)
library(doParallel)
library(mvtnorm)
cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)
parallel::stopCluster(cl)

sporozoite_all <- list()
sporozoite_record <- matrix(0,120,720)

for (j in 1:1000) {
  result <- foreach(i=1:120,.packages=c("rootSolve")) %dopar% {
    sporozoite_record <- sporozoite_prev(1000,round(game_male_record[j,i*24],0), round(game_female_record[j,i*24],0),0.87,2.4,0.085,13.6, 2.7, 17.8, 2.2)
    sporozoite_record
  }
  for (i in 1:120) {
    sporozoite_record[i,] <- result[[i]]
  }
  sporozoite_all[[j]] <- sporozoite_record
  print(j)
}

prob_sporozoite <- matrix(0,1000,120)
for (i in 1:1000) {
  prob_sporozoite[i,] <- sporozoite_all[[i]][,15*24]
}


total_median <- vector()
total_upper <- vector()
total_lower <- vector()

ring_median <- vector()
ring_upper <- vector()
ring_lower <- vector()
male_median <- vector()
male_upper <- vector()
male_lower <- vector()
female_median <- vector()
female_upper <- vector()
female_lower <- vector()
prob_median <- vector()
prob_upper <- vector()
prob_lower <- vector()
for (i in 1:2880) {
  total_median[i] <- median(total_record[,i])
  total_upper[i] <- quantile(total_record[,i],0.975)
  total_lower[i] <- quantile(total_record[,i],0.025)
  ring_median[i] <- median(ring_record[,i])
  ring_upper[i] <- quantile(ring_record[,i],0.975)
  ring_lower[i] <- quantile(ring_record[,i],0.025)
  male_median[i] <- median(game_male_record[,i])
  male_upper[i] <- quantile(game_male_record[,i],0.975)
  male_lower[i] <- quantile(game_male_record[,i],0.025)
  female_median[i] <- median(game_female_record[,i])
  female_upper[i] <- quantile(game_female_record[,i],0.975)
  female_lower[i] <- quantile(game_female_record[,i],0.025)
}
for (i in 1:120) {
  prob_median[i] <- median(prob_sporozoite[,i])
  prob_upper[i] <- quantile(prob_sporozoite[,i],0.975)
  prob_lower[i] <- quantile(prob_sporozoite[,i],0.025)
}

asymp_host <- matrix(0,12,2880)
asymp_host[1,] <- total_median
asymp_host[2,] <- total_upper
asymp_host[3,] <- total_lower
asymp_host[4,] <- male_median
asymp_host[5,] <- male_upper
asymp_host[6,] <- male_lower
asymp_host[7,] <- female_median
asymp_host[8,] <- female_upper
asymp_host[9,] <- female_lower
asymp_host[10,] <- ring_median
asymp_host[11,] <- ring_upper
asymp_host[12,] <- ring_lower
prob_trans <- matrix(0,3,120)
prob_trans[1,] <- prob_median
prob_trans[2,] <- prob_upper
prob_trans[3,] <- prob_lower
write.csv(asymp_host,"asymp4.csv",row.names = FALSE)
write.csv(prob_trans,"prob4.csv",row.names = FALSE)

asymp_4 <- read.csv("asymp4.csv")
asymp_5 <- read.csv("asymp5.csv")
asymp_6 <- read.csv("asymp6.csv")
asymp_7 <- read.csv("asymp7.csv")
prob_4 <- read.csv("prob4.csv")
prob_5 <- read.csv("prob5.csv")
prob_6 <- read.csv("prob6.csv")
prob_7 <- read.csv("prob7.csv")

###Figure 8
par(mar = c(5.1, 5.1, 4.1, 1.1))
plot(1:2880/24,asymp_4[1,],type="l",xlab = "Day post infection (Days)",ylab="Parasitemia /ml blood",log="y",ylim=c(1,1e+9),col = rgb(192, 192, 192, maxColorValue=255),,cex.axis=2,cex.lab=2,cex.main=2,yaxt="n")
axis(2,at=c(10^0,10^3,10^6,10^9),label=expression(10^0,10^3,10^6,10^9),cex.axis=2)
polygon(c(1:2880/24, rev(1:2880/24)), c(asymp_4[2,], rev(asymp_4[3,])),
        col = rgb(192, 192, 192, maxColorValue=255, alpha=100), border = NA)  
lines(1:2880/24,asymp_5[1,] ,col = rgb(128, 128, 128, maxColorValue=255))
polygon(c(1:2880/24, rev(1:2880/24)), c(asymp_5[2,] , rev(asymp_5[3,] )),
        col = rgb(128, 128, 128, maxColorValue=255, alpha=100), border = NA)
lines(1:2880/24,asymp_6[1,] ,col = rgb(64, 64, 64, maxColorValue=255))
polygon(c(1:2880/24, rev(1:2880/24)), c(asymp_6[2,] , rev(asymp_6[3,] )),
        col = rgb(64, 64, 64, maxColorValue=255, alpha=100), border = NA)
lines(1:2880/24,asymp_7[1,] )
polygon(c(1:2880/24, rev(1:2880/24)), c(asymp_7[2,] , rev(asymp_7[3,] )),
        col = rgb(0, 0, 0, maxColorValue=255, alpha=100), border = NA)
legend("bottomright",col=c(rgb(0, 0, 0, maxColorValue=255),rgb(64, 64, 64, maxColorValue=255),rgb(128, 128, 128, maxColorValue=255),rgb(192, 192, 192, maxColorValue=255)),legend=c(expression(10^7),expression(10^6),expression(10^5),expression(10^4)),lty=1,lwd=3,bty="n",cex=1.5)



plot(1:120,prob_4[1,],type="l",ylim=c(0,1),xlim=c(-0.5,120),xlab = "Day post infection (Days)",ylab = expression(p[i]),col = rgb(255,153,204, maxColorValue=255),cex.axis=2,cex.lab=2,cex.main=2,yaxt="n")
axis(2,at=c(0,0.25,0.5,0.75,1),cex.axis=2)
polygon(c(1:120, rev(1:120)), c(prob_4[2,], rev(prob_4[3,])),
        col = rgb(255,153,204, maxColorValue=255, alpha=150), border = NA)  
lines(1:120,prob_5[1,],col = rgb(255,51,153, maxColorValue=255))
polygon(polygon(c(1:120, rev(1:120)), c(prob_5[2,], rev(prob_5[3,])),
                col = rgb(255,51,153, maxColorValue=255, alpha=150), border = NA)  )
lines(1:120,prob_6[1,],col = rgb(204,0,102, maxColorValue=255))
polygon(polygon(c(1:120, rev(1:120)), c(prob_6[2,], rev(prob_6[3,])),
                col = rgb(204,0,102, maxColorValue=255, alpha=150), border = NA)  )
lines(1:120,prob_7[1,],col = rgb(102,0,51, maxColorValue=255))
polygon(polygon(c(1:120, rev(1:120)), c(prob_7[2,], rev(prob_7[3,])),
                col = rgb(102,0,51, maxColorValue=255, alpha=150), border = NA)  )
legend("topleft",col=c(rgb(102,0,51, maxColorValue=255),rgb(204,0,102, maxColorValue=255),rgb(255,51,153, maxColorValue=255),rgb(255,153,204, maxColorValue=255)),legend=c(expression(10^7),expression(10^6),expression(10^5),expression(10^4)),lty=1,lwd=3,bty="n",cex=1.5)



plot(c(10^4,10^5,10^6,10^7),c(prob_4[1,18],prob_5[1,18],prob_6[1,18],prob_7[1,18]),log="x",xlab="Parasitemia /ml",ylab=expression(p[i]),ylim=c(0,1),xaxt="n",cex.axis=2,cex.lab=2,cex.main=2,col=rgb(255,153,204, maxColorValue=255),pch=19,yaxt="n")
axis(2,at=c(0,0.25,0.5,0.75,1),cex.axis=2)
axis(1,at=c(10^4,10^5,10^6,10^7) ,lab=c(expression(10^4),expression(10^5),expression(10^6),expression(10^7)),las=TRUE,cex.axis=2)
lines(c(10^4,10^5,10^6,10^7),c(prob_4[1,18],prob_5[1,18],prob_6[1,18],prob_7[1,18]),col=rgb(255,153,204, maxColorValue=255))
lines(c(10^4,10^5,10^6,10^7),c(prob_4[2,18],prob_5[2,18],prob_6[2,18],prob_7[2,18]),lty=2,col=rgb(255,153,204, maxColorValue=255, alpha=150))
lines(c(10^4,10^5,10^6,10^7),c(prob_4[3,18],prob_5[3,18],prob_6[3,18],prob_7[3,18]),lty=2,col=rgb(255,153,204, maxColorValue=255, alpha=150))

points(c(10^4,10^5,10^6,10^7),c(prob_4[1,20],prob_5[1,20],prob_6[1,20],prob_7[1,20]),col=rgb(204,0,102, maxColorValue=255),pch=19)
lines(c(10^4,10^5,10^6,10^7),c(prob_4[1,20],prob_5[1,20],prob_6[1,20],prob_7[1,20]),col=rgb(204,0,102, maxColorValue=255, alpha=150))
lines(c(10^4,10^5,10^6,10^7),c(prob_4[2,20],prob_5[2,20],prob_6[2,20],prob_7[2,20]),lty=2,col=rgb(204,0,102, maxColorValue=255, alpha=150))
lines(c(10^4,10^5,10^6,10^7),c(prob_4[3,20],prob_5[3,20],prob_6[3,20],prob_7[3,20]),lty=2,col=rgb(204,0,102, maxColorValue=255, alpha=150))

points(c(10^4,10^5,10^6,10^7),c(prob_4[1,40],prob_5[1,40],prob_6[1,40],prob_7[1,40]),col=rgb(102,0,51, maxColorValue=255),pch=19)
lines(c(10^4,10^5,10^6,10^7),c(prob_4[1,40],prob_5[1,40],prob_6[1,40],prob_7[1,40]),col=rgb(102,0,51, maxColorValue=255, alpha=150))
lines(c(10^4,10^5,10^6,10^7),c(prob_4[2,40],prob_5[2,40],prob_6[2,40],prob_7[2,40]),lty=2,col=rgb(102,0,51, maxColorValue=255, alpha=150))
lines(c(10^4,10^5,10^6,10^7),c(prob_4[3,40],prob_5[3,40],prob_6[3,40],prob_7[3,40]),lty=2,col=rgb(102,0,51, maxColorValue=255, alpha=150))
legend("topleft",col=c(rgb(102,0,51, maxColorValue=255),rgb(204,0,102, maxColorValue=255),rgb(255,153,204, maxColorValue=255)),legend=c(expression(t==40),expression(t==20),expression(t==18)),lty=1,lwd=3,bty="n",cex=1.5)

