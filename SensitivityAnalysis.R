transmission_model <- function(N, Gm, Gf, V, rho, r, alpha_go, beta_go, alpha_os, beta_os){
  a <- 0.0050300 
  b <- 0.1326802
  c <- 0.3272420
  record <- matrix(0, N, 7)
  M <- rbinom(N,Gm,V/1000)
  F <- rbinom(N,Gf,V/1000)
  n <- vector()
  t <- vector()
  for (i in 1:N) {
    n[i] <- rbinom(1,round(min(rho*M[i],F[i]),0),r)
    p <- runif(1,0,1)
    t[i] <- 1/b*log(1+log(p)/(-a/b))
  }
  record[,c(1,4,5,6)] <- c(n,t,M,F)
  for (i in 1:N) {
    if (n[i]!=0){
      oocy <- rgamma(n[i],alpha_go,beta_go)
      sporo <- rgamma(n[i],alpha_os,beta_os)
      p <- runif(1,0,1)
      t <- 1/b*log(1+log(p)/(-a/b*exp(c)))
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
library(lhs)
library(rootSolve)
library(truncnorm)
library(foreach)
library(doParallel)
cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)
parallel::stopCluster(cl)

benchmark <- c(0.08,7,5.5,71,0.006,0.03,0.0004,0.007,0.02,0.2,0.87,2.4,0.081,13.6, 2.7, 17.8, 2.2)



#PRCC transmission probability
benchmark_mos <- c(4, 0.87, 2.40, 0.081,13.6, 2.7, 17.8, 2.2)
set.seed(1)
A <- randomLHS(1000, 4)
lhs_para <- matrix(0,1000,4)
for (i in 1:4) {
  lhs_para[,i] <- qunif(A[,i], min = 0.5*benchmark_mos[i], max = 1.5*benchmark_mos[i])
}

trans_prob <- vector()
sporozoite_record <- matrix(0,1000,100)
for (i in 1:1000){
  G <- 10^5
  GM <- round(G/(lhs_para[i,1]+1),0)
  GF <- round(G-GM,0)
  result <- foreach(k=1:100,.packages=c("rootSolve")) %dopar% {
    set.seed(k)
    sporozoite_record <- sporozoite_prev(1000,GM,GF,lhs_para[i,2],lhs_para[i,3],lhs_para[i,4],13.6, 2.7, 17.8, 2.2)
    sporozoite_record[15*24]
  }
  for (k in 1:100) {
    trans_prob[k] <- result[[k]]
  }
  sporozoite_record[i,] <- trans_prob
  print(i)
}

for (i in 1:1000){
  trans_prob[i] <- median(sporozoite_record[i,])
}


dataset <- cbind(lhs_para,trans_prob)

rank_trans_data <- matrix(0,1000,5)
for (i in 1:5){
  rank_trans_data[,i] <- rank(dataset[,i])
}
rank_trans_data <- as.data.frame(rank_trans_data)
names(rank_trans_data) <- c(expression(G[F]/G[M]),expression(V),expression(rho),expression(p[f]),"tp")
               #,expression(alpha[ZO]),expression(beta[ZO]),expression(alpha[OS]),expression(beta[OS])      

prcc_result <- vector()
for (i in 1:4){
  lm_p <- lm(tp~. , data = rank_trans_data[,-c(i)])
  lm_para <- lm(rank_trans_data[,i]~. ,data = rank_trans_data[,-c(i,5)])
  prcc_result[i] <- cor(lm_p$residuals,lm_para$residuals)
}

col_names <-  c(expression(G[F]/G[M]),expression(V),expression(rho),expression(p[f]),"tp")

prcc_result_order <- prcc_result[order(prcc_result), drop=TRUE]
col_names_order <- col_names[order(prcc_result), drop=TRUE]
par(mar = c(4, 3.5, 2, 1.1))
dotchart(prcc_result_order, labels=col_names_order,pch=19, xlim=c(-1,1),
         xlab = "PRCC values",main=expression(G==10^3~ml^-1),cex=2)
abline(v=0,lty=2)

load("~/Desktop/Within-vector/figure/PRCC/mos_trans_10^5_100_1000_GF_GM_75.RData")


#-0.7425985  0.3879726  0.8935809  0.9692221
#-0.7564083  0.8179426  0.9081534  0.9403906
#-0.8375764  0.8895180  0.9128998  0.9241975

#-0.6453472  0.3238268  0.8706070  0.9575849
#-0.5926387  0.7338752  0.8598597  0.9061015
#-0.6365572  0.7959825  0.8579182  0.8763471
#PRCC human
benchmark <- c(0.08,7,5.5,71,0.006,0.03,0.0004,0.007,0.02,0.2)
set.seed(1)
A <- randomLHS(1000, 9)
lhs_para <- matrix(0,1000,9)
for (i in 1:9) {
  lhs_para[,i] <- qunif(A[,i], min = 0.25*benchmark[i], max = 1.75*benchmark[i])
}
game_record <- matrix(0,1000,960)
for (i in 1:1000) {
  simulation <- try(sto_model((lhs_para[i,1])*4000,lhs_para[i,2],lhs_para[i,3],rep(lhs_para[i,4],2880),matrix(0,2880,3),lhs_para[i,5],lhs_para[i,6],lhs_para[i,7],lhs_para[i,8],lhs_para[i,8],lhs_para[i,9],benchmark[10],al=42,as=25,960),silent = TRUE)
  game_record[i,] <- as.numeric(simulation[,88]+simulation[,89])/4000
}

tG <- 10^4
time_game <- vector()
for (i in 1:1000) {
  time_game[i] <- min(which(game_record[i,]>tG))
}
dataset <- cbind(lhs_para,time_game)

rank_trans_data <- matrix(0,1000,10)
for (i in 1:10){
  rank_trans_data[,i] <- rank(dataset[,i])
}
rank_trans_data <- as.data.frame(rank_trans_data)
names(rank_trans_data) <- c(expression(P[init]),expression(mu),expression(sigma),expression(r[p]),expression(f),expression(delta[p]),expression(delta[g]),expression(delta[gc]),expression(m),"G")


prcc_result <- vector()
for (i in 1:9){
  lm_p <- lm(G~. , data = rank_trans_data[,-c(i)])
  lm_para <- lm(rank_trans_data[,i]~. ,data = rank_trans_data[,-c(i,10)])
  prcc_result[i] <- cor(lm_p$residuals,lm_para$residuals)
}

col_names <-  c(expression(P[init]),expression(mu),expression(sigma),expression(r[p]),expression(f),expression(delta[p]),expression(delta[g]),expression(delta[gc]),expression(m),"G")
prcc_result_order <- prcc_result[order(prcc_result), drop=TRUE]
col_names_order <- col_names[order(prcc_result), drop=TRUE]
par(mar = c(4, 3.5, 2, 1.1))
dotchart(prcc_result_order, labels=col_names_order,pch=19, xlim=c(-1,1),
         xlab = "PRCC values",cex=1.5)
abline(v=0,lty=2)

#sensitivity mosquito
cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)
parallel::stopCluster(cl)

trans_prob <- vector()
sensitivity_GM_G <- matrix(0,100,5)
sensitivity_GM_G_all <- list()
set.seed(2)
GF_GM <- c(2,4,8,12,16)
G <- c(10^3,10^4,10^5)
for (l in 1:3) {
  for (i in 1:5) {
    GM <- round(G[l]/(GF_GM[i]+1),0)
    GF <- round(G[l]-GM,0)
    result <- foreach(j=1:100,.packages=c("rootSolve")) %dopar% {
      set.seed(j)
      sporozoite_record <- sporozoite_prev(1000,GM, GF,benchmark_mos[2],benchmark_mos[3],benchmark_mos[4], benchmark_mos[5], benchmark_mos[6], benchmark_mos[7],benchmark_mos[8])
      sporozoite_record[15*24]
    }
    for (k in 1:100) {
      trans_prob[k] <- result[[k]]
    }
    sensitivity_GM_G[,i] <- trans_prob
  }
  sensitivity_GM_G_all[[l]] <- sensitivity_GM_G
  print(l)
}

data_GM_G <- matrix(0,500,2)
for (i in 1:5) {
  data_GM_G[((i-1)*100+1):(i*100),1] <- sensitivity_GM_G_all[[1]][,i]
  data_GM_G[((i-1)*100+1):(i*100),2] <- rep(GF_GM[i],100)
}
data_GM_G <- as.data.frame(data_GM_G)
names(data_GM_G) <- c("Prob","GM_G")

p=ggplot(data_GM_G, aes(x=GM_G, y=Prob, group=as.factor(GM_G))) +
  geom_boxplot()+xlab(expression(G[F]/G[M])) + ylab(expression(p[i]))+ theme(
    axis.title.x = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text.y = element_text(size = 30),panel.background = element_rect(colour="white" ,fill = "white",linetype = "solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "grey"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "grey"),
    axis.line = element_line(size = 0.5, linetype = "solid",
                             colour = "black"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) + scale_y_continuous(limits = c(0,0.02)) +scale_x_continuous(breaks = c(2,4,6,8,10,12,14,16))
p
ggsave("GM_G1.png", plot = p, width = 8, height = 6)

###V
trans_prob <- vector() 
sensitivity_V <- matrix(0,100,6)
sensitivity_V_all <- list()
set.seed(2)
V <- c(0.1,0.2,0.4,0.8,1.2,1.6)
G <- c(10^3,10^4,10^5)
for (l in 1:3) {
  for (i in 1:6) {
    GM <- round(G[l]/(benchmark_mos[1]+1),0)
    GF <- round(G[l]-GM,0)
    result <- foreach(j=1:100,.packages=c("rootSolve")) %dopar% {
      set.seed(j)
      sporozoite_record <- sporozoite_prev(1000,GM, GF,V[i],benchmark_mos[3],benchmark_mos[4], benchmark_mos[5], benchmark_mos[6], benchmark_mos[7],benchmark_mos[8])
      sporozoite_record[15*24]
    }
    for (k in 1:100) {
      trans_prob[k] <- result[[k]]
    }
    sensitivity_V[,i] <- trans_prob
  }
  sensitivity_V_all[[l]] <- sensitivity_V
  print(l)
}

data_V <- matrix(0,600,2)
for (i in 1:6) {
  data_V[((i-1)*100+1):(i*100),1] <- sensitivity_V_all[[1]][,i]
  data_V[((i-1)*100+1):(i*100),2] <- rep(V[i],100)
}
data_V <- as.data.frame(data_V)
names(data_V) <- c("Prob","V")


p <- ggplot(data_V, aes(x=V, y=Prob, group=as.factor(V))) +
  geom_boxplot()+xlab(expression(V~(mu~L^-1))) + ylab(expression(p[i]))+ theme(
    axis.title.x = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text.y = element_text(size = 30),panel.background = element_rect(colour="white" ,fill = "white",linetype = "solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "grey"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "grey"),
    axis.line = element_line(size = 0.5, linetype = "solid",
                             colour = "black"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) + scale_y_continuous(limits = c(0,0.03)) +scale_x_continuous(breaks=c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6))
p
ggsave("V1.png", plot = p, width = 8, height = 6)

###rho
trans_prob <- vector()
sensitivity_rho <- matrix(0,100,8)
sensitivity_rho_all <- list()
set.seed(2)
rho <- 1:8
G <- c(10^3,10^4,10^5)
for (l in 1:3) {
  for (i in 1:8) {
    GM <- round(G[l]/(benchmark_mos[1]+1),0)
    GF <- round(G[l]-GM,0)
    result <- foreach(j=1:100,.packages=c("rootSolve")) %dopar% {
      set.seed(j)
      sporozoite_record <- sporozoite_prev(1000,GM, GF,benchmark_mos[2],rho[i],benchmark_mos[4], benchmark_mos[5], benchmark_mos[6], benchmark_mos[7],benchmark_mos[8])
      sporozoite_record[15*24]
    }
    for (k in 1:100) {
      trans_prob[k] <- result[[k]]
    }
    sensitivity_rho[,i] <- trans_prob
  }
  sensitivity_rho_all[[l]] <- sensitivity_rho
  print(l)
}

data_rho <- matrix(0,800,2)
for (i in 1:8) {
  data_rho[((i-1)*100+1):(i*100),1] <- sensitivity_rho_all[[1]][,i]
  data_rho[((i-1)*100+1):(i*100),2] <- rep(rho[i],100)
}
data_rho <- as.data.frame(data_rho)
names(data_rho) <- c("Prob","rho")

p=ggplot(data_rho, aes(x=rho, y=Prob, group=as.factor(rho))) +
  geom_boxplot()+xlab(expression(rho)) + ylab(expression(p[i]))+ theme(
    axis.title.x = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text.y = element_text(size = 30),panel.background = element_rect(colour="white" ,fill = "white",linetype = "solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "grey"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "grey"),
    axis.line = element_line(size = 0.5, linetype = "solid",
                             colour = "black"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) + scale_y_continuous() 
p
ggsave("rho1.png", plot = p, width = 8, height = 6)

###p_f
trans_prob <- vector()
sensitivity_pf <- matrix(0,100,6)
sensitivity_pf_all <- list()
set.seed(2)
pf <- c(0.01,0.02,0.04,0.08,0.12,0.16) 
G <- c(10^3,10^4,10^5)
for (l in 1:3) {
  for (i in 1:6) {
    GM <- round(G[l]/(benchmark_mos[1]+1),0)
    GF <- round(G[l]-GM,0)
    result <- foreach(j=1:100,.packages=c("rootSolve")) %dopar% {
      set.seed(j)
      sporozoite_record <- sporozoite_prev(1000,GM, GF,benchmark_mos[2],benchmark_mos[3],pf[i], benchmark_mos[5], benchmark_mos[6], benchmark_mos[7],benchmark_mos[8])
      sporozoite_record[15*24]
    }
    for (k in 1:100) {
      trans_prob[k] <- result[[k]]
    }
    sensitivity_pf[,i] <- trans_prob
  }
  sensitivity_pf_all[[l]] <- sensitivity_pf
  print(l)
}

data_pf <- matrix(0,600,2)
for (i in 1:6) {
  data_pf[((i-1)*100+1):(i*100),1] <- sensitivity_pf_all[[1]][,i]
  data_pf[((i-1)*100+1):(i*100),2] <- rep(pf[i],100)
}
data_pf <- as.data.frame(data_pf)
names(data_pf) <- c("Prob","pf")

p=ggplot(data_pf, aes(x=pf, y=Prob, group=as.factor(pf))) +
  geom_boxplot()+xlab(expression(p[f])) + ylab(expression(p[i]))+ theme(
    axis.title.x = element_text(size = 30),
    axis.text.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text.y = element_text(size = 30),panel.background = element_rect(colour="white" ,fill = "white",linetype = "solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "grey"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "grey"),
    axis.line = element_line(size = 0.5, linetype = "solid",
                             colour = "black"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) +scale_y_continuous() + scale_x_continuous(breaks = c(0,0.04,0.08,0.12,0.16)) 
p
ggsave("pf1.png", plot = p, width = 8, height = 6)

###sensitivity human

rp <- c(30,50,70,90,110)
sensitivity_rp <- list()
for (i in 1:length(rp)) {
  game_record <- matrix(0,100,960)
  for (j in 1:100) {
    set.seed(j)
    simulation <- try(sto_model((benchmark[1])*4000,benchmark[2],benchmark[3],rep(rp[i],2880),matrix(0,2880,3),benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],al=42,as=25,960),silent = TRUE)
    game_record[j,] <- as.numeric(simulation[,88]+simulation[,89])/4000
  }
  sensitivity_rp[[i]] <- game_record
  print(i)
}

tG <- c(10^3,10^4,10^5)
time_game_rp_all <- list()
for (k in 1:3) {
  time_game_rp <- matrix(0,5,100)
  for (i in 1:5) {
    for (j in 1:100) {
      time_game_rp[i,j] <- min(which(sensitivity_rp[[i]][j,]>tG[k]))
    }
  }
  time_game_rp_all[[k]] <- time_game_rp
}


median_time_game_rp <- matrix(0,3,5)
upper_time_game_rp <- matrix(0,3,5)
lower_time_game_rp <- matrix(0,3,5)
for (k in 1:3) {
  for (i in 1:5) {
    median_time_game_rp[k,i] <- median(time_game_rp_all[[k]][i,])
    upper_time_game_rp[k,i] <- quantile(time_game_rp_all[[k]][i,],0.975)
    lower_time_game_rp[k,i] <- quantile(time_game_rp_all[[k]][i,],0.025)
  }
}
par(mar = c(5.1, 5.3, 4.1, 1.0))
plot(rp,median_time_game_rp[2,]/24,xlab=expression(r[p]),ylab=expression(t[G]),xaxt="n",cex.axis=2,cex.lab=2,cex.main=2,pch=19,yaxt="n")
axis(2,at=c(14,16,18,20),cex.axis=2)
axis(1,at=c(30,50,70,90,110),las=TRUE,cex.axis=2)
lines(rp,upper_time_game_rp[2,]/24,lty=2)
lines(rp,lower_time_game_rp[2,]/24,lty=2)

game_time_rp <- matrix(0,5,960)
game_median_rp <- matrix(0,5,960)
game_upper_rp <-  matrix(0,5,960)
game_lower_rp <-  matrix(0,5,960)
for (i in 1:5) {
  for (j in 1:960) {
    game_median_rp[i,j] <- median(sensitivity_rp[[i]][,j])
    game_upper_rp[i,j] <- quantile(sensitivity_rp[[i]][,j],0.975)
    game_lower_rp[i,j] <- quantile(sensitivity_rp[[i]][,j],0.025)
  }
}

par(mar = c(5.1, 5.1, 4.1, 1.1))
plot(1:960/24,game_median_rp[1,],type="l",xlab = "Day post infection (Days)",ylab="Gametocytemia /ml blood",log="y",col=rgb(139, 0, 0,maxColorValue = 255),ylim=c(1,1e+9),cex.axis=2,cex.lab=2,cex.main=2,yaxt="n",xlim=c(4,20))
axis(2,at=c(10^0,10^3,10^6,10^9),label=expression(10^0,10^3,10^6,10^9),cex.axis=2)
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_rp[1,], rev(game_lower_rp[1,])),col=rgb(139, 0, 0,maxColorValue = 255,alpha=50), border = NA)  
lines(1:960/24,game_median_rp[2,] ,col = rgb(255, 0, 0,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_rp[2,], rev(game_lower_rp[2,])),col=rgb(255, 0, 0,maxColorValue = 255,alpha=50), border = NA)  

lines(1:960/24,game_median_rp[3,] ,col = rgb(255, 165, 0,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_rp[3,], rev(game_lower_rp[3,])),col=rgb(255, 165, 0,maxColorValue = 255,alpha=50), border = NA)  

lines(1:960/24,game_median_rp[4,],col = rgb(0, 255, 0,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_rp[4,], rev(game_lower_rp[4,])),col=rgb(0, 255, 0,maxColorValue = 255,alpha=50), border = NA)  

lines(1:960/24,game_median_rp[5,],col = rgb(0, 191, 255,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_rp[5,], rev(game_lower_rp[5,])),col=rgb(0, 191, 255,maxColorValue = 255,alpha=50), border = NA)  

legend("topleft",col=c(rgb(139, 0, 0,maxColorValue = 255),rgb(255, 0, 0,maxColorValue = 255),rgb(255, 165, 0,maxColorValue = 255),rgb(0, 255, 0,maxColorValue = 255),rgb(0, 191, 255,maxColorValue = 255)),legend=c(expression(r[p]==30(PMF==8.5)),expression(r[p]==50(PMF==14.1)),expression(r[p]==70(PMF==19.7)),expression(r[p]==90(PMF==25.4)),expression(r[p]==110(PMF==31.0))),lty=1,lwd=1,bty="n",cex=1.5)




deltap <- c(0.01,0.02,0.03,0.04,0.05)
sensitivity_deltap <- list()
for (i in 1:length(deltap)) {
  game_record <- matrix(0,100,960)
  for (j in 1:100) {
    set.seed(j)
    simulation <- try(sto_model((benchmark[1])*4000,benchmark[2],benchmark[3],rep(benchmark[4],2880),matrix(0,2880,3),benchmark[5],deltap[i],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],al=42,as=25,960),silent = TRUE)
    game_record[j,] <- as.numeric(simulation[,88]+simulation[,89])/4000
  }
  sensitivity_deltap[[i]] <- game_record
  print(i)
}

tG <- c(10^3,10^4,10^5)
time_game_deltap_all <- list()
for (k in 1:3) {
  time_game_deltap <- matrix(0,5,100)
  for (i in 1:5) {
    for (j in 1:100) {
      time_game_deltap[i,j] <- min(which(sensitivity_deltap[[i]][j,]>tG[k]))
    }
  }
  time_game_deltap_all[[k]] <- time_game_deltap
}


median_time_game_deltap <- matrix(0,3,5)
upper_time_game_deltap <- matrix(0,3,5)
lower_time_game_deltap <- matrix(0,3,5)
for (k in 1:3) {
  for (i in 1:5) {
    median_time_game_deltap[k,i] <- median(time_game_deltap_all[[k]][i,])
    upper_time_game_deltap[k,i] <- quantile(time_game_deltap_all[[k]][i,],0.975)
    lower_time_game_deltap[k,i] <- quantile(time_game_deltap_all[[k]][i,],0.025)
  }
}

plot(deltap,median_time_game_deltap[2,]/24,xlab=expression(delta[p](h^-1)),ylab=expression(t[G]),xaxt="n",cex.axis=2,cex.lab=2,cex.main=2,pch=19,yaxt="n")
axis(2,at=c(12,14,16,18,20),cex.axis=2)
axis(1,at=deltap,las=TRUE,cex.axis=2)
lines(deltap,upper_time_game_deltap[2,]/24,lty=2)
lines(deltap,lower_time_game_deltap[2,]/24,lty=2)

game_time_deltap <- matrix(0,5,960)
game_median_deltap <- matrix(0,5,960)
game_upper_deltap <-  matrix(0,5,960)
game_lower_deltap <-  matrix(0,5,960)
for (i in 1:5) {
  for (j in 1:960) {
    game_median_deltap[i,j] <- median(sensitivity_deltap[[i]][,j])
    game_upper_deltap[i,j] <- quantile(sensitivity_deltap[[i]][,j],0.975)
    game_lower_deltap[i,j] <- quantile(sensitivity_deltap[[i]][,j],0.025)
  }
}


par(mar = c(5.1, 5.1, 4.1, 1.1))
plot(1:960/24,game_median_deltap[1,],type="l",xlab = "Day post infection (Days)",ylab="Gametocytemia /ml blood",log="y",col=rgb(139, 0, 0,maxColorValue = 255),ylim=c(1,1e+9),cex.axis=2,cex.lab=2,cex.main=2,yaxt="n",xlim=c(4,20))
axis(2,at=c(10^0,10^3,10^6,10^9),label=expression(10^0,10^3,10^6,10^9),cex.axis=2)
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_deltap[1,], rev(game_lower_deltap[1,])),col=rgb(139, 0, 0,maxColorValue = 255,alpha=50), border = NA)  
lines(1:960/24,game_median_deltap[2,] ,col = rgb(255, 0, 0,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_deltap[2,], rev(game_lower_deltap[2,])),col=rgb(255, 0, 0,maxColorValue = 255,alpha=50), border = NA)  

lines(1:960/24,game_median_deltap[3,] ,col = rgb(255, 165, 0,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_deltap[3,], rev(game_lower_deltap[3,])),col=rgb(255, 165, 0,maxColorValue = 255,alpha=50), border = NA)  

lines(1:960/24,game_median_deltap[4,],col = rgb(0, 255, 0,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_deltap[4,], rev(game_lower_deltap[4,])),col=rgb(0, 255, 0,maxColorValue = 255,alpha=50), border = NA)  

lines(1:960/24,game_median_deltap[5,],col = rgb(0, 191, 255,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_deltap[5,], rev(game_lower_deltap[5,])),col=rgb(0, 191, 255,maxColorValue = 255,alpha=50), border = NA)  

legend("topleft",col=c(rgb(139, 0, 0,maxColorValue = 255),rgb(255, 0, 0,maxColorValue = 255),rgb(255, 165, 0,maxColorValue = 255),rgb(0, 255, 0,maxColorValue = 255),rgb(0, 191, 255,maxColorValue = 255)),legend=c(expression(delta[p]==0.01(PMF==46.4)),expression(delta[p]==0.02(PMF==30.5)),expression(delta[p]==0.03(PMF==20.0)),expression(delta[p]==0.04(PMF==13.2)),expression(delta[p]==0.05(PMF==8.6))),lty=1,lwd=1,bty="n",cex=1.5)



m <- c(0.01,0.02,0.03,0.04,0.05)
sensitivity_m <- list()
for (i in 1:length(m)) {
  game_record <- matrix(0,100,960)
  for (j in 1:100) {
    set.seed(j)
    simulation <- try(sto_model((benchmark[1])*4000,benchmark[2],benchmark[3],rep(benchmark[4],2880),matrix(0,2880,3),benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],m[i],benchmark[10],al=42,as=25,960),silent = TRUE)
    game_record[j,] <- as.numeric(simulation[,88]+simulation[,89])/4000
  }
  sensitivity_m[[i]] <- game_record
  print(i)
}

tG <- c(10^3,10^4,10^5)
time_game_m_all <- list()
for (k in 1:3) {
  time_game_m <- matrix(0,5,100)
  for (i in 1:5) {
    for (j in 1:100) {
      time_game_m[i,j] <- min(which(sensitivity_m[[i]][j,]>tG[k]))
    }
  }
  time_game_m_all[[k]] <- time_game_m
}


median_time_game_m <- matrix(0,3,5)
upper_time_game_m <- matrix(0,3,5)
lower_time_game_m <- matrix(0,3,5)
for (k in 1:3) {
  for (i in 1:5) {
    median_time_game_m[k,i] <- median(time_game_m_all[[k]][i,])
    upper_time_game_m[k,i] <- quantile(time_game_m_all[[k]][i,],0.975)
    lower_time_game_m[k,i] <- quantile(time_game_m_all[[k]][i,],0.025)
  }
}

plot(m,median_time_game_m[2,]/24,xlab=expression(m(h^-1)),ylab=expression(t[G]),xaxt="n",cex.axis=2,cex.lab=2,cex.main=2,pch=19,yaxt="n",ylim=c(13,17))
axis(2,at=c(13,15,17),cex.axis=2)
axis(1,at=c(0.01,0.02,0.03,0.04,0.05,0.06),las=TRUE,cex.axis=2)
lines(m,upper_time_game_m[2,]/24,lty=2)
lines(m,lower_time_game_m[2,]/24,lty=2)

game_time_m <- matrix(0,5,960)
game_median_m <- matrix(0,5,960)
game_upper_m <-  matrix(0,5,960)
game_lower_m <-  matrix(0,5,960)
for (i in 1:5) {
  for (j in 1:960) {
    game_median_m[i,j] <- median(sensitivity_m[[i]][,j])
    game_upper_m[i,j] <- quantile(sensitivity_m[[i]][,j],0.975)
    game_lower_m[i,j] <- quantile(sensitivity_m[[i]][,j],0.025)
  }
}


par(mar = c(5.1, 5.1, 4.1, 1.1))
plot(1:960/24,game_median_m[1,],type="l",xlab = "Day post infection (Days)",ylab="Gametocytemia /ml blood",log="y",col=rgb(139, 0, 0,maxColorValue = 255),ylim=c(1,1e+9),cex.axis=2,cex.lab=2,cex.main=2,yaxt="n",xlim=c(4,20))
axis(2,at=c(10^0,10^3,10^6,10^9),label=expression(10^0,10^3,10^6,10^9),cex.axis=2)
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_m[1,], rev(game_lower_m[1,])),col=rgb(139, 0, 0,maxColorValue = 255,alpha=50), border = NA)  
lines(1:960/24,game_median_m[2,] ,col = rgb(255, 0, 0,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_m[2,], rev(game_lower_m[2,])),col=rgb(255, 0, 0,maxColorValue = 255,alpha=50), border = NA)  

lines(1:960/24,game_median_m[3,] ,col = rgb(255, 165, 0,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_m[3,], rev(game_lower_m[3,])),col=rgb(255, 165, 0,maxColorValue = 255,alpha=50), border = NA)  

lines(1:960/24,game_median_m[4,],col = rgb(0, 255, 0,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_m[4,], rev(game_lower_m[4,])),col=rgb(0, 255, 0,maxColorValue = 255,alpha=50), border = NA)  

lines(1:960/24,game_median_m[5,],col = rgb(0, 191, 255,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_m[5,], rev(game_lower_m[5,])),col=rgb(0, 191, 255,maxColorValue = 255,alpha=50), border = NA)  


legend("topleft",col=c(rgb(139, 0, 0,maxColorValue = 255),rgb(255, 0, 0,maxColorValue = 255),rgb(255, 165, 0,maxColorValue = 255),rgb(0, 255, 0,maxColorValue = 255),rgb(0, 191, 255,maxColorValue = 255)),legend=c(expression(m==0.01),expression(m==0.02),expression(m==0.03),expression(m==0.04),expression(m==0.05)),lty=1,lwd=1.5,bty="n",cex=1.5)



f <- c(0.0025,0.005,0.01,0.02,0.1,0.2)
sensitivity_f <- list()
for (i in 1:length(f)) {
  game_record <- matrix(0,100,960)
  for (j in 1:100) {
    set.seed(j)
    simulation <- try(sto_model((benchmark[1])*4000,benchmark[2],benchmark[3],rep(benchmark[4],2880),matrix(0,2880,3),f[i],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],al=42,as=25,960),silent = TRUE)
    game_record[j,] <- as.numeric(simulation[,88]+simulation[,89])/4000
  }
  sensitivity_f[[i]] <- game_record
  print(i)
}

tG <- c(10^3,10^4,10^5)
time_game_f_all <- list()
for (k in 1:3) {
  time_game_f <- matrix(0,6,100)
  for (i in 1:6) {
    for (j in 1:100) {
      time_game_f[i,j] <- min(which(sensitivity_f[[i]][j,]>tG[k]))
    }
  }
  time_game_f_all[[k]] <- time_game_f
}


median_time_game_f <- matrix(0,3,6)
upper_time_game_f <- matrix(0,3,6)
lower_time_game_f <- matrix(0,3,6)
for (k in 1:3) {
  for (i in 1:6) {
    median_time_game_f[k,i] <- median(time_game_f_all[[k]][i,])
    upper_time_game_f[k,i] <- quantile(time_game_f_all[[k]][i,],0.975)
    lower_time_game_f[k,i] <- quantile(time_game_f_all[[k]][i,],0.025)
  }
}

plot(f,median_time_game_f[2,]/24,xlab=expression(f),ylab=expression(t[G]),xaxt="n",cex.axis=2,cex.lab=2,cex.main=2,pch=19,yaxt="n",ylim=c(13,17))
axis(2,at=c(13,15,17),cex.axis=2)
axis(1,at=c(0,0.05,0.1,0.15,0.2),las=TRUE,cex.axis=2)
lines(f,upper_time_game_f[2,]/24,lty=2)
lines(f,lower_time_game_f[2,]/24,lty=2)


game_time_f <- matrix(0,6,960)
game_median_f <- matrix(0,6,960)
game_upper_f <-  matrix(0,6,960)
game_lower_f <-  matrix(0,6,960)
for (i in 1:6) {
  for (j in 1:960) {
    game_median_f[i,j] <- median(sensitivity_f[[i]][,j])
    game_upper_f[i,j] <- quantile(sensitivity_f[[i]][,j],0.975)
    game_lower_f[i,j] <- quantile(sensitivity_f[[i]][,j],0.025)
  }
}


par(mar = c(5.1, 5.1, 4.1, 1.1))
plot(1:960/24,game_median_f[1,],type="l",xlab = "Day post infection (Days)",ylab="Gametocytemia /ml blood",col=rgb(139, 0, 0,maxColorValue = 255),log="y",ylim=c(1,1e+9),cex.axis=2,cex.lab=2,cex.main=2,yaxt="n",xlim=c(7,20))
axis(2,at=c(10^0,10^3,10^6,10^9),label=expression(10^0,10^3,10^6,10^9),cex.axis=2)
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_f[1,], rev(game_lower_f[1,])),col=rgb(139, 0, 0,maxColorValue = 255,alpha=50), border = NA)  
lines(1:960/24,game_median_f[2,] ,col = rgb(255, 0, 0,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_f[2,], rev(game_lower_f[2,])),col=rgb(255, 0, 0,maxColorValue = 255,alpha=50), border = NA)  

lines(1:960/24,game_median_f[3,] ,col = rgb(255, 165, 0,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_f[3,], rev(game_lower_f[3,])),col=rgb(255, 165, 0,maxColorValue = 255,alpha=50), border = NA)  

lines(1:960/24,game_median_f[4,],col = rgb(0, 255, 0,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_f[4,], rev(game_lower_f[4,])),col=rgb(0, 255, 0,maxColorValue = 255,alpha=50), border = NA)  

lines(1:960/24,game_median_f[5,],col = rgb(0, 191, 255,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_f[5,], rev(game_lower_f[5,])),col=rgb(0, 191, 255,maxColorValue = 255,alpha=50), border = NA)  

lines(1:960/24,game_median_f[6,],col = rgb(216, 191, 216,maxColorValue = 255))
polygon(c(1:960/24, rev(1:960/24)), c(game_upper_f[6,], rev(game_lower_f[6,])),col=rgb(216, 191, 216,maxColorValue = 255,alpha=50), border = NA)  

legend("topleft",col=c(col=rgb(139, 0, 0,maxColorValue = 255),rgb(255, 0, 0,maxColorValue = 255),rgb(255, 165, 0,maxColorValue = 255),rgb(0, 255, 0,maxColorValue = 255),rgb(0, 191, 255,maxColorValue = 255),rgb(216, 191, 216,maxColorValue = 255)
                       ),legend=c(expression(f==0.0025),expression(f==0.005),expression(f==0.01),expression(f==0.02),expression(f==0.1),expression(f==0.2)),lty=1,lwd=1.5,bty="n",cex=1.5)


###subpatent patient
game_male_record <- matrix(0,1000,2880)
game_female_record <- matrix(0,1000,2880)
ring_record <- matrix(0,1000,2880)
total_record <- matrix(0,1000,2880)

asymp_time <- c(195, 234, 248, 286)
#10^4 195 -67.38 10^5 234  10^6 248  10^7 286


simulation <- try(sto_model((benchmark[1])*4000,benchmark[2],benchmark[3],c(rep(benchmark[4],234),rep(benchmark[4]-67.38,2880)),matrix(0,2880,3),benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],al=42,as=25,2880),silent = TRUE)
ring <- ifelse(rowSums(simulation[,1:25])==0,1,rowSums(simulation[,1:25]))/4000

total <- ifelse(rowSums(simulation[,c(1:25,60:83,88,89)])==0,1,rowSums(simulation[,c(1:25,60:83,88,89)]))/4000
plot(1:2880/24,ring,type="l",log = 'y',xlab = "Time post infection (Days)",ylab = "Parasitemia/ml",ylim=c(1,1e+8),xlim = c(1,120),cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
lines(1:2880/24,simulation[,88]/4000)
abline(h=1e+05)


j=4
for (i in 1:1000) {
  simulation <- try(sto_model((benchmark[1])*4000,benchmark[2],benchmark[3],c(rep(benchmark[4],asymp_time[j]),rep(benchmark[4]-67.38,2880)),matrix(0,2880,3),benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],al=42,as=25,2880),silent = TRUE)
  ring_record[i,] <- ifelse(rowSums(simulation[,1:25])==0,1,rowSums(simulation[,1:25]))/4000
  total_record[i,] <- ifelse(rowSums(simulation[,c(1:25,60:83,88,89)])==0,1,rowSums(simulation[,c(1:25,60:83,88,89)]))/4000
  game_male_record[i,] <- as.numeric(simulation[,88])/4000
  game_female_record[i,] <- as.numeric(simulation[,89])/4000
  print(i)
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

asymp10 <- matrix(0,12,2880)
asymp10[1,] <- total_median
asymp10[2,] <- total_upper
asymp10[3,] <- total_lower
asymp10[4,] <- male_median
asymp10[5,] <- male_upper
asymp10[6,] <- male_lower
asymp10[7,] <- female_median
asymp10[8,] <- female_upper
asymp10[9,] <- female_lower
asymp10[10,] <- ring_median
asymp10[11,] <- ring_upper
asymp10[12,] <- ring_lower
prob7 <- matrix(0,3,120)
prob7[1,] <- prob_median
prob7[2,] <- prob_upper
prob7[3,] <- prob_lower
write.csv(asymp10,"asymp10_7.csv",row.names = FALSE)
write.csv(prob7,"prob10_7.csv",row.names = FALSE)

asymp_4 <- read.csv("asymp10_4.csv")
asymp_5 <- read.csv("asymp10_5.csv")
asymp_6 <- read.csv("asymp10_6.csv")
asymp_7 <- read.csv("asymp10_7.csv")
prob_4 <- read.csv("prob10_4.csv")
prob_5 <- read.csv("prob10_5.csv")
prob_6 <- read.csv("prob10_6.csv")
prob_7 <- read.csv("prob10_7.csv")

par(mar = c(5.1, 5.1, 4.1, 1.1))
plot(1:2880/24,asymp_4[1,],type="l",xlab = "Day post infection (Days)",ylab="Parasitemia /ml blood",log="y",ylim=c(1,1e+9),col = rgb(192, 192, 192, maxColorValue=255),,cex.axis=2,cex.lab=2,cex.main=2,yaxt="n",xlim=c(-5,120))
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






plot(1:120,prob_4[1,],type="l",ylim=c(0,1),xlim=c(-5,120),xlab = "Day post infection (Days)",ylab = expression(p[i]),col = rgb(255,153,204, maxColorValue=255),cex.axis=2,cex.lab=2,cex.main=2,yaxt="n")
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

