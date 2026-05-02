library(lhs)
library(rootSolve)
library(truncnorm)
source("Models.R")
benchmark <- c(0.08,7,8,71,0.006,0.03,0.0004,0.007,0.02,0.2,3.4, 0.8, 0.029,13.8, 3.9, 19.3, 2.1)

#PRCC transmission probability
set.seed(1)
A <- randomLHS(1000, 13)
lhs_para <- matrix(0,1000,13)
for (i in 1:13) {
  lhs_para[,i] <- qunif(A[,i], min = 0.5*benchmark[i], max = 1.5*benchmark[i])
}

t_total <- 1440
game_record_male <- matrix(0,1000,1440)
game_record_female <- matrix(0,1000,1440)
asexual_record <- matrix(0,1000,1440)

for (i in 1:1000) {
  simulation <- try(sto_model((lhs_para[i,1])*4000,lhs_para[i,2],lhs_para[i,3],lhs_para[i,4],matrix(0,2880,3),lhs_para[i,5],lhs_para[i,6],lhs_para[i,7],lhs_para[i,8],lhs_para[i,8],lhs_para[i,9],lhs_para[i,10],al=42,as=25,t_total,10^6.5),silent = TRUE)
  asexual_record[i,] <- rowSums(simulation[,1:24])/4000
  game_record_male[i,] <- as.numeric(simulation[,88])/4000
  game_record_female[i,] <- as.numeric(simulation[,89])/4000
  print(i)
}


library(foreach)
library(doParallel)
library(mvtnorm)
cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)
parallel::stopCluster(cl)

time_record <- numeric(1000)

for (k in 1:10) {
  result <- foreach(i=1:100,.packages=c("tmvtnorm")) %dopar% {
    prob <- vector()
    d <- 0
    j <- 0
    while (d<0.0005){
      j <- j+1
      if (j>1440){
        break
      }
      if (game_record_male[(k-1)*100+i,j]>0){
        oocyst_record <- oocyst_prev(100000,round(game_record_male[(k-1)*100+i,j]),round(game_record_female[(k-1)*100+i,j]),lhs_para[(k-1)*100+i,11],lhs_para[(k-1)*100+i,12],lhs_para[(k-1)*100+i,13],benchmark[14],benchmark[15],benchmark[16],benchmark[17])
        d <- oocyst_record[8.5*24]
        prob[j] <- d
      }else{
        d <- 0
        prob[j] <- d
      }
    }
   j
  }
  for (i in 1:100) {
    time_record [(k-1)*100+i] <- result[[i]]
  }
  print(k)
}

prob_record <- vector()
for (i in 1:1000) {
  oocyst_record <- oocyst_prev(100000,round(game_record_male[i,t_total]),round(game_record_female[i,t_total]),lhs_para[i,11],lhs_para[i,12],lhs_para[i,13],benchmark[14],benchmark[15],benchmark[16],benchmark[17])
  prob_record[i] <- oocyst_record[8.5*24]
  print(i)
}

##t_p <- time_record ##t_{0.001}
t_p <- prob_record ##p_{60}

dataset <- cbind(lhs_para,t_p)
rank_trans_data <- matrix(0,1000,14)
for (i in 1:14){
  rank_trans_data[,i] <- rank(dataset[,i])
}
rank_trans_data <- as.data.frame(rank_trans_data)
names(rank_trans_data) <- c(expression(P[init]),expression(mu),expression(sigma),expression(r[p]),expression(f),expression(delta[p]),expression(delta[g]),expression(delta[gc]),expression(m),expression(f[g]),expression(V),expression(rho),expression(p[f]),"t_p")

name_para <- c(expression(P[init]),expression(mu),expression(sigma),expression(r[p]),expression(f),expression(delta[p]),expression(delta[g]),expression(delta[gc]),expression(m),expression(f[g]),expression(V),expression(rho),expression(p[f]),"t_p")

prcc_result <- vector()
for (i in 1:13){
  lm_p <- lm(t_p~. , data = rank_trans_data[,-c(i)])
  lm_para <- lm(rank_trans_data[,i]~. ,data = rank_trans_data[,-c(i,14)])
  prcc_result[i] <- cor(lm_p$residuals,lm_para$residuals)
}
prcc_result

col_names <-  c(expression(P[init]),expression(mu),expression(sigma),expression(r[p]),expression(f),expression(delta[p]),expression(delta[g]),expression(delta[gc]),expression(m),expression(f[g]),expression(V),expression(rho),expression(p[f]),"t_p")
human_order <- order(prcc_result[1:10])
mos_order <- order(prcc_result[11:13])+10
prcc_result_order <- prcc_result[c(mos_order,human_order), drop=TRUE]

##record the order of t_{0.001} in high transmission settings
order_high <- c(mos_order,human_order)
prcc_result_order<- prcc_result[order_high, drop=TRUE]
par(mar = c(4, 4.5, 2, 2))
bp <- barplot(prcc_result_order, names.arg = col_names[order_high],
        horiz = TRUE,          
        xlim = c(-1, 1),
        las = 1,              
        col = "grey",
        xlab = "PRCC",cex.axis = 2, width = 0.3,
        space = 0.6,yaxt="n",cex.lab=2)

axis(2, at = bp, labels = col_names[order_high], las = 1,cex.axis=2)

### marginal relationship
k=1  #totally 13 parameters, the order is the same as in col_names
x <- lhs_para[,k]
y <- t_p
x_seq <- seq(min(x), max(x), length.out = 500)
get_local_stats <- function(x0, x, y, frac) {
  d <- abs(x - x0)
  k <- max(50, ceiling(frac * length(x)))
  idx <- order(d)[1:k]
  c(
    med  = median(y[idx], na.rm = TRUE),
    low  = quantile(y[idx], 0.025, na.rm = TRUE),
    high = quantile(y[idx], 0.95, na.rm = TRUE)
  )
}
stats <- t(sapply(x_seq, get_local_stats, x = x, y = y, frac = 0.2))
fit_med  <- loess(stats[, 1]  ~ x_seq, span = 0.5)
fit_low  <- loess(stats[, 2]  ~ x_seq, span = 0.5)
fit_high <- loess(stats[, 3] ~ x_seq, span = 0.5)
med_s  <- predict(fit_med)
low_s  <- predict(fit_low)
high_s <- predict(fit_high)


par(mar = c(5.1, 5, 2, 1.5))
plot(x, y, col = adjustcolor("grey40", 0.3), pch = 16,xlim=c(benchmark[k]*0.5,benchmark[k]*1.5),ylim=c(0,1),xaxt="n",yaxt="n",ylab=expression(p[60]),xlab=name_para[k],cex.lab=2,cex.axis=2,,bty="l")
benchmark[k]*0.5
axis(1,at=c(0.02,0.03,0.04),cex.axis=2)
axis(2,at=c(0,0.25,0.5,0.75,1),cex.axis=2)
polygon(c(x_seq, rev(x_seq)),
        c(high_s, rev(low_s)),
        col = adjustcolor("black", 0.2), border = NA)
lines(x_seq, med_s, lwd = 1.5)

