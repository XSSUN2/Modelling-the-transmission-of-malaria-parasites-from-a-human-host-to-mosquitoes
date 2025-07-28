transmission_model <- function(N, Gm, Gf, V,rho, r, alpha_go, beta_go, alpha_os, beta_os){
  a <- 0.0050300 
  b <- 0.1326802
  c <- 0.3272420
  record <- matrix(0, N, 7)
  M <- rbinom(N,Gm,V/1000)
  F <- rbinom(N,Gf,V/1000)
  n <- vector()
  t <- vector()
  for (i in 1:N) {
    n[i] <- rbinom(1,round(min(rho*M[i],F[i])),r)
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


oocyst_prev <- function(N, Gm, Gf, V, rho, r, alpha_go, beta_go, alpha_os, beta_os){
  a <- transmission_model(N, Gm, Gf, V, rho, r, alpha_go, beta_go, alpha_os, beta_os)
  survival <- matrix(0,N,240)
  for (i in 1:N) {
    if (ceiling(a[i,4]*24)<240){
      survival[i,(ceiling(a[i,4]*24)):240] <- 1
    }
  }
  oocyst <- matrix(0,N,240)
  for (i in 1:N) {
    if (ceiling(a[i,2]*24)<240 & a[i,2]!=0 ){
      oocyst[i, (ceiling(a[i,2]*24)):min((ceiling(a[i,7]*24)),240)] <- 1
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
cl <- parallel::makeCluster(12)
doParallel::registerDoParallel(cl)
parallel::stopCluster(cl)
library(rootSolve)

time.step=20
n=5000
d_total <- matrix(0,time.step,n)
weights <- matrix(1/n,time.step,n) 
para_rec <- matrix(0,n,3)
times <- matrix(0,time.step,n)
para <- list()
epsilon <- rep(0,time.step)
epsilon[1] <- 0.3699087


for (t in 11:13){
  if(t==1){
    rec <- matrix(0,n,4)
    result <- foreach(i=1:n,.packages=c("tmvtnorm")) %dopar% {
      d <- epsilon[t]+1
      time <- 0
      while(d > epsilon[t]){
        set.seed(n*time+i)
        proposed.mu <- c(runif(1,0,1)*runif(1,1,8),runif(1,0,1))
        oocyst_record <- vector()
        for (j in 1:33){
          simulation <- oocyst_prev(30,round(as.numeric(data[j,2]),0),round(as.numeric(data[j,3]),0),0.87,proposed.mu[1],proposed.mu[2],13.7, 2.7, 17.8, 2.2)
          oocyst_record[j] <- simulation[8.5*24]
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
          oocyst_record[j] <- simulation[8.5*24]
        }
        time <- time+1
        d <- sqrt(sum((oocyst_record-as.numeric(data[,4])/as.numeric(data[,5]))^2))
      }
      mu.weights.denominator=0
      for (k in 1:n){
        mu.weights.denominator <- mu.weights.denominator+weights[t-1,k]*dtmvnorm(proposed.mu,as.numeric(para[[t-1]][k,]),cov_matr,c(0,0) ,c(8,1))
      }
      prior_rho <- ifelse(proposed.mu[1]<1,log(8)/7,log(8/proposed.mu)/7)
      mu.weights.numerator = dunif(proposed.mu[2],0 ,1)*prior_rho
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
  print(t)
  save.image("~/Desktop/Untitled.RData")
}

result <- as.data.frame(para[[13]])
names(result) <- c("rho","p_f")

#hist
z<-1:80/10-0.05
f_Z <- ifelse(z >= 0 & z <= 1, (1/7) * log(8),  (1/7) * log(8/z))
integral_value <- sum(f_Z) 
f_Z_normalized <- f_Z / integral_value
prior <- f_Z_normalized*5000
prior_df <- as.data.frame(cbind(z,prior))
names(prior_df) <- c("x","y")

par(mar=c(5,5, 2,2))
hist(para[[12]][,1],breaks=80,main="",xlim=c(0,8),ylim=c(0,400),xlab=expression(rho),cex.axis=2,cex.lab=2)
abline(v=0.933,col="red")
lines(z, prior,lty=2, lwd = 1)
legend("topright",lty=c(1,2),col=c("red",1),legend=c("Estimate","Prior"),cex=2,bty="n",lwd=2)


prior_pf_x <- 0:50/50
prior_pf_y <- rep(50,51)
prior_pf <- as.data.frame(cbind(prior_pf_x,prior_pf_y))
names(prior_pf) <- c("x","y")

hist(para[[12]][,2],breaks=100,main="",xlab=expression(p[f]),xlim=c(0,0.5),ylim=c(0,300),cex.axis=2,cex.lab=2,yaxt="n")
axis(2,at=c(0,100,200,300),cex.axis=2)
lines(0:1, rep(50,2),lty=2, lwd = 1)
abline(v=0.0954,col="red")
legend("topright",lty=c(2,2),col=c("red",1),legend=c("Estimate","Prior"),cex=2,bty="n",lwd=2)

##HDPI
hdi(para[[12]][,1], credMass = 0.95)  # credMass = 0.95 表示95%置信区间
median(para[[12]][,2])


library(MASS)
library(ggplot2)
library(grid)
library(gridExtra)

kde <- kde2d(para[[12]][,1],para[[12]][,2], n = 100)
dx <- diff(kde$x[1:2])
dy <- diff(kde$y[1:2])
cell_area <- dx * dy
frequency_matrix <- kde$z * cell_area * n

kde_df <- expand.grid(x = kde$x, y = kde$y)
kde_df$z <- as.vector(frequency_matrix)
kde_df[which.max(kde$z),]

blue_to_yellow <- c("#18244B", "#1D3B72", "#2A5B9D", "#4B7BB7", "#7697C4",
                    "#9DB4D3", "#C8D4E1", "#E5EBD5", "#F5F2AE", "#FCE77F")




p_main <- ggplot(kde_df, aes(x = x, y = y, fill = z)) +
  geom_raster(interpolate = TRUE) +
  geom_vline(xintercept = 0.933, color = "red", size = 0.5,linetype="dashed")+
  geom_hline(yintercept = 0.0954, color = "red", size = 0.5,linetype="dashed")+
  scale_fill_gradientn(
    colours = c("white", blue_to_yellow),
    name = "Frequency"
  ) +
  labs(
    title = "Posterior Density Heatmap",
    x = expression(rho),
    y = expression(p[f])
  ) +
  theme_minimal()+theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 24),       # axis tick labels
    axis.title = element_text(size = 24),      # axis titles
    plot.title = element_blank(),  # title
    #legend.position = "none"
  )+xlim(0, 8)+ylim(0,0.5)+  # Set x-axis range from 0 to 6
  geom_point(aes(x = 0.933, y = 0.0954), size = 2,col="red")

p_main
#kde_df[which.max(kde_df$z),]

df <- as.data.frame(para[[12]])
colnames(df) <- c("rho", "pf")


# Top histogram with prior
p_top <- ggplot(df, aes(x = rho)) +
  geom_histogram(bins = 80, fill = "lightblue", color = "white") +
  geom_line(data = prior_df, aes(x = x, y = y), color = "black", size = 0.5,linetype="dashed") +
  geom_vline(xintercept = 0.933, color = "red", size = 0.5,linetype="dashed")+
  coord_cartesian(xlim = c(0, 8)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.title.y = element_blank(),
    axis.text.y  = element_text(color = NA),  
    axis.ticks.y = element_line(color = NA),
    plot.margin  = margin(0, 5, 0, 5),         
    legend.position = "none"
  )
p_top
# Right histogram (flipped)
p_right <- ggplot(df, aes(x = pf)) +
  geom_histogram(bins = 100, fill = "lightblue", color = "white") +
  geom_line(data = prior_pf, aes(x = x, y = y), color = "black", size = 0.5,linetype="dashed") +
  geom_vline(xintercept = 0.0954, color = "red", size = 0.5,linetype="dashed")+
  coord_flip(xlim=c(0,0.5)) + 
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),   
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),   
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(5, 0, 5, 0)            
  )
# Convert to grobs
g_top   <- ggplotGrob(p_top)
g_main  <- ggplotGrob(p_main)
g_right <- ggplotGrob(p_right)
g_top$widths   <- g_main$widths
g_right$heights <- g_main$heights

# 组合三图
grid.newpage()
grid.layout <- grid.layout(nrow = 2, ncol = 2,
                           widths  = unit(c(4, 1), "null"),
                           heights = unit(c(1, 4), "null"))

pushViewport(viewport(layout = grid.layout))
vplayout <- function(row, col) viewport(layout.pos.row = row, layout.pos.col = col)

pushViewport(vplayout(1, 1)); grid.draw(g_top);    upViewport()
pushViewport(vplayout(2, 1)); grid.draw(g_main);   upViewport()
pushViewport(vplayout(2, 2)); grid.draw(g_right);  upViewport()

