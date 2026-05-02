###Simulation of asymptomatic infection
source("Models.R")
t_total <- 2880
par(mar = c(5.1, 7, 2, 1.1))

para <- 10^6.5  #4.5 for low transmission settings 6.5 for high transmission settings
benchmark <- c(0.08,7,8,71,0.006,0.03,0.0004,0.007,0.02,0.2,3.4, 0.8, 0.029,13.8, 3.9, 19.3, 2.1)
high_trans <- try(sto_model((benchmark[1])*4000,benchmark[2],benchmark[3],benchmark[4],matrix(0,2880,3),benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],al=42,as=25,t_total,10^6.5),silent = TRUE)
plot(1:2880/24,rowSums(high_trans[,1:24])/4000,type="l",xlab = "Time post infection (Days)",ylab="Parasitemia /ml blood",log="y",ylim=c(10,1e+8),col = "#B22222",cex.axis=2,cex.lab=2,cex.main=2,yaxt="n",xlim=c(0,80),bty="l",xaxt="n")
axis(2,at=c(10^2,10^4,10^6,10^8),label=expression(10^2,10^4,10^6,10^8),cex.axis=2)
axis(1,at=c(0,20,40,60,80),cex.axis=2)

low_trans <- try(sto_model((benchmark[1])*4000,benchmark[2],benchmark[3],benchmark[4],matrix(0,2880,3),benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],al=42,as=25,t_total,10^4.5),silent = TRUE)
lines(1:2880/24,rowSums(low_trans[,1:24])/4000 ,col = "#0072B2")
legend("topleft",col=c("#B22222","#0072B2"),legend=c("High","Low"),lty=1,lwd=1.5,bty="n",cex=1.5)

high_prob <- vector()
low_prob <- vector()
high_male <- as.numeric(high_trans[,88])/4000
high_female <- as.numeric(high_trans[,89])/4000
low_male <- as.numeric(low_trans[,88])/4000
low_female <- as.numeric(low_trans[,89])/4000
for (i in 1:120) {
  oocyst_record <- oocyst_prev(100000,round(high_male[i*24]),round(high_female[i*24]),benchmark[11],benchmark[12],benchmark[13],benchmark[14],benchmark[15],benchmark[16],benchmark[17])
  high_prob[i] <- oocyst_record[8.5*24]
  oocyst_record <- oocyst_prev(100000,round(low_male[i*24]),round(low_female[i*24]),benchmark[11],benchmark[12],benchmark[13],benchmark[14],benchmark[15],benchmark[16],benchmark[17])
  low_prob[i] <- oocyst_record[8.5*24]
  print(i)
}

plot(1:120,high_prob,type="l",log="y",xlim=c(0,80),ylim=c(0.0005,1),xlab = "Time post infection (Days)",ylab = "Human-to-mosquito\n transmission probability",col = "#B22222",cex.axis=2,cex.lab=2,cex.main=2,yaxt="n",bty="l",xaxt="n")
axis(2,at=c(0.001,0.01,0.1,1),cex.axis=1.8)
axis(1,at=c(0,20,40,60,80),cex.axis=2)
lines(1:120,low_prob,col ="#0072B2")


legend(60,0.15,col=c("#B22222","#0072B2"),legend=c("High","Low"),lty=1,lwd=2,bty="n",cex=1.5)

segments(x0=-5, y0=0.001,
         x1=18.4, y1=0.001,
         lty=2,
         lwd=1,
         col="black")
segments(x0=18.4, y0=0.00001,
         x1=18.4, y1=0.001,
         lty=2,
         lwd=1,
         col="black")

segments(x0=-5, y0=low_prob[60],
         x1=60, y1=low_prob[60],
         lty=2,
         lwd=1,
         col="black")
segments(x0=60, y0=0.00001,
         x1=60, y1=low_prob[60],
         lty=2,
         lwd=1,
         col="black")



segments(x0=12.8, y0=0.00001,
         x1=12.8, y1=0.001,
         lty=2,
         lwd=1,
         col="black")

segments(x0=-5, y0=high_prob[60],
         x1=60, y1=high_prob[60],
         lty=2,
         lwd=1,
         col="black")
segments(x0=60, y0=0.00001,
         x1=60, y1=high_prob[60],
         lty=2,
         lwd=1,
         col="black")



