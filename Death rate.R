##Appendix S1

Infected_data <- read.csv("Infected_mosquito.csv")
Uninfected_data <- read.csv("Uninfected_mosquito.csv")
data <- cbind(rbind(Infected_data,Uninfected_data),c(rep(1,17),rep(0,2)))
ft <- -log(data[,2])/data[,1]
data <- cbind(data,ft)
names(data) <- c("Time","Percentage","Infected","Transform")

Exp_nls <- nls(Transform ~ a*exp(b*Time+c*Infected),data = data ,start=list(a=0.1,b=0.1,c=0.1))
summary(Exp_nls)

plot(data[1:17,1],exp(-0.0043140*data[1:17,1]*exp(0.0874910*data[1:17,1]+0.3604109)),type="l",xlab="Time post blood feed (Days)",ylab="Mosquito survival",ylim=c(0,1),yaxt="n",cex.axis=1.2,cex.lab=1.2,cex.main=1.2,col="red")
axis(2,at=c(0,0.2,0.4,0.6,0.8,1) ,lab=paste0(c(0,0.2,0.4,0.6,0.8,1) * 100,"%"),las=TRUE,cex.axis=1.2)
points(data[1:17,1],data[1:17,2],col="red")
lines(data[1:17,1],exp(-0.0043140*data[1:17,1]*exp(0.0874910*data[1:17,1])),col="blue")
points(data[18:19,1],data[18:19,2],col="blue")
legend("topright",col=c("blue","red"),lty=1,legend=c("Uninfected","Infected"))
