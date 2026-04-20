#Within-human model
sto_model <- function(N,mu,sd,rp,kd,f,detp,detg,detgm,detgf,m,fg,al=42,as=25,t,thres){
  
  pinit <- rtruncnorm(max(round(N,0),1),1,al,mu,sd) 
  pinit <- round(pinit,digits = 0)
  rec <- matrix(0,t,(2*al+5))
  for (i in 1:al) { #record initial number of parasites
    rec[1,i] <- length(which(pinit==i)) 
  }
  
  rec_f <- vector()
  rec_f[1] <- 0
  rec_fg <- vector()
  rec_fg[1] <- 0
  
  para_hist <- vector()
  para_hist[1] <- sum(rec[1,1:42])/4000
  sw1 <- FALSE
  sw2 <- FALSE
  sw3 <- FALSE
  t_sw1 <- NA_integer_
  t_sw2 <- NA_integer_
  
  for (i in 2:t) { 
    for (j in 1:al) { #first compute asexual parasites and in this part drug is effective
      if (j==1){ #age = 1 hour
        # number of parasites in state 48 last time * replication rate - number of death (detp--death rate ,kd--killing effect )
        # Here kd will be explained in the PQP file.
        death <- rpois(1,(detp+kd[i])*rec[i-1,al])
        
        p.a.t <- (rec[i-1,al]-death)*rp
        p.a.t <- max(0,round(p.a.t,0))
        rec[i,j] <- p.a.t
      }
      else{if (j == as+1){
        death <- rpois(1,(detp+kd[i])*rec[i-1,as])
        
        p.a.t <- rec[i-1,as]-death
        p.a.t <- max(0,p.a.t)
        
        rec_f[i] <- rbinom(1,p.a.t,f)
        rec[i,j] <- p.a.t-rec_f[i]
      }
        else{
          death <- rpois(1,(detp+kd[i])*rec[i-1,j-1])
          
          p.a.t <- rec[i-1,j-1]-death
          p.a.t <- max(0,p.a.t)
          
          rec[i,j] <- p.a.t
          
        }
      }
    }
    for (j in (al+1):(2*al-1)){ # 47 states for sexual development
      #column al+1 is Pg(as+1), al+2 is Pg(as+2), ... 
      if (j==(2*al-as+1)){ # age 1 of the sexually committed parasites
        # replication 
        p.a.t <- (rec[i-1,2*al-as]-rpois(1,detp*rec[i-1,2*al-as]))*rp
        p.a.t <- max(0,round(p.a.t,0))
        rec[i,j] <- p.a.t
      }else{if (j==(al+1)){ 
        #fraction
        p.a.t <- rec_f[i]-rpois(1,detp*rec_f[i])
        p.a.t <- max(0,p.a.t)
        rec[i,j] <- p.a.t
      }
        else{
          p.a.t <- rec[i-1,j-1]-rpois(1,detp*rec[i-1,j-1])
          p.a.t <- max(0,p.a.t)
          rec[i,j] <- p.a.t
        }
      }
    }
    # j in (2*al):(2*al+4) G1 to G5
    #random a number to be parasites going from g1 to g2, compare with number in g1
    #to make sure g1 has enough parasites to go to g2.
    #the maturation rate is m
    g1tog2 <- min(rec[i-1,2*al],rpois(1,m*rec[i-1,2*al]))
    #all the parasites in (Pg(as-1) - nature death (death rate is detp)) will go to state g1 
    #+ parasites in g1 last time - nature death (death rate is detg)
    #- parasites goes to g2
    #= number in g1 state
    rec[i,2*al] <- max(0,rec[i-1,2*al]+rec[i-1,2*al-1]-rpois(1,detg*rec[i-1,2*al])-rpois(1,detp*rec[i-1,2*al-1])-g1tog2)
    
    #next g2, g3, g4 
    g2tog3 <- min(rec[i-1,2*al+1],rpois(1,m*rec[i-1,2*al+1]))
    rec[i,2*al+1] <- max(0,rec[i-1,2*al+1]+g1tog2-g2tog3-rpois(1,detg*rec[i-1,2*al+1]))
    
    g3tog4 <- min(rec[i-1,2*al+2],rpois(1,m*rec[i-1,2*al+2]))
    rec[i,2*al+2] <- max(0,rec[i-1,2*al+2]+g2tog3-g3tog4-rpois(1,detg*rec[i-1,2*al+2]))
    
    g4tog5 <- min(rec[i-1,2*al+3],rpois(1,m*rec[i-1,2*al+3]))
    rec[i,2*al+3] <- max(0,rec[i-1,2*al+3]+g3tog4-g4tog5-rpois(1,detg*rec[i-1,2*al+3]))
    
    rec_fg[i] <- rbinom(1, g4tog5,fg)
    rec[i,2*al+4] <- max(0,rec[i-1,2*al+4]+rec_fg[i]-rpois(1,detgm*(rec[i-1,2*al+4]+rec_fg[i])))
    rec[i,2*al+5] <- max(0,rec[i-1,2*al+5]+g4tog5-rec_fg[i]-rpois(1,detgf*(rec[i-1,2*al+5]+g4tog5-rec_fg[i])))
    if (anyNA(rec[i,])==TRUE){
      rec=rec[-(i:t),]
      break
    }
    para_hist[i] <- sum(rec[i, 1:42])/4000
    
    
    if (!sw1 && para_hist[i] >= (thres / (rp / 4))) {
      PMF_target <- thres / para_hist[i]
      rp <- PMF_target / (1 - detp)^al / (1 - f)
      sw1 <- TRUE
      t_sw1 <- i
    }
    
    
    if (sw1 && !sw2 && !is.na(t_sw1) && i >= t_sw1 + al) {
      PMF_target <- thres / para_hist[i]
      rp <- PMF_target / (1 - detp)^al / (1 - f)
      sw2 <- TRUE
      t_sw2 <- i
    }
    
    
    if (sw1 && sw2 && !sw3 && !is.na(t_sw2) && i >= t_sw2 + al) {
      PMF_target <- 1
      rp <- PMF_target / (1 - detp)^al / (1 - f)
      sw3 <- TRUE
    }
    
  }
  
  rec <- data.frame(rec)
  return(rec)
}

##mosquito feeding model and simulating the direct feeding assays
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

###Calculating the probability of human-to-mosquito transmission
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