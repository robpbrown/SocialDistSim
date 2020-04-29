#!/usr/local/bin/Rscript

source ("seir_kissler.r")

startstate = c(1,0,0,0,0,0,0)
startt = 70
xd=seq(as.Date("2020-1-1"),length=100,by="3 months")

# Policy parameters controls the effects of SD
# Assuming SD reduces transmission by 0.6
out = evolve (startstate=startstate, startt = startt, stopt = 1000, policyparameters = list(distancingeffect = 0.6, hi = 37.5/1e4, low = 10/1e4), deltat=1/10)
pdf ("social-distancing.60.pdf",width=8,height=8)
par(mfrow=c(2,1))
plot (out$state[,3]*10000,ylab="Infections per 10K",type="l",col="red",xlab="Days",xaxt="n", ylim=c(0,60))
int=seq(1,nrow(out$states),by=900)
axis(1,at=int,labels=xd[1:length(int)])
plot (out$state[,6]*10000,ylab="Critical cases per 10K",type="l",col="blue",xlab="Days",xaxt="n",ylim=c(0,1.2))
int=seq(1,nrow(out$states),by=900)
axis(1,at=int,labels=xd[1:length(int)])
dev.off()

# Assuming SD reduces transmission by 0.9
out = evolve (startstate=startstate, startt = startt, stopt = 1000, policyparameters = list(distancingeffect = 0.9, hi = 37.5/1e4, low = 10/1e4), deltat=1/10)
pdf ("social-distancing.90.pdf",width=8,height=8)
par(mfrow=c(2,1))
plot (out$state[,3]*10000,ylab="Infections per 10K",type="l",col="red",xlab="Days",xaxt="n", ylim=c(0,60))
int=seq(1,nrow(out$states),by=900)
axis(1,at=int,labels=xd[1:length(int)])
plot (out$state[,6]*10000,ylab="Critical cases per 10K",type="l",col="blue",xlab="Days",xaxt="n",ylim=c(0,1.2))
int=seq(1,nrow(out$states),by=900)
axis(1,at=int,labels=xd[1:length(int)])

dev.off()

