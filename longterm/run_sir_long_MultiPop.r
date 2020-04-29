#!/usr/local/bin/Rscript
rm(list=ls())
source ("seir_kissler_MultiPop.r")
library(ggplot2)
library(reshape2)

#startstate = c(1,0,0,0,0,0,0)
#startstate=matrix(c(.3,.7,.2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=3,ncol=7)
startt = 70
xd=seq(as.Date("2020-1-1"),length=100,by="3 months")

# Policy parameters controls the effects of SD
# Assuming SD reduces transmission by 0.6
#out = evolve (startstate=startstate, startt = startt, stopt = 1000, policyparameters =list(distancingeffect = 0.6, hi = 37.5/1e4, low = 10/1e4, Qdur=14, NQdur=7, stagger=7), deltat=1/10)

startstate=matrix(c(0.33,0.33,0.34,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=3,ncol=7)
outModel06Q = evolve (startstate=startstate, startt = startt, stopt = 1000, policyparameters =list(distancingeffect = sqrt(0.6), hi = 37.5/1e4, low = 10/1e4, Qdur=14, NQdur=7, stagger=7), deltat=1/10)
outModel08Q = evolve (startstate=startstate, startt = startt, stopt = 1000, policyparameters =list(distancingeffect = sqrt(0.8), hi = 37.5/1e4, low = 10/1e4, Qdur=14, NQdur=7, stagger=7), deltat=1/10)
outModel10Q = evolve (startstate=startstate, startt = startt, stopt = 1000, policyparameters =list(distancingeffect = sqrt(1.0), hi = 37.5/1e4, low = 10/1e4, Qdur=14, NQdur=7, stagger=7), deltat=1/10)

startstate=matrix(c(1.00,0.0,0.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nrow=3,ncol=7)
outModel06QONEPOP = evolve (startstate=startstate, startt = startt, stopt = 1000, policyparameters =list(distancingeffect = sqrt(0.6), hi = 37.5/1e4, low = 10/1e4, Qdur=14, NQdur=7, stagger=7), deltat=1/10)
outModel10QONEPOP = evolve (startstate=startstate, startt = startt, stopt = 1000, policyparameters =list(distancingeffect = sqrt(1.0), hi = 37.5/1e4, low = 10/1e4, Qdur=14, NQdur=7, stagger=7), deltat=1/10)





####### plot outModel06Q ###########
dfplot=data.frame(outModel06Q$states)
dfplot$times=outModel06Q$times
rownames(dfplot)=c(1:nrow(dfplot))
dfplot=melt(dfplot,id.vars = 'times')
dfplot$vartype=substring(dfplot$variable, first = 1, last = 2)
head(dfplot)
dftmp=dfplot
head(dftmp)
#dfplotCUM=aggregate( dftmp$value, dftmp[,c(1,3)], FUN = sum )

dffilt=dftmp[which(dftmp$vartype=='I.'),]
pdf ("Var.I.06.MultiPop.pdf",width=8,height=8)
p=ggplot(dffilt, aes(x=times, y=value, color=variable))+geom_point()
print(p)
dev.off()

dffilt=dftmp[which(dftmp$vartype=='S.'),]
pdf ("Var.S.06.MultiPop.pdf",width=8,height=8)
p=ggplot(dffilt, aes(x=times, y=value, color=variable))+geom_point()
print(p)
dev.off()

dffilt=dftmp[which(dftmp$vartype=='R.'),]
pdf ("Var.R.06.MultiPop.pdf",width=8,height=8)
p=ggplot(dffilt, aes(x=times, y=value, color=variable))+geom_point()
print(p)
dev.off()

### THis shows where the I are

tmp=dfplot[,c(1,3,4)]
head(tmp)
dfplotCUM=aggregate( tmp$value, tmp[,c(1,3)], FUN = sum )
head(dfplotCUM)
head(dfplot)
head(dfplotCUM)
tmp=merge(dfplot,dfplotCUM,by=c('times','vartype'))
head(tmp)
colnames(tmp)[5]='TotalVar'
tmp$PropOf=tmp$value/tmp$TotalVar
head(tmp)

dffilt=tmp[which(tmp$vartype=='I.'),]
head(dffilt)
pdf ("../Graphs/Var.Proportion.I.06.MultiPop.pdf",width=8,height=8)
p=ggplot(dffilt, aes(x=times, y=PropOf, color=variable))+geom_point()
print(p)
dev.off()



ALLPLOT=data.frame(times=numeric(),vartype=character(),x=numeric(),Model=character())

####### plot outModel06Q ###########
dfplot=data.frame(outModel06Q$states)
dfplot$times=outModel06Q$times
rownames(dfplot)=c(1:nrow(dfplot))
dfplot=melt(dfplot,id.vars = 'times')
dfplot$vartype=substring(dfplot$variable, first = 1, last = 2)
head(dfplot)
dftmp=dfplot[,c(1,3,4)]
head(dftmp)
dfplotCUM=aggregate( dftmp$value, dftmp[,c(1,3)], FUN = sum )
pdf ("social-distancing.06.MultiPop.pdf",width=8,height=8)
p=ggplot(dfplotCUM, aes(x=times, y=x, color=vartype))+geom_point()
print(p)
dev.off()

outline=data.frame(times=dfplotCUM$times,vartype=dfplotCUM$vartype,x=dfplotCUM$x, Model='MultiPop06')
ALLPLOT=rbind(ALLPLOT,outline)

####### plot outModel08Q ###########
dfplot=data.frame(outModel08Q$states)
dfplot$times=outModel08Q$times
rownames(dfplot)=c(1:nrow(dfplot))
dfplot=melt(dfplot,id.vars = 'times')
dfplot$vartype=substring(dfplot$variable, first = 1, last = 2)
head(dfplot)
dftmp=dfplot[,c(1,3,4)]
head(dftmp)
dfplotCUM=aggregate( dftmp$value, dftmp[,c(1,3)], FUN = sum )
pdf ("social-distancing.08.MultiPop.pdf",width=8,height=8)
p=ggplot(dfplotCUM, aes(x=times, y=x, color=vartype))+geom_point()
print(p)
dev.off()

outline=data.frame(times=dfplotCUM$times,vartype=dfplotCUM$vartype,x=dfplotCUM$x, Model='MultiPop08')
ALLPLOT=rbind(ALLPLOT,outline)



####### plot outModel10Q ###########
dfplot=data.frame(outModel10Q$states)
dfplot$times=outModel10Q$times
rownames(dfplot)=c(1:nrow(dfplot))
dfplot=melt(dfplot,id.vars = 'times')
dfplot$vartype=substring(dfplot$variable, first = 1, last = 2)
head(dfplot)
dftmp=dfplot[,c(1,3,4)]
head(dftmp)
dfplotCUM=aggregate( dftmp$value, dftmp[,c(1,3)], FUN = sum )
pdf ("social-distancing.10.MultiPop.pdf",width=8,height=8)
p=ggplot(dfplotCUM, aes(x=times, y=x, color=vartype))+geom_point()
print(p)
dev.off()

outline=data.frame(times=dfplotCUM$times,vartype=dfplotCUM$vartype,x=dfplotCUM$x, Model='MultiPop10')
ALLPLOT=rbind(ALLPLOT,outline)

####### plot outModel06QONEPOP ###########
dfplot=data.frame(outModel06QONEPOP$states)
dfplot$times=outModel06QONEPOP$times
rownames(dfplot)=c(1:nrow(dfplot))
dfplot=melt(dfplot,id.vars = 'times')
dfplot$vartype=substring(dfplot$variable, first = 1, last = 2)
head(dfplot)
dftmp=dfplot[,c(1,3,4)]
head(dftmp)
dfplotCUM=aggregate( dftmp$value, dftmp[,c(1,3)], FUN = sum )
pdf ("social-distancing.06.ONEPOP.pdf",width=8,height=8)
p=ggplot(dfplotCUM, aes(x=times, y=x, color=vartype))+geom_point()
print(p)
dev.off()

outline=data.frame(times=dfplotCUM$times,vartype=dfplotCUM$vartype,x=dfplotCUM$x, Model='OnePop06')
ALLPLOT=rbind(ALLPLOT,outline)

####### plot outModel10QONEPOP ###########
dfplot=data.frame(outModel10Q$states)
dfplot$times=outModel10Q$times
rownames(dfplot)=c(1:nrow(dfplot))
dfplot=melt(dfplot,id.vars = 'times')
dfplot$vartype=substring(dfplot$variable, first = 1, last = 2)
head(dfplot)
dftmp=dfplot[,c(1,3,4)]
head(dftmp)
dfplotCUM=aggregate( dftmp$value, dftmp[,c(1,3)], FUN = sum )
pdf ("social-distancing.10.ONEPOP.pdf",width=8,height=8)
p=ggplot(dfplotCUM, aes(x=times, y=x, color=vartype))+geom_point()
print(p)
dev.off()

outline=data.frame(times=dfplotCUM$times,vartype=dfplotCUM$vartype,x=dfplotCUM$x, Model='OnePop10')
ALLPLOT=rbind(ALLPLOT,outline)

ALLPLOT$x10000=ALLPLOT$x*10000
for (var in unique(ALLPLOT$vartype)){
ALLPLOTI=ALLPLOT[which(ALLPLOT$vartype==var),]
pdf(paste("../Graphs/ValuesAcrossModels",var,".pdf",sep=''),width=8,height=8)
p=ggplot(ALLPLOTI, aes(x=times, y=x10000, color=Model))+geom_point()+ggtitle(var)
print(p)
dev.off()
}

