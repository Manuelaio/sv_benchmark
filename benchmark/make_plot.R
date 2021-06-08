data=read.table("all_result.txt", header=T, sep="\t", dec=".")

del = data[grepl("*DEL", data$SV_TYPE), ]
#del=del[-c(5),]
dup = data[grepl("*DUP", data$SV_TYPE), ]

del$Precision= (del$Precision)/100
del$Recall=(del$Recall)/100
dup$Precision= (dup$Precision)/100
dup$Recall=(dup$Recall)/100


#par(mfrow=c(1,2))

#Mycol<-c("darkorange", "dodgerblue")
Mycol<-c("purple","darkorange","dodgerblue", "chartreuse","magenta")
#Mycol<-c("purple","darkorange","dodgerblue", "magenta")
xlab="Precision"
ylab<-"Recall"
for (i in unique(del$SIZE)) { 
plot(del[del$SIZE==i,]$Precision,del[del$SIZE==i,]$Recall,xlim=c(0,1),ylim=c(0,1),col=Mycol,pch=19,main="Precision-Recall  DELETIONS",cex=1.5,xlab=xlab,ylab=ylab)
legend("topleft", legend=c("DELLY", "ERDS", "lumpy", "Manta", "SVABA") ,
       col=Mycol, pch=19, box.lwd=par(1))

P<-seq(0.05,2,by=0.001)
FVec<-seq(0.1,0.9,by=0.1)
for (i in 1:9)
{
  
  F<-FVec[i]
  R<-F*P/(2*P-F)
  if (i!=9 & i!=7)
  {
    lines(P,R,col="gray80")
  }
  if (i==9)
  {
    lines(P[700:1000],R[700:1000],col="gray80")
  }
  if (i==7)
  {
    lines(P[400:1000],R[400:1000],col="gray80")
  }
}
text(0.05,0.1,"f=0.1",col="gray80")
text(0.15,0.2,"f=0.2",col="gray80")
text(0.25,0.3,"f=0.3",col="gray80")
text(0.35,0.4,"f=0.4",col="gray80")
text(0.45,0.5,"f=0.5",col="gray80")
text(0.55,0.6,"f=0.6",col="gray80")
text(0.65,0.7,"f=0.7",col="gray80")
text(0.75,0.8,"f=0.8",col="gray80")
text(0.85,0.9,"f=0.9",col="gray80")}
#############################
#######DEL#####à

for (i in unique(dup$SIZE)) { 
plot(dup[dup$SIZE==i,]$Precision,dup[dup$SIZE==i,]$Recall,xlim=c(0,1),ylim=c(0,1),col=Mycol,pch=19,main="Precision-Recall DUPLICATIONS",cex=1.5,xlab=xlab,ylab=ylab)
legend("bottomright", legend=c("DELLY", "ERDS", "lumpy", "Manta", "SVABA"),
       col=Mycol, pch=19, box.lwd=par(1))

P<-seq(0.05,2,by=0.001)
FVec<-seq(0.1,0.9,by=0.1)
for (i in 1:9)
{
  
  F<-FVec[i]
  R<-F*P/(2*P-F)
  if (i!=9 & i!=7)
  {
    lines(P,R,col="gray80")
  }
  if (i==9)
  {
    lines(P[700:1000],R[700:1000],col="gray80")
  }
  if (i==7)
  {
    lines(P[400:1000],R[400:1000],col="gray80")
  }
}
text(0.05,0.1,"f=0.1",col="gray80")
text(0.15,0.2,"f=0.2",col="gray80")
text(0.25,0.3,"f=0.3",col="gray80")
text(0.35,0.4,"f=0.4",col="gray80")
text(0.45,0.5,"f=0.5",col="gray80")
text(0.55,0.6,"f=0.6",col="gray80")
text(0.65,0.7,"f=0.7",col="gray80")
text(0.75,0.8,"f=0.8",col="gray80")
text(0.85,0.9,"f=0.9",col="gray80")}



library(ggplot2)
library(tidyr)
library(rlist)


del_reshape= del %>% gather (Qual, tot, TP,FP,FN)

for (i in unique(del_reshape$SIZE)) { 
print(ggplot(del_reshape[del_reshape$SIZE==i,], aes(x=Caller, y= tot, fill= Qual)) +  geom_bar(stat = "identity",position = "dodge") + ggtitle(i))
}





dup_reshape= dup %>% gather (Qual, tot, TP,FP,FN)

for (i in unique(dup_reshape$SIZE)) { 
  print(ggplot(dup_reshape[dup_reshape$SIZE==i,], aes(x=Caller, y= tot, fill= Qual)) +  geom_bar(stat = "identity",position = "dodge") + ggtitle("Dup ",i))
}




