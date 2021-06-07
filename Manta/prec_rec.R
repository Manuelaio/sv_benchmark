library(rlist)

#####Rscript 
args = commandArgs(trailingOnly=TRUE)
vcf<-  args[1]
goldset<- args[2]

####PRECISION AND RECALL FUNCTION#############

PR_del <- function(data, data2){
  file.DEL= subset(data, type=="DEL")
  intersect.del= file.DEL[file.DEL$intersect!=0,]
  True.p= intersect.del[which(intersect.del$type==intersect.del$type_T),]
  TP_no_dup=True.p[!duplicated(True.p[1:4]),]
  discordanti= intersect.del[which(intersect.del$type!=intersect.del$type_T),]
  discordanti.na.omit= na.omit(discordanti)
  truset.del= subset(data2, V4=="DEL")
  TP= nrow(TP_no_dup)
  FN= (nrow(truset.del)) - TP - (nrow(discordanti.na.omit))
  FP= sum((nrow(file.DEL[file.DEL$intersect==0,])), (nrow(discordanti.na.omit)))
  PRECISIONE= ((TP)/(TP+FP))*100
  
  RECALL=  ((TP)/(TP+ abs(FN)))*100
  f1 = 2 * PRECISIONE * RECALL / (PRECISIONE + RECALL)
  df <- data.frame(PRECISIONE,RECALL,f1, TP, FN, FP, length(file.DEL$type))
  colnames(df) <- c("Precision","Recall","FCoef","TP", "FN", "FP", "number")
  return(df)
  
}

PR_dup <- function(data, data2){
  file.DUP= subset(data, type=="DUP")
  intersect.dup= file.DUP[file.DUP$intersect!=0,]
  True.p=intersect.dup[intersect.dup$type_T=="DUP" | intersect.dup$type_T=="INS",]
  discordanti_dup= intersect.dup[intersect.dup$type_T=='DEL',]
  truset.dup= subset(data2, V4=="DUP")
  TP= nrow(True.p[!duplicated(True.p[1:4]),])
  FN= (nrow(truset.dup)) - TP - (nrow(discordanti_dup))
  FP= sum((nrow(file.DUP[file.DUP$intersect==0,])), (nrow(discordanti_dup)))
  PRECISIONE= ((TP)/(TP+FP))*100
  RECALL=  ((TP)/(TP+ abs(FN)))*100
  f1 = 2 * PRECISIONE * RECALL / (PRECISIONE + RECALL)
  df <- data.frame(PRECISIONE,RECALL,f1, TP, FN, FP, length(file.DUP$type))
  colnames(df) <- c("Precision","Recall","FCoef","TP", "FN", "FP", "number")
  return(df)
  return(df)
}

###################################################################


Manta= read.table(vcf, stringsAsFactors = F)

trueset=read.table(goldset, stringsAsFactors = F)
# sv SS 50 - 100 bp
# sv S 100 - 1 kb
# sv M 1 - 50 Kb
# sv XM 50 KB - 100 Kb 
# sv L > 100Kb 


colnames(Manta)= c("chr", "start", "end", "type","manta_len", "geno","chr_T", "start_T", "end_T", "type_T", "source", "length_T", "intersect")

Manta$manta_len= as.numeric(Manta$manta_len)
Manta$length_T= as.numeric(Manta$length_T)



############################### DELETION BY SIZE ######################

SS= subset(Manta, manta_len>50 & manta_len<=100)
S= subset(Manta, manta_len>100 & manta_len<=1000)
M= subset(Manta, manta_len>1000 & manta_len<=50000)
XM= subset(Manta, manta_len>50000 & manta_len<=100000)
L=subset(Manta, manta_len>100000)
fino_50KB= subset(Manta, manta_len>50 & manta_len<=50000)
#XLS=subset(all, manta_len <50 | manta_len>1000000)

truest.del.SS= trueset[trueset$V6>= 50 & trueset$V6<=100,]
truest.del.S= trueset[trueset$V6> 100 & trueset$V6<=1000,]
truest.del.M= trueset[trueset$V6> 1000 & trueset$V6<=50000,]
truest.del.XM= trueset[trueset$V6> 50000 & trueset$V6<=100000,]
truest.del.L= trueset[trueset$V6> 100000 ,]
trueste.XLS=trueset[trueset$V6<50 | trueset$V6>1000000,]
truset_fino_50= trueset[trueset$V6> 50 & trueset$V6<=50000,]


sink('analysis-output.txt')



cat("=============================\n")

cat("Precision and recall  deletion\n")
cat("=============================\n")
cat("all\n")
PR_del(Manta, trueset)
cat("=============================\n")
cat("SS\n")
PR_del(SS, truest.del.SS)
cat("=============================\n")
cat("S\n")
PR_del(S, truest.del.S)
cat("=============================\n")
cat("M\n")
PR_del(M, truest.del.M)
cat("=============================\n")
cat("XM\n")
PR_del(XM, truest.del.XM)
cat("=============================\n")
cat("L\n")
PR_del(L, truest.del.L)
cat("=========from 50 bp to 50 KB=============\n")
PR_del(fino_50KB, truset_fino_50)


cat("=============================\n")
cat("Precision and recall  duplication\n")
cat("=============================\n")
cat("all\n")
PR_dup(Manta, trueset)
cat("=============================\n")
cat("SS\n")
PR_dup(SS, truest.del.SS)
cat("=============================\n")
cat("S\n")
PR_dup(S, truest.del.S)
cat("=============================\n")
cat("M\n")
PR_dup(M, truest.del.M)
cat("=============================\n")
cat("XM\n")

PR_dup(XM, truest.del.XM)
cat("=============================\n")
cat("L\n")
PR_dup(L, truest.del.L)

cat("=============================\n")

cat("=========from 50 bp to 50 KB=============\n")
PR_dup(fino_50KB, truset_fino_50)


cat("==========End================\n")




# Stop writing to the file
sink()
