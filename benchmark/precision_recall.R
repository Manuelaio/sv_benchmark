library(rlist)

#####Rscript 
args = commandArgs(trailingOnly=TRUE)
vcf<-  args[1]
goldset<- args[2]
output<- args[3]

if (length(args)==0) {
  stop("Two arguments must be supplied: caller bed file intersected with truset and goldset /work/emanuela.iovino/intersect_SV/new_analysis//DGV_paiper/NA12878_dgv_long_read.bed.n", call.=FALSE)
} 

if (length(args[2])==0) {
  stop("Gold set must be supplied:/work/emanuela.iovino/intersect_SV/new_analysis//DGV_paiper/NA12878_dgv_long_read.bed.n", call.=FALSE)
} 

####PRECISION AND RECALL FUNCTION#############

PR_del <- function(data, data2){
  file.DEL= subset(data, type=="DEL")
  intersect.del= file.DEL[file.DEL$intersect!=0,]
  True.p= intersect.del[which(intersect.del$type==intersect.del$type_T),]
  TP_no_dup=True.p[!duplicated(True.p[1:4]),]
  seg_dup_tp= nrow(!is.na(TP_no_dup[!is.na(TP_no_dup$sg),]))
  discordanti= intersect.del[which(intersect.del$type!=intersect.del$type_T),]
  discordanti.na.omit= na.omit(discordanti)
  truset.del= subset(data2, V4=="DEL")
  colnames(truset.del)= c("chr_T", "start_T", "end_T", "type_T", "source", "length_T", "SegG","RepR")
  fn= dplyr::anti_join(truset.del,TP_no_dup, all=T)
  seg_dup_fn=nrow(!is.na(fn[!is.na(fn$SegG),]))
  TP= nrow(TP_no_dup)
  FN= (nrow(fn))- (nrow(discordanti.na.omit))
  FP= sum((nrow(file.DEL[file.DEL$intersect==0,])), (nrow(discordanti.na.omit)))
  PRECISIONE= ((TP)/(TP+FP))*100
  RECALL=  ((TP)/(TP+ abs(FN)))*100
  f1 = 2 * PRECISIONE * RECALL / (PRECISIONE + RECALL)
  df <- data.frame(PRECISIONE,RECALL,f1, TP, FN, FP,seg_dup_tp,seg_dup_fn )
  colnames(df) <- c("Precision","Recall","FCoef","TP", "FN", "FP", "number_SD_TP", "number_SD_FN")
  return(df)
  
}

PR_dup <- function(data, data2){
  file.DUP= subset(data, type=="DUP")
  intersect.dup= file.DUP[file.DUP$intersect!=0,]
  True.p=intersect.dup[intersect.dup$type_T=="DUP" | intersect.dup$type_T=="INS",]
  seg_dup_tp= nrow(!is.na(True.p[!is.na(True.p$sg),]))
  discordanti_dup= intersect.dup[intersect.dup$type_T=='DEL',]
  truset.dup= subset(data2, V4=="DUP" )
  colnames(truset.dup)= c("chr_T", "start_T", "end_T", "type_T", "source", "length_T", "SegG","RepR")
  TP= nrow(True.p[!duplicated(True.p[1:4]),])
  TP_no_dup= True.p[!duplicated(True.p[1:4]),]
  fn= dplyr::anti_join(truset.dup,TP_no_dup, all=T)
  seg_dup_fn=nrow(!is.na(fn[!is.na(fn$SegG),]))
  #FN= (nrow(truset.dup)) - TP -  (nrow(discordanti_dup))
  FN=nrow(fn)
  FP= sum((nrow(file.DUP[file.DUP$intersect==0,])), (nrow(discordanti_dup)))
  PRECISIONE= ((TP)/(TP+FP))*100
  RECALL=  ((TP)/(TP+ abs(FN)))*100
  f1 = 2 * PRECISIONE * RECALL / (PRECISIONE + RECALL)
  df <- data.frame(PRECISIONE,RECALL,f1, TP, FN, FP,seg_dup_tp,seg_dup_fn )
  colnames(df) <- c("Precision","Recall","FCoef","TP", "FN", "FP", "number_SD_TP", "number_SD_FN")
  return(df)
  return(df)
}


###################################################################


caller= read.table(vcf, stringsAsFactors = F)

trueset=read.table(goldset, stringsAsFactors = F)
# sv SS 50 - 100 bp
# sv S 100 - 1 kb
# sv M 1 - 100 Kb
# sv L 100 -1 Mb


colnames(caller)= c("chr", "start", "end", "type","len", "geno","sg","rr","chr_T", "start_T", "end_T", "type_T", "source", "length_T","SG","RR" ,"intersect")

caller$len= as.numeric(caller$len)
caller$length_T= as.numeric(caller$length_T)
caller$start_T=as.numeric(caller$start_T)



############################### DELETION BY SIZE ######################

SS= subset(caller, len>50 & len<=100)
S= subset(caller, len>100 & len<=1000)
M= subset(caller, len>1000 & len<=50000)
XM= subset(caller, len>50000 & len<=100000)
L=subset(caller, len>100000)
fino_50KB= subset(caller, len>50 & len<=50000)
#XLS=subset(all, len <50 | len>1000000)

truest.del.SS= trueset[trueset$V6>= 50 & trueset$V6<=100,]
truest.del.S= trueset[trueset$V6> 100 & trueset$V6<=1000,]
truest.del.M= trueset[trueset$V6> 1000 & trueset$V6<=50000,]
truest.del.XM= trueset[trueset$V6> 50000 & trueset$V6<=100000,]
truest.del.L= trueset[trueset$V6> 100000 ,]
trueste.XLS=trueset[trueset$V6<50 | trueset$V6>1000000,]
truset_fino_50= trueset[trueset$V6> 50 & trueset$V6<=50000,]

##################################################################
caller2= caller[is.na(caller$sg),]
trueset2=trueset[is.na(trueset$V7),]
SS_no= subset(caller2, len>50 & len<=100)
S_no= subset(caller2, len>100 & len<=1000)
M_no= subset(caller2, len>1000 & len<=50000)
XM_no= subset(caller2, len>50000 & len<=100000)
L_no=subset(caller2, len>100000)
fino_50KB_no= subset(caller2, len>50 & len<=50000)
#XLS=subset(all, len <50 | len>1000000)

truest.del.SS_no= trueset2[trueset2$V6>= 50 & trueset2$V6<=100,]
truest.del.S_no= trueset2[trueset2$V6> 100 & trueset2$V6<=1000,]
truest.del.M_no= trueset2[trueset2$V6> 1000 & trueset2$V6<=50000,]
truest.del.XM_no= trueset2[trueset2$V6> 50000 & trueset2$V6<=100000,]
truest.del.L_no= trueset2[trueset2$V6> 100000 ,]
trueste.XLS_no=trueset2[trueset2$V6<50 | trueset2$V6>1000000,]
truset_fino_50_no= trueset2[trueset2$V6> 50 & trueset2$V6<=50000,]


####

sink(output)



cat("=============================\n")

cat("Precision and recall  deletion\n")
cat("=============================\n")
cat("all\n")
PR_del(caller, trueset)
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
cat("=========until50=============\n")
PR_del(fino_50KB, truset_fino_50)


cat("=============================\n")
cat("Precision and recall  duplication\n")
cat("=============================\n")
cat("all\n")
PR_dup(caller, trueset)
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

cat("=========until50=============\n")

PR_dup(fino_50KB, truset_fino_50)


cat("==========End================\n")

cat("Precision and recall  deletion no segmental duplication \n")
cat("=============================\n")
cat("all\n")
PR_del(caller2, trueset2)
cat("=============================\n")
cat("SS\n")
PR_del(SS_no, truest.del.SS_no)
cat("=============================\n")
cat("S\n")
PR_del(S_no, truest.del.S_no)
cat("=============================\n")
cat("M\n")
PR_del(M_no, truest.del.M_no)
cat("=============================\n")
cat("XM\n")
PR_del(XM_no, truest.del.XM_no)
cat("=============================\n")
cat("L\n")
PR_del(L_no, truest.del.L_no)
cat("=========until50=============\n")
PR_del(fino_50KB_no, truset_fino_50_no)


cat("=============================\n")
cat("Precision and recall  duplication\n")
cat("=============================\n")
cat("all\n")
PR_dup(caller2, trueset2)
cat("=============================\n")
cat("SS\n")
PR_dup(SS_no, truest.del.SS_no)
cat("=============================\n")
cat("S\n")
PR_dup(S_no, truest.del.S_no)
cat("=============================\n")
cat("M\n")
PR_dup(M_no, truest.del.M_no)
cat("=============================\n")
cat("XM\n")

PR_dup(XM_no, truest.del.XM_no)
cat("=============================\n")
cat("L\n")
PR_dup(L_no, truest.del.L_no)

cat("=============================\n")

cat("=========until50=============\n")

PR_dup(fino_50KB_no, truset_fino_50_no)


cat("==========End================\n")





# Stop writing to the file
sink()
