#########################################################################
########## packages to load #######
library(data.table)
library(stringr)
library(plyr)
########## prepare genotype and phenotype files for GWAS ###############
G<-as.matrix(fread("g_est.txt",header=FALSE,sep=","))
LYCID<-read.table("indiv_ids.txt",header=TRUE,sep=",")
head(LYCID)
ids<-read.table("ids.txt",header=FALSE)
P<-read.csv("phenotypes.csv")
head(P)
G[1:5,1:5]
dim(G)

########################################################
Gdata<-data.frame(ind=LYCID$ind,pop=ids,G)
Gdata[1:5,1:5]
Gdata$ind[Gdata$ind=="VIC-20-09-10b"]<-"VIC-20-09-01"
Gdata$ind[Gdata$ind=="VIC-20-09-10a"]<-"VIC-20-09-10"
Gdata$ind[Gdata$ind=="BST-21-07-08"]<-"BST-21-07-09"

P$indID<-as.numeric(P$individual.ID)
P$iFID<-as.numeric(P$Female.ID)

P$individual.IDnew[!is.na(P$indID)]<-sprintf("%02d",as.numeric(P$individual.ID[!is.na(P$indID)]))
P$individual.IDnew[is.na(P$indID)]<-paste(sprintf("%02d",as.numeric(str_extract(P$individual.ID[is.na(P$indID)],"[1-9]+"))),
      str_extract(P$individual.ID[is.na(P$indID)],"[aA-zZ]+"),sep="")

Gdata$ind[order(Gdata$ind)]


P$Female.IDnew<-sprintf("%02d",P$Female.ID)

P$ind<-paste(P$Butterfly.population,P$Year,P$Female.IDnew,P$individual.IDnew,sep="-")
head(P$ind)
P$ind[order(P$ind)]

### IDoverlap<-Gdata[Gdata$ind %in% P$ind,c("ind","V1")]
### BSTID<-IDoverlap[IDoverlap$V1=="BST",]
### length(BSTID[,1])

BSTIDG<-Gdata[Gdata$V1=="BST",c("ind","V1")]
head(BSTIDG)
length(BSTIDG[,1])

BSTIDP<-P$ind[P$Butterfly.population=="BST"]
head(P$ind)
length(P$ind[P$Butterfly.population=="BST"])

### BSTcomp<-data.frame(PID=BSTIDP[order(BSTIDP)],GID=c(BSTIDG$ind[order(BSTIDG$ind)],rep(NA,39)))
#### write.csv(BSTcomp,"BSTcomp.csv")
Gdatanew<-Gdata[Gdata$ind %in% P$ind,]
length(Gdatanew[,1])
Gdatanew[1:5,1:5]

Pnew<-P[P$ind %in% Gdatanew$ind,]
head(Pnew)
length(Pnew[,1])
Pnew<-Pnew[order(Pnew$ind),]
head(Pnew)
head(Gdatanew[,1:3])
Pnew$dessication<-(Pnew$weight.before-Pnew$weight.after)/Pnew$weight.before
Pnew$growth<-(Pnew$weight.14-Pnew$weight.7)/Pnew$weight.7
min<-c()
sec<-c()
for(i in 1:length(Pnew[,1])){
  min[i]<-as.numeric(strsplit(Pnew$heat.time, ":")[[i]][1])
  sec[i]<-as.numeric(strsplit(Pnew$heat.time, ":")[[i]][2])
}

Pnew$heat.resist<-min*60+sec

BST_P<-Pnew[Pnew$Butterfly.population=="BST",]
length(BST_P[,1])
unique(BST_P$Rearing.temperature)

BST_P<-BST_P[,c("ind","Rearing.temperature","weight.7","weight.14","growth","dessication","heat.resist")]
head(BST_P)
############################################################################################

BSTh_P<-BST_P[BST_P$Rearing.temperature=="36/30",]
length(BSTh_P[,1])
BSTl_P<-BST_P[BST_P$Rearing.temperature=="26/20",]
length(BSTl_P[,1])
head(BSTh_P)

write.table(BSTh_P$weight.7,"BSTh_weight7.txt",col.names=F,row.names = F)
write.table(BSTh_P$weight.14,"BSTh_weight14.txt",col.names=F,row.names = F)
write.table(BSTh_P$growth,"BSTh_growth.txt",col.names=F,row.names = F)
write.table(BSTh_P$dessication,"BSTh_dessication.txt",col.names=F,row.names = F)
write.table(BSTh_P$heat.resist,"BSTh_heatresist.txt",col.names=F,row.names = F)

write.table(BSTl_P$weight.7,"BSTl_weight7.txt",col.names=F,row.names = F)
write.table(BSTl_P$weight.14,"BSTl_weight14.txt",col.names=F,row.names = F)
write.table(BSTl_P$growth,"BSTl_growth.txt",col.names=F,row.names = F)
write.table(BSTl_P$dessication,"BSTl_dessication.txt",col.names=F,row.names = F)
write.table(BSTl_P$heat.resist,"BSTl_heatresist.txt",col.names=F,row.names = F)

BSTh_G<-Gdatanew[Gdatanew$ind %in% BSTh_P$ind,]
length(BSTh_G[,1])
BSTl_G<-Gdatanew[Gdatanew$ind %in% BSTl_P$ind,]
length(BSTl_G[,1])

snps<-read.table("snps.txt",header=FALSE,sep=" ")
head(snps)
snps$loci<-paste("Scaffold",paste(snps$V1,snps$V2,sep="."),sep="_")
head(snps)

BSTh_Gnew<-BSTh_G[,-c(1,2)]
BSTh_Gnew[1:5,1:5]
BSTh_Gnew<-t(BSTh_Gnew)
BSTh_Gnew[1:5,1:5]
BSTh_Gnew2<-data.frame(V1=snps$loci,V2=rep("T",length(snps$loci)),V3=rep("A",length(snps$loci)),BSTh_Gnew)
BSTh_Gnew2[1:5,1:5]
write.table(BSTh_Gnew2,"BSTh_G.txt",col.names=F,row.names=F)

BSTl_Gnew<-BSTl_G[,-c(1,2)]
BSTl_Gnew[1:5,1:5]
BSTl_Gnew<-t(BSTl_Gnew)
BSTl_Gnew[1:5,1:5]
BSTl_Gnew2<-data.frame(V1=snps$loci,V2=rep("T",length(snps$loci)),V3=rep("A",length(snps$loci)),BSTl_Gnew)
BSTl_Gnew2[1:5,1:5]
write.table(BSTl_Gnew2,"BSTl_G.txt",col.names=F,row.names=F)

