#################################################
library(ggplot2)
library(plyr)
library(cowplot)
#################################################
###### estimation PVE and polygenic scores ######
### create a function to calculate model average ###
modelaverage<-function(dataset){
  datn<-fread(dataset[1],header=TRUE)[,c(1:4)]
  datn$para<-0
  datn$gamma<-0
  for(i in 1:10){
    dat<-fread(dataset[i],header=TRUE)
    dat$para<-dat$alpha+dat$beta*dat$gamma
    datn$para<-dat$para+datn$para
    datn$gamma<-dat$gamma+datn$gamma
  }
  datn$para<-datn$para/10
  datn$gamma<-datn$gamma/10
  datn
}


####### read data file #######
BSTh_dessication<-list.files(full.names=TRUE,pattern="BSTh_dessication_ch")
BSTh_growth<-list.files(full.names=TRUE,pattern="BSTh_growth_ch")
BSTh_heatresist<-list.files(full.names=TRUE,pattern="BSTh_heatresist_ch")
BSTh_weight7<-list.files(full.names=TRUE,pattern="BSTh_weight7_ch")
BSTh_weight14<-list.files(full.names=TRUE,pattern="BSTh_weight14_ch")

pos<-function(data,vect){
  data1<-data
  split_list1<-strsplit(vect, "_")
  # Extract integer and decimal parts
  split_list2<- sapply(split_list1, function(x)(x[2]))
  split_list3<-strsplit(split_list2, "\\.")
  data1$scaffold<-sapply(split_list3, function(x)(x[1]))
  data1$position<-sapply(split_list3, function(x)(x[2]))
  data1
}

BSTh_desc_avg<-modelaverage(BSTh_dessication)
BSTh_grow_avg<-modelaverage(BSTh_growth)
BSTh_heatres_avg<-modelaverage(BSTh_heatresist)
BSTh_wg7_avg<-modelaverage(BSTh_weight7)
BSTh_wg14_avg<-modelaverage(BSTh_weight14)


BSTh_desc_avg<-pos(BSTh_desc_avg,BSTh_desc_avg$rs)
BSTh_grow_avg<-pos(BSTh_grow_avg,BSTh_grow_avg$rs)
BSTh_heatres_avg<-pos(BSTh_heatres_avg,BSTh_heatres_avg$rs)
BSTh_wg7_avg<-pos(BSTh_wg7_avg,BSTh_wg7_avg$rs)
BSTh_wg14_avg<-pos(BSTh_wg14_avg,BSTh_wg14_avg$rs)

BSTl_dessication<-list.files(full.names=TRUE,pattern="BSTl_dessication_ch")
BSTl_growth<-list.files(full.names=TRUE,pattern="BSTl_growth_ch")
BSTl_heatresist<-list.files(full.names=TRUE,pattern="BSTl_heatresist_ch")
BSTl_weight7<-list.files(full.names=TRUE,pattern="BSTl_weight7_ch")
BSTl_weight14<-list.files(full.names=TRUE,pattern="BSTl_weight14_ch")


BSTl_desc_avg<-modelaverage(BSTl_dessication)
BSTl_grow_avg<-modelaverage(BSTl_growth)
BSTl_heatres_avg<-modelaverage(BSTl_heatresist)
BSTl_wg7_avg<-modelaverage(BSTl_weight7)
BSTl_wg14_avg<-modelaverage(BSTl_weight14)

BSTl_desc_avg<-pos(BSTl_desc_avg,BSTl_desc_avg$rs)
BSTl_grow_avg<-pos(BSTl_grow_avg,BSTl_grow_avg$rs)
BSTl_heatres_avg<-pos(BSTl_heatres_avg,BSTl_heatres_avg$rs)
BSTl_wg7_avg<-pos(BSTl_wg7_avg,BSTl_wg7_avg$rs)
BSTl_wg14_avg<-pos(BSTl_wg14_avg,BSTl_wg14_avg$rs)

head(BSTh_desc_avg)
BSTh_desc_avg<-data.frame(BSTh_desc_avg[,c(2,5:8)])
BSTh_grow_avg<-data.frame(BSTh_grow_avg[,c(2,5:8)])
BSTh_heatres_avg<-data.frame(BSTh_heatres_avg[,c(2,5:8)])
BSTh_wg7_avg<-data.frame(BSTh_wg7_avg[,c(2,5:8)])
BSTh_wg14_avg<-data.frame(BSTh_wg14_avg[,c(2,5:8)])

BSTl_desc_avg<-data.frame(BSTl_desc_avg[,c(2,5:8)])
BSTl_grow_avg<-data.frame(BSTl_grow_avg[,c(2,5:8)])
BSTl_heatres_avg<-data.frame(BSTl_heatres_avg[,c(2,5:8)])
BSTl_wg7_avg<-data.frame(BSTl_wg7_avg[,c(2,5:8)])
BSTl_wg14_avg<-data.frame(BSTl_wg14_avg[,c(2,5:8)])

### number of loci ####
head(BSTh_desc_avg)
length(BSTh_desc_avg$para)
head(BSTl_heatres_avg)
### sum of effect size ####
sum(BSTh_desc_avg$para)
sum(BSTh_grow_avg$para)
sum(BSTh_heatres_avg$para)
sum(BSTh_wg7_avg$para)
sum(BSTh_wg14_avg$para)

sum(BSTh_desc_avg$para)
sum(BSTl_grow_avg$para)
sum(BSTl_heatres_avg$para)
sum(BSTl_wg7_avg$para)
sum(BSTl_wg14_avg$para)

#### manhattan plot ###
head(BSTh_desc_avg)
lgpos<-read.table("lgs.txt",header=TRUE)
colnames(lgpos)<-c("scaffold","LG")
head(lgpos)
lgpos$LG<-factor(lgpos$LG,levels=c("1","2","3","4","5","6","7",
                                              "8","9","10","11","12","13",
                                              "14","15","16","17","18","19",
                                              "20","21","22","Z","NA"))
BSTh_desc_avgnew<-merge(BSTh_desc_avg,lgpos,by="scaffold",all.x=TRUE)
head(BSTh_desc_avgnew)
length(BSTh_desc_avgnew[,1])
BSTh_desc_avgnew<-BSTh_desc_avgnew[order(BSTh_desc_avgnew$LG,BSTh_desc_avgnew$scaffold,
                                  BSTh_desc_avgnew$position),]
BSTh_grow_avgnew<-merge(BSTh_grow_avg,lgpos,by="scaffold",all.x=TRUE)
BSTh_grow_avgnew<-BSTh_grow_avgnew[order(BSTh_grow_avgnew$LG,BSTh_grow_avgnew$scaffold,
                                         BSTh_grow_avgnew$position),]

BSTh_heatres_avgnew<-merge(BSTh_heatres_avg,lgpos,by="scaffold")
BSTh_heatres_avgnew<-BSTh_heatres_avgnew[order(BSTh_heatres_avgnew$LG,BSTh_heatres_avgnew$scaffold,
                                         BSTh_heatres_avgnew$position),]
BSTh_wg7_avgnew<-merge(BSTh_wg7_avg,lgpos,by="scaffold")
BSTh_wg7_avgnew<-BSTh_wg7_avgnew[order(BSTh_wg7_avgnew$LG,BSTh_wg7_avgnew$scaffold,
                                               BSTh_wg7_avgnew$position),]

BSTh_wg14_avgnew<-merge(BSTh_wg14_avg,lgpos,by="scaffold")
BSTh_wg14_avgnew<-BSTh_wg14_avgnew[order(BSTh_wg14_avgnew$LG,BSTh_wg14_avgnew$scaffold,
                                       BSTh_wg14_avgnew$position),]

BSTl_desc_avgnew<-merge(BSTl_desc_avg,lgpos,by="scaffold")
head(BSTh_desc_avgnew)
BSTl_desc_avgnew<-BSTl_desc_avgnew[order(BSTl_desc_avgnew$LG,BSTl_desc_avgnew$scaffold,
                                         BSTl_desc_avgnew$position),]
BSTl_grow_avgnew<-merge(BSTl_grow_avg,lgpos,by="scaffold")
BSTl_grow_avgnew<-BSTl_grow_avgnew[order(BSTl_grow_avgnew$LG,BSTl_grow_avgnew$scaffold,
                                         BSTl_grow_avgnew$position),]
BSTl_heatres_avgnew<-merge(BSTl_heatres_avg,lgpos,by="scaffold")
BSTl_heatres_avgnew<-BSTl_heatres_avgnew[order(BSTl_heatres_avgnew$LG,BSTl_heatres_avgnew$scaffold,
                                               BSTl_heatres_avgnew$position),]
BSTl_wg7_avgnew<-merge(BSTl_wg7_avg,lgpos,by="scaffold")
BSTl_wg7_avgnew<-BSTl_wg7_avgnew[order(BSTl_wg7_avgnew$LG,BSTl_wg7_avgnew$scaffold,
                                       BSTl_wg7_avgnew$position),]
BSTl_wg14_avgnew<-merge(BSTl_wg14_avg,lgpos,by="scaffold")
BSTl_wg14_avgnew<-BSTl_wg14_avgnew[order(BSTl_wg14_avgnew$LG,BSTl_wg14_avgnew$scaffold,
                                       BSTl_wg14_avgnew$position),]
head(BSTh_desc_avgnew)


# Number of levels in the factor
num_levels <- length(levels(BSTh_desc_avgnew$LG))

# Create a dataframe for the rectangles
rect_data <- data.frame(
  xmin = seq(0.5, by = 1, length.out = ceiling(num_levels)),
  xmax = seq(1.5, by = 1, length.out = ceiling(num_levels)),
  ymin = -Inf, # Adjust if you have a minimum y value
  ymax = Inf   # Adjust if you have a maximum y value
)

# Adjust if the number of levels is odd
if (num_levels %% 2 != 0) {
  rect_data <- rbind(rect_data, data.frame(xmin = num_levels, xmax = num_levels, ymin = -Inf, ymax = Inf))
}
# Create alternating colors
rect_data$fill <- ifelse(seq_along(rect_data$xmin) %% 2 == 0, "grey", "white")
#######################################################################
summary(BSTh_desc_avgnew$gamma)
BSTh_desc_avgnew$trait<-"desccication"
BSTh_grow_avgnew$trait<-"growth rate"
BSTh_heatres_avgnew$trait<-"heat resistance"
BSTh_wg7_avgnew$trait<-"7 day weight"
BSTh_wg14_avgnew$trait<-"14 day weight"

BSTh_avg<-rbind(BSTh_desc_avgnew,BSTh_grow_avgnew,BSTh_heatres_avgnew,
                BSTh_wg7_avgnew,BSTh_wg14_avgnew)
BSTh_avg<-BSTh_avg[order(BSTh_avg$LG,BSTh_avg$scaffold,BSTh_avg$position),]  
BSTh_avgnew<-BSTh_avg[BSTh_avg$gamma>0.01,]

BSTl_desc_avgnew$trait<-"desccication"
BSTl_grow_avgnew$trait<-"growth rate"
BSTl_heatres_avgnew$trait<-"heat resistance"
BSTl_wg7_avgnew$trait<-"7 day weight"
BSTl_wg14_avgnew$trait<-"14 day weight"

BSTl_avg<-rbind(BSTl_desc_avgnew,BSTl_grow_avgnew,BSTl_heatres_avgnew,
                BSTl_wg7_avgnew,BSTl_wg14_avgnew)

BSTl_avg<-BSTl_avg[order(BSTl_avg$LG,BSTl_avg$scaffold,BSTl_avg$position),]  
BSTl_avgnew<-BSTl_avg[BSTl_avg$gamma>0.01,]

BSTh_mat_plot<-ggplot()+
  geom_point(data=BSTh_avgnew,aes(x=LG,y=gamma))+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1))+
  geom_rect(data = rect_data, aes(xmin = as.numeric(as.character(xmin)), xmax = as.numeric(as.character(xmax)), ymin = ymin, ymax = ymax, fill = fill), inherit.aes = FALSE) +
  scale_fill_identity()+
  geom_point(data=BSTh_avgnew,aes(x=LG,y=gamma,col=trait,shape=trait),size=4)+
  geom_hline(yintercept =c(0.1,0.5),linetype="dotted")+
  scale_color_manual(values=c("#FFDB6D",   "#D16103",  "#52854C", "#4E84C4", "#293352"))+
  scale_shape_manual(values=15:20)+
  theme_bw()+labs(y="Inclusion probability",x="Chromosome",title=
                    "Caterpillars rearing in high temperature")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title=element_text(size=26),
        legend.text=element_text(size=24),
        plot.title =element_text(size = 30,hjust=0.5),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text = element_text(colour = 'black',size=20))

BSTl_mat_plot<-ggplot()+
  geom_point(data=BSTl_avgnew,aes(x=LG,y=gamma))+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1))+
  geom_rect(data = rect_data, aes(xmin = as.numeric(as.character(xmin)), xmax = as.numeric(as.character(xmax)), ymin = ymin, ymax = ymax, fill = fill), inherit.aes = FALSE) +
  scale_fill_identity()+
  geom_point(data=BSTl_avgnew,aes(x=LG,y=gamma,col=trait,shape=trait),size=4)+
  geom_hline(yintercept =c(0.1,0.5),linetype="dotted")+
  scale_color_manual(values=c("#FFDB6D",   "#D16103",  "#52854C", "#4E84C4", "#293352"))+
  scale_shape_manual(values=15:20)+
  theme_bw()+labs(y="Inclusion probability",x="Chromosome",title=
                    "Caterpillars rearing in benigh temperature")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title=element_text(size=26),
        legend.text=element_text(size=24),
        plot.title =element_text(size = 30,hjust=0.5),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28),
        axis.text = element_text(colour = 'black',size=20))

plot_grid(BSTh_mat_plot,BSTl_mat_plot,labels=c("A","B"),nrow=2,label_size=40)


#############################################################################
postH_dessi<-read.table("post_output_BSThigh_dessication.txt",header=FALSE)
length(postH_dessi[,1])
postL_dessi<-read.table("post_output_BSTlow_dessication.txt",header=FALSE)
length(postL_dessi[,1])

postH_heatr<-read.table("post_output_BSThigh_heatresist.txt",header=FALSE)
length(postH_heatr[,1])
postL_heatr<-read.table("post_output_BSTlow_heatresist.txt",header=FALSE)
length(postL_heatr[,1])

postH_growth<-read.table("post_output_BSThigh_growth.txt",header=FALSE)
length(postH_growth[,1])
postL_growth<-read.table("post_output_BSTlow_growth.txt",header=FALSE)
length(postL_growth[,1])

postH_weight7<-read.table("post_output_BSThigh_weight7.txt",header=FALSE)
length(postH_weight7[,1])
postL_weight7<-read.table("post_output_BSTlow_weight7.txt",header=FALSE)
length(postL_weight7[,1])

postH_weight14<-read.table("post_output_BSThigh_weight14.txt",header=FALSE)
length(postH_weight14[,1])
postL_weight14<-read.table("post_output_BSTlow_weight14.txt",header=FALSE)
length(postL_weight14[,1])

H_pve<-rbind(postH_dessi[1,1],postH_heatr[1,1],postH_growth[1,1],postH_weight7[1,1],postH_weight14[1,1])
H_lb<-rbind(postH_dessi[1,2],postH_heatr[1,2],postH_growth[1,2],postH_weight7[1,2],postH_weight14[1,2])
H_ub<-rbind(postH_dessi[1,3],postH_heatr[1,3],postH_growth[1,3],postH_weight7[1,3],postH_weight14[1,3])
H_pge<-rbind(postH_dessi[1,1]*postH_dessi[1,4], postH_heatr[1,1]*postH_heatr[1,4],postH_growth[1,1]*postH_growth[1,4],
             postH_weight7[1,1]*postH_weight7[1,4],postH_weight14[1,1]*postH_weight14[1,4])

L_pve<-rbind(postL_dessi[1,1],postL_heatr[1,1],postL_growth[1,1],postL_weight7[1,1],postL_weight14[1,1])
L_lb<-rbind(postL_dessi[1,2],postL_heatr[1,2],postL_growth[1,2],postL_weight7[1,2],postL_weight14[1,2])
L_ub<-rbind(postL_dessi[1,3],postL_heatr[1,3],postL_growth[1,3],postL_weight7[1,3],postL_weight14[1,3])
L_pge<-rbind(postL_dessi[1,1]*postL_dessi[1,4], postL_heatr[1,1]*postL_heatr[1,4],postL_growth[1,1]*postL_growth[1,4],
             postL_weight7[1,1]*postL_weight7[1,4],postL_weight14[1,1]*postL_weight14[1,4])


################ compute polygenic scores ##################
BSThigh_G<-fread("BSTh_G.txt",header=FALSE)
BSTlow_G<-fread("BSTl_G.txt",header=FALSE)
dim(BSThigh_G)
dim(BSTlow_G)
BST_G<-cbind(BSThigh_G,BSTlow_G[,-c(1:3)])
dim(BST_G)
BST_Gnew<-t(BST_G[,-c(1:3)])
hN<-dim(BST_Gnew)[1] ### number of individuals ###
hL<-dim(BST_Gnew)[2] ### number of loci ###
dim(BST_Gnew)
#### center G
for(i in 1:hL){
  BST_Gnew[,i]<-BST_Gnew[,i]-mean(BST_Gnew[,i])
}
BST_Gnew<-t(BST_Gnew)
BST_Gnew[1:6,1:6]
dim(BST_Gnew)

BSThigh_dessG<-BST_Gnew[BST_G$V1 %in% BSTh_desc_avg$rs,]
dim(BSThigh_dessG)
BSThigh_dessG<-as.matrix(t(BSThigh_dessG))

BSThigh_growG<-BST_Gnew[BST_G$V1 %in% BSTh_grow_avg$rs,]
dim(BSThigh_growG)
BSThigh_growG<-as.matrix(t(BSThigh_growG))

BSThigh_heatG<-BST_Gnew[BST_G$V1 %in% BSTh_heatres_avg$rs,]
dim(BSThigh_heatG)
BSThigh_heatG<-as.matrix(t(BSThigh_heatG))

BSThigh_wg7G<-BST_Gnew[BST_G$V1 %in% BSTh_wg7_avg$rs,]
dim(BSThigh_wg7G)
BSThigh_wg7G<-as.matrix(t(BSThigh_wg7G))

BSThigh_wg14G<-BST_Gnew[BST_G$V1 %in% BSTh_wg14_avg$rs,]
dim(BSThigh_wg14G)
BSThigh_wg14G<-as.matrix(t(BSThigh_wg14G))
BSThigh_wg14G[1:6,1:6]

BSTlow_dessG<-BST_Gnew[BST_G$V1 %in% BSTl_desc_avg$rs,]
dim(BSTlow_dessG)
BSTlow_dessG<-as.matrix(t(BSTlow_dessG))

BSTlow_growG<-BST_Gnew[BST_G$V1 %in% BSTl_grow_avg$rs,]
dim(BSTlow_growG)
BSTlow_growG<-as.matrix(t(BSTlow_growG))

BSTlow_heatG<-BST_Gnew[BST_G$V1 %in% BSTl_heatres_avg$rs,]
dim(BSTlow_heatG)
BSTlow_heatG<-as.matrix(t(BSTlow_heatG))

BSTlow_wg7G<-BST_Gnew[BST_G$V1 %in% BSTl_wg7_avg$rs,]
dim(BSTlow_wg7G)
BSTlow_wg7G<-as.matrix(t(BSTlow_wg7G))

BSTlow_wg14G<-BST_Gnew[BST_G$V1 %in% BSTl_wg14_avg$rs,]
dim(BSTlow_wg14G)
BSTlow_wg14G<-as.matrix(t(BSTlow_wg14G))

psH_decc<-BSThigh_dessG %*% BSTh_desc_avg$para
psH_grow<-BSThigh_growG %*% BSTh_grow_avg$para
psH_heat<-BSThigh_heatG %*% BSTh_heatres_avg$para
psH_wg7<-BSThigh_wg7G %*% BSTh_wg7_avg$para
psH_wg14<-BSThigh_wg14G %*% BSTh_wg14_avg$para

psL_decc<-BSTlow_dessG %*% BSTl_desc_avg$para
psL_grow<-BSTlow_growG %*% BSTl_grow_avg$para
psL_heat<-BSTlow_heatG %*% BSTl_heatres_avg$para
psL_wg7<-BSTlow_wg7G %*% BSTl_wg7_avg$para
psL_wg14<-BSTlow_wg14G %*% BSTl_wg14_avg$para

psH<-cbind(psH_decc,psH_grow,psH_heat,psH_wg7,psH_wg14)
psL<-cbind(psL_decc,psL_grow,psL_heat,psL_wg7,psL_wg14)
head(psH)
dim(psH)
head(psL)

cor.traits<-data.frame(matrix(nrow=5,ncol=2))
colnames(cor.traits)<-c("trait","cor")
cor.traits$trait<-traits
cor.traits[1,2]<-cor(psH[,1],psL[,1])
cor.traits[2,2]<-cor(psH[,2],psL[,2])
cor.traits[3,2]<-cor(psH[,3],psL[,3])
cor.traits[4,2]<-cor(psH[,4],psL[,4])
cor.traits[5,2]<-cor(psH[,5],psL[,5])

write.csv(cor.traits,"cor.traits")
par(mar=c(4.5,4.5,2.5,2.5),mfrow=c(2,1))

corH<-data.frame(matrix(nrow=5,ncol=5))
colnames(corH)<-traits
rownames(corH)<-traits
for(i in 1:5){
  for (j in 1:5){
    corH[i,j]<-cor(psH[1:195,i],psH[1:195,j])
  }
}

corL<-data.frame(matrix(nrow=5,ncol=5))
colnames(corL)<-traits
rownames(corL)<-traits
for(i in 1:5){
  for (j in 1:5){
    corL[i,j]<-cor(psL[196:474,i],psL[196:474,j])
  }
}

write.csv(corH,"cortraits_high.csv")
write.csv(corL,"cortraits_low.csv")
## panel A ##
cs<-brewer.pal(n=5,"Dark2")
rownames(H_pve)<-traits
H_pve<-data.frame(H_pve)
dotchart(H_pve$H_pve,labels=traits,pch=NA,cex.lab=cl,xlim=c(0,1),xlab="Variance explained",ylab="")

for(j in 1:5){
    segments(H_lb[j,1],j,H_ub[j,1],j,col=cs[j])
   points(H_pve[j,1],j,pch=19,col=cs[j])
  }

title(main="(a) Caterpillar performance heritability under high temperature",cex.main=cm)

H_h<-data.frame(H_pve,H_lb,H_ub)
## panel B ##

L_pve<-data.frame(L_pve)
dotchart(L_pve$L_pve,labels=traits,pch=NA,cex.lab=cl,xlim=c(0,1),xlab="Variance explained",ylab="")

for(j in 1:5){
  segments(L_lb[j,1],j,L_ub[j,1],j,col=cs[j])
  points(L_pve[j,1],j,pch=19,col=cs[j])
}

title(main="(b) Caterpillar performance heritability under benign temperature",cex.main=cm)
L_h<-data.frame(L_pve,L_lb,L_ub)
### save as 10*10 ###
write.csv(H_h,"Hightemp_heritability.csv")
write.csv(L_h,"Lowtemp_heritability.csv")
