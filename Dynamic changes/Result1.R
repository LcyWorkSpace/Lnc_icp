#Figure 1B 整体相关性情况
library(ggplot2)
a<-list.files("cancer_nets\\",full.names=T)
cn<-gsub(".*\\\\TCGA-(.*?)\\..*", "\\1", a)
i<-1
t1<-lapply(a,function(x){
  net1<-read.table(x,sep="\t",header=T)
  okcol<-gsub("[.]"," ",colnames(net1))
  s1<-apply(net1,1,function(h){ 
    return(h[c(1,2,match(h[grep("stage",okcol)],okcol))]) 
  })
  s1<-t(s1)
  da<-data.frame(cancer=cn[i],cor=as.numeric(s1[,3]))
  i<<-i+1
  da$group<-ifelse(da$cor>0,"+","-")
  return(da)
})
t2<-do.call(rbind,t1)
t2$cor<-abs(t2$cor)
t2$group<-factor(t2$group,levels=c("+","-"))

#所有相关性概览
ggplot(t2)+
  geom_boxplot(aes(x=cancer,y=cor,color=group),lwd=0.3)+
  scale_y_continuous(
    name = "Positive correlation",
    sec.axis = sec_axis(~.*(-1), name="Negative correlation")) +
  theme_light()+
  scale_color_manual(values = c("-" = "#87b390", "+" = "#e3604e"))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        axis.text = element_text(size=10),
        axis.title.x = element_blank())



#气泡饼图-lnc、gene比例
t1<-lapply(a,function(x){
  net1<-read.table(x,sep="\t",header=T)
  cc<-tapply(1:dim(net1)[1],net1$stage,function(y){
    lnc<-length(unique(net1[y,1]))
    gene<-length(unique(net1[y,2]))
    return(c(lnc,gene))
  })
  do.call(cbind.data.frame,cc)
})
library(plyr)
s1<-rbind.fill(lapply(t1,function(x){x[1,]}))
lnc<-s1[,c(4,1:3)]
lnc[is.na(lnc)]<-0
rownames(lnc)<-cn
s2<-rbind.fill(lapply(t1,function(x){x[2,]}))
gene<-s2[,c(4,1:3)]
gene[is.na(gene)]<-0
rownames(gene)<-cn

library(ggplot2)
library(scatterpie)
library(RColorBrewer)
zuob<-c()
for(a in 1:15){
  x<-1+(a-1)*6
  for(y in 1:4){
    zuob<<-rbind(zuob,c(x,1+(y-1)*5))
  }
}
lncc<-unlist(sapply(15:1,function(x){lnc[x,]}))
stgs<-gsub("Stg","stage",colnames(lnc))
genee<-unlist(sapply(15:1,function(x){gene[x,]}))
r<-log((lncc+genee+1),10)
data<-data.frame(x=zuob[,2],y=zuob[,1],lnc=lncc,gene=genee,r)
data<-data[data$r!=0,]

ggplot() + 
  geom_scatterpie(aes(x=x, y=y, r=r),data=data,cols=c("lnc","gene"), color=NA) + 
  coord_equal()+
  scale_x_continuous(breaks = c(1,6,11,16), expand = c(0,0),labels = stgs) +
  scale_y_continuous(expand = c(0, 0), breaks = unique(data$y),labels = gsub(" ","",rev(rownames(lnc))), sec.axis = dup_axis())+
  theme_bw()+
  theme(
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 0.6, vjust = 0.7,size=15),
    axis.text.y = element_text(size=13),
    legend.title = element_text(size=13),
    legend.text = element_text(size=13)
  )+
  scale_fill_manual(values=c("tomato","lightblue"))+
  geom_text(aes(zuob[,2],zuob[,1], label = lncc+genee))


#气泡饼图-正负相关对数
t1<-lapply(a,function(x){
  net1<-read.table(x,sep="\t",header=T)
  colnames(net1)<-sub("[.]"," ",colnames(net1))
  cc<-tapply(1:dim(net1)[1],net1$stage,function(y){y})
  i<<-0
  mm<-lapply(cc,function(h){
    i<<-i+1
    val<-net1[h,match(names(cc)[i],colnames(net1))]
    zc<-length(which(val>0))
    fc<-length(which(val<0))
    return(c(zc,fc))
  })
  do.call(cbind.data.frame,mm)
})

s1<-rbind.fill(lapply(t1,function(x){x[1,]}))
zc<-s1[,c(4,1:3)]
zc[is.na(zc)]<-0
rownames(zc)<-cn
s2<-rbind.fill(lapply(t1,function(x){x[2,]}))
fc<-s2[,c(4,1:3)]
fc[is.na(fc)]<-0
rownames(fc)<-cn
aa<-unlist(sapply(15:1,function(x){zc[x,]}))
bb<-unlist(sapply(15:1,function(x){fc[x,]}))
r<-log((aa+bb+1),10)
data<-data.frame(x=zuob[,2],y=zuob[,1],positive=aa,negative=bb,r)

ggplot() + 
  geom_scatterpie(aes(x=x, y=y, r=r),data=data,cols=c("positive","negative"), color=NA) + 
  coord_equal()+
  scale_x_continuous(breaks = c(1,6,11,16), expand = c(0,0),labels = stgs) +
  scale_y_continuous(expand = c(0, 0), breaks = unique(data$y),labels = gsub(" ","",rev(rownames(zc))), sec.axis = dup_axis())+
  theme_bw()+
  theme(
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 0.6, vjust = 0.7,size=15),
    axis.text.y = element_text(size=13),
    legend.title = element_text(size=13),
    legend.text = element_text(size=13)
  )+
  scale_fill_manual(values=c("tomato2","darkseagreen"))+
  geom_text(aes(zuob[,2],zuob[,1], label = aa+bb)) 



#Figure 1C
i<<-0
t1<-lapply(a,function(x){
  i<<-i+1
  net1<-read.table(x,sep="\t",header=T)
  data.frame(net1[,c("lnc_name","gene_name","stage")],cancer=cn[i])
})
all_net<-do.call(rbind.data.frame,t1)
early<-all_net[all_net$stage %in% c("Stg I","Stg II"),]
all_early_lnc<-unique(early$lnc_name)

adv<-all_net[all_net$stage %in% c("Stg III","Stg IV"),]
all_advance_lnc<-unique(adv$lnc_name)


library(RISmed)
options(timeout=60)
all_lnc<-sample1[which(sample1[,3]=="lncRNA"),2]         #sample1是tcga的随意一个样本

num2<-c()
i<-0
bbb<-sapply(all_lnc,function(x){
  i<<-i+1
  search_topic <- paste("(",x,")"," AND ","(early stage)",sep="")
  eee<-EUtilsSummary(search_topic,db="pubmed", retmax=10000,mindate=2000, maxdate=2023)@count
  num2[i]<<-eee
  return(NULL)
})
saveRDS(num2,"num2.rds")


num3<-c()
i<-1
ccc<-sapply(all_lnc,function(x){
  i<<-i+1  
  search_topic <- paste("(",x,")"," AND ","(metastasis)",sep="")
  eee<-EUtilsSummary(search_topic,db="pubmed", retmax=10000,mindate=2000, maxdate=2023)@count
  num3[i]<<-eee
  return(NULL)
})
saveRDS(num3,"num3.rds")



early_stage<-readRDS("num2.rds")
names(early_stage)<-all_lnc
length(which(early_stage!=0))
metastasis<-readRDS("num3.rds")
names(metastasis)<-all_lnc
length(which(metastasis!=0))

length(intersect(all_early_lnc,names(early_stage)[which(early_stage!=0)]))
phyper(97-1, 361, 16901-361, 780,lower.tail = F)           #p=1.153515e-47

length(intersect(all_advance_lnc,names(metastasis)[which(metastasis!=0)]))
phyper(241-1, 1114, 16901-1114, 944,lower.tail = F)       #p=3.842642e-82

#与lnc2cancer交
stg4<-sample1[match(unique(as.character(unlist(s44))),gg),2]
length(intersect(stg4,yzg))
phyper(122-1, 1051, 16901-1051, 764,lower.tail = F)

hh<-unique(c(stg3,stg4))
length(intersect(hh,yzg))
phyper(137-1, 1051, 16901-1051, 944,lower.tail = F)


#Figure 1D
cancer_name<-str_extract(a,"TCGA-(\\w+).txt",group=1)
i<-0
incancer<-lapply(a,function(x){
  i<<-i+1
  b<-read.table(x,sep="\t",header=T)
  d<-data.frame(gene=c(b[,"lnc_id"],b[,"gene_id"]),cancer=cancer_name[i])
  unique(d)
})
all_gene_cancer<-do.call(rbind.data.frame,incancer)
mm<-tapply(all_gene_cancer$cancer,all_gene_cancer$gene,function(x){length(x)})

gs<-data.frame(type=c("1","2-5","6-10",">10"),
               num=c(2750,1105,72,17))

ggplot(gs, aes(x = '', y = num, fill = type)) +  
  geom_bar(stat = "identity", color=NA, width = 1) +  
  scale_fill_manual(values = c("#e86032","#71caf2","#83a0bf","#f6c949")) +  
  coord_polar("y", start = 0) +  
  theme_void() +  
  theme(legend.position = "right") +   
  geom_text(aes(x=-0.7,label='')) 


#Figure 1E
instage<-lapply(a,function(x){
  b<-read.table(x,sep="\t",header=T)
  d<-rbind(as.matrix(b[,c("lnc_id","stage")]),as.matrix(b[,c("gene_id","stage")]))
  d<-unique(d)
  m<-tapply(d[,2],d[,1],function(k){length(k)})
  as.data.frame(table(m)/sum(table(m)))
})
all_cancer_instage<-Reduce(function(x,y){
  merge.data.frame(x,y,by="m",all=T)
},instage[2:15],init=instage[[1]])
library(stringr)
colnames(all_cancer_instage)[-1]<-str_extract(a,"TCGA-(\\w+).txt",group=1)
library(reshape2)
data<-melt(all_cancer_instage)
data<-na.omit(data)
library(ggplot2)
library(ggbreak) 

ggplot(data, aes(x = variable, y = value, fill = as.factor(m))) +
  geom_bar(stat = "identity", position = "stack",width=0.6) +
  theme_void() +
  scale_y_break(c(0.2,0.9))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())



tcga<-read.table("tcga.txt",sep="\t",header=T)            #tcga是tcga中的任意一个样本，主要取ensg——id与基因名和基因类型的对应关系
g10<-tcga[match(names(which(mm>10)),tcga$gene_id),2]
paste0(g10,collapse=",")

dd<-lapply(a,function(x){
  b<-read.table(x,sep="\t",header=T)
  d<-rbind(as.matrix(b[,c("lnc_id","stage")]),as.matrix(b[,c("gene_id","stage")]))
  d<-unique(d)
  m<-tapply(d[,2],d[,1],function(k){length(k)})
  tcga[match(names(m)[which(m==4)],tcga$gene_id),2]
})
ff<-unique(unlist(dd))

paste0(intersect(g10,ff),collapse=",")
paste0(setdiff(g10,ff),collapse=",")
paste0(setdiff(ff,g10),collapse=",")

tcga[match(intersect(g10,ff),tcga$gene_name),2:3]

