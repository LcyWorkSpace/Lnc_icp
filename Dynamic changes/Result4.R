library(IOBR)
tcga<-read.table("D:\\tcga.txt",sep="\t",header=T)

#拿表达谱
a<-list.files("D:\\整合完的表达谱",full.name=T)
b<-list.files("D:\\整合完的表达谱")
bb<-sapply(b,function(x){unlist(strsplit(x,"[.]"))[1]})

#看需要哪些癌
aa<-list.files("D:\\新互作边")
aaa<-sapply(aa,function(x){unlist(strsplit(x,"[.]"))[1]})
aaa<-aaa[-10]
needbdp<-a[match(aaa,bb)]
timer_names<-tolower(gsub("TCGA-","",aaa))

#从初始表达谱里提取用到的样本
d<-list.files("D:\\处理过的表达谱",full.name=T)
dd<-list.files("D:\\处理过的表达谱")
ddd<-sapply(dd,function(x){unlist(strsplit(x,"[.]"))[1]})
needdd<-d[match(aaa,ddd)]



test<-lapply(1:length(needbdp),function(x){
  
  #读表达谱+去重
  bdp_file<-needbdp[x]
  bdp<-read.table(bdp_file,sep="\t",header=T,check=F)
  kk<-tcga[match(rownames(bdp),tcga[,1]),2]
  rms<-kk[which(!duplicated(kk))]
  nbdp<-bdp[which(!duplicated(kk)),]
  rownames(nbdp)<-rms
  
  #取需要的样本
  sample_file<-needdd[x]
  sam_file<-readRDS(sample_file)
  nbdp<-nbdp[,colnames(nbdp) %in% colnames(sam_file)]
  
  xcell<-deconvo_tme(eset = nbdp, method = "xcell",arrays=F)
  return(xcell)
})
saveRDS(test,"D:\\免疫浸润结果.rds")


#每个stage里做相关
tcga<-read.table("D:\\tcga.txt",sep="\t",header=T)
ensg<-tcga$gene_id

res<-readRDS("D:\\免疫浸润结果.rds")
names(res)<-aaa

#把xcell的细胞类型筛选一下
del_names<-c("aDC_xCell","Adipocytes_xCell","Astrocytes_xCell","Chondrocytes_xCell","Erythrocytes_xCell","HSC_xCell","Mesangial_cells_xCell","mv_Endothelial_cells_xCell",
             "Preadipocytes_xCell","Smooth_muscle_xCell","Endothelial_cells_xCell","Melanocytes_xCell","Myocytes_xCell","Pericytes_xCell","ImmuneScore_xCell",
             "Keratinocytes_xCell","MPP_xCell","Plasma_cells_xCell","Sebocytes_xCell","StromaScore_xCell","CMP_xCell","Hepatocytes_xCell","ly_Endothelial_cells_xCell",
             "MEP_xCell","MSC_xCell","Neurons_xCell","Osteoblast_xCell","Platelets_xCell","Skeletal_muscle_xCell","MicroenvironmentScore_xCell")

#取出xcell中剩下的细胞类型
xcells<-lapply(res,function(x){
  xcell<-x[[3]]
  xcell[,!colnames(xcell) %in% del_names]
})
names(xcells)<-aaa


b<-list.files("D:\\表达谱配套stage\\",full.names=T)
stgs<-lapply(b,function(x){read.table(file=x,sep="\n")})
names(stgs)<-b


d<-list.files("D:\\处理过的表达谱",full.name=T)
dd<-list.files("D:\\处理过的表达谱")
ddd<-sapply(dd,function(x){unlist(strsplit(x,"[.]"))[1]})
needbdp<-d[match(aaa,ddd)]

s<-list.files("D:\\所有癌症四个阶段网络\\",full.name=T)
af<-lapply(c(6,8,10,12),function(x){
  a<-read.table(s[x],sep="\t")
  a
})
afs<-do.call(rbind,af)

a<-split(afs,~V3+V4,drop=T)
i<<-1
test<-lapply(a,function(x){
  bdp<-readRDS(d[match(x[1,3],ddd)])
  ss<-stgs[[grep(x[1,3],names(stgs))]]
  sams<-bdp[,which(ss[,1]==x[1,4])]
  
  #在每个stage样本中计算gene,lnc,细胞的相关性
  imm<-xcells[[match(x[1,3],names(xcells))]]
  cells<-imm[match(colnames(sams),imm$ID),-1]
  
  #取lnc表达谱
  lnc<-ensg[tcga[,2] %in% x[,1]]
  lnc_bdp<-sams[rownames(sams) %in% lnc,]
  
  #取gene表达谱
  gene<-ensg[tcga[,2] %in% x[,2]]
  gene_bdp<-sams[rownames(sams) %in% gene,]
  
  res_lnc<-apply(cells,2,function(m){
    s1<-apply(lnc_bdp,1,function(n){
      cor_res<-cor.test(m,n,method="spearman")
      if(cor_res$p.value<0.05){return(cor_res$estimate)} else {return(0)}
    },simplify = F)
    do.call(rbind,s1)
  })
  
  
  res_gene<-apply(cells,2,function(m){
    s1<-apply(gene_bdp,1,function(n){
      cor_res<-cor.test(m,n,method="spearman")
      if(cor_res$p.value<0.05){return(cor_res$estimate)} else {return(0)}
    },simplify = F)
    do.call(rbind,s1)
  })
  
  if(dim(x)[1]==1){return(list(lnc=res_lnc,gene=res_gene))}
  
  
  rownames(res_lnc)<-paste(tcga[match(rownames(lnc_bdp),ensg),2],rownames(lnc_bdp),sep=";")
  qh<-sapply(1:dim(res_lnc)[1],function(h){
    d<-as.numeric(res_lnc[h,])
    if(length(which(d==0))/dim(res_lnc)[2]==1){return(0)} else{return(1)}
  })
  
  
  rownames(res_gene)<-paste(tcga[match(rownames(gene_bdp),ensg),2],rownames(gene_bdp),sep=";")
  qh2<-sapply(1:dim(res_gene)[1],function(h){
    d<-as.numeric(res_gene[h,])
    if(length(which(d==0))/dim(res_gene)[2]==1){return(0)} else{return(1)}
  })
  print(names(a)[i])
  i<<-i+1
  return(list(lnc=res_lnc[which(qh==1),],gene=res_gene[which(qh2==1),]))
})


#因为BRCA在stg3只有一个对，所以和其他的格式不一样，单独调整一下它
tz<-matrix(test[["TCGA-BRCA.Stg III"]][["lnc"]],nrow=1)
colnames(tz)<-names(test[["TCGA-BRCA.Stg III"]][["lnc"]])
rownames(tz)<-"EGOT;ENSG00000235947"
tz2<-matrix(test[["TCGA-BRCA.Stg III"]][["gene"]],nrow=1)
colnames(tz2)<-names(test[["TCGA-BRCA.Stg III"]][["gene"]])
rownames(tz2)<-"PBX3;ENSG00000167081"
test[["TCGA-BRCA.Stg III"]][["lnc"]]<-tz
test[["TCGA-BRCA.Stg III"]][["gene"]]<-tz2

library(reshape2)

i<<-0
all_res<-lapply(test,function(x){
  i<<-i+1
  s1<-lapply(x,function(y){
    melt(y)
  })
  data.frame(do.call(rbind,s1),cancer=names(test)[i])
})


t1<-do.call(rbind,all_res)
t1<-t1[t1$value!=0,]
saveRDS(t1,"D:\\基因和免疫细胞相关所有结果.rds")


t2<-split(t1,~Var1+Var2,drop=T)
t6<-sapply(t2,function(x){dim(x)[1]})



#Figure 4C####
t8<-tapply(t1$value,t1$Var1,function(x){length(which(x!=0))})
ggplot()+geom_density(aes(t8),fill="grey")+labs(x="Number of cancer progression",y="Density of gene/lncRNA count")+theme_bw()
da<-sort(t8,decreasing=T)[1:20]
data<-data.frame(gene=factor(sub(";.*","",names(da)),levels = sub(";.*","",names(da))),da)
ggplot(data)+geom_bar(aes(x=gene,y=da),stat="identity",fill="grey")+labs(x="",y="Number of cancer progression")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))


#Figure 4D####
#展示MALAT1,FENDRR,FUS,CRNDE在不同阶段，不同细胞类型，在不同癌症中相关性的正负
tq<-t1[t1[,1] %in% c("MALAT1;ENSG00000251562","FENDRR;ENSG00000268388","FUS;ENSG00000089280","CRNDE;ENSG00000245694"),]
tq[,1]<-sapply(as.character(tq[,1]),function(x){unlist(strsplit(x,";"))[1]})
cf<-sapply(as.character(tq[,4]),function(x){unlist(strsplit(x,"[.]"))})
cc<-t(cf)
newdf<-data.frame(tq[,1:3],cc)


pp<-data.frame(x=1:dim(newdf)[1],y=newdf$value)
kk<-ggplot(pp, aes(x = x,y=y)) +
  geom_point(aes(color=y)) +  # 绘制蓝色的点
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
print(kk)
ggdata<-ggplot_build(kk)$data[[1]]
ggcol<-ggdata[,c(3,1)]
ggcol<-unique(ggcol)

#分别生成
malat1<-newdf[newdf[,1]=="MALAT1",]
malat1_plot<-split(malat1,~Var2+X2,drop=F)
fendrr<-newdf[newdf[,1]=="FENDRR",]
fendrr_plot<-split(fendrr,~Var2+X2,drop=F)
fus<-newdf[newdf[,1]=="FUS",]
fus_plot<-split(fus,~Var2+X2,drop=F)
crnde<-newdf[newdf[,1]=="CRNDE",]
crnde_plot<-split(crnde,~Var2+X2,drop=F)
all_plot<-c(malat1_plot,fendrr_plot,fus_plot,crnde_plot)

void_plot<-function(){
  ggplot() +
    theme_void()
}

non_void_plot<-function(x){
  sx<-sort(unique(x$X1))
  ymin<-0
  ymax<-1
  testp<-tapply(1:dim(x)[1],x$X1,function(cancer){
    z1<-sprintf('geom_rect(aes(xmin=0,xmax=3,ymin=%s,ymax=%s),fill="%s")',ymin,ymax,coll[match(sx[ymax],coll[,1]),2])
    z2<-sprintf('geom_rect(aes(xmin=4,xmax=5,ymin=%s,ymax=%s),fill="%s")',ymin,ymax,ggcol[match(x[cancer,3],ggcol[,1]),2])
    ymin<<-ymin+1
    ymax<<-ymax+1
    paste0(z1,"+",z2)
  })
  paste0('ggplot()+',paste0(testp,collapse = "+"),'+coord_polar("y", start = 0)+theme_void()')
}


coll<-rbind(
  c("TCGA-BRCA","#725783"),
  c("TCGA-COAD","#2cb5c0"),
  c("TCGA-HNSC","#33a65c"),
  c("TCGA-KIRC","#c9d6e3"),
  c("TCGA-KIRP","#f7a694"),
  c("TCGA-LIHC","#b9d7bd"),
  c("TCGA-LUAD","#e5eaaa"),
  c("TCGA-PAAD","#ff7171"),
  c("TCGA-SKCM","#ef5845"),
  c("TCGA-STAD","#7873c0"),
  c("TCGA-THCA","#fff2d3"),
  c("TCGA-UCEC","#bab0ac"),
  c("TCGA-ACC","#fabfd2"),
  c("TCGA-BLCA","#d383ab"),
  c("TCGA-OV","#c37171"))

tttt<-lapply(all_plot,function(da){
  if(dim(da)[1]!=0){
    eval(parse(text=non_void_plot(da)))
  }else{
    void_plot()
  }
})

wrap_plots(tttt, nrow = 37,byrow=F)


#Figure 4E####
t2[which(t6==24)]
ggplot(t2[[which(t6==24)]],aes(x=cancer,y=value))+
  geom_bar(stat="identity")+
  theme_bw()+  
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))+
  labs(x="Cancer progress",y="cor")


col<-rbind(
  c("BRCA","#725783"),
  c("COAD","#2cb5c0"),
  c("HNSC","#33a65c"),
  c("KIRC","#c9d6e3"),
  c("KIRP","#f7a694"),
  c("LIHC","#b9d7bd"),
  c("LUAD","#e5eaaa"),
  c("PAAD","#ff7171"),
  c("SKCM","#ef5845"),
  c("STAD","#7873c0"),
  c("THCA","#fff2d3"),
  c("UCEC","#bab0ac"),
  c("ACC","#fabfd2"),
  c("BLCA","#d383ab"),
  c("OV","#c37171"))


col1<-sapply(1:dim(col)[1],function(y){
  x<-col[y,1]
  if(length(grep(x,t2[which(t6==22)][[1]][,4]))!=0){
    return(data.frame(can_pos=grep(x,t2[which(t6==22)][[1]][,4]),col=col[y,2]))  
  }
})
col1<-do.call(rbind,col1)[order(do.call(rbind,col1)[,1]),2]

ggplot()+
  geom_bar(data=t2[which(t6==22)][[1]],aes(x=cancer,y=value),stat="identity",fill=col1)+
  theme_bw()+  
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))+
  labs(x="Cancer progress",y="cor")






#Figure 4F####
t9<-t2[which(t6>=10)]
t10<-lapply(t9,function(x){
  return(c(as.character(x[1,1]),as.character(x[1,2]),length(which(x$value>0)),length(which(x$value<0))))
})
t11<-as.data.frame(do.call(rbind,t10))
rownames(t11)<-NULL
t11$V3<-as.numeric(t11$V3)
t11$V4<-as.numeric(t11$V4)
rr<-log((t11$V3+t11$V4),10)*0.4
t11<-data.frame(t11,x=match(t11$V1,unique(t11$V1)),y=match(t11$V2,unique(t11$V2)),rr)
colnames(t11)[3:4]<-c("Positive","Negetive")

library(ggplot2)
library(scatterpie)
ggplot() + geom_scatterpie(data=t11, aes(x=x, y=y,r=rr), cols=c("Positive","Negetive"),color=NA) + coord_fixed()+
  geom_text(data=t11,aes(x=x,y=y, label = Positive+Negetive))+
  scale_fill_manual(values=c("tomato2","lightblue"))+
  scale_x_continuous(breaks = 1:length(unique(t11$V1)), expand = c(0,0),labels = gsub(";.*","",unique(t11$V1))) +
  scale_y_continuous(breaks = 1:length(unique(t11$V2)), expand = c(0,0),labels = gsub("_xCell","",unique(t11$V2)))+
  theme_bw()+
  theme(
    axis.text.x = element_text(size=13,angle=90,vjust = 0.5,hjust=1),
    axis.text.y = element_text(size=13),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank()
  )


