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

#把样本按stage顺序排一下
ss<-list.files("D:\\表达谱配套stage",full.name=T)
ssd<-list.files("D:\\表达谱配套stage")
need_stage<-ss[match(aaa,gsub(".txt","",ssd))]

samp_stage<-lapply(1:length(res),function(x){
  sample_file<-needdd[x]
  sam_file<-readRDS(sample_file)
  stgs<-read.table(need_stage[x],sep="\n")
  a<-data.frame(id=colnames(sam_file),stage=stgs[,1])
  a[order(a$stage),]
})
saveRDS(samp_stage,"D:\\正在使用的样本和stage")


#在早-晚差异的细胞类型及其相关的lnc
samp_stage<-readRDS("D:\\正在使用的样本和stage")      #这个根据样本的stage排序了 不是最开始的

aa<-list.files("D:\\新互作边")
aaa<-sapply(aa,function(x){unlist(strsplit(x,"[.]"))[1]})
aaa<-aaa[-10]

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
lnc<-tapply(afs[,1],afs[,3],function(x){unique(x)})



#生成p值
#早晚中
xcell_in_EA_P<-lapply(1:length(xcells),function(each_cancer){
  if(each_cancer!=11){
    
    sample_and_stage<-samp_stage[[each_cancer]]
    cancer_immune<-xcells[[each_cancer]]
    order_sample<-cancer_immune[match(sample_and_stage[,1],cancer_immune$ID),]
    
    #在早-晚中的p值
    r1<-apply(order_sample[,-1],2,function(z){
      eora<-rep("advance_stage",dim(sample_and_stage)[1])
      eora[which(sample_and_stage[,2] %in% c("Stg I","Stg II"))]<-"early_stage"
      b<-wilcox.test(z~eora)$p.value
      return(b)
    })
    r1
  }
})
p1<-do.call(cbind,xcell_in_EA_P)
colnames(p1)<-gsub("TCGA-","",aaa[-11])
rownames(p1)<-gsub("_xCell","",rownames(p1))
p1<-(-log10(p1))
saveRDS(p1,"D:\\p_of_EA_xcell.rds")


#转移非转移中
clean_site<-readRDS("D:\\整理完的转移位点.rds")
has_site<-clean_site[c(1:3,5:8)]
dres <- xcells[match(names(has_site),aaa)]
bdp_with_met<-d[match(names(has_site),ddd)]
lnc_with_met<-lnc[match(names(has_site),names(lnc))]



#转移非转移单p值
met_with_imm<-lapply(1:length(dres),function(each_cancer){
  imm<-dres[[each_cancer]]
  met_site<-rep("non",dim(imm)[1])
  met_site[imm$ID %in% has_site[[each_cancer]][,1]]<-"met"
  p1<-apply(imm[,-1],2,function(x){wilcox.test(x~met_site)$p.value})
  p1
})
p2<-do.call(cbind,met_with_imm)
colnames(p2)<-gsub("TCGA-","",names(has_site))
rownames(p2)<-gsub("_xCell","",rownames(p2))
p2<-(-log10(p2))
saveRDS(p2,"D:\\p_of_met_xcell.rds")



tcga<-read.table("D:\\tcga.txt",sep="\t",header=T)
ensg<-tcga$gene_id

a<-readRDS("D:\\上传的代码\\基因和免疫细胞相关所有结果.rds")
b<-sapply(a$cancer,function(x){
  unlist(strsplit(x,"[.]"))
})
cc<-sapply(as.character(a$Var1),function(x){
  unlist(strsplit(x,";"))
})
a<-data.frame(a,t(b),t(cc),type=sub(".*\\.([^.]+)\\..*", "\\1", rownames(a)))
a$Var2<-as.character(a$Var2)
colnames(a)[5:8]<-c("can","stg","gene_name","gene_id")


lg<-tapply(1:dim(a)[1],a$can,function(x){
  cancer<-a[x,]
  s1<-tapply(1:dim(cancer)[1],cancer$Var2,function(y){
    cancer_cell<-cancer[y,]
    ucc<-unique(cancer_cell[,c("gene_name","type")])
    cd<-table(ucc$type)
    data.frame(type=names(cd),num=as.numeric(cd))
  })
  kk<-Reduce(function(m,n){
    merge(m,n,by="type",all=T)
  },s1[2:length(s1)],init=s1[[1]])
  colnames(kk)[-1]<-names(s1)
  kk
})

cell_type<-sort(unique(a$Var2))
cancer<-sort(unique(a$can))
bq<-lapply(lg,function(x){
  rownames(x)<-x[,1]
  if(dim(x)[2]!=38){
    mat<-matrix(0,nrow=2,ncol=38-dim(x)[2])
    colnames(mat)<-setdiff(cancer,colnames(x))
    x<-data.frame(x,mat,check.names=F)
  }
  x[,match(cell_type,colnames(x))]
  x[is.na(x)]<-0
  x[,-1]
})

g<-lapply(bq,function(x){
  x[which(rownames(x)=="gene"),]
})
gene1<-do.call(rbind.data.frame,g)
rownames(gene1)<-sub("TCGA-","",rownames(gene1))

l<-lapply(bq,function(x){
  x[which(rownames(x)=="lnc"),]
})
lnc1<-do.call(rbind.data.frame,l)
rownames(lnc1)<-sub("TCGA-","",rownames(lnc1))


#读入那几个数据并对齐
lnc1<-readRDS("D:\\lnc1.rds")
p1<-readRDS("D:\\p_of_EA_xcell.rds")
gene1<-readRDS("D:\\gene1.rds")
p2<-readRDS("D:\\p_of_met_xcell.rds")


n_p<-matrix(0,37,15)
n_p[,1:7]<-p2
new_order<-c(match(colnames(p2),rownames(lnc1)),setdiff(1:15,match(colnames(p2),rownames(lnc1))))
lnc1<-lnc1[new_order,]
lnc1[which(lnc1!=0)]<-log(lnc1[which(lnc1!=0)],5)
gene1[which(gene1!=0)]<-log(gene1[which(gene1!=0)],5)
gene1<-gene1[new_order,]


#读入那几个数据并对齐
lnc1<-readRDS("D:\\lnc1.rds")
p1<-readRDS("D:\\p_of_EA_xcell.rds")
gene1<-readRDS("D:\\gene1.rds")
p2<-readRDS("D:\\p_of_met_xcell.rds")



p11<-data.frame(p1,PAAD=0)
p11<-p11[,match(sub("TCGA-","",rownames(lnc1)),colnames(p11))]



library(ggplot2)
library(ggforce)

circle4Fun <- function(center = c(0,0), diameter = 1, npoints = 100, n_circle = 1){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  lab <- rep(paste(n_circle,1:4,sep="-"),each=npoints/4)
  circle_pos<-data.frame(x = c(xx,rep(center[1],4)), y = c(yy,rep(center[2],4)),n_circle,each_piece=c(rep(1:4,each=npoints/4),1:4),bq=c(lab,paste(n_circle,1:4,sep="-")))
  return(circle_pos)
}


#生成所有坐标
pos<-lapply(1:37,function(y){
  s1<-lapply(1:15,function(x){
    circle4Fun(c(1+2*(y-1),1+2*(x-1)),2,npoints = 100,n_circle = paste(y,x,sep="-"))
  })
  do.call(rbind.data.frame,s1)
})
all_pos<-do.call(rbind.data.frame,pos)


#把p1填上去
pp1<-apply(all_pos,1,function(x){
  n<-as.numeric(unlist(strsplit(x[3],"-")))
  if(x[4]==2){
    return(p11[n[1],n[2]]) 
  } 
  if(x[4]==1){
    return(n_p[n[1],n[2]]) 
  } 
  if(x[4]==3){
    return(lnc1[n[2],n[1]]) 
  } 
  if(x[4]==4){
    return(gene1[n[2],n[1]]) 
  }  
})
plot_data<-data.frame(all_pos,p1=pp1)

library("ggnewscale")
library("RColorBrewer")

q<-tapply(1:dim(plot_data)[1],plot_data$n_circle,function(x){
  step1<-plot_data[x,]
  data.frame(step1[which(step1$each_piece %in% c(1,2)),],step1[which(step1$each_piece %in% c(3,4)),])
})
lookok<-do.call(rbind.data.frame,q)

bkx<-sapply(1:37,function(x){1+2*(x-1)})
bky<-sapply(1:15,function(x){1+2*(x-1)})


ggplot(lookok) + 
  scale_x_continuous(breaks=bkx,labels=sub("_xCell","",colnames(lnc1)),position = "top")+
  scale_y_continuous(breaks=bky,labels=rownames(lnc1))+labs(x="",y="")+
  geom_shape(aes(x,y,fill = p1, group = bq), expand = unit(-0.1, 'mm'))+
  scale_fill_gradientn(colors = c("white",colorRampPalette(c("#fbead8","#df5d47","#d8261c"))(20)))+
  new_scale("fill") +
  geom_shape(aes(x.1,y.1,fill = p1.1, group = bq.1), expand = unit(-0.1, 'mm'))+
  theme_bw()+
  scale_fill_gradientn(colors = c("white",colorRampPalette(c("#e6e6fa","#b8b8c9","#7d7d88"))(20)))+
  theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=1),
        panel.grid = element_blank())+coord_fixed()


#Figure 4B
tcga<-read.table("D:\\tcga.txt",sep="\t",header=T)
ensg<-tcga$gene_id
#计算与所有细胞类型相关的lnc个数
xcell_with_lnc2<-lapply(1:length(xcells),function(each_cancer){
  
  cancer_immune<-as.data.frame(xcells[[each_cancer]])
  
  #取lnc表达谱
  bdp<-readRDS(needbdp[each_cancer])
  ne<-ensg[tcga[,2] %in% lnc[[each_cancer]]]
  lnc_bdp<-bdp[rownames(bdp) %in% ne,]
  duiq<-lnc_bdp[,match(cancer_immune[,1],colnames(lnc_bdp))]
  
  s2<-apply(cancer_immune[,-1],2,function(x){
    s1<-apply(duiq,1,function(y){
      cor_res<-cor.test(x,y,method="spearman")
      if(cor_res$p.value<0.05){return(cor_res$estimate)} else {return(0)}
    },simplify = F)
    do.call(rbind,s1)
  })
  rownames(s2)<-paste(tcga[match(rownames(lnc_bdp),ensg),2],rownames(lnc_bdp),sep=";")
  qh<-sapply(1:dim(s2)[1],function(h){
    d<-as.numeric(s2[h,])
    if(length(which(d==0))/dim(s2)[2]==1){return(0)} else{return(1)}
  })
  return(s2[which(qh==1),])                                      #返回具体的基因
  #apply(s2[which(qh==1),],2,function(zz){length(which(zz!=0))})   #统计个数
})
names(xcell_with_lnc2)<-aaa


#展示一些相关性很强的
a<-xcell_with_lnc2[["TCGA-PAAD"]]
h<-which(a>0.8,arr.ind = T)
cell<-colnames(a)[h[2]]
gene<-rownames(a)[h[1]]
cancer_immune<-as.data.frame(xcells$`TCGA-PAAD`)
lnc_bdp<-readRDS(needbdp[grep("paad",needbdp,ignore.case = T)])
duiq<-lnc_bdp[which(rownames(lnc_bdp)=="ENSG00000251442"),]
pd<-data.frame(lnc=as.numeric(duiq[1,]),cell=cancer_immune[match(colnames(lnc_bdp),cancer_immune$ID),match(cell,colnames(cancer_immune))])
ggplot(pd,aes(x=lnc,y=cell))+geom_point()+geom_smooth(method="lm")+theme_bw()


a<-xcell_with_lnc2[["TCGA-KIRP"]]
h<-which(a>0.7,arr.ind = T)
cell<-colnames(a)[h[2]]
gene<-rownames(a)[h[1]]
cancer_immune<-as.data.frame(xcells$`TCGA-KIRP`)
lnc_bdp<-readRDS(needbdp[grep("kirp",needbdp,ignore.case = T)])
duiq<-lnc_bdp[which(rownames(lnc_bdp)=="ENSG00000224137"),]
pd<-data.frame(lnc=as.numeric(duiq[1,]),cell=cancer_immune[match(colnames(lnc_bdp),cancer_immune$ID),match(cell,colnames(cancer_immune))])
ggplot(pd,aes(x=lnc,y=cell))+geom_point()+geom_smooth(method="lm")+theme_bw()



gene<-tapply(afs[,2],afs[,3],function(x){unique(x)})
xcell_with_gene2<-lapply(1:length(xcells),function(each_cancer){
  
  cancer_immune<-as.data.frame(xcells[[each_cancer]])
  
  #取gene表达谱
  bdp<-readRDS(needbdp[each_cancer])
  ne<-ensg[tcga[,2] %in% gene[[each_cancer]]]
  gene_bdp<-bdp[rownames(bdp) %in% ne,]
  duiq<-gene_bdp[,match(cancer_immune[,1],colnames(gene_bdp))]
  
  s2<-apply(cancer_immune[,-1],2,function(x){
    s1<-apply(duiq,1,function(y){
      cor_res<-cor.test(x,y,method="spearman")
      if(cor_res$p.value<0.05){return(cor_res$estimate)} else {return(0)}
    },simplify = F)
    do.call(rbind,s1)
  })
  rownames(s2)<-paste(tcga[match(rownames(gene_bdp),ensg),2],rownames(gene_bdp),sep=";")
  qh<-sapply(1:dim(s2)[1],function(h){
    d<-as.numeric(s2[h,])
    if(length(which(d==0))/dim(s2)[2]==1){return(0)} else{return(1)}
  })
  return(s2[which(qh==1),])
  #apply(s2[which(qh==1),],2,function(zz){length(which(zz!=0))})
})
names(xcell_with_gene2)<-aaa

#展示一些相关性很强的
a<-xcell_with_gene2[["TCGA-ACC"]]
h<-which(a>0.8,arr.ind = T)
cell<-colnames(a)[h[1,2]]
gene<-rownames(a)[h[1,1]]
cancer_immune<-as.data.frame(xcells$`TCGA-ACC`)
gene_bdp<-readRDS(needbdp[grep("acc",needbdp,ignore.case = T)])
duiq<-gene_bdp[which(rownames(gene_bdp)=="ENSG00000138180"),]
pd<-data.frame(gene=as.numeric(duiq[1,]),cell=cancer_immune[match(colnames(gene_bdp),cancer_immune$ID),match(cell,colnames(cancer_immune))])
ggplot(pd,aes(x=gene,y=cell))+geom_point()+geom_smooth(method="lm")+theme_bw()


a<-xcell_with_gene2[["TCGA-SKCM"]]
h<-which(a>0.8,arr.ind = T)
cell<-colnames(a)[h[4,2]]
gene<-rownames(a)[h[4,1]]
cancer_immune<-as.data.frame(xcells$`TCGA-SKCM`)
gene_bdp<-readRDS(needbdp[grep("skcm",needbdp,ignore.case = T)])
duiq<-gene_bdp[which(rownames(gene_bdp)=="ENSG00000066336"),]
pd<-data.frame(gene=as.numeric(duiq[1,]),cell=cancer_immune[match(colnames(gene_bdp),cancer_immune$ID),match(cell,colnames(cancer_immune))])
ggplot(pd,aes(x=gene,y=cell))+geom_point()+geom_smooth(method="lm")+theme_bw()


