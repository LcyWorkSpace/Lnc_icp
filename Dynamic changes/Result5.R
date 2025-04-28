#处理单细胞样本
library(Matrix)
Biermann <- readMM("Exp_data_UMIcounts.mtx")
genes <- read.table("Genes.txt",sep="\t")
cellsid <- read.csv("Cells.csv")
Biermann<-Biermann[!duplicated(genes[,1]), ]

Biermann <- CreateSeuratObject(Biermann)
colnames(Biermann) <- cellsid$cell_name
rownames(Biermann) <- genes[!duplicated(genes[,1]),1]
Biermann$cell_type <- cellsid$cell_type
Biermann$cell_subtype <- cellsid$cell_subtype
Biermann$orig.ident <- cellsid$sample
Biermann$source <- cellsid$source

# 取出脑转移的样本
needdata<-subset(Biermann,source=="Brain_Metastasis")

library(dplyr)
allen_reference <- SCTransform(needdata, ncells = 15000, verbose = FALSE)
aa<-allen_reference %>% RunPCA(verbose = FALSE) 
ElbowPlot(aa, ndims = 50)
aa %>% RunUMAP(dims = 1:30) -> allen_reference
DimPlot(allen_reference, group.by = "cell_subtype", label = TRUE)


library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library('glmGamPoi')
library(ggplot2)
library(cowplot)


#处理空间样本
a<-Load10X_Spatial(
  data.dir="GSE250636/GSM7983366",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = T
)
brain <- SCTransform(a, assay = "Spatial", verbose = FALSE)
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

# the annotation is stored in the 'subclass' column of object metadata
anchors <- FindTransferAnchors(reference = allen_reference, query = brain, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$cell_subtype, prediction.assay = TRUE,
                                  weight.reduction = brain[["pca"]], dims = 1:30)
brain[["predictions"]] <- predictions.assay
DefaultAssay(brain) <- "predictions"
saveRDS(brain,"GSM7983366.rds")


#Figure 5AB-left
library(SeuratObject)
dq<-lapply(ss[1:9],function(x){
  readRDS(x)
})

SpatialDimPlot(dq[[1]], crop = TRUE, label = TRUE,pt.size.factor =2.5)

SpatialDimPlot(dq[[2]], crop = TRUE, label = TRUE, pt.size.factor = 1.8) +
  scale_fill_manual(values = c("#e5c494","#B3DE69","#53b400","#f8766d","#00b6eb","#FFD92F","#fb61d7"))

SpatialDimPlot(dq[[3]], crop = TRUE, label = TRUE, pt.size.factor = 1.6) 

SpatialDimPlot(dq[[4]], crop = TRUE, label = TRUE, pt.size.factor = 1.6) 

SpatialDimPlot(dq[[5]], crop = TRUE, label = TRUE, pt.size.factor = 1.6) 

SpatialDimPlot(dq[[6]], crop = TRUE, label = TRUE, pt.size.factor = 1.6)+
  scale_fill_manual(values = c("#f8766d","#e5c494","#53b400","#00b6eb","#B3DE69"))

SpatialDimPlot(dq[[7]], crop = TRUE, label = TRUE, pt.size.factor = 1.6)+
  scale_fill_manual(values = c("#f8766d","#e5c494","#53b400","#00b6eb","#FFD92F"))

SpatialDimPlot(dq[[8]], crop = TRUE, label = TRUE, pt.size.factor = 1.6)

SpatialDimPlot(dq[[9]], crop = TRUE, label = TRUE, pt.size.factor = 1.6)




#Figure 5AB-right
library(ggplot2)
library(ggforce)
library("ggnewscale")
circle2Fun <- function(center = c(0,0), diameter = 0.2, npoints = 50, n_circle = 1){
  r = diameter / 2
  tt <- seq(pi/2,5*pi/2,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  opp_point<-c(center[1] + r * cos(3*pi/2),center[2] + r * sin(3*pi/2))
  lab <- rep(paste(n_circle,1:2,sep="-"),each=npoints/2)
  circle_pos<-data.frame(x = c(xx,opp_point[1],rep(center[1],2),opp_point[1]), y = c(yy,opp_point[2],rep(center[2],2),opp_point[2]),n_circle,each_piece=c(rep(1:2,each=npoints/2),c(1,1,2,2)),bq=c(lab,paste(n_circle,rep(1:2,each=2),sep="-")))
  return(circle_pos)
}


skcm<-read.table("cancer_nets\\TCGA-SKCM.txt",sep="\t",header=T)
stage4<-skcm[skcm$stage=="Stg IV",]
needgene<-unique(c(stage4$lnc_name,stage4$gene_name))

ss<-list.files("cancer_nets\\",full.names=T)
j<-0
kkk<-lapply(ss[1:9],function(ff){
  j<<-j+1
  brain<-readRDS(ff)
  #pdf_names<-sprintf("D:\\sam%s.pdf",j)
  #pdf(pdf_names)
  
  gene_scores<-brain@assays[["SCT"]]@data
  needgene<-needgene[needgene %in% rownames(gene_scores)]
  data<-gene_scores[match(needgene,rownames(gene_scores)),]
  inskcm<-skcm[skcm$lnc_name %in% needgene & skcm$gene_name %in% needgene,]
  spatial_cluster<-brain@meta.data[["seurat_clusters"]]
  

  d1<-brain@images[["slice1"]]@coordinates
  d1<-as.data.frame(d1)
  d1$col<-(d1$col)/1.5  

  
  
  i<-0
  test<-apply(inskcm[,c("lnc_name","gene_name")],1,function(y){
    i<<-i+1
    #print(i)
    y<-as.character(y)
    a<-data[which(rownames(data)==y[1]),]
    b<-data[which(rownames(data)==y[2]),]
    w1<-length(intersect(which(a!=0),which(b!=0)))
    w2<-length(intersect(which(a==0),which(b!=0)))
    w3<-length(intersect(which(a!=0),which(b==0)))
    w4<-length(intersect(which(a==0),which(b==0)))
    bb <- matrix(c(w1,w2,w3,w4), nrow = 2, byrow = TRUE)
    p<-fisher.test(bb)$p.value
    if(!is.na(p) && p<0.05){
      eval(parse(text=sprintf("pos<-data.frame(d1[,2:3],%s=a,%s=b)",gsub("-","_",y[1]),gsub("-","_",y[2]))))
      jj<-pos[pos[,3]!=0 & pos[,4]!=0,]
      if(dim(jj)[1]!=0){
        pp<-apply(jj,1,function(mm){
          data.frame(circle2Fun(center=mm[1:2],n_circle = paste(mm[1],mm[2],sep="-"),diameter = 1.1, npoints = 130),mm[3],mm[4])
        })
        all_pos<-do.call(rbind.data.frame,pp)
        colnames(all_pos)[6:7]<-gsub("-","_",y)
        all_info<-paste0("cor in tcga:",round(inskcm$Stg.IV[i],2),";p:",round(p,2),";or:",round((w1+w4)/(w2+w3),2))
        eval(parse(text=sprintf('p1<-ggplot() + 
          geom_shape(data=all_pos[all_pos$each_piece==1,],aes(y,-x,fill = %s, group = bq))+
          geom_point(data=pos[pos[,3]!=0 & pos[,4]==0,],aes(col,-row,fill = %s),shape=21,size=2.2,color="white")+
          scale_fill_gradient2(low="white",high="red")+
          new_scale("fill") +
          geom_point(data=pos[pos[,3]==0 & pos[,4]!=0,],aes(col,-row,fill = %s),shape=21,size=2.2,color="white")+
          geom_shape(data=all_pos[all_pos$each_piece==2,],aes(y,-x,fill = %s, group = bq))+
          scale_fill_gradient2(low="white",high="blue")+
          theme_void()+
          coord_fixed()+
          ggtitle(label=all_info)',gsub("-","_",y[1]),gsub("-","_",y[1]),gsub("-","_",y[2]),gsub("-","_",y[2]))))
        print(p1)
      } else {
        all_info<-paste0("cor in tcga:",round(inskcm$Stg.IV[i],2),";p:",round(p,2),";or:",round((w1+w4)/(w2+w3),2))
        eval(parse(text=sprintf('p1<-ggplot() + 
          geom_point(data=pos[pos[,3]!=0 & pos[,4]==0,],aes(col,-row,fill = %s),shape=21,size=2.2,color="white")+
          scale_fill_gradient2(low="white",high="red")+
          new_scale("fill") +
          geom_point(data=pos[pos[,3]==0 & pos[,4]!=0,],aes(col,-row,fill = %s),shape=21,size=2.2,color="white")+
          scale_fill_gradient2(low="white",high="blue")+
          theme_void()+
          coord_fixed()+
          ggtitle(label=all_info)',gsub("-","_",y[1]),gsub("-","_",y[1]),gsub("-","_",y[2]),gsub("-","_",y[2]))))
        print(p1)
      }
      return(c(y,all_info))
    }
      
  })
  dev.off()
  return(test)
})




brain<-readRDS(ss[1])
d1<-brain@images[["slice1"]]@coordinates
d1<-as.data.frame(d1)
gene_scores<-a@assays[["SCT"]]@data
d1<-data.frame(d1,SNHG1=gene_scores[rownames(gene_scores)=="SNHG1",],IRF4=gene_scores[rownames(gene_scores)=="IRF4",])
d1<-data.frame(d1,spatial_cluster)
table(spatial_cluster)
#生成所有坐标
pos<-lapply(1:dim(d1[d1$SNHG1!=0 & d1$IRF4!=0,])[1],function(x){
  x<-as.numeric(d1[d1$SNHG1!=0 & d1$IRF4!=0,][x,])
  
  data.frame(circle2Fun(center=x[2:3],n_circle = paste(x[2],x[3],sep="-"),diameter = 0.7, npoints = 100),SNHG1=x[6],IRF4=x[7])
})
all_pos<-do.call(rbind,pos)


ggplot() + 
  geom_shape(data=all_pos[all_pos$each_piece==1,],aes(y,-x,fill = SNHG1, group = bq))+
  geom_point(data=d1[d1$SNHG1!=0 & d1$IRF4==0,],aes(col,-row,fill = SNHG1),shape=21,size=2,color="white")+
  scale_fill_gradient2(low="white",high="red")+
  new_scale("fill") +
  geom_point(data=d1[d1$SNHG1==0 & d1$IRF4!=0,],aes(col,-row,fill = IRF4),shape=21,size=2,color="white")+
  geom_shape(data=all_pos[all_pos$each_piece==2,],aes(y,-x,fill = IRF4, group = bq))+
  scale_fill_gradient2(low="white",high="blue")+
  theme_void()+
  coord_fixed()



#Figure 5C
cell_scores<-brain@assays[["predictions"]]@data[-25,]
spatial_cluster<-brain@meta.data[["seurat_clusters"]]
r1<-apply(cell_scores,1,function(x){
  table(spatial_cluster[which(x!=0)])
})


cells<-r1[,colSums(r1)!=0]
cells<-cells[,-5]
cell_cluster<-melt(data.frame(cells,cluster=paste0("Cluster ",rownames(cells))))
cell_cluster<-cell_cluster[cell_cluster$value!=0,]
test<-apply(cell_cluster,1,function(x){
  a<-which(as.character(spatial_cluster)==substring(x[1],9,9))
  ce<-unlist(strsplit(as.character(x[2]),"[.]"))[1]
  nd<-cell_scores[grep(ce,rownames(cell_scores)),a]
  mean(nd[nd!=0])
})
pdd<-data.frame(cell_cluster,mean_value=test)
tapply(1:dim(pdd)[1],pdd$variable,function(x){
  pdd[x,4]<<-scale(pdd[x,4])
  return()
})

ggplot(pdd,aes(x=cluster,y=variable))+
  geom_point(aes(size=value,color=mean_value))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90))+
  scale_color_gradient(low="#ffdbd2",high="#e24f39")




#Figure 5D
library(SeuratObject)
library(Seurat)
brain<-readRDS("cancer_nets\\GSM7983360.rds")

SpatialFeaturePlot(brain, features = c("Plasma cells","Tumor cells","CD4+ T cells","Oligodendrocytes","Endothelial cells","Astrocytes","Neurons"), 
                   pt.size.factor = 1.6, ncol = 4, crop = TRUE)




#Figure 5E
brain<-readRDS(ss[1])
d1<-brain@images[["slice1"]]@coordinates
d1<-as.data.frame(d1)
gene_scores<-brain@assays[["SCT"]]@data
cell_scores<-brain@assays[["predictions"]]@data
d1<-data.frame(d1,PHLDA3=gene_scores[rownames(gene_scores)=="PHLDA3",],Neurons=cell_scores[rownames(cell_scores)=="Neurons",])
d1<-data.frame(d1,spatial_cluster)
d1$col<-(d1$col)/1.5
table(spatial_cluster)
#生成所有坐标


pos<-lapply(1:dim(d1[d1$PHLDA3!=0 & d1$Neurons!=0,])[1],function(x){
  x<-as.numeric(d1[d1$PHLDA3!=0 & d1$Neurons!=0,][x,])
  
  data.frame(circle2Fun(center=x[2:3],n_circle = paste(x[2],x[3],sep="-"),diameter = 1, npoints = 100),PHLDA3=x[6],Neurons=x[7])
})
all_pos<-do.call(rbind,pos)


ggplot() + 
  geom_shape(data=all_pos[all_pos$each_piece==1,],aes(y,-x,fill = PHLDA3, group = bq))+
  geom_point(data=d1[d1$PHLDA3!=0 & d1$Neurons==0,],aes(col,-row,fill = PHLDA3),shape=21,size=2,color="white")+
  scale_fill_gradient2(low="white",high="red")+
  new_scale("fill") +
  geom_point(data=d1[d1$PHLDA3==0 & d1$Neurons!=0,],aes(col,-row,fill = Neurons),shape=21,size=2,color="white")+
  geom_shape(data=all_pos[all_pos$each_piece==2,],aes(y,-x,fill = Neurons, group = bq))+
  scale_fill_gradient2(low="white",high="blue")+
  theme_void()+
  coord_fixed()



SpatialFeaturePlot(brain, features = tt2[25:31],alpha = c(0.1, 1),ncol=4)



d1<-brain@images[["slice1"]]@coordinates
d1<-as.data.frame(d1)
d1<-data.frame(d1,MLXIP=gene_scores[rownames(gene_scores)=="MLXIP",],CD4=cell_scores[rownames(cell_scores)=="CD4+ T cells",])
d1<-data.frame(d1,spatial_cluster)
d1$col<-(d1$col)/1.5
table(spatial_cluster)
#生成所有坐标

pos<-lapply(1:dim(d1[d1$MLXIP!=0 & d1$CD4!=0,])[1],function(x){
  x<-as.numeric(d1[d1$MLXIP!=0 & d1$CD4!=0,][x,])
  
  data.frame(circle2Fun(center=x[2:3],n_circle = paste(x[2],x[3],sep="-"),diameter = 1, npoints = 100),MLXIP=x[6],CD4=x[7])
})
all_pos<-do.call(rbind,pos)


ggplot() + 
  geom_shape(data=all_pos[all_pos$each_piece==1,],aes(y,-x,fill = MLXIP, group = bq))+
  geom_point(data=d1[d1$MLXIP!=0 & d1$CD4==0,],aes(col,-row,fill = MLXIP),shape=21,size=2,color="white")+
  scale_fill_gradient2(low="white",high="red")+
  new_scale("fill") +
  geom_point(data=d1[d1$MLXIP==0 & d1$CD4!=0,],aes(col,-row,fill = CD4),shape=21,size=2,color="white")+
  geom_shape(data=all_pos[all_pos$each_piece==2,],aes(y,-x,fill = CD4, group = bq))+
  scale_fill_gradient2(low="white",high="blue")+
  theme_void()+
  coord_fixed()

