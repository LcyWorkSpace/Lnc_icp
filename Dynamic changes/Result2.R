#Figure 2A  癌症在各个stage中在9个cluster中的分布
nodes<-read.csv("all_nodes_in_cluster.csv")             #网络所有节点        
a<-list.files("cancer_nets\\",full.names=T)
library(stringr)
cancer_name<-str_extract(a,"TCGA-(\\w+).txt",group=1)
bb<-lapply(a,function(x){
  b<-read.table(x,sep="\t",header=T)
  d<-rbind(as.matrix(b[,c("lnc_id","stage","lnc_name")]),as.matrix(b[,c("gene_id","stage","gene_name")]))
  d<-unique(d)
  k<-tapply(d[,3],d[,2],function(y){
    label_c<-rep(NA,9)
    cluster<-nodes[match(y,nodes[,2]),1]
    ff<-table(cluster)
    label_c[match(as.numeric(names(ff)),1:9)]<-ff
    label_c
  })
  do.call(cbind.data.frame,k)
})
stages<-lapply(c("Stg I","Stg II","Stg III","Stg IV"),function(y){
  s1<-lapply(bb,function(m){
    if(!is.null(m[[y]])){
      m[[y]]
    } else {return(rep(NA,9))}
  })
  all_one_stage<-do.call(cbind.data.frame,s1)
  colnames(all_one_stage)<-cancer_name
  all_one_stage
})
library("pheatmap")
all_cluster_stage<-do.call(rbind,stages)
all_cluster_stage<-t(all_cluster_stage)
f1<-all_cluster_stage
f1[is.na(f1)]<-""

pheatmap(log2(all_cluster_stage),cluster_cols = F,cluster_rows = F,display_numbers = f1,
         gaps_col = c(9,18,27),na_col = "#ffffff",legend = F,fontsize = 12)


#各个cluster的节点个数
bap<-data.frame(table(nodes$glayCluster))
library(ggplot2)
ggplot(bap,aes(x=Var1,y=log2(Freq)))+
  geom_bar(stat="identity",fill="#f19085")+
  theme_void()+
  geom_text(stat = "identity", aes(label = Freq), vjust = -0.5)



#Figure 2C
rr<-readRDS("cluster_function.rds")
pcc<-lapply(rr, function(x){
  if(dim(x)[1]!=0){
    if(length(grep("REACTOME",x$ID))!=0){
      a<-match(gsub("REACTOME_","",x$ID),ReactomePathways2$V2)
      aa<-match(ReactomePathways2[a,1],p2t[,1])
      bb<-match(ReactomePathways2[a,1],p2t[,3])
      aa[is.na(aa)]<-bb[!is.na(bb)]
      r1<-na.omit(p2t[aa,c(2,4)])
      return(data.frame(r1))
    } else {
      data.frame(TopicName=x$ID)
    }
  }
})

m<-Reduce(function(x,y){
  if(length(y)!=0){
    a<-table(unlist(y$TopicName))
    b<-data.frame(pathway=names(a),num=as.numeric(a))
    merge(x,b,by="pathway",all=T)
  } else {merge(x,data.frame(pathway=character(),num=numeric()),by="pathway",all=T)}
},pcc,init=data.frame(pathway=character()))
colnames(m)[-1]<-paste0("cluster",1:9)
rownames(m)<-gsub(' \\(Homo sapiens\\)',"",m[,1])

df_sorted <- m[do.call(order, m[,-1]), ]
f1<-df_sorted[,-1]
f1[is.na(f1)]<-""

library(pheatmap)
pheatmap(log2(df_sorted[,-1]),cluster_rows = F,cluster_cols = F,display_numbers = f1,
         na_col = "#ffffff",legend = F,fontsize = 12)


#Figure 2D
a<-list.files("cancer_nets\\",full.names=T)
kegg<-download_KEGG("hsa", keggType = "KEGG", keyType = "kegg")

cc<-lapply(a,function(x){
  b<-read.table(x,sep="\t",header=T)
  d<-rbind(as.matrix(b[,c("lnc_id","stage","lnc_name")]),as.matrix(b[,c("gene_id","stage","gene_name")]))
  d<-unique(d)
  cluster<-nodes[match(d[,3],nodes[,2]),1]
  m<-data.frame(d,cluster)
  stage_cluster<-apply(m,1,function(h){paste0(h[c(2,4)],collapse="_")})
  rr<-tapply(m$lnc_id,stage_cluster,function(kk){
    if(length(kk)>=3){
      gene1<-bitr(kk,fromType="ENSEMBL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
      js<-kegg[["KEGGPATHID2EXTID"]][which(kegg[["KEGGPATHID2EXTID"]][,1]=="hsa03015"),2]
      ensg<-gene1[gene1[,2] %in% js,1]
      unique(m[m$lnc_id %in% ensg,"lnc_name"])
    } else {return(NULL)}
  })
})
gg<-nodes[match(unique(unlist(cc)),nodes$name),]    #通路交集基因属于哪些cluster



#Figure 2E
all_plot<-readRDS("all_plot.rds")

void_plot<-function(){
  ggplot() +
    theme_void()
}

non_void_plot<-function(x){
  sx<-sort(unique(x$cluster))
  ymin<-0
  ymax<-1
  testp<-lapply(sx,function(cluster){
    cc<-paste0("c",cluster)
    z1<-sprintf('geom_rect(aes(xmin=0,xmax=3,ymin=%s,ymax=%s),fill="%s")',ymin,ymax,coll[match(cc,coll[,1]),2])
    ymin<<-ymin+1
    ymax<<-ymax+1
    z1
  })
  paste0('ggplot()+',paste0(testp,collapse = "+"),'+coord_polar("y", start = 0)+theme_void()')
}

coll<-rbind(
  c("c1","#8DD3C7"),
  c("c2","#FFFFB3"),
  c("c3","#BEBADA"),
  c("c4","#FB8072"),
  c("c5","#80B1D3"),
  c("c6","#FDB462"),
  c("c7","#B3DE69"),
  c("c8","#FCCDE5"),
  c("c9","#D9D9D9"))

tttt<-lapply(all_plot,function(da){
  if(dim(da)[1]!=0){
    eval(parse(text=non_void_plot(da)))
  }else{
    void_plot()
  }
})

library(patchwork)
wrap_plots(tttt, nrow = 4,byrow=F)


#Figure 2F
ff<-readRDS("figure_2f.rds")
ggplot(ff,aes(clu,stg))+
  geom_rect(aes(xmin=0,xmax=10,ymin=0,ymax=3.5),fill="#fcdfe8")+ 
  geom_rect(aes(xmin=0,xmax=10,ymin=3.5,ymax=6.4),fill="#e9c1d5")+
  geom_rect(aes(xmin=0,xmax=10,ymin=6.4,ymax=8.5),fill="#b8abc1")+
  geom_rect(aes(xmin=0,xmax=10,ymin=8.5,ymax=11.5),fill="#95dadf")+
  geom_rect(aes(xmin=0,xmax=10,ymin=11.5,ymax=13.5),fill="#99d2ad")+
  geom_rect(aes(xmin=0,xmax=10,ymin=13.5,ymax=17.5),fill="#e4eaf1")+
  geom_rect(aes(xmin=0,xmax=10,ymin=17.5,ymax=21.5),fill="#fbd2c9")+
  geom_rect(aes(xmin=0,xmax=10,ymin=21.5,ymax=24.5),fill="#dcebde")+
  geom_rect(aes(xmin=0,xmax=10,ymin=24.5,ymax=28.5),fill="#f2f4d4")+
  geom_rect(aes(xmin=0,xmax=10,ymin=28.5,ymax=29.5),fill="#e1b8b8")+
  geom_rect(aes(xmin=0,xmax=10,ymin=29.5,ymax=31.5),fill="#ffb8b8")+
  geom_rect(aes(xmin=0,xmax=10,ymin=31.5,ymax=34.5),fill="#f7aba2")+
  geom_rect(aes(xmin=0,xmax=10,ymin=34.5,ymax=37.5),fill="#bbb9df")+
  geom_rect(aes(xmin=0,xmax=10,ymin=37.5,ymax=41.5),fill="#fff8e9")+
  geom_rect(aes(xmin=0,xmax=10,ymin=41.5,ymax=45.5),fill="#dcd7d5")+
  geom_point(aes(size=Freq,col=Freq))+
  theme_bw()+
  theme(axis.text.y=element_text(size=rel(1),hjust = 1))+
  theme(axis.ticks.length =unit(0,"cm"))+
  scale_colour_gradient(low="moccasin",high="orangered3")+
  theme(panel.grid = element_blank(),
        axis.title = element_blank())



