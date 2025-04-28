#Figure 3A
#统计每个癌1个stage共表达，2个stage共表达...
a<-list.files("cancer_nets\\",full.names=T)
bl<-lapply(a,function(x){
  b<-read.table(x,sep="\t",header=T)
  stg_col<-b[,grep("Stg",colnames(b))]
  f0<-apply(stg_col,1,function(h){
    length(which(h!=0))
  })
  d<-as.data.frame(table(f0))
  colnames(d)[1]<-"Category"
  d
})

library(patchwork)
library(ggplot2)

create_pie_chart <- function(data) {
  ggplot(data, aes(x = '', y = Freq, fill = Category)) +  
    geom_bar(stat = "identity", color=NA, width = 1) +  
    scale_fill_manual(values = c("#e86032","#71caf2","#83a0bf","#f6c949")) +  
    coord_polar("y", start = 0) +  
    theme_void() +  
    theme(legend.position = "right") +   
    geom_text(aes(x=-0.5,y=5,label=''))
}

plots <- lapply(bl, function(df) create_pie_chart(df))
wrap_plots(plots, nrow = 1) + 
  plot_layout(guides = "collect")  


#识别共表达模式
a<-list.files("cancer_nets\\",full.names=T)
library(stringr)
cancer_name<-str_extract(a,"TCGA-(\\w+).txt",group=1)

def_val<-function(k){
  if(k>0){s_v<-"+"}
  if(k<0){s_v<-"-"}
  return(s_v)
}

trans_type<-function(x,y){
  x1<-x[which.max(abs(x))]
  y1<-y[which.max(abs(y))]
  s_x<-def_val(x1)
  s_y<-def_val(y1)
  return(paste0(s_x," -> ",s_y))
}

dd<-lapply(a,function(m){
  b<-read.table(m,sep="\t",header=T)
  stgs<-b[,grep("Stg",colnames(b))]
  test<-sapply(1:dim(stgs)[1],function(x){
    h<-stgs[x,]
    early<-h[colnames(h) %in% c("Stg.I","Stg.II")]
    adv<-h[colnames(h) %in% c("Stg.III","Stg.IV")]
    if(length(which(h==0))==dim(stgs)[2]-1){
      return(colnames(stgs)[which(h!=0)])
    } else if(max(abs(early))==0){
      return("advanced specific")
    } else if((length(adv)==0) || (max(abs(adv))==0)){
      return("early specific")
    } else if(!all(early>=0) && !all(early<=0)){
      return("early有正有负")
    } else if(!all(adv>=0) && !all(adv<=0)){
      return("adv有正有负")
    } else return(trans_type(early,adv))
  })
  data.frame(b,test)
})

s1<-lapply(dd,function(x){
  final_label<-x$test
  final_label[final_label %in% c("Stg.I","Stg.II")]<-"early specific"
  final_label[final_label %in% c("Stg.IV","Stg.III")]<-"advanced specific"
  data.frame(x,final_label)
})

s2<-lapply(s1,function(x){
  zf<-which(x$final_label=="adv有正有负")
  if(length(zf)!=0){
    x<-x[-zf,]
    return(x)
  } else {return(x)}
})

#Figure 3C
kk<-Reduce(function(x,y){
  merge(x,as.data.frame(table(y$final_label)),by="Var1",all=T)
},s2,init=data.frame(Var1 = character(), Freq = numeric()))
kk<-kk[,-2]
colnames(kk)[-1]<-cancer_name
rownames(kk)<-kk[,1]
library(pheatmap)
f1<-kk[,-1]
f1[is.na(f1)]<-""

pheatmap(log2(kk[,-1]),cluster_cols = F,cluster_rows = F,display_numbers = f1,
         na_col = "#ffffff",legend = F,fontsize = 12)



#Figure 3D
#各个癌症中各模式中的中心lnc/gene
st<-lapply(s1,function(x){
  aa<-data.frame(gene=c(x$lnc_name,x$gene_name),mod=rep(x$final_label,times=2))
  hh<-tapply(aa$gene,aa$mod,function(y){
    k<-table(y)
    data.frame(nn=names(k),num=as.numeric(k))
  })
  ff<-Reduce(function(m,n){
    merge(m,n,by="nn",all=T)
  },hh,init=data.frame(nn=character()))
  colnames(ff)[-1]<-names(hh)
  ff
})

i<-0
st2<-lapply(st,function(x){
  i<<-i+1
  x1<-data.frame(x[,-1])
  x1[is.na(x1)]<-0
  aa<-apply(x1,2,function(y){
    which(y>=10)
  })
  bb<-x[unique(unlist(aa)),]
  colnames(bb)[-1]<-paste0(cancer_name[i],"_",colnames(bb)[-1])
  if(dim(bb)[2]!=2){
    all_sum<-colSums(bb[,-1],na.rm = T)
    d10<-which(all_sum>=10)
    if(length(d10)>0){
      return(bb[,c(1,d10+1)])
    }
  } else {return(bb)}
})

all_gene10<-Reduce(function(x,y){
  if(!is.null(y)){
    x<-merge(x,y,by="nn",all=T)
  } else {x<-x}
},st2,init=data.frame(nn=character()))

library(ggplot2)
library(reshape2)
library(RColorBrewer)
cb<-melt(all_gene10)
cb$cancer<-sapply(as.character(cb$variable),function(x){unlist(strsplit(x,"_"))[1]})
df_sorted <- all_gene10[do.call(order, all_gene10[,-1]), ]
cb$nn<-factor(cb$nn,levels = rev(df_sorted[,1]))

ggplot(cb,aes(x=variable,y=nn))+
  geom_point(aes(size=log2(value+1),color=log2(value+1)))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90),
        strip.text = element_text(angle = 90))+
  scale_color_gradient2(low="#4575b4",high="#e14b36")+
  facet_grid(~ cancer,drop=T,scales = "free_x",space = "free_x")


tcga<-read.table("D:\\tcga.txt",sep="\t",header=T)
dim(tcga[tcga$gene_name %in% df_sorted[,1],])
tcga[match(df_sorted[,1],tcga$gene_name),c("gene_name","gene_type")]



#泛癌中同一个对在不同癌中的正负不同
a<-list.files("cancer_nets\\",full.names=T)
library(stringr)
cancer_name<-str_extract(a,"TCGA-(\\w+).txt",group=1)
i<-0
dd<-lapply(a,function(m){
  i<<-i+1
  b<-read.table(m,sep="\t",header=T)
  stg<-grep("stage",colnames(b))
  stgs<-b[,grep("Stg",colnames(b))]
  b$stage<-gsub(" ",".",b$stage)
  aa<-apply(b,1,function(x){
    stage_cor<-x[match(x[stg],colnames(b))]
    stage<-x[stg]
    return(data.frame(rev(x)[2],rev(x)[1],stage,stage_cor))
  },simplify = F)
  bb<-do.call(rbind.data.frame,aa)
  colnames(bb)[1:2]<-c("lnc_name","gene_name")
  rownames(bb)<-NULL
  data.frame(bb,cancer=cancer_name[i])
})
all_edges<-do.call(rbind.data.frame,dd)
bind_names<-apply(all_edges,1,function(x){paste0(x[1:2],collapse=";")})
t1<-tapply(1:dim(all_edges)[1],bind_names,function(x){all_edges[x,]})
t2<-sapply(t1,function(x){dim(x)[1]})
t1[which(t2==4)]

#Figure 3F
#统计每个lnc-gene对参与了多少cancer
ttt<-data.frame(p=names(t2),num=as.numeric(t2))
pair<-sapply(ttt$p,function(x){unlist(strsplit(x,";"))})
bbb<-data.frame(t(pair),num=ttt$num)
tj<-tapply(bbb[,3],bbb[,1],function(x){
  data.frame(table(x))
})
me<-Reduce(function(x,y){
  merge(x,y,by="x",all=T)
},tj,init=data.frame(x=factor()))
colnames(me)[-1]<-names(tj)


ggplot(me,aes(x=x,y=log2(MALAT1+1)))+
  geom_bar(aes(fill=x),stat = "identity")+
  geom_text(stat = "identity", aes(label = MALAT1), vjust = -0.5)+
  theme_void()



ggplot(me,aes(x=x,y=log2(CRNDE+1)))+
  geom_bar(aes(fill=x),stat = "identity")+
  geom_text(stat = "identity", aes(label = CRNDE), vjust = -0.5)+
  theme_void()



#Figure 3G
ma<-c("MALAT1;IL11RA",
      "MALAT1;ALKBH5",
      "MALAT1;DCST2",
      "MALAT1;MGRN1",
      "MALAT1;PHLDA3",
      "MALAT1;SPG7")

cr<-c("CRNDE;BEX1",
      "CRNDE;CIAPIN1",
      "CRNDE;COL1A2",
      "CRNDE;COL5A1",
      "CRNDE;HEY1",
      "CRNDE;MGP")

tt1<-lapply(c(ma,cr),function(x){t1[[x]]})
p1<-do.call(rbind.data.frame,tt1)
p1$stage_cor<-as.numeric(p1$stage_cor)
library(ggplot2)

ggplot(p1,aes(x=cancer,y=stage))+
  geom_point(aes(color=stage_cor,size=abs(stage_cor)))+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90),
        strip.text = element_text(angle = 90))+
  scale_color_gradient2(low="#4575b4",high="#e14b36")+
  facet_grid(~ gene_name,drop=T,scales = "free_x",space = "free_x")



#附图
a<-list.files("cancer_nets\\",full.names=T)
g2<-lapply(a,function(x){
  b<-read.table(x,sep="\t",header=T)
  stg<-b[,grep("Stg",colnames(b))]
  mm<-apply(stg,1,function(y){
    paste0(gsub("Stg.","",colnames(stg))[which(y!=0)],collapse = "-")
  })
  as.data.frame(table(mm))
})


hb<-Reduce(function(x,y){
  merge(x,g2[[y]],by="mm",all=T)
},2:15,init=g2[[1]])
rownames(hb)<-hb[,1]
hb<-hb[,-1]
colnames(hb)<-stringr::str_extract(a,"(?<=TCGA-)[A-Z]+")
library(pheatmap)
pheatmap(log2(hb+1),cluster_cols = F,cluster_rows = F,display_numbers = hb,na_col = "#ffffff")


