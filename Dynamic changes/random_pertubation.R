library("rio")

clinlic_files<-list.files(path="临床tsv\\",full.names =T,recursive=T)
clin<-import_list(clinlic_files)
#筛选有stage的样本   每个癌型的stage列不同，这里匹配定义一下
col_of_stage<-lapply(clin,function(x){grep("Stage",x)})
col_of_stage[[8]]<-28
col_of_stage[[10]]<-28
col_of_stage[[13]]<-28
col_of_stage[[28]]<-28
col_of_stage[[33]]<-28

#只获取需要计算的癌种
s<-list.files("新互作边\\")
ns<-gsub(".txt","",s)

sample_And_Stage<-lapply(jsd,function(aa){
  x<-clin[[aa]]
  m1<-unique(x[,c(2,col_of_stage[[aa]])])
  m2<-m1[grep("Stage",m1[,2]),]
  m2
})

#单独处理乳腺癌
rxa<-sample_And_Stage[[3]]
sample_And_Stage[[3]]<-rxa[-grep("X",rxa[,2]),]
names(sample_And_Stage)<-names(clin)[jsd]

for(i in ns){               
  lujing1<-paste0("整合完的表达谱\\",i,'.txt')
  aa<-read.table(lujing1,header = T,row.names = 1,check=F)
  w<-grep(i,names(sample_And_Stage))
  
  
  #临床信息里有stage的样本和表达谱里的样本有一些对不上，这里重新取交集
  bothIn<-intersect(sample_And_Stage[[w]][,1],colnames(aa))
  ins<-match(bothIn,sample_And_Stage[[w]][,1])
  inaa<-match(bothIn,colnames(aa))
  profile_with_stage<-aa[,inaa]
  rm(aa)
  ss<-sample_And_Stage[[w]][ins,2]
  st1<-gsub("[A-Ea-e]","",ss)
  stagess<-gsub("[0-9]","",st1)
  stagess<-gsub("Stg IS","Stg I",stagess)
  
  als<-table(stagess)
  remain_stage<-names(als)[als>=10]
  if(length(remain_stage)>1){
    rema<-stagess %in% remain_stage
    new_stage<-stagess[rema]
    new_profile<-profile_with_stage[,rema]
    rownames(new_profile)<-sapply(rownames(new_profile),function(k){unlist(strsplit(k,"[.]"))[1]})
    rm(profile_with_stage)
    lujing5<-paste0("新互作边\\",i,".txt")
    a<-read.table(lujing5,sep="\t",header=T)
    ngs<-unique(c(a[,1],a[,2]))
    npf<-new_profile[rownames(new_profile) %in% ngs,]
    lujing2<-paste0("D:\\处理过的表达谱\\",i,".rds")
    saveRDS(npf,file=lujing2)
    lujing3<-paste0("D:\\表达谱配套stage\\",i,".txt")
    write.table(new_stage,file=lujing3,row.names = F,col.names = F,sep="\n",quote=F)
    
  }
}



#计算相关性
s<-list.files("D:\\新互作边\\")
ns<-gsub(".txt","",s)
set.seed(2023)

for(i in ns){ 
  
  lujing2<-paste0("D:\\处理过的表达谱\\",i,".rds")
  b<-readRDS(lujing2)
  nb<-rownames(b)
  lujing3<-paste0("D:\\表达谱配套stage\\",i,".txt")
  sg<-read.table(file=lujing3,sep="\n")
  lujing4<-paste0("D:\\新互作边\\",i,".txt")
  bia<-read.table(file=lujing4,sep="\t",header=T)
  
  #计算真实网络的评分
  p_in_lie<-grep("p.value",colnames(bia))
  pvals<-bia[,p_in_lie]
  cors<-bia[,p_in_lie-1]
  lie<-dim(pvals)[2]
  nns<-gsub("..estimate","",colnames(cors))
  nns<-gsub("[.]"," ",nns)
  okdata<-sapply(1:dim(pvals)[1],function(z1){
    val<-rep(0,lie)
    xz<-which(pvals[z1,]<0.05)
    if(length(xz)!=0){
      val<-cors[z1,]
      val[setdiff(1:lie,xz)]<-0
    }
    return(val)
  })

  zong_sc<-apply(okdata,2,function(val){
    zong_score<-0
    val<-abs(unlist(val))
    if(length(which(val!=0))!=0){
      aa<-val/max(val)
      specific_score<-(lie-sum(aa))/(lie-1)
      strong_score<-max(val)/0.75
      zong_score<-specific_score*strong_score
    }
    return(zong_score)
  })

  nl<-apply(okdata,2,function(val){
    lx<-"none"
    if(length(which(val!=0))!=0){
      aa<-abs(unlist(val))
      lx<-nns[which.max(aa)]
    }
    return(lx)
  })

  bian<-bia[,1:2]

  #随机循环开始
  
  sjzf<-sapply(1:1000,function(sjcs){
 
  sg1<-sg[sample(dim(sg)[1]),1]
  stageee<-tapply(1:length(sg1),sg1,function(stg){stg})
  mm<-apply(bian,1,function(bian){
    fg<-which(nb==bian[1])
    yh<-which(nb==bian[2])
    if(length(fg)!=0 & length(yh)!=0){
      sapply(1:length(stageee),function(kk){
        xx<-stageee[[kk]]
        d<-cor.test(as.numeric(b[fg,xx]),as.numeric(b[yh,xx]),method="pearson")
        return(c(d$estimate,d$p.value))
      })
      }
  })
  rownames(mm)<-paste(rep(names(stageee),each=2),c("-estimate","-p.value"))
  #计算随机情况的评分
  mm<-t(mm)
  p_in_lie2<-grep("p.value",colnames(mm))
  pvals2<-mm[,p_in_lie2]
  cors2<-mm[,p_in_lie2-1]
  lie2<-dim(pvals2)[2]
  nns<-gsub(" -estimate","",colnames(cors2))
  okdata2<-sapply(1:dim(pvals2)[1],function(z1){
    val<-rep(0,lie2)
    xz<-which(pvals2[z1,]<0.05)
    if(length(xz)!=0){
      val<-cors2[z1,]
      val[setdiff(1:lie2,xz)]<-0
    }
    return(val)
  })
  
  
  zong_sc<-apply(okdata2,2,function(val){
    zong_score<-0
    val<-abs(unlist(val))
    if(length(which(val!=0))!=0){
      aa<-val/max(val)
      specific_score<-(lie2-sum(aa))/(lie2-1)
      strong_score<-max(val)/0.75
      zong_score<-specific_score*strong_score
    }
      return(zong_score)
  })
  
  nl<-apply(okdata2,2,function(val){
    lx<-"none"
    if(length(which(val!=0))!=0){
      aa<-abs(val)
      lx<-nns[which.max(aa)]
    }
    return(lx)
  })
  
  return(data.frame(zong_sc,nl))
  })
  lujing<-paste0("D:\\随机矩阵\\",i,".rds")
  saveRDS(sjzf,file=lujing)
    
}


#统计随机结果，输出显著对

library(stringr)
tcga<-read.table("D:\\tcga.txt",sep="\t",header=T)

s<-list.files("D:\\新互作边\\")
ns<-gsub(".txt","",s)

test<-lapply(ns,function(i){
  
  lujing4<-paste0("D:\\新互作边\\",i,".txt")
  bia<-read.table(file=lujing4,sep="\t",header=T)
  #计算真实网络的评分
  p_in_lie<-grep("p.value",colnames(bia))
  pvals<-bia[,p_in_lie]
  cors<-bia[,p_in_lie-1]
  lie<-dim(pvals)[2]
  nns<-gsub("..estimate","",colnames(cors))
  nns<-gsub("[.]"," ",nns)
  okdata<-sapply(1:dim(pvals)[1],function(z1){
    val<-rep(0,lie)
    xz<-which(pvals[z1,]<0.05)
    if(length(xz)!=0){
      val<-cors[z1,]
      val[setdiff(1:lie,xz)]<-0
    }
    return(val)
  })
  
  zong_sc<-apply(okdata,2,function(val){
    zong_score<-0
    val<-abs(unlist(val))
    if(length(which(val!=0))!=0){
      aa<-val/max(val)
      specific_score<-(lie-sum(aa))/(lie-1)
      strong_score<-max(val)/0.75
      zong_score<-specific_score*strong_score
    }
    return(zong_score)
  })
  
  nl<-apply(okdata,2,function(val){
    lx<-"none"
    if(length(which(val!=0))!=0){
      aa<-abs(unlist(val))
      lx<-nns[which.max(aa)]
    }
    return(lx)
  })
  
  #开始与随机结果比较
  lujing_s<-paste0("D:\\随机矩阵\\",i,".rds")
  sdfg<-readRDS(lujing_s)
  sjpf<-do.call(cbind,sdfg[1,])
  sjsl<-do.call(cbind,sdfg[2,])
  bj<-sapply(1:length(nl),function(x){
    if(nl[x]!="none"){
      w<-which(sjsl[x,]==nl[x])
      return(length(which(sjpf[x,w]>zong_sc[x]))/1000)
    }
  })
  
  dd<-sapply(1:length(bj),function(s){
    if(!is.null(bj[[s]])){
      if(bj[[s]]<0.05 & zong_sc[s]>0.5){
        return(s)
      }
    }
  })
  
  ud<-unlist(dd)
  nets<-bia[ud,1:2]
  lb<-nl[ud]
  lncn<-tcga[match(bia[ud,1],tcga[,1]),2]
  genen<-tcga[match(bia[ud,2],tcga[,1]),2]
  ccrs<-matrix(unlist(t(okdata)[ud,]),ncol=dim(cors)[2])
  
  data<-data.frame(bia[ud,1:2],ccrs,stage=nl[ud],zong_sc[ud],lncn,genen)
  colnames(data)<-c("lnc_id","gene_id",str_extract(colnames(cors),"(.+)\\.\\.",group=1),"stage","score","lnc_name","gene_name")
  output_path<-paste0("D:\\重新生成网络\\",i,".txt")
  write.table(data,file=output_path,sep="\t",quote=F,row.names = F,col.names = T)
  return()
  
})


