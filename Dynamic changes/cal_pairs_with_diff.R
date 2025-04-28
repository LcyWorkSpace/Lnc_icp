library("rio")

#整理样本的stage
clinlic_files<-list.files(path="临床tsv\\",full.names =T,recursive=T)
clin<-import_list(clinlic_files)
#筛选有stage的样本   每个癌型的stage列不同，这里匹配定义一下
col_of_stage<-lapply(clin,function(x){grep("Stage",x,ignore.case = T)})
col_of_stage[[8]]<-28
col_of_stage[[10]]<-28
col_of_stage[[13]]<-28
col_of_stage[[28]]<-28
col_of_stage[[33]]<-28

#暂时去掉stage不明晰的数据集
jsd<-setdiff(1:33,c(9,14,15,22,23,25))

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

cancers<-read.table("D:\\cancers.txt")       #癌症名称

#阶段差异
for(i in cancers[jsd,1]){               
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
  #去掉样本小于10的阶段
  remain_stage<-names(als)[als>=10]
  if(length(remain_stage)>1){
    rema<-stagess %in% remain_stage
    new_stage<-stagess[rema]
    new_profile<-profile_with_stage[,rema]
    rm(profile_with_stage)
    number_of_stage<-length(unique(new_stage))
    calc_fangcha<-apply(new_profile,1,function(x){
      h<-log2(x)
      zhengtai<-tapply(h,new_stage,function(xiaozu){shapiro.test(xiaozu)$p.value})
      if(length(which(zhengtai>0.05))==number_of_stage){
        qici<-bartlett.test(h~new_stage)$p.value
        if(qici>0.05){
          return(1)
        } else{return(0)}
      } else{return(0)}
    })
    #查看分布，大于70%的基因符合正态+方差齐就方差分析
    if(sum(calc_fangcha)>=dim(new_profile)[1]*0.7){
      pvals<-apply(new_profile,1,function(n){
        h<-log2(n)
        e<-aov(h~new_stage)
        p_of_fangcha<-summary(e)[[1]][1,5]
        p_of_fangcha
      })
    } else{
      pvals<-apply(new_profile,1,function(n){
        kruskal.test(n~new_stage)$p.value
      })
    }
  
    b<-sapply(rownames(new_profile),function(x){unlist(strsplit(x,"[.]"))[1]})
    a<-data.frame(genename=b,pval=pvals)
    lujing2<-paste0("所有基因阶段差异检验结果\\",i,".txt")
    write.table(a,file=lujing2,row.names = T,col.names = T,sep="\t",quote=F)
  }
}



#计算相关性
library("qvalue")
okNet<-readRDS("提取的lnc-gene互作总网.rds")
all_genes<-unique(c(okNet[,2],okNet[,3]))
cg<-list.files(path="所有基因阶段差异检验结果\\",full.names =T,recursive=T)

for(i in cancers[jsd,1]){               
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
    nn<-sapply(rownames(new_profile),function(x){unlist(strsplit(x,"[.]"))[1]})
    rm(profile_with_stage)

    cyf<-cg[grep(i,cg)]
    a<-read.table(file=cyf,sep="\t",header=T)
    b<-a[a[,1] %in% all_genes,]
    
    jiaozheng<-qvalue(b[,2])
    final_pval<-jiaozheng$qvalues
    asd<-which(final_pval<0.05)
    length(asd)
    
    s1<-okNet[,3] %in% b[asd,1]
    s2<-okNet[,2] %in% b[asd,1]
    final_net<-okNet[s1 | s2,2:3]

    #计算差异对的spearman
    lujing2<-paste0("新互作边\\",i,".txt")
    stageee<-tapply(1:length(new_stage),new_stage,function(stg){stg})
    cat("lnc","gene",paste(rep(names(stageee),each=2),c("$estimate","$p.value")),sep="\t",file=lujing2)
    cat("\n",file=lujing2,append = T)
    mm<-apply(final_net,1,function(bian){
      a<-which(nn==bian[1])
      b<-which(nn==bian[2])
      if(length(a)!=0 & length(b)!=0){
        cat(bian[1],bian[2],sep="\t",file=lujing2,append = T)
        for(kk in 1:length(stageee)){
          xx<-stageee[[kk]]
          d<-cor.test(as.numeric(new_profile[a,xx]),as.numeric(new_profile[b,xx]),method="pearson")
          cat("\t",file=lujing2,append = T)  
          cat(as.numeric(d$estimate),d$p.value,file=lujing2,sep="\t",append = T)
        }
        cat("\n",file=lujing2,append = T)
      }
    })
  }
}





