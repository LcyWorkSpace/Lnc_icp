#批量提取TPM值

library("rio")

cancers<-c("TCGA-ACC","TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-CHOL","TCGA-COAD","TCGA-DLBC","TCGA-ESCA","TCGA-GBM","TCGA-HNSC","TCGA-KICH",
           "TCGA-KIRC","TCGA-KIRP","TCGA-LAML","TCGA-LGG","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-MESO","TCGA-OV","TCGA-PAAD","TCGA-PCPG",
           "TCGA-PRAD","TCGA-READ","TCGA-SARC","TCGA-SKCM","TCGA-STAD","TCGA-TGCT","TCGA-THCA","TCGA-THYM","TCGA-UCEC","TCGA-UCS","TCGA-UVM")

for(i in cancers){
  
  lujing1<-paste0("D:\\TCGA数据\\表达谱\\",i,'\\')
  acc_files<-list.files(path=lujing1,full.names =T,recursive=T)
  
  acc<-import_list(acc_files)
  ttt<-sapply(acc,function(x){x[,7]})
  rownames(ttt)<-acc[[1]][,1]
  rm(acc)
  
  lujing2<-paste0("D:\\samplesheets\\",i,'.tsv')
  acc_samplesheets<-read.table(file=lujing2,sep="\t",header=T)
  
  sampleID<-sapply(colnames(ttt),function(n){
    l<-grep(n,acc_samplesheets[,2])
    acc_samplesheets[l,7]
  })
  
  colnames(ttt)<-sampleID
  ttt<-ttt[-(1:4),]                    #整合完所有样本
  
  #去掉正常样本
  b<-sapply(colnames(ttt),function(x){substr(x,14,15)})
  tum<-ttt[,which(as.numeric(b)<=10)]
  
  #去掉缺失>30%的基因
  qq<-apply(tum,1,function(x){ 
    length(which(x==0))/dim(tum)[2]
  })
  tum<-tum[-which(qq>0.3),]

  
  #有的样本来自同一个体，如"TCGA-BL-A0C8-01B" "TCGA-BL-A0C8-01A" "TCGA-BL-A0C8-01A"，求几列均值
  samm<-sapply(colnames(tum),function(x){substr(x,1,12)})
  test1<-apply(tum,1,function(hh){
    tapply(1:length(samm),samm,function(yyy){mean(hh[yyy])})
  })
  test1<-t(test1)
  
  #用行均值补缺失值
  buwanqueshi<-apply(test1,1,function(x){
    x[x==0]<-mean(x[x!=0])
    x
  })
  buwanqueshi<-t(buwanqueshi)
  
  lujing3<-paste0("D:\\整合完的表达谱\\",i,'.txt')
  write.table(buwanqueshi,file=lujing3,sep="\t",quote = F)

}








