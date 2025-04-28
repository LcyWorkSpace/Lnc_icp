#整理生存信息
library("rio")
a<-list.files("D:\\临床tsv\\",full.names=T,recursive=T)
b<-import_list(a)
test<-lapply(b,function(x){
  dd1<-x[,c("case_submitter_id","vital_status","days_to_death","days_to_last_follow_up")]
  dd2<-dd1[dd1$vital_status!="Not Reported",]
  q<-intersect(which(dd2$days_to_death=="'--"),which(dd2$days_to_last_follow_up=="'--"))
  if(length(q)!=0){
    dd<-dd2[-q,]
  } else {dd<-dd2}
  os_time<-dd$days_to_death
  os_time[which(os_time=="'--")]<-(dd$days_to_last_follow_up)[which(os_time=="'--")]
  os<-rep(1,dim(dd)[1])
  os[which(dd$vital_status=="Alive")]<-0
  pp<-unique(data.frame(id=dd$case_submitter_id,os,os_time=as.numeric(os_time)))
})
names(test)<-sapply(names(test),function(x){unlist(strsplit(x,"[.]"))[3]})
saveRDS(test,"D:\\处理完的TCGA生存.rds")


#Figure 6A
tcga<-read.table("D:\\tcga.txt",sep="\t",header=T)
ensg<-tcga$gene_id

s<-list.files("D:\\所有癌症四个阶段网络\\",full.name=T)
af<-lapply(c(6,8,10,12),function(x){
  a<-read.table(s[x],sep="\t")
  a
})

a<-list.files("D:\\处理过的表达谱",full.name=T)
b<-list.files("D:\\处理过的表达谱")
bb<-sapply(b,function(x){unlist(strsplit(x,"[.]"))[1]})

ss<-list.files("D:\\表达谱配套stage",full.name=T)

ttest<-function(x,y){
  x<-log2(x)
  y<-log2(y)
  qixing<-bartlett.test(c(x,y)~c(rep("tumor",length(x)),rep("normal",length(y))))$p.value
  if(qixing>0.05){
    pval<-t.test(x,y,var.equal=T,paired=F)$p.value}
  else{
    pval<-t.test(x,y,var.equal=F,paired=F)$p.value}
  pval
}

# 找到lnc+gene
lncs<-tapply(1:dim(af[[4]])[1],af[[4]][,3],function(x){unique(c(af[[4]][x,1],af[[4]][x,2]))})
#找到表达谱
needbdp<-a[match(names(lncs),bb)]
needstage<-ss[match(names(lncs),bb)]
markers<-sapply(1:length(needbdp),function(cancertype){
  bdp<-readRDS(needbdp[cancertype])
  needgene<-bdp[rownames(bdp) %in% ensg[tcga$gene_name %in% lncs[[cancertype]]],]
  st<-read.table(needstage[cancertype],sep="\n")
  dq<-which(st[,1] %in% "Stg IV")
  vals<-apply(needgene,1,function(hang){
    tp<-ttest(as.numeric(hang[dq]),as.numeric(hang[-dq]))
    fc<-mean(hang[dq])/mean(hang[-dq])
    if(tp<0.05){return(fc)}
  })
  vals<-unlist(vals)
  low<-vals[which(vals<=1/1.5)]
  high<-vals[which(vals>=1.5)]
  return(c(low,high))
})
names(markers)<-names(lncs)
saveRDS(markers,"D:\\markers.rds")



library("GSVA")
library(survival)
library("survminer")
sc<-readRDS("D:\\处理完的TCGA生存.rds")


#循环
i<<-1
new_marker<-markers[which(lengths(markers)!=0)]
ttt<-lapply(1:length(new_marker),function(x){
  bdp<-readRDS(needbdp[grep(names(new_marker)[x],needbdp)])
  genelist<-list(low_expression=names(new_marker[[x]])[which(new_marker[[x]]<1)],
                 high_expression=names(new_marker[[x]])[which(new_marker[[x]]>1)])
  test<-gsva(expr=as.matrix(log2(bdp)),genelist,kcdf='Gaussian',method="plage")
  sc_data<-sc[[match(names(new_marker)[x],names(sc))]]
  dysc<<-sc_data[match(colnames(bdp),sc_data[,1]),]
  tit<-paste(names(new_marker)[x],rownames(test),sep="_")
  fz<-apply(test,1,function(y){
    f<<-ifelse(y>median(y),"high","low")                         
    b<-survfit(Surv(os_time,os)~f,data=cbind(dysc,f))
    kk<-ggsurvplot(b,pval=T,title=tit[i])
    i<<-i+1
    print(kk)
  })
  i<<-1
  return()
})




#mlr3模型
library(mlr3pipelines)
library(mlr3verse)
library(mlr3proba)
library(survival)
library(mlr3tuningspaces)
library(mlr3learners)
library(mlr3)
library(mlr3mbo)
library(data.table)
learners_space[["classif.nnet"]] <- ", size = 2"     
lnc <- readRDS("Melanoma-PRJEB23709.rds")
lnc[,1:399]<-as.data.frame(apply(lnc[,-c("response")],2,function(x){log2(x+1)}))
table(lnc$response)
# N  R 
# 33 40 
wrapper_learner <- c("classif.AdaBoostM1", "classif.C50", "classif.IBk", "classif.J48", 
                     "classif.JRip", "classif.LMT", "classif.OneR", "classif.PART", 
                     "classif.abess", "classif.bart", "classif.cforest", "classif.ctree", 
                     "classif.cv_glmnet", "classif.earth", "classif.fnn", "classif.gamboost", 
                     "classif.gausspr", "classif.gbm", "classif.glmboost", "classif.imbalanced_rfsrc", 
                     "classif.kknn", "classif.ksvm", "classif.lda", "classif.liblinear", 
                     "classif.lightgbm", "classif.log_reg", "classif.mob", "classif.multinom", 
                     "classif.naive_bayes", "classif.nnet", "classif.randomForest", 
                     "classif.ranger", "classif.rfsrc", "classif.rpart", "classif.svm", "classif.xgboost",
                     "classif.glmnet", "classif.lssvm")

set.seed(225)
#先留出一部分样本单纯用于验证 留20%
task_valid <- as_task_classif(lnc, target = "response", id = "holdout")
task_valid$set_col_roles("response",c("target","stratum"))
rsmp_holdout = rsmp("holdout", ratio = 0.8)
rsmp_holdout$instantiate(task_valid)
# 查看分配到训练集和测试集的索引
train_set = rsmp_holdout$train_set(1)  # 获取训练集的索引
test_set = rsmp_holdout$test_set(1)    # 获取测试集的索引

#把剩下的80%用于交叉验证
#构建任务
task_train <- as_task_classif(lnc[rsmp_holdout$train_set(1), ], target = "response", id = "mela")
task_train$set_col_roles("response",c("target","stratum"))
rsmp_cv5=rsmp("cv", folds=5)
rsmp_cv5$instantiate(task_train)
#开始遍历每一个分类器
res_list<-lapply(wrapper_learner,function(x){
  eee <- tryCatch({
    instance = fsi(
      task = task_train,
      learner = eval(parse(text = paste0("lrn(x", learners_space[[x]], ", predict_type = 'prob')"))),
      resampling = rsmp_cv5,
      measure = msr("classif.auc"),      
      terminator = trm("evals", n_evals = 100),
      store_benchmark_result = T,
      store_models = T
    )
    fselector = fs("genetic_search")
    fselector$optimize(instance)
    111
  },error = function(e){"error"})
  if(eee!="error" && !is.null(instance$result)) {
    best_auc<-as.data.table(instance$result)[,"classif.auc"]
    if(TRUE %in% (best_auc>=0.75)){      #如果最佳模型auc都不能大于0.75，就不要这个模型
      scores100<-instance$archive$benchmark_result$aggregate(msrs(c("classif.auc","classif.acc")))
      goodmodels_pos<-which(scores100$classif.auc>=0.75)
      features100<-instance$archive$data  
      each_result<-lapply(goodmodels_pos,function(goodmodels){
        learner_id<-scores100$learner_id[goodmodels]
        auc<-scores100$classif.auc[goodmodels]
        acc<-scores100$classif.acc[goodmodels]
        all_in<-unlist(features100$x_domain[goodmodels])
        feature<-names(all_in)[which(all_in)]
        fv<-paste0(feature,collapse=",")
        return(data.frame(learner_id=learner_id,auc=auc,acc=acc,feature=fv,n_feature=length(feature)))
      })
      do.call(rbind.data.frame,each_result)
    }
  }
})
step1<-do.call(rbind.data.frame,res_list)
step1<-readRDS("D:\\225_step1.rds")
step1<-readRDS("D:\\step1.rds")

#处理独立验证数据集
dat<-"Melanoma-GSE115821.rds"
test1 <- readRDS(dat)
test1$response<-as.factor(test1$response)
test1[,non_response_col]<-as.data.frame(apply(test1[,-c("response")],2,function(x){log2(x+1)}))

#预测预留出的测试集+独立验证集
pre<-function(learner_row){
  feature<-sapply(step1[learner_row,"feature"],function(x){unlist(strsplit(x,","))})
  tsk_valid=as_task_classif(lnc[train_set,match(c(feature,"response"),colnames(lnc)),with=F], target = "response", id = "valid")
  lrn_best=lrn(step1$learner_id[learner_row],predict_type = 'prob')
  lrn_best$train(tsk_valid)#训练模型
  data<-lnc[test_set,match(c(feature,"response"),colnames(lnc)),with=F]
  data$response<-as.factor(data$response)
  prediction=lrn_best$predict_newdata(newdata=data,task=tsk_valid)#预测
  pp1<-prediction$score(msrs(c("classif.acc","classif.auc")))
  pred2=lrn_best$predict_newdata(newdata=test1,task=tsk_valid)#预测
  pp2<-pred2$score(msrs(c("classif.acc","classif.auc")))
  c(pp1,pp2)
}
test<-sapply(1:dim(step1)[1],function(x){
  pre(x)
})
test<-t(test)
ff<-data.frame(step1,test)


#Figure 6C
library(ggplot2)

ggplot()+
  geom_histogram(aes(x=step1[which(step1$auc<step1[57,2]),2],y=after_stat(density)),binwidth = 0.01, color='black',fill='#bfefff',linewidth=0.5)+
  geom_density(aes(x=step1[which(step1$auc<step1[57,2]),2],y=after_stat(density)),color='black', fill='#333d3e',alpha=0.25,linewidth=0.5)+
  theme_bw()+
  labs(x="AUC", y='Count of classifier')+
  theme_classic(base_size = 14)+
  theme(axis.text = element_text(color = 'black'))


#Figure 6D
plotdata<-data.frame(type=c("Train","Test","Validation"),auc=as.numeric(ff[57,c(2,7,9)]))

ggplot(plotdata,aes(x=type,y=auc))+
  geom_bar(fill="lightblue1",stat="identity")+
  labs(y="AUC", x="")+
  theme_classic(base_size = 14)+
  theme(text=element_text(angle=45))




library(data.table)
#治疗前治疗后改变的基因
data<-readRDS("Melanoma-PRJEB23709.Response.Rds")
clin<-read.table("Melanoma-PRJEB23709.Response.tsv",sep="\t",header=T)
cl<-readRDS("Melanoma-PRJEB23709.rds")
features_bdp<-data[match(gsub("_","-",colnames(cl)[-400]),data$GENE_SYMBOL),]


library(GSVA)
d1<-as.matrix(data[,-1])
rownames(d1)<-data$GENE_SYMBOL
tt1<-gsva(d1,list(gsub("_","-",colnames(cl)[-400])),method="plage")
boxplot(as.numeric(tt1)~clin$response)

df<-data.frame(score=as.numeric(tt1),group=paste(clin$Treatment,clin$response_NR,sep="_"),sam=clin$patient_name,therapy=clin$Therapy)


library(ggplot2)
library(dplyr)

df$group<-factor(df$group,levels = c("PRE_N","EDT_N","PRE_R","EDT_R"))


#Figure 6E
library(ggpubr)

ggplot(df, aes(x = group, y = score, group = sam)) +
  geom_boxplot(aes(group = group),width=0.5,outlier.size = 0.6,fill="#ee9f90") +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs( y = "Plage Score") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank())+  
  geom_signif(                       
    comparisons=list(c("PRE_N","EDT_N"),
                     c("PRE_R","EDT_R")), 
    step_increase = 0,
    test="wilcox.test",                    
    map_signif_level=F               
  )




#Boruta方法
library(Boruta)
data<-readRDS("Melanoma-PRJEB23709.rds")
set.seed(225)
data<-as.data.frame(data)
data$response<-as.factor(data$response)
test<-Boruta(response~.,data=data,doTrace=0)
import<-names(test$finalDecision[which(test$finalDecision=="Confirmed")])
v22<-read.table("D:\\gencode参考基因组\\id_name_v22.txt",sep="\t")
v22[,1]<-sapply(v22[,1],function(x){unlist(strsplit(x,"[.]"))[1]})
ensg<-v22[match(gsub("_","-",import),v22[,2]),1]

skcm<-read.table("cancer_nets\\TCGA-SKCM.txt",sep="\t",header=T)
need_variable<-skcm[(skcm[,1] %in% ensg) | (skcm[,2] %in% ensg),]



#Figure 6F-left
plot(test)
allat<-test[["ImpHistory"]]
a<-allat[,colnames(allat) %in% getSelectedAttributes(test)]
b<-melt(a)
b$Var2<-factor(b$Var2,levels=colnames(a)[order(apply(a,2,median))])

ggplot(b,aes(y=Var2,x=value))+
  geom_boxplot(fill="#c1deaa",width=0.5,size=0.5,outlier.size = 0.6)+
  theme_bw()+
  theme(axis.title=element_blank())


#Figure 6F-right
fea<-data[,colnames(data) %in% getSelectedAttributes(test)]
mean_fea<-apply(fea,2,function(x){
  x<-scale(x)
  tapply(x,data$response,mean)
})
p<-apply(fea,2,function(x){
  wilcox.test(x~data$response)$p.value
})
mean_fea<-mean_fea[,order(apply(a,2,median))]
mean_fea<-t(mean_fea)
library(pheatmap)

pheatmap(mean_fea,cluster_cols = F,cluster_rows = F)



#Figure 6G
library(ggplot2)
pair<-apply(need_variable,1,function(x){paste0(x[9:10],collapse=";")})
plot_data<-data.frame(need_variable[,3:6],pair=pair)
library(reshape2)
pd<-melt(plot_data)
pd<-pd[pd$value!=0,]
pd$pair<-factor(pd$pair,levels=pair[order(need_variable$stage)])

ggplot(pd,aes(x=variable,y=pair))+
  geom_point(aes(size=abs(value),color=value))+
  theme_bw()+
  scale_color_gradient2(low="blue",high="red",mid="white",midpoint = 0)+
  theme(axis.title = element_blank())



