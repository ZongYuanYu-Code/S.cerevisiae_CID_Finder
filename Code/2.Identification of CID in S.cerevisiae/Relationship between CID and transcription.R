#########################################
####################################


library(ConsensusClusterPlus)
library(ggsignif)
library(ggplot2)
library(ggthemes)
library(cowplot)
gene_site<-read.table("test_gene_site.bed")
gene_site$site1<-apply(gene_site,1,function(x){round(median(x[2]:x[3]))})
gene_site$site2<-gene_site$site1+1
gene_site2<-gene_site[,c(1,6,7,4,5)]
dat_plot1<-glouble_test_1
dat_plot1<-dat_plot1[!dat_plot1$cv==0,]
dat_plot1$len<-(dat_plot1$V3-dat_plot1$V2)/1000
summary(dat_plot1$len)
hist(dat_plot1$len,breaks =50,xlab = "Length from CID (Kb)",col="blue",border="black")

dat_plot1$glouble_name<-paste(dat_plot1$V1,dat_plot1$V4,sep = "_")

glouble_gene<-bt.intersect(a=dat_plot1,b=gene_site,wo=dat_plot1)
glouble_exp<-aggregate(glouble_gene$V15,by=list(glouble_gene$V10),FUN=max)
dat_plot1$exp<-glouble_exp$x[match(dat_plot1$glouble_name,glouble_exp$Group.1)]
dat_plot1<-dat_plot1[!is.na(dat_plot1$exp),]

dat_plot1$V1<-factor(dat_plot1$V1,levels = paste("chr",1:16,sep = ""))
q1<-ggplot(dat_plot1,aes(V1,ig,fill=V1))+
  stat_boxplot(geom ='errorbar', width = 0.2)+
  geom_boxplot(width = 0.5)+
  labs(x="",y="Contact preference score")+
  theme_bw()+
  theme(legend.position = "none")+
  ylim (0, 1)

chr_glouble<-as.data.frame(table(dat_plot1$V1))
q2<-ggplot(chr_glouble, mapping = aes(x = factor(Var1), y = Freq)) + 
  geom_bar(stat = 'identity', position = 'stack',fill="#65C3A4",width = 0.6)+
  # scale_fill_manual(values = c("#ABDDA5","#65C3A4"))+
  theme_bw()+
  labs(x="",y="The number of CIDs")
plot_grid(q1,q2,ncol = 1)


feature_matrix<-as.matrix(dat_plot1[,c(5:8)])
row.names(feature_matrix)<-dat_plot1$glouble_name
feature_matrix<-apply(feature_matrix,2,function(x){log2(x)})
d = sweep(t(feature_matrix),1, apply(t(feature_matrix),1,median,na.rm=T))

pheatmap::pheatmap(dd,show_rownames = T,show_colnames = F,cluster_cols = T,
                   cluster_rows = F,clustering_distance_cols="euclidean",clustering_method="ward.D")
aa<-ConsensusClusterPlus(d, maxK = 10, reps=100, pItem=0.7, 
                         pFeature=1, clusterAlg="km",
                         title="hic",
                         distance="euclidean",verbose=F,plot = "pdf")
cluster2<-aa[[3]][["consensusClass"]]
dat_plot1$cluster<-as.character(cluster2[match(dat_plot1$glouble_name,names(cluster2))])

p1<-ggplot(dat_plot1,aes(cluster,log2(exp+1),fill=cluster))+
  stat_boxplot(geom ='errorbar', width = 0.2)+
  geom_boxplot(width = 0.5,outlier.shape = NA)+
  scale_fill_manual(values =c('#72D0BB','#1DB7C5','#0290BE','#075FA6'))+
  geom_signif(comparisons = list(c("1","2"),c("1","3"),c("2","3")),
              step_increase=0.1,map_signif_level = T,test = wilcox.test,color="black")+
  labs(x="",y="log2(TPM+1)")+
  theme_few()
p2<-ggplot(dat_plot1,aes(cluster,dis,fill=cluster))+
  stat_boxplot(geom ='errorbar', width = 0.2)+
  geom_boxplot(width = 0.5,outlier.shape = NA)+
  scale_fill_manual(values =c('#72D0BB','#1DB7C5','#0290BE','#075FA6'))+
  geom_signif(comparisons = list(c("1","2"),c("1","3"),c("2","3")),
              step_increase=0.1,map_signif_level = T,test = t.test,color="black")+
  labs(x="",y="Distance")+
  theme_few()
p3<-ggplot(dat_plot1,aes(cluster,ig,fill=cluster))+
  stat_boxplot(geom ='errorbar', width = 0.2)+
  geom_boxplot(width = 0.5)+
  scale_fill_manual(values =c('#72D0BB','#1DB7C5','#0290BE','#075FA6'))+
  geom_signif(comparisons = list(c("1","2"),c("1","3"),c("2","3")),
              step_increase=0.1,map_signif_level = T,test = t.test,color="black")+
  labs(x="",y="Contact preference score")+
  theme_few()
p4<-ggplot(dat_plot1,aes(cluster,cv,fill=cluster))+
  stat_boxplot(geom ='errorbar', width = 0.2)+
  geom_boxplot(width = 0.5,outlier.shape = NA)+
  scale_fill_manual(values =c('#72D0BB','#1DB7C5','#0290BE','#075FA6'))+
  geom_signif(comparisons = list(c("1","2"),c("1","3"),c("2","3")),
              step_increase=0.1,map_signif_level = T,test = t.test,color="black")+
  labs(x="",y="Coefficient of variation")+
  theme_few()
p5<-ggplot(dat_plot1,aes(cluster,stren,fill=cluster))+
  stat_boxplot(geom ='errorbar', width = 0.2)+
  geom_boxplot(width = 0.5,outlier.shape = NA)+
  scale_fill_manual(values =c('#72D0BB','#1DB7C5','#0290BE','#075FA6'))+
  geom_signif(comparisons = list(c("1","2"),c("1","3"),c("2","3")),
              step_increase=0.1,map_signif_level = T,test = t.test,color="black")+
  labs(x="",y="Contact score")+
  theme_few()
plot_grid(p1,p2,p3,p4,p5,ncol = 3)


aa<-cor.test(log2(dat_plot1[,5]),log2(dat_plot1$exp+1),method ="spearman")
plot(log2(dat_plot1[,5]),log2(dat_plot1$exp+1),xlab = "Contact preference score (log)",
     ylab = "log2(TPM+1)",pch=16,lwd=2,col="#185ABD")
text(x = -4,y = 14,
     labels = paste("rho = ",signif(aa[["estimate"]][["rho"]],3),";","P-value = ",signif(aa[["p.value"]],3),sep = ""),
     cex = 1)
text(x = -6,y = 12,
     labels = paste("P-value = ",signif(aa[["p.value"]],3),sep = ""),
     cex = 1)

##################Construct regression model##############
library(nnet)
library(pROC)
logit_dat<-dat_plot1[,c(5:8)]
mult.model<-multinom(cluster~.,data=logit_dat)
pre_logistic<-predict(mult.model,newdata = logit_dat)
dat_plot1$cluster<-pre_logistic
table(logit_dat$cluster,pre_logistic)

test_dat=dat_plot1[,c(5:8,15)]
log1<-glm(cluster2~.,family=binomial(link='logit'),data=logit_dat)
summary(log1)
log2<-step(log1)
summary(log2)
pred<-predict(log2,newdata = logit_dat)
library(ROCR)
pred2<-prediction(pred,logit_dat$cluster2)
performance(pred2,'auc')@y.values
perf<-performance(pred2,'tpr','fpr')
plot(perf,col=2,type="l",lwd=2)
f=function(x){ y=x;return(y)}
curve(f(x),0,1,col=4,lwd=2,lty=2,add=T)
roc1<-roc(logit_dat$cluster2,pred,plot=TRUE, print.thres=TRUE, print.auc=TRUE,levels = c(0,1),direction = "<")
