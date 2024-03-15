library(dplyr)
library(pheatmap)
library(ggsci)
library(clusterProfiler)
library(ggplot2)
library(RColorBrewer)
library("org.Sc.sgd.db")
library(reshape2)
library(ggpubr)
library(ggthemes)
#################"AZF1","SFP1","SFL1"
colnames(dat_exp)<-unlist(lapply(strsplit(colnames(dat_exp),split = "_"),function(x){x[1]}))
aa<-melt(as.matrix(dat_exp[c("AZF1","SFP1","SFL1"),]))
aa$Var2<-as.character(aa$Var2)
aa$Var2[aa$Var2=="R1"]<-"T20"
aa$Var2[aa$Var2=="R2"]<-"T58"
aa$Var2[aa$Var2=="R6"]<-"T48"
aa$Var2<-factor(aa$Var2,levels = c("T58","T20","T48"))
p1<-ggplot(aa,aes(Var2,log2(value),fill=Var2))+
  facet_wrap(~Var1)+
  stat_boxplot(geom ='errorbar', width = 0.15,position=position_dodge(0.8))+
  geom_boxplot(width = 0.5,outlier.shape = NA,position=position_dodge(0.8))+
  stat_compare_means(method = "t.test",comparisons = list(c("T20","T58"),c("T20","T48"),c("T58","T48")),aes(label=..p.signif..))+
  labs(x="",y="TPM(log2)")+
  scale_fill_manual(values = c("#79BF87","#F96262","#72A0CC"))+
  theme_few()+
  theme(legend.position = "none")