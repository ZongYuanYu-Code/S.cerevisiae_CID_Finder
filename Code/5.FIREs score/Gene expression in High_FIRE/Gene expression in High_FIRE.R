rm(list=ls())
library(tidyr)
library(bedtoolsr)
options(bedtools.path = "/home/yzy/anaconda3/envs/bedtools2.30.0/bin/")

gene_site <- read.table("test_H5_gene_site.bed")
window <- read.table("test_combined_2kb_10.txt",header = T)
window <- window[,c(1,3,4,5,6,7,8,9,2)]
med_fire <- median(window$FIRE)
colnames(gene_site) <- c("chr","start","end","gene","expr")
gene_order <- gene_site[order(gene_site$expr,na.last = TRUE, decreasing = T),]
high_gene <- gene_order[gene_order$expr >= median(gene_order$expr),]
##############################################################################
select_PC <- window[window$PC>=0.7,]
select_fire <- select_PC[order(select_PC$FIRE,na.last = T,decreasing = T),]
high_fire <- select_fire[select_fire$FIRE>=median(select_fire$FIRE),]
low_fire <- select_fire[select_fire$FIRE<median(select_fire$FIRE),]
high_fire$type <- "FIRE"
low_fire$type <- "Non-FIRE"
fire_all <- rbind(high_fire,low_fire)
glouble_gene <- bt.intersect(a=fire_all,b=gene_site,wo=T)
glouble_gene$V17 <- paste0(glouble_gene$V14,"_",glouble_gene$V10)
glouble_exp <- aggregate(glouble_gene$V15,by=list(glouble_gene$V17),FUN=mean)
for (i in 1:nrow(glouble_exp)){
  glouble_exp[i,3] <- strsplit(glouble_exp$Group.1,"_")[[i]][1]
  glouble_exp[i,4] <- strsplit(glouble_exp$Group.1,"_")[[i]][2]
}
glouble_exp <- glouble_exp[,-1]
glouble_exp <- glouble_exp[,c(2,3,1)]
colnames(glouble_exp) <- c("gene","type","value")
glouble_exp$value <- log2(glouble_exp$value+1)
p <- ggplot(glouble_exp,aes(type,value,fill=type))+
  stat_boxplot(geom ='errorbar', width = 0.1)+
  geom_boxplot(width = 0.5,outlier.shape = NA)+
  # geom_violin(trim = FALSE)+
  scale_fill_manual(values = c("#DEB887","#A4D3EE"))+
  geom_signif(comparisons = list(c("FIRE","Non-FIRE")),
              step_increase=0.1,map_signif_level = T,test = wilcox.test,color="black")+
  labs(x="",y="log2(TPM+1)",title ="Wild")+
  theme_few()+
  theme(legend.position = "none")
ggsave("Gene expression in High_FIRE.pdf",p,width = 3, height = 4.5,dpi = 300)

