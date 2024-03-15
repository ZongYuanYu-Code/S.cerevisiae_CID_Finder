

boundary_hDEGs<-function(glouble_dat,sample_gene_bed,sample_id){
  boundary_dat<-data.frame(chr=rep(glouble_test_1$V1,2),site=c(glouble_test_1$V2,glouble_test_1$V3-2000),
                           site2=c(glouble_test_1$V2+2000,glouble_test_1$V3))
  boundary_dat<-boundary_dat[order(boundary_dat$chr),]
  boundary_dat<-boundary_dat[!duplicated(boundary_dat),]
  
  #####################The number of highly expressed genes upstream and downstream in boundary region###############
  all_bin_site<-read.table("BY4741_2000_abs.bed")
  all_bin_site<-all_bin_site[-grep("chrM",all_bin_site$V1),]
  colnames(boundary_dat)<-c("V1","V2","V3")
  boundary_dat<-bt.intersect(a=boundary_dat,b=all_bin_site,wo=boundary_dat)
  boundary_dis<-vector()
  for (i in 1:dim(all_bin_site)[1]) {
    dat_chr<-boundary_dat[boundary_dat$V1%in%all_bin_site[i,1],]
    aa<-dat_chr$V2-all_bin_site[i,2]
    if(length(which(aa>=0))>0){
      min_a<-min(aa[aa>=0])
    } else {min_a= Inf}
    if(length(which(aa<=0))>0){
      min_b<-max(aa[aa<=0])
    } else {min_b= -Inf}
    
    if(abs(min_a)>abs(min_b)){
      boundary_dis<-c(boundary_dis,min_b)
    } else {
      boundary_dis<-c(boundary_dis,min_a)
    }
    
  }
  all_bin_site$boundary_dis<-boundary_dis/1000
  all_bin_site$boundary_dis<-floor(all_bin_site$boundary_dis)
  
  gene_site<-read.table(sample_gene_bed)
  heg<-gene_site[gene_site$V5>=summary(gene_site$V5)[3],]
  
  heg<-heg[,c(1,2,3,4,5)]
  boundary_test_1<-boundary_dat
  boundary_test_1_up1<-all_bin_site[all_bin_site$boundary_dis==2,]
  boundary_test_1_up2<-all_bin_site[all_bin_site$boundary_dis==4,]
  boundary_test_1_up3<-all_bin_site[all_bin_site$boundary_dis==6,]
  boundary_test_1_up4<-all_bin_site[all_bin_site$boundary_dis==8,]
  boundary_test_1_up5<-all_bin_site[all_bin_site$boundary_dis==10,]
  boundary_test_1_up6<-all_bin_site[all_bin_site$boundary_dis==12,]
  boundary_test_1_up7<-all_bin_site[all_bin_site$boundary_dis==14,]
  boundary_test_1_up8<-all_bin_site[all_bin_site$boundary_dis==16,]
  boundary_test_1_down1<-all_bin_site[all_bin_site$boundary_dis== -2,]
  boundary_test_1_down2<-all_bin_site[all_bin_site$boundary_dis== -4,]
  boundary_test_1_down3<-all_bin_site[all_bin_site$boundary_dis== -6,]
  boundary_test_1_down4<-all_bin_site[all_bin_site$boundary_dis== -8,]
  boundary_test_1_down5<-all_bin_site[all_bin_site$boundary_dis== -10,]
  boundary_test_1_down6<-all_bin_site[all_bin_site$boundary_dis== -12,]
  boundary_test_1_down7<-all_bin_site[all_bin_site$boundary_dis== -14,]
  boundary_test_1_down8<-all_bin_site[all_bin_site$boundary_dis== -16,]
  bon_gene<-bt.intersect(a=boundary_test_1,b=heg,v=boundary_test_1)
  bon_gene_up1<-bt.intersect(a=boundary_test_1_up1,b=heg,v=boundary_test_1_up1)
  bon_gene_up2<-bt.intersect(a=boundary_test_1_up2,b=heg,v=boundary_test_1_up2)
  bon_gene_up3<-bt.intersect(a=boundary_test_1_up3,b=heg,v=boundary_test_1_up3)
  bon_gene_up4<-bt.intersect(a=boundary_test_1_up4,b=heg,v=boundary_test_1_up4)
  bon_gene_up5<-bt.intersect(a=boundary_test_1_up5,b=heg,v=boundary_test_1_up5)
  bon_gene_up6<-bt.intersect(a=boundary_test_1_up6,b=heg,v=boundary_test_1_up6)
  bon_gene_up7<-bt.intersect(a=boundary_test_1_up7,b=heg,v=boundary_test_1_up7)
  bon_gene_up8<-bt.intersect(a=boundary_test_1_up8,b=heg,v=boundary_test_1_up8)
  bon_gene_down1<-bt.intersect(a=boundary_test_1_down1,b=heg,v=boundary_test_1_down1)
  bon_gene_down2<-bt.intersect(a=boundary_test_1_down2,b=heg,v=boundary_test_1_down2)
  bon_gene_down3<-bt.intersect(a=boundary_test_1_down3,b=heg,v=boundary_test_1_down3)
  bon_gene_down4<-bt.intersect(a=boundary_test_1_down4,b=heg,v=boundary_test_1_down4)
  bon_gene_down5<-bt.intersect(a=boundary_test_1_down5,b=heg,v=boundary_test_1_down5)
  bon_gene_down6<-bt.intersect(a=boundary_test_1_down6,b=heg,v=boundary_test_1_down6)
  bon_gene_down7<-bt.intersect(a=boundary_test_1_down7,b=heg,v=boundary_test_1_down7)
  bon_gene_down8<-bt.intersect(a=boundary_test_1_down8,b=heg,v=boundary_test_1_down8)
  a<-c(dim(bon_gene_down8)[1]/dim(boundary_test_1_down8)[1],
       dim(bon_gene_down7)[1]/dim(boundary_test_1_down7)[1],
       dim(bon_gene_down6)[1]/dim(boundary_test_1_down6)[1],
       dim(bon_gene_down5)[1]/dim(boundary_test_1_down5)[1],
       dim(bon_gene_down4)[1]/dim(boundary_test_1_down4)[1],
       dim(bon_gene_down3)[1]/dim(boundary_test_1_down3)[1],
       dim(bon_gene_down2)[1]/dim(boundary_test_1_down2)[1],
       dim(bon_gene_down1)[1]/dim(boundary_test_1_down1)[1],
       dim(bon_gene)[1]/dim(boundary_test_1)[1],
       dim(bon_gene_up1)[1]/dim(boundary_test_1_up1)[1],
       dim(bon_gene_up2)[1]/dim(boundary_test_1_up2)[1],
       dim(bon_gene_up3)[1]/dim(boundary_test_1_up3)[1],
       dim(bon_gene_up4)[1]/dim(boundary_test_1_up4)[1],
       dim(bon_gene_up5)[1]/dim(boundary_test_1_up5)[1],
       dim(bon_gene_up6)[1]/dim(boundary_test_1_up6)[1],
       dim(bon_gene_up7)[1]/dim(boundary_test_1_up7)[1],
       dim(bon_gene_up8)[1]/dim(boundary_test_1_up8)[1]
       
  )
  
  non_boundary<-all_bin_site[all_bin_site$boundary_dis!=0,]
  non_bon_gene<-bt.intersect(a=non_boundary,b=heg,v=non_boundary)
  fisher_a<-matrix(c(dim(boundary_test_1)[1]-dim(bon_gene)[1],dim(bon_gene)[1],
                     dim(non_boundary)[1]-dim(non_bon_gene)[1],dim(non_bon_gene)[1]),ncol = 2)
  a_p<-fisher.test(fisher_a)[["p.value"]]
  pdf(paste("HEGs in boundary ","of ",sample_id,".pdf",sep = ""),width = 5.16,height = 5.65)
  plot(seq(-16, 16, 2),1-a,type = "o",xlab = "Distance from CID boundary (Kb)",
       ylab = "The number of HEGs",xlim = c(-16, 16), ylim = c(0, 1),pch=16,xaxt = "n",lwd=2,col="blue")
  axis(1,seq(-16, 16, 2), seq(-16, 16, 2))
  text(x = -2,y = 0.9,
       labels = paste("P-value = ",signif(a_p,3),sep = ""),
       cex = 1)
  dev.off()
  ##################Association gene expression #######
  library(ggsignif)
  library(ggthemes)
  dat3<-data.frame(chr=rep(glouble_test_1$V1,2),site=c(glouble_test_1$V2,glouble_test_1$V3-2000),
                   site2=c(glouble_test_1$V2+2000,glouble_test_1$V3))
  dat3$site[dat3$site<0]<-0
  glouble_gene<-bt.intersect(a=dat3,b=gene_site,wo=dat3)
  glouble_exp<-aggregate(glouble_gene$V8,by=list(glouble_gene$V7),FUN=mean)
  non_glouble_gene<-gene_site[!gene_site$V4%in%glouble_gene$V7,]
  dat_plot3<-data.frame(type=c(rep("Boundary",length(log2(unique(glouble_exp$x)+1))),
                               rep("Non-boundary",length(log2(non_glouble_gene$V5+1)))),
                        value=c(log2(unique(glouble_exp$x)+1),log2(non_glouble_gene$V5+1)))
  q1<-ggplot(dat_plot3,aes(type,value,fill=type))+
    geom_violin()+
    stat_boxplot(geom ='errorbar', width = 0.1)+
    geom_boxplot(width = 0.5,outlier.shape = NA)+
    # geom_violin(trim = FALSE)+
    scale_fill_manual(values = c("#2E89BE","#D63C4E"))+
    geom_signif(comparisons = list(c("Boundary","Non-boundary")),
                step_increase=0.1,map_signif_level = T,test = wilcox.test,color="black")+
    labs(x="",y="log2(TPM+1)",title ="Wild")+
    theme_few()+
    theme(legend.position = "none")
  return(q1)
}
