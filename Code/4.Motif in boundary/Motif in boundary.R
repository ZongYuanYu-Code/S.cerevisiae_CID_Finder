

heg<-read.table("test_motif.tsv",sep = "\t",header = T)
boundary_dat<-data.frame(chr=rep(glouble_h5_1$V1,2),site=c(glouble_h5_1$V2,glouble_h5_1$V3-2000),
                         site2=c(glouble_h5_1$V2+2000,glouble_h5_1$V3))
boundary_dat<-boundary_dat[order(boundary_dat$chr),]
boundary_dat<-boundary_dat[!duplicated(boundary_dat),]

all_bin_site<-as.data.frame(h5_2k[["IDX"]])
all_bin_site<-all_bin_site[-grep("chrM",all_bin_site$V1),]
all_bin_site$V1<-factor(all_bin_site$V1,levels = paste("chr",1:16,sep = ""))
all_bin_site<-all_bin_site[order(all_bin_site$V1),]
all_bin_site$ratio<-direction_h5_1$ratio

colnames(boundary_dat)<-c("V1","V2","V3")
boundary_dat<-bt.intersect(a=boundary_dat,b=all_bin_site,wo=boundary_dat)
boundary_dat2<-boundary_dat[abs(boundary_dat$V8)> -1,1:3]
boundary_dis<-vector()
for (i in 1:dim(all_bin_site)[1]) {
  dat_chr<-boundary_dat2[boundary_dat2$V1%in%all_bin_site[i,1],]
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

boundary_h5_1<-boundary_dat[,1:3]
boundary_h5_1_up1<-all_bin_site[all_bin_site$boundary_dis==2,]
boundary_h5_1_up2<-all_bin_site[all_bin_site$boundary_dis==4,]
boundary_h5_1_up3<-all_bin_site[all_bin_site$boundary_dis==6,]
boundary_h5_1_up4<-all_bin_site[all_bin_site$boundary_dis==8,]
boundary_h5_1_up5<-all_bin_site[all_bin_site$boundary_dis==10,]
boundary_h5_1_up6<-all_bin_site[all_bin_site$boundary_dis==12,]
boundary_h5_1_up7<-all_bin_site[all_bin_site$boundary_dis==14,]
boundary_h5_1_up8<-all_bin_site[all_bin_site$boundary_dis==16,]
boundary_h5_1_down1<-all_bin_site[all_bin_site$boundary_dis== -2,]
boundary_h5_1_down2<-all_bin_site[all_bin_site$boundary_dis== -4,]
boundary_h5_1_down3<-all_bin_site[all_bin_site$boundary_dis== -6,]
boundary_h5_1_down4<-all_bin_site[all_bin_site$boundary_dis== -8,]
boundary_h5_1_down5<-all_bin_site[all_bin_site$boundary_dis== -10,]
boundary_h5_1_down6<-all_bin_site[all_bin_site$boundary_dis== -12,]
boundary_h5_1_down7<-all_bin_site[all_bin_site$boundary_dis== -14,]
boundary_h5_1_down8<-all_bin_site[all_bin_site$boundary_dis== -16,]
bon_gene<-bt.intersect(a=boundary_h5_1,b=heg,v=boundary_h5_1)
bon_gene_up1<-bt.intersect(a=boundary_h5_1_up1,b=heg,v=boundary_h5_1_up1)
bon_gene_up2<-bt.intersect(a=boundary_h5_1_up2,b=heg,v=boundary_h5_1_up2)
bon_gene_up3<-bt.intersect(a=boundary_h5_1_up3,b=heg,v=boundary_h5_1_up3)
bon_gene_up4<-bt.intersect(a=boundary_h5_1_up4,b=heg,v=boundary_h5_1_up4)
bon_gene_up5<-bt.intersect(a=boundary_h5_1_up5,b=heg,v=boundary_h5_1_up5)
bon_gene_up6<-bt.intersect(a=boundary_h5_1_up6,b=heg,v=boundary_h5_1_up6)
bon_gene_up7<-bt.intersect(a=boundary_h5_1_up7,b=heg,v=boundary_h5_1_up7)
bon_gene_up8<-bt.intersect(a=boundary_h5_1_up8,b=heg,v=boundary_h5_1_up8)
bon_gene_down1<-bt.intersect(a=boundary_h5_1_down1,b=heg,v=boundary_h5_1_down1)
bon_gene_down2<-bt.intersect(a=boundary_h5_1_down2,b=heg,v=boundary_h5_1_down2)
bon_gene_down3<-bt.intersect(a=boundary_h5_1_down3,b=heg,v=boundary_h5_1_down3)
bon_gene_down4<-bt.intersect(a=boundary_h5_1_down4,b=heg,v=boundary_h5_1_down4)
bon_gene_down5<-bt.intersect(a=boundary_h5_1_down5,b=heg,v=boundary_h5_1_down5)
bon_gene_down6<-bt.intersect(a=boundary_h5_1_down6,b=heg,v=boundary_h5_1_down6)
bon_gene_down7<-bt.intersect(a=boundary_h5_1_down7,b=heg,v=boundary_h5_1_down7)
bon_gene_down8<-bt.intersect(a=boundary_h5_1_down8,b=heg,v=boundary_h5_1_down8)
a<-c(dim(bon_gene_down8)[1]/dim(boundary_h5_1_down8)[1],
     dim(bon_gene_down7)[1]/dim(boundary_h5_1_down7)[1],
     dim(bon_gene_down6)[1]/dim(boundary_h5_1_down6)[1],
     dim(bon_gene_down5)[1]/dim(boundary_h5_1_down5)[1],
     dim(bon_gene_down4)[1]/dim(boundary_h5_1_down4)[1],
     dim(bon_gene_down3)[1]/dim(boundary_h5_1_down3)[1],
     dim(bon_gene_down2)[1]/dim(boundary_h5_1_down2)[1],
     dim(bon_gene_down1)[1]/dim(boundary_h5_1_down1)[1],
     dim(bon_gene)[1]/dim(boundary_h5_1)[1],
     dim(bon_gene_up1)[1]/dim(boundary_h5_1_up1)[1],
     dim(bon_gene_up2)[1]/dim(boundary_h5_1_up2)[1],
     dim(bon_gene_up3)[1]/dim(boundary_h5_1_up3)[1],
     dim(bon_gene_up4)[1]/dim(boundary_h5_1_up4)[1],
     dim(bon_gene_up5)[1]/dim(boundary_h5_1_up5)[1],
     dim(bon_gene_up6)[1]/dim(boundary_h5_1_up6)[1],
     dim(bon_gene_up7)[1]/dim(boundary_h5_1_up7)[1],
     dim(bon_gene_up8)[1]/dim(boundary_h5_1_up8)[1]
     
)
non_boundary<-all_bin_site[all_bin_site$boundary_dis!=0,]
non_bon_gene<-bt.intersect(a=non_boundary,b=heg,v=non_boundary)
fisher_a<-matrix(c(dim(boundary_h5_1)[1]-dim(bon_gene)[1],dim(bon_gene)[1],
                   dim(non_boundary)[1]-dim(non_bon_gene)[1],dim(non_bon_gene)[1]),ncol = 2)
a_p<-fisher.test(fisher_a)[["p.value"]]

pdf("motif_number_plot.pdf",width = 5,height = 4.6)
plot(seq(-16, 16, 2),1-a,type = "o",xlab = "Distance from CID boundary (Kb)",
     ylab = "The number of motif",xlim = c(-16, 16), ylim = c(0, 0.8),pch=16,xaxt = "n",lwd=2,col="blue")
axis(1,seq(-16, 16, 2), seq(-16, 16, 2))
text(x = -2,y = 0.4,
     labels = paste("P-value = ",signif(a_p,3),sep = ""),
     cex = 1)
dev.off()
