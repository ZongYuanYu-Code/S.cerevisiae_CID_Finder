library(Matrix)
library(rhdf5)
library(ggplot2)
library(ggsci)
library(ggsignif)
library(RColorBrewer)
library(ggthemes)
library(bedtoolsr)
options(bedtools.path = "/home/zhang/bedtools2/bin/")
file1<-list.files(pattern = "_10000.h5")
load("Site_TEL_CEN.rdata")
NaN_replay<-function(bb){
  nan_replay=unlist(lapply(which(bb=="NaN"),
                           function(x){site2=names(sort(abs(which(bb!="NaN")-x)))[1:2];
                           return(mean(as.vector(bb[site2])))}))
  bb[names(nan_replay)]=as.vector(nan_replay)
  return(bb)
}
title_name<-c("T20","T58","J05","J27","Wild","T48")
all_p<-list()
for (i in 1:length(file1)) {
  scaleFUN <- function(x){sprintf("%.3f", x)}
  dat_h5<-h5read(file1[i],"matrix")
  dat2_h5<-as.data.frame(h5read(file1[i],"intervals"))
  dat2_h5$name<-paste(dat2_h5$chr_list,dat2_h5$start_list,dat2_h5$end_list,sep = "_")
  dat2_h5<-dat2_h5[,c(1,4,2,3,5)]
  counts <- dat_h5$data
  col_name <- dat2_h5$name
  indices <- dat_h5$indices
  indptr <- dat_h5$indptr
  shape <- dat_h5$shape
  mat <- sparseMatrix(i = as.vector(indices), p = as.vector(indptr),x = as.numeric(counts[]), dims =shape[],dimnames=list(as.vector(col_name),as.vector(col_name)), index1=FALSE,symmetric=T)
  mat<-as.matrix(mat)
  mat16<-mat[-grep("chrM",colnames(mat)),-grep("chrM",colnames(mat))]
  for (j in paste("chr",1:16,"_",sep = "")) {
    mat16[grep(j,colnames(mat)),grep(j,colnames(mat))]<-0
  }
  cen_bin<-bt.intersect(a=cen_chr,b=dat2_h5,wo=cen_chr)
  tel_bin<-bt.intersect(a=tel_chr,b=dat2_h5,wo=tel_chr)
  aa<-colSums(mat16)
  cen_infor<-vector()
  for (k in cen_bin$V7) {
    aa2<-aa
    if(aa[cen_bin$V12[cen_bin$V7%in%k]]==0){
      aa2[cen_bin$V12[cen_bin$V7%in%k]]<-1
      aa2<-aa2[!aa2==0]
      site<-which(names(aa2)%in%cen_bin$V12[cen_bin$V7%in%k])
      site_infor<-names(aa2)[c((site-2):(site-1),(site+1):(site+3))]
    } else {
      aa2<-aa2[!aa2==0]
      site<-which(names(aa2)%in%cen_bin$V12[cen_bin$V7%in%k])
      site_infor<-names(aa2)[(site-2):(site+2)]
    }
    cen_infor<-c(cen_infor,site_infor)
  }
  mat16_cen<-mat16[cen_infor,]
  bb<-apply(mat16_cen,2,function(x){x=x[x!=0];return(mean(x))})
  bb<-NaN_replay(bb)
  
  dat_plot<-data.frame(name=names(bb),value=as.vector(bb))
  dat_plot_bin<-data.frame(name=cen_bin$V7,min=which(dat_plot$name%in%cen_bin$V12)-2,
                           max=which(dat_plot$name%in%cen_bin$V12)+2)
  dat_plot$type<-dat2_h5$chr_list[-grep("chrM",colnames(mat))]
  dat_plot$x<-1:dim(dat_plot)[1]
  p1<-ggplot(dat_plot, aes(x=x, y=value)) +   
    geom_rect(data=dat_plot_bin,inherit.aes=FALSE,aes(xmin=min, xmax=max, ymin=-Inf, ymax=Inf,group=name,fill='#F9D6E4'),alpha = .5)+
    # geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1) +
    geom_line(color="black",size=.5) +
    # geom_point()+
    theme_few()+
    xlab("The length of Genome (10Kb for one bin)")+ylab("Contact score(cen)")+
    scale_x_discrete(expand = c(0, 0))+
    theme(legend.position = "none")+
    # ylim(0,100)+
    ggtitle(paste("Inter-chromosomal interaction of ",title_name[i],sep = ""))+
    scale_y_continuous(labels = scaleFUN)+
    geom_vline(xintercept=unlist(lapply(1:16,function(x){sum(table(dat_plot$type)[paste("chr",1:16,sep = "")][1:x])})),size=0.1)
  
  tel_left_infor<-vector()
  left_site<-vector()
  for (k in unique(tel_bin$V7[grep("L",tel_bin$V7)])) {
    aa2<-aa
    left_site<-c(left_site,tel_bin$V12[tel_bin$V7%in%k][1])
    if(sum(aa[tel_bin$V12[tel_bin$V7%in%k]])==0){
      aa2[tel_bin$V12[tel_bin$V7%in%k]]<-1
      aa2<-aa2[!aa2==0]
      site<-which(names(aa2)%in%tel_bin$V12[tel_bin$V7%in%k])[1]
      site_infor<-names(aa2)[(site+1):(site+5)]
    } else {
      aa2<-aa2[!aa2==0]
      site<-which(names(aa2)%in%tel_bin$V12[tel_bin$V7%in%k])[1]
      site_infor<-names(aa2)[site:(site+4)]
    }
    tel_left_infor<-c(tel_left_infor,site_infor)
  }
  
  tel_right_infor<-vector()
  right_site<-vector()
  for (k in unique(tel_bin$V7[grep("R",tel_bin$V7)])) {
    aa2<-aa
    right_site<-c(right_site,rev(tel_bin$V12[tel_bin$V7%in%k])[1])
    if(sum(aa[tel_bin$V12[tel_bin$V7%in%k]])==0){
      aa2[tel_bin$V12[tel_bin$V7%in%k]]<-1
      aa2<-aa2[!aa2==0]
      site<-which(names(aa2)%in%tel_bin$V12[tel_bin$V7%in%k])[1]
      site_infor<-names(aa2)[(site-5):(site-1)]
    } else {
      aa2<-aa2[!aa2==0]
      site<-which(names(aa2)%in%tel_bin$V12[tel_bin$V7%in%k])[1]
      site_infor<-names(aa2)[(site-4):site]
    }
    tel_right_infor<-c(tel_right_infor,site_infor)
  }
  
  
  mat16_tel<-mat16[c(tel_left_infor,tel_left_infor),]
  bb<-apply(mat16_tel,2,function(x){x=x[x!=0];return(mean(x))})
  bb<-NaN_replay(bb)
  
  dat_plot2<-data.frame(name=names(bb),value=as.vector(bb))
  dat_plot_bin2<-data.frame(name=c(unique(tel_bin$V7[grep("L",tel_bin$V7)]),unique(tel_bin$V7[grep("R",tel_bin$V7)])),
                            min=c(which(dat_plot2$name%in%left_site)-1,which(dat_plot2$name%in%right_site)-4),
                           max=c(which(dat_plot2$name%in%left_site)+4,which(dat_plot2$name%in%right_site)+1))
  dat_plot2$type<-dat2_h5$chr_list[-grep("chrM",colnames(mat))]
  dat_plot2$x<-1:dim(dat_plot2)[1]
  p2<-ggplot(dat_plot2, aes(x=x, y=value,label=round(value, digits = 3))) +   
    geom_rect(data=dat_plot_bin2,inherit.aes=FALSE,aes(xmin=min, xmax=max, ymin=-Inf, ymax=Inf,group=name,fill='#F9D6E4'),alpha = .5)+
    # geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1) +
    geom_line(color="black",size=.5) +
    # geom_point()+
    theme_few()+
    xlab("The length of Genome (1Kb for one bin)")+ylab("Contact score(tel)")+
    scale_x_discrete(expand = c(0, 0))+
    theme(legend.position = "none")+
    scale_y_continuous(labels = scaleFUN)+
    geom_vline(xintercept=unlist(lapply(1:16,function(x){sum(table(dat_plot$type)[paste("chr",1:16,sep = "")][1:x])})),size=0.1)
  all_p<-c(all_p,list(p1),list(p2))

}
library(cowplot)
plot_grid(p1,p2,ncol = 1)