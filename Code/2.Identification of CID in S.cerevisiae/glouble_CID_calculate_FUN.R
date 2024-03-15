glouble_define<-function(chr_name,contact_matrix,infer_cut,ig_cutoff,bin_num){
  all_chr_glouble<-data.frame()
  all_hic_bin_stream<-data.frame()
  ig_cutoff<-log2(ig_cutoff)
  infer_cut<-infer_cut
  # dis_cutoff<-0.2
  # cv_cutoff<-0.4
  chr4_name<-chr_name
  h5_2k<-contact_matrix
  for (i in chr4_name) {
    chr<-i
    chr4 <- select_subset(h5_2k,chr, 0, 100000000)
    chr4_contact <- chr4$z
    chr4_name<-h5_2k[["IDX"]][h5_2k[["IDX"]]$V1%in%chr,]
    chr4_name$site<-1:dim(chr4_name)[1]
    chr4_name<-as.data.frame(chr4_name)
    colnames(chr4_contact)<-1:dim(chr4_contact)[1]
    row.names(chr4_contact)<-1:dim(chr4_contact)[1]
    chr4_contact<-apply(chr4_contact,1,function(x){x/sum(x+1)})
    fc<-unlist(lapply(1:dim(chr4_contact)[1],function(x){
      up_vector<-rev(chr4_contact[max(1,x-bin_num):max(1,x-1),x])
      up_vector<-up_vector[1:min(length(up_vector),bin_num)]
      down_vector<-chr4_contact[min(dim(chr4_contact)[1],x+1):dim(chr4_contact)[1],x]
      down_vector<-down_vector[1:min(length(down_vector),bin_num)]
      site1<-which(up_vector!=0)
      site2<-which(down_vector!=0)
      union_site<-intersect(site1,site2)
      fc1<-log2(sum(up_vector[union_site]))-log2(sum(down_vector[union_site]))
      return(fc1)
    }))
    fc[fc=="NaN"]<-0
    # log2(circle_up/circle_down)
    # hic_bin_stream<-data.frame(name=colnames(chr4_contact),up_score=circle_up,down_score=circle_down)
    # hic_bin_stream<-hic_bin_stream[!(hic_bin_stream$up_score==0&hic_bin_stream$down_score==0),]
    hic_bin_stream<-data.frame(name=colnames(chr4_contact),ratio=fc)
    # hic_bin_stream$ratio<-log2(hic_bin_stream$up_score/hic_bin_stream$down_score)
    hic_bin_stream$type<-"up"
    hic_bin_stream$type[hic_bin_stream$ratio<0]<-"down"
    hic_bin_stream$a<-1:dim(hic_bin_stream)[1]
    hic_bin_stream$chr<-chr
    all_hic_bin_stream<-rbind(all_hic_bin_stream,hic_bin_stream)
    hic_bin_stream$boundary<-"True"
    hic_bin_stream$boundary[abs(hic_bin_stream$ratio)<log2(infer_cut)]<-"False"
    hic_bin_stream2<-hic_bin_stream[!abs(hic_bin_stream$ratio)<log2(infer_cut),]
    site<-vector()
    for (k in 1:dim(hic_bin_stream2)[1]) {
      if(k<dim(hic_bin_stream2)[1]){
        if(hic_bin_stream2$ratio[k]>0 & hic_bin_stream2$ratio[k+1]<0){
          site<-c(site,k)
        }
      }
    }
    site<-c(site,dim(hic_bin_stream2)[1])
    site2<-c(1,(site[1:(length(site)-1)]+1))
    hic_bin_stream2$glouble<-NA
    for (k in 1:length(site)) {
      hic_bin_stream2$glouble[site[k]:site2[k]]<-paste("glouble",k,sep = "_")
    }
    
    hic_bin_stream$boundary2<-0
    hic_bin_stream$glouble<-NA
    hic_bin_stream3<-hic_bin_stream2[hic_bin_stream2$boundary=="True",]
    for (k in unique(hic_bin_stream3$glouble)) {
      test_dat<-hic_bin_stream3[hic_bin_stream3$glouble%in%k,]
      if(dim(test_dat)[1]>1){
        if(max(test_dat$ratio)>0&min(test_dat$ratio)<0){
          hic_bin_stream$boundary2[hic_bin_stream$name%in%test_dat$name[1]]<-1
          hic_bin_stream$boundary2[hic_bin_stream$name%in%test_dat$name[dim(test_dat)[1]]]<-2
        }
      }
    }
    
    left_bon<-which(hic_bin_stream$boundary2==1)
    right_bon<-which(hic_bin_stream$boundary2==2)
    ####################Calculate the index in each glouble####################
    glouble_site<-vector()
    for (k in 1:length(left_bon)) {
      hic_bin_stream$glouble[left_bon[k]:right_bon[k]]<-paste("glouble",k,sep = "-")
      test_dat2<-chr4_name[chr4_name$site%in%hic_bin_stream$name[left_bon[k]:right_bon[k]],]
      glouble_site<-rbind(glouble_site,c(unique(test_dat2$V1),min(test_dat2$V2),max(test_dat2$V3),paste("glouble",k,sep = "-")))
    }
    glouble_site<-as.data.frame(glouble_site)
    glouble_site$V2<-as.numeric(glouble_site$V2)
    glouble_site$V3<-as.numeric(glouble_site$V3)
    glouble<-unique(hic_bin_stream$glouble[grep("glouble",hic_bin_stream$glouble)])
    dis_hic<-1-cor(chr4_contact)
    dis_hic[is.na(dis_hic)]<-0
    glouble_index<-glouble_calculate(dis_hic,hic_bin_stream,chr4_contact)
    glouble_site<-cbind(glouble_site,glouble_index[match(glouble_site$V4,glouble_index$glouble),2:5])
    
    
    ig_site<-which(glouble_site$ig<ig_cutoff)
    ig_site2<-vector()
    glouble_site_non<-glouble_site[glouble_site$ig<ig_cutoff,]
    if(dim(glouble_site_non)[1]>0){
      for (j in 1:(length(ig_site)-1)) {
        ig_site2<-c(ig_site2,ig_site[j+1]-ig_site[j])
      }
      non_site<-c(which(ig_site2>1),length(ig_site))
      non_site2<-c(1,which(ig_site2>1)+1)
      new_glouble_site<-vector()
      for (j in 1:length(non_site)) {
        test_dat4<-glouble_site_non[non_site2[j]:non_site[j],1:3]
        new_glouble_site<-rbind(new_glouble_site,c(unique(test_dat4$V1),min(test_dat4$V2),max(test_dat4$V3),paste("new_glouble",j,sep = "-")))
      }
      new_glouble_site<-as.data.frame(new_glouble_site)
      new_glouble_site$V2<-as.numeric(new_glouble_site$V2)
      new_glouble_site$V3<-as.numeric(new_glouble_site$V3)
      new_hic_glouble<-bt.intersect(a=new_glouble_site,b=chr4_name,wo=new_glouble_site)
      new_hic_glouble<-new_hic_glouble[new_hic_glouble$V9%in%as.numeric(colnames(chr4_contact)),]
      new_hic_bin<-hic_bin_stream[hic_bin_stream$name%in%new_hic_glouble$V9,1:7]
      new_hic_bin$glouble<-new_hic_glouble$V4[match(new_hic_bin$name,new_hic_glouble$V9)]
      new_glouble_index<-glouble_calculate(dis_hic,new_hic_bin,chr4_contact)
      new_glouble_site<-cbind(new_glouble_site,new_glouble_index[match(new_glouble_site$V4,new_glouble_index$glouble),2:5])
      all_glouble_site<-rbind(glouble_site[glouble_site$ig>=ig_cutoff,],new_glouble_site)
    } else {
      all_glouble_site<-glouble_site
    }
    # all_glouble_site2<-all_glouble_site[all_glouble_site$dis<0.25&all_glouble_site$cv<cv_cutoff,]
    all_chr_glouble<-rbind(all_chr_glouble,all_glouble_site)
    
  }
  return(list(all_chr_glouble,all_hic_bin_stream))
}
##################################################################################
glouble_calculate<-function(dis_hic,hic_bin,chr_contact){
  all_ig<-vector()
  all_dis<-vector()
  all_cv<-vector()
  all_stren<-vector()
  hic_bin_stream<-hic_bin
  dis_hic<-dis_hic
  chr4_contact<-chr_contact
  glouble<-unique(hic_bin_stream$glouble[grep("glouble",hic_bin_stream$glouble)])
  for (j in glouble) {
    test_dat3<-hic_bin_stream[hic_bin_stream$glouble%in%j,]
    test_dat3$ratio[test_dat3$ratio%in%c(Inf)]<-1
    test_dat3$ratio[test_dat3$ratio%in%c(-Inf)]<- -1
    ig<-mean(c(max(test_dat3$ratio[test_dat3$ratio>0]),abs(min(test_dat3$ratio[test_dat3$ratio<0]))))
    all_ig<-c(all_ig,ig)
    dis_glouble<-dis_hic[as.character(test_dat3$name),as.character(test_dat3$name)]
    dis<-as.vector(dis_glouble[lower.tri (dis_glouble, diag = FALSE)])
    dis<-mean(dis[dis!=0])
    all_dis<-c(all_dis,dis)
    contact1<-as.data.frame(melt(chr4_contact[as.character(test_dat3$name),as.character(test_dat3$name)]))
    contact1<-contact1[!contact1$value==0,]
    # contact1[!lower.tri(contact1, diag = TRUE)] <- 0;
    contact1$a<-apply(contact1,1,function(x){max(x[1:2])});
    contact1$b<-apply(contact1,1,function(x){min(x[1:2])});
    contact1$ab<-paste(contact1$a,contact1$b,sep = "_")
    contact1<-contact1[!duplicated(contact1$ab),]
    contact1$len<-contact1$a-contact1$b
    stren<-mean(contact1$value)
    all_stren<-c(all_stren,stren)
    if(dim(contact1)[1]>1){
      contact_mean<-aggregate(contact1$value,by=list(contact1$len),FUN=mean)
      contact1$mean<-contact_mean$x[match(contact1$len,contact_mean$Group.1)]
      contact1$r<-contact1$value/contact1$mean
      cv<-sd(contact1$r)/mean(contact1$r)
      all_cv<-c(all_cv,cv)
    } else {
      cv<- 0
      all_cv<-c(all_cv,cv)
    }
  }
  mat1<-data.frame(glouble=glouble,ig=all_ig,dis=all_dis,cv=all_cv,stren=all_stren)
  return(mat1)
}
