


load("Gene_anno_information.rdata")
cen_chr<-gene_infor[gene_infor$type%in%"centromere",]
cen_chr<-cen_chr[cen_chr$width==10,1:6]
cen_chr$name<-paste("cen",1:16,sep = "")

tel_chr<-gene_infor[gene_infor$type%in%"telomere",]
tel_chr<-tel_chr[,1:6]
tel_chr$name<-"aa"
tel_chr$name[tel_chr$start==1]<-paste("tel",1:16,"L",sep = "")
tel_chr$name[tel_chr$start!=1]<-paste("tel",1:16,"R",sep = "")
save(cen_chr,tel_chr,file = "Site_TEL_CEN.rdata")

