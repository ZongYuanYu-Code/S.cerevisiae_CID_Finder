
library(RIdeogram)
load("Gene_anno_information.rdata")

##################Chromosome location ############
chr_infor<-gene_infor[gene_infor$type%in%"region",1:3]
chr_infor<-chr_infor[!chr_infor$seqnames%in%"MT",]
cen_chr<-gene_infor[gene_infor$type%in%"centromere",]
cen_chr<-cen_chr[cen_chr$width>100,1:4]
chr_infor$censtart<-cen_chr$start[match(chr_infor$seqnames,cen_chr$seqnames)]
chr_infor$censend<-cen_chr$end[match(chr_infor$seqnames,cen_chr$seqnames)]
colnames(chr_infor)<-c("Chr","Start","End","CE_start","CE_end")
chr_infor$Chr<-as.character(chr_infor$Chr)
chr_infor$CE_start<-chr_infor$CE_start-5000
chr_infor$CE_end<-chr_infor$CE_end+5000
write.table(chr_infor,"sacc_karyotype.txt",row.names = F,quote = F,sep = "\t")
#################Gene density###########
gene_density <- GFFex(input = "GCF_000146045.2_R64_genomic.gff", karyotype = "sacc_karyotype.txt", feature = "gene", window = 10000)
###################Genes are grouped according to expression#########
load("RNA_exp_TPM.rdata")
load("Gene_anno_information.rdata")
row.names(tpm)<-tpm$V2
tpm<-tpm[,-1]
zero_site<-apply(tpm,1,function(x){sum(x>0)})
dat_tpm<-tpm[zero_site>1,]
table(row.names(dat_tpm)%in%gene_infor$ID)
anno_gene<-gene_infor[match(row.names(dat_tpm),gene_infor$ID),]
row.names(dat_tpm)<-anno_gene$Name
gene_site<-anno_gene[,c(1:3,8,10)]
gene_site<-gene_site[!gene_site$seqnames%in%"MT",]
H5_exp<-dat_tpm[,13:15]
zero_site2<-apply(H5_exp,1,function(x){sum(x>0)})
gene_site2<-gene_site[gene_site$Name%in%row.names(H5_exp)[zero_site2>1],]
gene_exp<-H5_exp[gene_site2$Name,]
table(row.names(gene_exp)==gene_site2$Name)
gene_exp<-apply(gene_exp,1,function(x){max(x)})
gene_site2$exp<-gene_exp
gene_site2<-gene_site2[order(gene_site2$exp,decreasing = T),]
gene_500<-gene_site2[c(1:250,(dim(gene_site2)[1]-249):dim(gene_site2)[1]),c(1:3,5)]
colnames(gene_500)<-c("Chr","Start","End","Type")
gene_500$Chr<-as.character(gene_500$Chr)
gene_500$Type<-c(rep("High-expression",250),rep("Low-expression",250))
gene_500$Shape<-"triangle"
gene_500$color<-"D78697"
gene_500$color[gene_500$Type%in%"Low-expression"]<-"ABDDA5"

gene_500<-gene_500[,c(4,5,1:3,6)]

exp_disturb<-gene_site2[,c(1:3,6)]
exp_disturb$seqnames<-as.character(exp_disturb$seqnames)
colnames(exp_disturb)<-colnames(gene_density)
exp_disturb$Value<-log2(exp_disturb$Value+1)
exp_disturb2<-exp_disturb
exp_disturb2<-site_fun(exp_disturb2,chr_infor,window = 10000)
exp_disturb2$Color<-"009292"
###################Mapping the distribution of chromosomes#########
ideogram(karyotype = chr_infor,label = gene_500, label_type = "marker")
ideogram(karyotype = chr_infor,label = gene_500, label_type = "marker",overlaid = gene_density,colorset1 = c("#2E89BC","#65C3A6","#AADEA3", "#FDAF61","#EE8175", "#E00C29"))######
convertSVG("chromosome.svg", device = "pdf")
ideogram(karyotype = chr_infor, label = exp_disturb2, label_type = "line")
convertSVG("chromosome.svg", device = "png")
ideogram(karyotype = chr_infor,overlaid = exp_disturb,label = exp_disturb2, label_type = "line",colorset1 = brewer.pal(9,'YlOrBr'))
convertSVG("chromosome.svg", device = "pdf")

#####################################
site_fun<-function (input, karyotype, window = 1e+06) 
{
  gff <- input
  karyotype <- karyotype
  gff <- subset(gff, gff$Chr %in% karyotype$Chr)
  list_chr <- vector("list", length(names(table(gff$Chr))))
  names(list_chr) <- names(table(gff$Chr))
  for (i in 1:(length(list_chr))) {
    list_chr[[i]] <- as.data.frame(cut(subset(gff,gff$Chr == names(list_chr[i]))$Start, 
                                       breaks = c(seq(0,subset(karyotype, karyotype$Chr == names(list_chr[i]))[1,3], window), subset(karyotype, karyotype$Chr == names(list_chr[i]))[1, 3]))
    )
    list_chr[[i]] <- tidyr::separate(list_chr[[i]], 1, into = c("Start", 
                                                                "End"), sep = ",")
    list_chr[[i]]$Start <- gsub("\\(", "", list_chr[[i]]$Start)
    list_chr[[i]]$End <- gsub("\\]", "", list_chr[[i]]$End)
    list_chr[[i]]$Start <- as.numeric(list_chr[[i]]$Start)
    list_chr[[i]]$End <- as.numeric(list_chr[[i]]$End)
    list_chr[[i]]$Start <- list_chr[[i]]$Start + 1
    list_chr[[i]]$Chr <- names(list_chr[i])
    list_chr[[i]] <- cbind(list_chr[[i]][, 3], list_chr[[i]][,1:2])
    list_chr[[i]]$Value<-subset(gff,gff$Chr == names(list_chr[i]))$Value
    colnames(list_chr[[i]]) <- c("Chr", "Start", "End","Value")
    list_chr[[i]]$Chr <- as.character(list_chr[[i]]$Chr)
    list_cop<-list_chr[[i]]
    list_cop$chr_name<-paste(list_cop$Chr,list_cop$Start,list_cop$End,sep = "_")
    list_cop2<-aggregate(list_cop[,4],by=list(list_cop$chr_name),FUN = max)
    list_cop3<-list_cop[match(list_cop2$Group.1,list_cop$chr_name),1:3]
    list_cop3$Value<-list_cop2$x
    list_chr[[i]]<-list_cop3
    list_chr[[i]]<-list_chr[[i]][order(list_chr[[i]]$Start),]
    
  }
  l <- data.frame()
  for (i in 1:(length(list_chr))) {
    df.now <- list_chr[[i]]
    l <- rbind(l, df.now)
  }
  data.frame(l)
}