options(scipen=200)
library(GENOVA)
source("glouble_calculate.R")
library(GENOVA)
library(bedtoolsr)
options(bedtools.path = "/home/zhang/bedtools2/bin/")
library(reshape2)
library(ggvenn)
library(ggsignif)
library(ggthemes)

test_2k<-load_contacts("test_2000.cool", 
                     balancing =T, 
                     sample_name = "H5", 
                     colour = "red")

glouble_test_1<-glouble_define(paste("chr",1:16,sep = ""),test_2k,infer_cut = 1,ig_cutoff = 0,10)
direction_test_1<-glouble_test_1[[2]]
glouble_test_1<-glouble_test_1[[1]]


all_bin_site<-as.data.frame(test_2k[["IDX"]])
all_bin_site<-all_bin_site[-grep("chrM",all_bin_site$V1),]
all_bin_site$V1<-factor(all_bin_site$V1,levels = paste("chr",1:16,sep = ""))
all_bin_site<-all_bin_site[order(all_bin_site$V1),]
all_bin_site$ratio<-direction_test_1$ratio
direction_test_dat<-cbind(all_bin_site[,1:3],direction_test_1$ratio)
write.table(direction_test_dat,"direction_test_dat.bedgraph",sep = "\t",quote = F,col.names = F,row.names = F)
myBegraphToBigwig("direction_test_dat.bedgraph","direction_test_dat.bw",tmp.seqlength = as.vector(table(all_bin_site$V1)),tmp.levels = paste("chr",1:16,sep = ""))
###################
library(RColorBrewer)
library(dplyr)
library(BiocGenerics)
library(GenomeInfoDb)
library(cowplot)
library(patchwork)
library(Biostrings)
p_hic<-pyramid(exp = test_2k,
               chrom = 'chr5',
               start = 10000,
               end=500000,
               colour = scale_fill_gradientn(colours = colorRampPalette(c("white",brewer.pal(9,"YlOrRd")[6:9]))(100),
                                             limits=c(0,5000)),
               crop_y=c(0,150000)
)
p_segment<-dps_segment_plot(glouble_test_1,chromosome="chr5",region=c(10000,500000),tmp.ylimits = c(0,0.1))
p_DSB<- myBigwigTrack(region = as("chr5:10000-500000","GRanges"),
                      bigwig = "R5_1Aligne.bw",
                      smooth = 1000,
                      lognorm = F,
                      type = "coverage",
                      y_label = "RNA",
                      fontsize=10,
                      track.color="#FE0084",
                      tmp.ylimits=c(0,200),
                      max.downsample = 3000,
                      downsample.rate = 0.1,
                      tmp.seed=44)
p_atac<- myBigwigTrack(region = as("chr5:10000-500000","GRanges"),
                       bigwig = "A5-1.last.bw",
                       smooth = 600,
                       lognorm = F,
                       type = "coverage",
                       y_label = "ATAC",
                       fontsize=10,
                       track.color="#4776FE",
                       tmp.ylimits=c(0,50),
                       max.downsample = 3000,
                       downsample.rate = 0.1,
                       tmp.seed=41)

p_ig<- myBigwigTrack(region = as("chr5:10000-500000","GRanges"),
                     bigwig = "direction_test_dat.bw",
                     smooth = 600,
                     lognorm = F,
                     type = "line",
                     y_label = "DPS",
                     fontsize=10,
                     track.color="#00A800",
                     tmp.ylimits=c(-1,1),
                     max.downsample = 3000,
                     downsample.rate = 0.1,
                     tmp.seed=41)
pdf("CID_plot_test_chr5.pdf",width = 8.6,height = 7.4)
p_hic +p_segment+p_ig+ p_DSB+ p_atac+
  plot_layout(ncol = 1,heights = c(1.5,0.1,1,1,1))
dev.off()
p_hic +p_segment+plot_layout(ncol = 1,heights = c(1,0.1))
