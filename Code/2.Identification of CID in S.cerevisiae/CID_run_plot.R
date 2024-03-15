
source("glouble_calculate.R")
library(GENOVA)
library(bedtoolsr)
options(bedtools.path = "/home/zhang/bedtools2/bin/")
library(reshape2)
library(ggvenn)
test_2k_1<- load_contacts(signal_path = "BY4741_2000_iced.matrix",
                      indices_path = "BY4741_2000_abs.bed",
                      sample_name = "Wild_1",
                      colour = "black")

glouble_test_1<-glouble_define(paste("chr",1:16,sep = ""),test_2k_1,infer_cut = 1,ig_cutoff = 0,10)
# glouble_test_1<-glouble_test_1[-grep("new_glouble",glouble_test_1$V4),]
direction_test_1<-glouble_test_1[[2]]
glouble_test_1<-glouble_test_1[[1]]




glouble_site2<-glouble_test_1
hic_matrixplot(exp1 = test_2k_1,
               # exp2 = Hap1_SCC4_10kb,
               chrom = 'chr5',
               start = 460000,
               end=480000,
               chip = "R5_1Aligned.bw",
               symmAnn = F,
               colour_bar=T,
               tads = list(glouble_site2[grep("new_glouble",glouble_site2$V4),],
                           glouble_site2[grep("^glouble",glouble_site2$V4),]), # see ATA
               tads.type = 'lower', # only plot in lower triangle
               tads.colour = c('red',"blue") # green TAD-borders
               # colour = scale_fill_gradientn(colours = rev(hcl.colors(100,palette = "YlGnBu")))
               # colour_lim = c(0,9000)
)
ATA_Hap1_WTcalls <- ATA(test_2k_1,
                        bed = glouble_site2,dist_thres=c(5000,60000))
visualise(ATA_Hap1_WTcalls)


