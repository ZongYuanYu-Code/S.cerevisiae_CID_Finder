library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(RColorBrewer)
source("utils.R")
kegmt<-read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")

geneList<-test_diff_gene$logFC
names(geneList)<-row.names(test_diff_gene)
geneList<-geneList[abs(geneList)>0.58]
geneList<-sort(geneList,decreasing = T)

kegg_go<-GSEA(geneList,TERM2GENE = kegmt,pvalueCutoff = 1)

kegg_kegg<-GSEA(geneList,TERM2GENE = kegmt,pvalueCutoff = 0.05)
mycol<-brewer.pal(n = 11, "RdYlGn")
gseaplot2(kegg_go,c(1,2),
          color=c("#25A9E0","#BE1E2D"),pvalue_table = T,
          ES_geom="line")
