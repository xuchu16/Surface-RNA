library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
args <- commandArgs()
filename = args[6]
f=read.csv(filename,sep = '\t',header=T)
gene = f[,1]
print(gene)
gene.df <- bitr(gene,
	fromType = "SYMBOL",
	toType = c("ENTREZID","SYMBOL"),
	OrgDb = org.Mm.eg.db)
print(gene.df)

GO_BP = enrichGO(gene = gene,
	OrgDb = org.Mm.eg.db,
	keyType = "SYMBOL",
	ont = "BP",
#	ont = "ALL",
	pAdjustMethod = "BH",
	pvalueCutoff = 0.05)
print(GO_BP)
dotplot(GO_BP,showCategory=10)
ggsave(paste(filename,'GOBP.pdf'),width=16, height=20)
write.table(GO_BP,file = paste(filename,"GOBP.xls"),sep = '\t')

GO_MF = enrichGO(gene = gene,
        OrgDb = org.Mm.eg.db,
        keyType = "SYMBOL",
        ont = "MF",
#       ont = "ALL",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05)
dotplot(GO_MF,showCategory=10)
ggsave(paste(filename,'GOMF.pdf'),width=16, height=20)
write.table(GO_MF,file = paste(filename,"GOMF.xls"),sep = '\t')

GO_CC = enrichGO(gene = gene,
        OrgDb = org.Mm.eg.db,
        keyType = "SYMBOL",
#        ont = "CC",
        ont = "ALL",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05)
dotplot(GO_CC,showCategory=10)
ggsave(paste(filename,'GOALL.pdf'),width=16, height=20)
write.table(GO_CC,file = paste(filename,"GOALL.xls"),sep = '\t')

kegg <- enrichKEGG(gene = gene.df$ENTREZID,
                    organism  = 'mmu',
                 pvalueCutoff  = 0.05)
dotplot(kegg,showCategory=10)
ggsave(paste(filename,'KEGG.pdf'),width=16, height=20)
write.table(kegg,file = paste(filename,"KEGG.xls"),sep='\t')
