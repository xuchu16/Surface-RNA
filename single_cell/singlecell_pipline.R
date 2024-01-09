#usr/bin/R
library(Seurat)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(ggsci)

themes_pop <- function() {
  p=theme(
    plot.title=element_text(size=24, face='bold', margin=margin(0,0,3,0), hjust=0.5),
    line = element_line(colour = "black", lineend = "round", linetype = "solid"),
    rect = element_rect(fill = "white",colour = "black",linetype = "solid"),
    legend.text = element_text(colour = 'midnightblue',size = 12, face = 'bold'),
    legend.title = element_text(colour = 'midnightblue', size = 24, face = 'bold'),
    axis.text.x= element_text(size=15,colour = "black",face = 'bold'),#angle = 90),
    axis.text.y= element_text(size=15,colour = "black",face = 'bold'),
    axis.title.x=element_text(size = 24,face='bold'),
    axis.title.y= element_text(size = 24,face = 'bold'),
    panel.grid = element_blank(),
    legend.key.size=unit(0.3,'cm'),
    panel.border = element_rect(fill = NA),#,color ='grey50'),
    # legend.key.height =unit(100,'cm'),
    legend.background = element_rect(colour = NA),
    legend.position = 'right')#  legend.position="bottom",
}
ref <- readRDS("pbmc_multimodal_2023.rds")
a = Read10X('filtered_feature_bc_matrix/')
data  <- CreateSeuratObject(counts = a, project = "Y1", min.cells = 3, min.features = 200)
data <- SCTransform(data, verbose = FALSE)
anchor <- FindTransferAnchors(reference = ref,
			      query = data,
			      reference.reduction = "spca",
			      normalization.method = "SCT",
			      dims = 1:50)
data <- MapQuery(anchorset = anchor,
		 query = data,
		 reference = ref,
		 refdata = list(celltype.l1 = "celltype.l1",
				celltype.l2 = "celltype.l2",
				celltype.l3 = "celltype.l3",
				predicted_ADT = "ADT"),
		 reference.reduction = "spca", "wnn.umap")


celltype1 = table(data@meta.data$predicted.celltype.l1)
write.table(celltype1,'Y4_celltype1.stat.xls',sep = '\t')
celltype2 = table(data@meta.data$predicted.celltype.l2)
write.table(celltype2,'Y4_celltype2.stat.xls',sep = '\t')
celltype3 = table(data@meta.data$predicted.celltype.l3)
write.table(celltype3,'Y4_celltype3.stat.xls',sep = '\t')


