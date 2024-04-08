library(tidyverse)
library(ggthemes)
library(ggrepel)
library(ggsci)
library(RColorBrewer)
library(wesanderson)
library(ggplot2)
themes_pop <- function() {
  p=theme(
    plot.title=element_text(size=24, face='bold', margin=margin(0,0,3,0), hjust=0.5),
    line = element_line(colour = "black", lineend = "round", linetype = "solid"),
    rect = element_rect(fill = "white",colour = "black",linetype = "solid"),
    legend.text = element_text(colour = 'midnightblue',size = 12, face = 'bold'),
    legend.title = element_text(colour = 'midnightblue', size = 24, face = 'bold'),
    axis.text.x= element_text(size=15,colour = "black"),
    axis.text.y= element_text(size=15,colour = "black"),
    axis.title.x=element_text(size = 24,face='bold'),
    axis.title.y= element_text(size = 24,face = 'bold'),
    panel.grid = element_blank(),
    legend.key.size=unit(0.3,'cm'),
    panel.border = element_rect(fill = NA),#,color ='grey50'),
    # legend.key.height =unit(100,'cm'),
    legend.background = element_rect(colour = NA),
    legend.position = 'right')#  legend.position="bottom",
}

d3 = c(
  '#1f77b4',
  '#ff7f0e',
  '#2ca02c',
  '#d62728',
  '#9467bd',
  '#8c564b',
  '#e377c2', 
  '#7f7f7f',
  '#bcbd22',
  '#17becf')
adjustdata <- function(data) {data<-cbind(rownames(data),data)}
pal <- wes_palette(30, name = "Zissou1", type = "continuous")
############read data PCA and distance

library('openxlsx')
human_surface = read.xlsx('/home/xuchu/data/RNASeq/surfaceRNA/human_surfaceRNA.TPM.20230809.xlsx',
      colNames = TRUE,rowNames = TRUE)
human_surface_zero = human_surface[which(rowSums(human_surface)==0),]
human_surface_zero = adjustdata(human_surface_zero)
write.table(human_surface_zero,'human_surface_zreo.xls',
            sep = '\t',row.names = FALSE, quote  = FALSE)
head(human_surface)
human_surface = read.csv('Human.blood.cell.surfece.all.xls',
                         header = T,row.names = 1,sep ='\t')
head(human_surface$gene_name_type)
head(human_surface$H1)
human_surface = human_surface[which(rowSums(human_surface)>0),]
data2.pca <- prcomp(human_surface, center=F, scale=F)
data.pca = data.frame(data2.pca$rotation)
PC1 = data.pca$PC1
PC2 = data.pca$PC2
sample = rownames(data.pca)
data=data.frame(PC1,PC3,sample)
data$Sample[which(data$sample == 'H1')] = 'Hematopoietic stem cell'
data$Sample[which(data$sample == 'H2')] = 'Hematopoietic stem cell'
data$Sample[which(data$sample == 'H2')] = 'Hematopoietic stem cell'
data$Sample[which(data$sample == 'T1')] = 'T cell'
data$Sample[which(data$sample == 'T2')] = 'T cell'
data$Sample[which(data$sample == 'T3')] = 'T cell'
data$Sample[which(data$sample == 'NK1')] = 'NK cell'
data$Sample[which(data$sample == 'NK2')] = 'NK cell'
data$Sample[which(data$sample == 'NK3')] = 'NK cell'
data$Sample[which(data$sample == 'B1')] = 'B cell'
data$Sample[which(data$sample == 'B2')] = 'B cell'
data$Sample[which(data$sample == 'B3')] = 'B cell'
data$Sample[which(data$sample == 'Mono1')] = 'Monocyte'
data$Sample[which(data$sample == 'Mono2')] = 'Monocyte'
data$Sample[which(data$sample == 'Mono3')] = 'Monocyte'
p = ggplot(data,aes(x=PC1, y=PC2,col=Sample
    ))+geom_point(size=3)+theme_bw()+themes_pop()+scale_color_jama()
p
ggsave('mouse_PCA_PC1_PC3.pdf',width=9,height=6)
#############heatmap
library(pheatmap)
M=cor(log(human_surface+1))
pheatmap(M,display_numbers = TRUE)
#############cluster plot
library(factoextra)
human = scale(human_surface)
human = human_surface
human = human[order(-human$H1),]
human = human[1:10000,]
human = t(human)
dd <- dist(human, method = "manhattan") 
hc <- hclust(dd, method = "ward.D2") 
fviz_dend(hc, k = 6, 
          cex = 0.8, 
          k_colors = c("#2E9FDF", "#00AFBB","cyan4", "#E7B800", "#FC4E07",'firebrick4'),
          color_labels_by_k = FALSE, 
          rect_border = c("#2E9FDF", "#00AFBB","cyan4", "#E7B800", "#FC4E07",'firebrick4'),
          rect = TRUE, 
          type ="phylogenic",
          rect_fill = TRUE)
################different expression
human_surface = read.xlsx('/home/xuchu/data/RNASeq/surfaceRNA/human_surfaceRNA.TPM.20230809.xlsx',
     colNames = TRUE,rowNames = TRUE)
human_surface = read.xlsx('../mouse_surfaceRNA.TPM.20230809.xlsx',
     colNames = TRUE,rowNames = TRUE)
human_surface = human_surface[which(rowSums(human_surface)>0),]
human_surface$hspc_mean = rowMeans(human_surface[, 1:3])
human_surface$other_mean = rowMeans(human_surface[, 4:12])
p = vector()
for (i in 1:102168){
        h = human_surface[i,][,1:3];
        o = human_surface[i,][,4:12];
        a = t.test(h,o);
        p=append(p,a$p.value);
}
human_surface$Log2Foldchange=log2((human_surface$hspc_mean+1)/(human_surface$other_mean+1))
human_surface$pvalue = p
human_surface_hspc_vs_other = human_surface

hspc_high = subset(human_surface_hspc_vs_other,
   human_surface_hspc_vs_other$hspc_mean>50)
write.table(hspc_high,'hspc_high_mouse.xls',sep = '\t')
###############point plot
library(limma)
human_surface = read.xlsx('../mouse_surfaceRNA.TPM.20230809.xlsx',
     colNames = TRUE,rowNames = TRUE)
write.table(human_surface,'mouse_surface.xls',sep = '\t')
human_surface = human_surface[which(rowSums(human_surface)>0),]

human_surface = data.frame(human_surface)
#human_surface=normalizeBetweenArrays(human_surface)
boxplot(human_surface,outline=FALSE, notch=T, las=2)
human_surface$hspc_mean = rowMeans(human_surface[, 1:3])
#human_surface$T_mean = rowMeans(human_surface[, 10:12])
hspc_high = subset(human_surface,
    human_surface$hspc_mean>50)
write.table(hspc_high,'mouse_HSPC_TPM_50.xls',sep = '\t')
boxplot(hspc_high,outline=FALSE, notch=T, las=2)
p = vector()
for (i in 1:190004){
  h = human_surface[i,][,1:3];
  o = human_surface[i,][,7:9];
  a = t.test(h,o);
  p=append(p,a$p.value);
}
human_surface$Log2Foldchange=log2((human_surface$hspc_mean)/(human_surface$T_mean))
human_surface$pvalue = p

human_surface_hspc_vs_T = human_surface
##############################vacanoplot

human_surface_hspc_vs_T = subset(human_surface_hspc_vs_T,human_surface_hspc_vs_T$hspc_mean > 1 
     & human_surface_hspc_vs_T$T_mean > 1)
human_surface_hspc_vs_T$Group = 'Not Significant'
human_surface_hspc_vs_T$Group[which(human_surface_hspc_vs_T$Log2Foldchange >= 1 
     & human_surface_hspc_vs_T$pvalue <0.05)] = 'Hspc Up-regulated'
human_surface_hspc_vs_T$Group[which(human_surface_hspc_vs_T$Log2Foldchange <= -1
     & human_surface_hspc_vs_T$pvalue <0.05)] = 'Hspc Down-regulated'

human_surface_hspc_vs_T$Logpvalue = -log(human_surface_hspc_vs_T$pvalue)
write.table(human_surface_hspc_vs_T,'human_surface_Hspc_vs_T_diff.xls',sep = '\t')

p = ggplot(human_surface_hspc_vs_T,aes(x=human_surface_hspc_vs_T$Log2Foldchange, y=human_surface_hspc_vs_T$Logpvalue,col=human_surface_hspc_vs_T$Group
))+geom_point(size = 1)+theme_bw()+themes_pop()
p = p+ylab('-log(Pvalue)')+xlab('Log2(FoldChange)')+
  guides(color=guide_legend(title = "Group"))+geom_vline(xintercept = 1,linetype = "dashed")+
  geom_hline(yintercept = -log(0.05),linetype = "dashed")+geom_vline(xintercept = -1,linetype = "dashed")
p = p+scale_color_manual(values = c("midnightblue","#FC4E07","#CCCCCC"))
p
hspc_vs_T_up = subset(human_surface_hspc_vs_T,human_surface_hspc_vs_T$Group == 'Hspc Up-regulated')
hspc_vs_T_down = subset(human_surface_hspc_vs_T,human_surface_hspc_vs_T$Group == 'Hspc Down-regulated')
write.table(hspc_vs_T_down,'human_hspc_vs_B_down.xls',sep = '\t')
write.table(hspc_vs_T_up,'human_hspc_vs_B_up.xls',sep = '\t')

hspc_diff = rbind(hspc_vs_T_up,hspc_vs_T_down)
for_heatmap = log(hspc_diff[1:12]+1)
pheatmap(hspc_diff)
data = data.frame(hspc_diff$H1,hspc_diff$H2,hspc_diff$H3,hspc_diff$T1,
                  hspc_diff$T2,hspc_diff$T3)
data = log(data+1)
pheatmap(data,cluster_rows =F,cluster_cols = T,show_rownames = F)

###################################pie plot
a = read.csv('./RNAtype/eAMOUR_RNAtype.xls',sep = '\t',header = F)
a = a[order(-a$V2),]lab.font
a$Types = factor(a$V1,levels=a$V1)
a = data.frame(Type = a$Types,Number = a$V2)
a = data.frame(a)
p<-ggpie(a,"Number",label= "Type",
   fill="Type",color=c("white","black"),lab.pos="out")
p=p+theme(legend.position = "right",legend.text = element_text(size=10))
p=pp+labs(fill="RNA Type")+scale_fill_d3()#c(pal[1],pal[2],pal[3],pal[4],pal[5],pal[6],pal[7],pal[8],pal[9],pal[10],pal[11],pal[12],pal[13]))  #lab.pos=c("out","in")
ggsave('eAMOUR_RNAtype.pdf',width = 12,height = 8)

####################################protein
f = read.csv('Y4.xls',header = T,sep = '\t')
data = data.frame(f$avg_con,f$avg,f$Log2Foldchange,f$name_2)
data$anno = 'Not Signicicant'
data$anno[which(data$f.Log2Foldchange>= 1)] = 'Binding'
data$s = 'Not Signicicant'
data$s[which(data$f.Log2Foldchange>= 1)] = 'Binding'
data$s[which(data$f.name_2 != 'NA')] = 'Histone'
data$s = factor(data$s,levels = c('Not Signicicant','Binding','Histone'))
data$size = 'other'
data$size[which(data$f.name_2 != 'NA')] = 'Histone'
data$size = factor(data$size,levels = c('other','Histone'))
  p = ggplot(data,aes(x=data$f.avg_con, y=data$f.avg,col=data$s,label = data$f.name_2,size = data$size
))+geom_point()+theme_bw()+themes_pop()
p = p+ylab('Y4 Average')+xlab('NC Y4 Average')+xlim(0,2000)+ylim(0,2000)
  guides(color=guide_legend(title = "Group"))#+geom_vline(xintercept = 1,linetype = "dashed")+
 # geom_hline(yintercept = -log(0.05),linetype = "dashed")+geom_vline(xintercept = -1,linetype = "dashed")
p = p+scale_color_manual(values = c("#CCCCCC","midnightblue","#FC4E07"))+
  geom_abline(intercept = 0, slope = 2,linetype = "dashed")
p = p+geom_text_repel()
#p
ggsave('Y4.pdf',width = 12,height = 10)
################################protein2
f = read.csv('Y3.xls',header = T,sep = '\t')
data = data.frame(f$Log2Foldchange,f$Pvalue,f$name)
data$s = 'Not Signicicant'
data$s[which(data$f.Log2Foldchange>= 1 & data$f.Pvalue<0.05)] = 'Binding'
data$s[which(data$f.name != 'NA')] = 'Histone'
data$s = factor(data$s,levels = c('Not Signicicant','Binding','Histone'))
data$size = 'other'
data$size[which(data$f.name != 'NA')] = 'Histone'
data$size = factor(data$size,levels = c('other','Histone'))
data$p2 = -log(data$f.Pvalue)
p = ggplot(data,aes(x=data$f.Log2Foldchange, y=data$p2,col=data$s,label = data$f.name,size = data$size
))+geom_point()+theme_bw()+themes_pop()
p = p+ylab('-log(p-value)')+xlab('Log2Foldchange')+
  guides(color=guide_legend(title = "Group"))+geom_vline(xintercept = 1,linetype = "dashed")+
geom_hline(yintercept = -log(0.05),linetype = "dashed")+geom_vline(xintercept = -1,linetype = "dashed")
p = p+scale_color_manual(values = c("#CCCCCC","midnightblue","#FC4E07"))#+
#  geom_abline(intercept = 0, slope = 2,linetype = "dashed")
p = p+geom_text_repel()
#p
ggsave('Y3_2.pdf',width = 12,height = 10)
###################################plot-Correlation of mouse and human
library(biomaRt)
listMarts()
human_surface = read.xlsx('../human_surfaceRNA.TPM.20230809.xlsx',
    colNames = TRUE,rowNames = TRUE)
human_surface = human_surface[which(rowSums(human_surface)>0),]
adjustdata <- function(data) {data<-cbind(rownames(data),data)}
human_surface<-adjustdata(human_surface)
write.table(human_surface,'human_surface.xls',sep = '\t',row.names = FALSE, quote  = FALSE)
a = read.csv('human_surface_new.xls',header = T,sep = '\t')
mouse <- useMart('ensembl',dataset = "mmusculus_gene_ensembl")
human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl")
name = a$gene_name
MtoH <- getLDS(attributes = "hgnc_symbol", # 要转换符号的属性，这里基因名（第3步是基因名）
               filters = "hgnc_symbol", #参数过滤
               mart = human, #需要转换的基因名的种属来源，也就是第2步的mouse
               values = name, #要转换的基因集
               attributesL = "mgi_symbol", #要同源转换的目标属性，这里还是转为基因名，也可加其他
               martL = mouse, #要同源转换的目标种属，也就是第2步的human
               uniqueRows = TRUE)
library("homologene")

o = homologene(name, inTax = 9606, outTax = 10090)
###
a = setwd('/home/xuchu/data/RNASeq/surfaceRNA/Correlation/')
a = read.csv('Human_Mouse_homologene_matrix.xls',header = T,sep = '\t')
human = a[3:32]

human_HSPC = rowMeans(data.frame(human$human.Hspc.rep1,human$human.Hspc.rep2,human$human.Hspc.rep3))
human_B = rowMeans(data.frame(human$human.B.rep1,human$human.B.rep2,human$human.B.rep3))
human_T = rowMeans(data.frame(human$human.T.rep1,human$human.T.rep2,human$human.T.rep3))
human_NK = rowMeans(data.frame(human$human.NK.rep1,human$human.NK.rep2,human$human.NK.rep3))
human_Mono = rowMeans(data.frame(human$human.Mono.rep1,human$human.Mono.rep2,human$human.Mono.rep3))

mouse_HSPC = rowMeans(data.frame(human$mouse.Hspc.rep1,human$mouse.Hspc.rep2,human$mouse.Hspc.rep3))
mouse_B = rowMeans(data.frame(human$mouse.B.rep1,human$mouse.B.rep2,human$mouse.B.rep3))
mouse_T = rowMeans(data.frame(human$mouse.T.rep1,human$mouse.T.rep2,human$mouse.T.rep3))
mouse_NK = rowMeans(data.frame(human$mouse.NK.rep1,human$mouse.NK.rep2,human$mouse.NK.rep3))
mouse_Mono = rowMeans(data.frame(human$mouse.Mono.rep1,human$mouse.Mono.rep2,human$mouse.Mono.rep3))

human2 = data.frame(human_HSPC,human_B,human_T,human_NK,human_Mono,mouse_HSPC,mouse_B,mouse_T,mouse_NK,mouse_Mono)
human2 = human2[which(rowSums(human2)>5),]
human2 = scale(human2)
human2 = t(human2)
dd <- dist(human2, method = "manhattan") 
hc <- hclust(dd, method = "ward.D2") 
fviz_dend(hc, k = 2, 
          cex = 0.8, 
          k_colors = c("#2E9FDF", "#FC4E07"),
          color_labels_by_k = TRUE, 
          rect_border = c("#2E9FDF", "#FC4E07"),
          rect = TRUE, 
          repel = TRUE,
          horiz =TRUE,
          type = "circular",
          rect_fill = TRUE)
ggsave('distance_4.pdf',width=10,height=10)

human2 = human[which(rowSums(human)>0),]
data2.pca <- prcomp(human2, center=F, scale=F)
data.pca = data.frame(data2.pca$rotation)
sample = rownames(data.pca)
PC1 = data.pca$PC1
PC2 = data.pca$PC2
data = data.frame(PC1,PC2,sample)
data$Sample = NA
data$Sample[which(data$sample == 'human.B.rep1' | 
                    data$sample == 'human.B.rep2'|
                    data$sample == 'human.B.rep3')] = 'Human B cell'
data$Sample[which(data$sample == 'human.T.rep1' | 
                    data$sample == 'human.T.rep2'|
                    data$sample == 'human.T.rep3')] = 'Human T cell'
data$Sample[which(data$sample == 'human.NK.rep1' | 
                    data$sample == 'human.NK.rep2'|
                    data$sample == 'human.NK.rep3')] = 'Human NK cell'
data$Sample[which(data$sample == 'human.Hspc.rep1' | 
                    data$sample == 'human.Hspc.rep2'|
                    data$sample == 'human.Hspc.rep3')] = 'Human Hspc cell'
data$Sample[which(data$sample == 'human.Mono.rep1' | 
                    data$sample == 'human.Mono.rep2'|
                    data$sample == 'human.Mono.rep3')] = 'Human Mono cell'

data$Sample[which(data$sample == 'mouse.Hspc.rep1' | 
                    data$sample == 'mouse.Hspc.rep2'|
                    data$sample == 'mouse.Hspc.rep3')] = 'mouse Hspc cell'
data$Sample[which(data$sample == 'mouse.B.rep1' | 
                    data$sample == 'mouse.B.rep2'|
                    data$sample == 'mouse.B.rep3')] = 'mouse B cell'
data$Sample[which(data$sample == 'mouse.T.rep1' | 
                    data$sample == 'mouse.T.rep2'|
                    data$sample == 'mouse.T.rep3')] = 'mouse T cell'
data$Sample[which(data$sample == 'mouse.NK.rep1' | 
                    data$sample == 'mouse.NK.rep2'|
                    data$sample == 'mouse.NK.rep3')] = 'mouse NK cell'
data$Sample[which(data$sample == 'mouse.Mono.rep1' | 
                    data$sample == 'mouse.Mono.rep2'|
                    data$sample == 'mouse.Mono.rep3')] = 'mouse Mono cell'

p = ggplot(data,aes(x=PC1, y=PC2,col=Sample
))+geom_point(size=3)+theme_bw()+themes_pop()+scale_color_d3()
p
ggsave('Mouse_human_PC1_PC2.pdf',width=9,height=7)
#############################################
setwd('/home/xuchu/data/RNASeq/surfaceRNA/293T/')
HEK293T = read.xlsx('./HEK293T_surfaceRNA.TPM.20230903.xlsx',
   colNames = TRUE,rowNames = TRUE)
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}
a = adjustdata (data)
write.table(a,file = "HEK293T_surfaceRNA.xls",sep = "\t",quote = FALSE,row.names = FALSE)
data = data.frame(HEK293T)
data$mean = rowMeans(data[, 1:2])
data = data[order(-data$mean),]
data = data[1:1000,]
par(mar= c(10,10,10,10))
p = ggplot(data,aes(x=Rep1, y=Rep2
))+geom_point(size=1)+theme_bw()+themes_pop()+scale_color_d3(
)+geom_smooth(method = 'lm', se = F, color = 'red')+
  stat_cor(data=data, method = "spearman"
)+xlab('HEK293T Rep1 (TPM)')+ylab('HEK293T Rep2 (TPM)')
p= p+theme(
#  plot.background = element_rect(fill = "orange", color = "black", size = 10),
  plot.title = element_text(hjust = 1, color = "red", face = "italic"),
  plot.margin = margin(t = 30, r = 30, b = 30, l = 30, unit = "pt")
)
p
###################
f = read.table('HEK293T.table.xls',header = T,sep ='\t')
data = data.frame(f[1:15,])
data$name = factor(data$gene_name,levels = data$gene_name)
data
p = ggplot(data,aes(x=name,y=mean
,fill = name))+geom_bar(stat = 'identity')+theme_bw()+themes_pop()
p = p+scale_fill_manual(values = pal)+ylab('TPM')+xlab('Gene name')
p
hspc_high = subset(data,
   data$mean>50)
write.table(hspc_high,'HEK293T_TPM50.xls',sep = '\t')
###################WGA
WGA = read.xlsx('./WGA_surfaceRNA.TPM.20230907.xlsx',
                    colNames = TRUE,rowNames = TRUE)
WGA = WGA[which(rowSums(WGA)>0),]
WGA$TWGA_mean = rowMeans(WGA[, 1:2])
WGA$WGA_mean = rowMeans(WGA[, 3:4])
p = vector()
for (i in 1:124720){
  h = WGA[i,][,1:2];
  o = WGA[i,][,3:4];
  a = t.test(h,o);
  p=append(p,a$p.value);
}
pv_adjust = p.adjust(p,method = 'BH')
WGA$Log2Foldchange=log2((WGA$TWGA_mean+1)/(WGA$WGA_mean+1))
WGA$pvalue = p
WGA$FDR = pv_adjust
data = adjustdata(WGA)
write.table(data,'TWGA_VS_WGA.xls',sep = '\t',row.names = FALSE, quote  = FALSE)
TWGA_down = subset(data,
   data$Log2Foldchange<-1&data$pvalue<0.05)
write.table(TWGA_down,'TWGA_VS_WGA_down.xls',sep = '\t',row.names = FALSE, quote  = FALSE)

########################
library(viridis)
library(Polychrome)
Polychrome::swatch(viridis(20))
#######################################TopTPM for surface RNA 
library('ggbreak')
a = read.csv('T-WGA_count_unique',sep = '\t',header = F)
a$TPM = a$V2/2
data = data.frame(a)
data = data[order(-data$TPM),]
data_unique = data[!duplicated(data$TPM),]
data_unique$V1[data_unique$V1=='ENSG00000286171']='Y5' 
data_unique$V1[data_unique$V1=='ENSG00000283907']='LINC01766'
data_unique = data_unique %>% 
   filter(!data_unique$V1 %in% c('RNA5-8SN3','RNA5S8',"RNA5S1",'5_8S_rRNA','RNA5S9'))
data_plot = data_unique[1:15,]
#data_plot$V1[data_plot$V1=='RNA5-8SN3']='5_8S_rRNA'

data_plot$name = factor(data_plot$V1,levels = data_plot$V1)
p = ggplot(data_plot,aes(x=name,y=TPM
   ,fill = name))+geom_bar(stat = 'identity')+theme_bw()+themes_pop()
p = p+scale_fill_manual(values = rev(pal))+ylab('Read count')+theme( axis.text.x= element_text(size=15,
    colour = "black",angle = 45,vjust = 0.5))+
  xlab('Gene name')+ggtitle('T-WGA_top_count')##+scale_y_break(c(40000, 180000), scales=0.4)
p
ggsave('T-WGA_TOP_count.pdf',width = 12,height = 8)
#########################################################
#########################################################
library('ggbreak')
a = read.csv('Mouse.NK.count_unique',sep = '\t',header = F)
a$TPM = a$V2/3
data = data.frame(a)
data = data[order(-data$TPM),]
data_unique = data[!duplicated(data$TPM), ]
data_unique$V1[data_unique$V1=='ENSG00000286171']='Y5' 
data_unique$V1[data_unique$V1=='ENSG00000283907']='LINC01766'
data_unique = data_unique %>% 
  filter(!data_unique$V1 %in% c('Rn18s-rs5','n-R5s100',"n-R5s106",'5_8S_rRNA'))
data_plot = data_unique[1:15,]
#data_plot$V1[data_plot$V1=='RNA5-8SN3']='5_8S_rRNA'
data_plot$name = factor(data_plot$V1,levels = data_plot$V1)
p = ggplot(data_plot,aes(x=name,y=TPM
  ,fill = name))+geom_bar(stat = 'identity')+theme_bw()+themes_pop()
p = p+scale_fill_manual(values = rev(pal))+ylab('Read count')+theme( axis.text.x= element_text(size=15,
  colour = "black",angle = 45,vjust = 0.5))+
  xlab('Gene name')+ggtitle('mouse_NK_top_count')##+scale_y_break(c(40000, 180000), scales=0.4)
p
ggsave('mouse_NK_TOP_count.pdf',width = 12,height = 8)
###########################################################
###########################################################
####
WGA = read.xlsx('./WGA_surfaceRNA.TPM.20230907.xlsx',
                colNames = TRUE,rowNames = TRUE)
WGA = WGA[which(rowSums(WGA)>0),]
WGA$TWGA_mean = rowMeans(WGA[, 1:2])
WGA$WGA_mean = rowMeans(WGA[, 3:4])
###########################################################
a = read.csv('HEK293T_rep1_rep2_stat',sep = '\t',header = F)
data = data.frame(a)
data = data[order(-data$V4),]
data$V1= factor(data$V1,levels = unique(data$V1))
p = ggplot(data,aes(x=data$V3,y=data$V4
 ,fill = data$V1))+geom_bar(stat = 'identity',position=position_dodge()
 )+theme_bw()+themes_pop()+scale_fill_discrete(name = "hahah")
p = p+ylab('Proportion')+xlab('Sample name')+scale_fill_d3("category20")+
  scale_y_continuous(labels = scales::percent_format(scale = 100))+#,limits=c(0, 1))
  labs(fill = 'RNA type')+ggtitle('293T TPM > 50 gene proportion')
p
ggsave('HEK293T_rep1_rep2_proportion.pdf',width = 12,height = 8)
##########################################################
a = read.csv('HEK293T_rep1_rep2_stat',sep = '\t',header = F)
data = data.frame(a)
data = data[order(-data$V2),]
data$V1= factor(data$V1,levels = unique(data$V1))
p = ggplot(data,aes(x=data$V3,y=data$V2
  ,fill = data$V1))+geom_bar(stat = 'identity',position=position_dodge()
  )+theme_bw()+themes_pop()+scale_fill_discrete(name = "hahah")
p = p+ylab('Number')+xlab('Sample name')+scale_fill_d3("category20")+
 # scale_y_continuous(labels = scales::percent_format(scale = 100))+#,limits=c(0, 1))
  labs(fill = 'RNA type')+ggtitle('293T TPM > 50 gene number')
p
ggsave('HEK293T_rep1_rep2_number.pdf',width = 12,height = 8)
##############################################################
#########
#########  RNA analysis from  feturecount
##############################################################
library('RNAontheBENCH')
a = read.csv('../shaoxin/eAMOUR_counts.txt',sep ='\t',header = T)
e5 = counts2fpkm(a$eAMOUR0905Aligned.sortedByCoord.out.bam,lengths=a$Length)
eAMOUR_rep1 = fpkm2tpm(e5)
e19 = counts2fpkm(a$eAMOUR0919Aligned.sortedByCoord.out.bam,lengths=a$Length)
eAMOUR_rep2 = fpkm2tpm(e19)
data = data.frame(a$Geneid,eAMOUR_rep1,eAMOUR_rep2)
write.table(data,'eAMOUR_TPM.xls',sep = '\t',row.names = FALSE, quote  = FALSE)
#############################################################
setwd('/home/xuchu/data/RNASeq/shaoxin/combine')
a = read.csv('overlap',header = T,sep = '\t')
data = a[2:7]
M=cor(log(data+1))
pheatmap(M,display_numbers = TRUE)
############################################################
datas = t(data)
dd <- dist(datas, method = "manhattan") 
hc <- hclust(dd, method = "ward.D2") 
fviz_dend(hc, k = 4, 
          cex = 0.8, 
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = FALSE, 
          rect_border = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          rect = TRUE, 
         # type ="circular",
          rect_fill = TRUE)
################################################################
a = read.csv('./count/eAMOUR_count_for_plot.xls',header = F,sep = '\t')
a$Mean = rowMeans(a[, 2:3])
a = a[order(-a$Mean),]
data = data.frame(a[1:15,])
data$name = factor(data$V4,levels = data$V4)
data
p = ggplot(data,aes(x=name,y=Mean
 ,fill = name))+geom_bar(stat = 'identity')+theme_bw()+themes_pop()
p = p+scale_fill_manual(values = rev(pal))+
   ylab('Reads count')+theme( axis.text.x= element_text(size=15,
   colour = "black",angle = 45,vjust = 0.5))+
  xlab('Gene name')+ggtitle('eAMOUR_top_Tcount')##+scale_y_break(c(40000, 180000), scales=0.4)
p
ggsave('eAMOUR_top_count.pdf',width = 12,height = 8)
###############################################################
a = read.csv('././../eAMOUR_TPM.table.xls',header = T,sep = '\t')
a$Mean = rowMeans(a[, 2:3])
hspc_high = subset(a,
      a$Mean>50)
write.table(hspc_high,'eAMOUR_TPM_50.xls',sep = '\t')
################################################################
SampleName = c('HEK293T rep1','HEK293T rep2',
  'WGA rep1','WGA rep2','eAMOUR rep1','eAMOUR rep2')
RNY1 = c(2919.98,2989.3,512.93,620.387,0,0)
RNY3 = c(2366.96,2265.39,1411.72,1227.5,0,0)
RNY5 = c(2827.86,2756.38,883.317,365.223,0,0)
data = data.frame(SampleName,RNY1,RNY3,RNY5)
a = read.csv('target.xls',header = F,sep = '\t')
p = ggplot(a,aes(x=a$V2,y=a$V1
   ,fill = a$V3))+geom_bar(stat = 'identity',position='dodge')+theme_bw()+themes_pop()
p = p+scale_fill_d3()+
  ylab('TPM')+theme(axis.text.x= element_text(size=15,
      colour = "black",angle = 45,vjust = 0.5))+
  xlab('Sample')+ggtitle('RNY TPM')+labs(fill = 'RNA Name')##+scale_y_break(c(40000, 180000), scales=0.4)
p
ggsave('RNY_TPM.pdf',width = 12,height = 8)
################################################################
a = read.csv('RNAtype_plot1.xls',sep = '\t',header = T)
data = data.frame(a)
data = data[order(-data$HEK293T),]
data$group= factor(data$group,levels = unique(data$group))
a = read.csv('RNAtype_plot2.xls',sep = '\t',header = T)
data2 = data.frame(a)
data2$groups = factor(data2$group,levels = data$group)
data2$Strategy = factor(data2$Strategy,levels = c('HEK293T','WGA','eAMOUR'))
p = ggplot(data2,aes(x=data2$Strategy,y=data2$percentage
                    ,fill = data2$groups))+geom_bar(stat = 'identity',position=position_dodge()
                    )+theme_bw()+themes_pop()+scale_fill_discrete(name = "hahah")
p = p+ylab('Percentage')+xlab('Strategy')+scale_fill_d3("category20")+
  theme(axis.text.x= element_text(size=15,
  colour = "black",angle = 45,vjust = 0.5))+
   scale_y_continuous(labels = scales::percent_format(scale = 100))+#,limits=c(0, 1))
  labs(fill = 'RNA type')+ggtitle('RNA type of three Strategy')
p
ggsave('RNAtype_three_strategy.pdf',width = 14,height = 8)
#################################################################
#######
f = read.csv('MTRNR2.xls',header = T,sep = '\t')
data = data.frame(f$avg_con,f$avg,f$Log2Foldchange,f$name,f$Pvalue)
data$anno = 'Not Signicicant'
data$anno[which(data$f.Log2Foldchange>= 1)] = 'Binding'
data$s = 'Not Signicicant'
data$s[which(data$f.Log2Foldchange>= 1)] = 'Binding'
data$s[which(data$f.name != 'NA' & data$f.Log2Foldchange>= 1 )] = 'Histone'
data$s = factor(data$s,levels = c('Not Signicicant','Binding','Histone'))
data$size = 'other'
data$size[which(data$f.name != 'NA')] = 'Histone'
data$size[which(data$f.Log2Foldchange < 1)] = 'other'
data$size = factor(data$size,levels = c('other','Histone'))
p = ggplot(data,aes(x=data$f.Log2Foldchange, y=-log(data$f.Pvalue),
                    col=data$s,label = data$f.name,size = data$size
))+geom_point()+theme_bw()+themes_pop()
p = p+ylab('R2+UV')+xlab('R2-UV')+#+xlim(0,25000)+ylim(0,25000)
 guides(color=guide_legend(title = "Group"))+geom_vline(xintercept = 1,linetype = "dashed")+
 geom_hline(yintercept = -log(0.05),linetype = "dashed")+geom_vline(xintercept = -1,linetype = "dashed")
p = p+scale_color_manual(values = c("#CCCCCC","midnightblue","#FC4E07"))
p = p+geom_text_repel()
p
ggsave('MTRNR2222.pdf',width = 12,height = 10)
##################################################################
#######

library('RNAontheBENCH')
a = read.csv('Moncyte_readcount.txt',sep ='\t',header = T)
human_monocyte1 = counts2fpkm(a$h.m.1Aligned.sortedByCoord.out.bam,lengths=a$Length)
human_monocyte_rep1 = fpkm2tpm(human_monocyte1)
human_monocyte2 = counts2fpkm(a$h.m.2Aligned.sortedByCoord.out.bam,lengths=a$Length)
human_monocyte_rep2 = fpkm2tpm(human_monocyte2)
human_monocyte3 = counts2fpkm(a$h.m.3Aligned.sortedByCoord.out.bam,lengths=a$Length)
human_monocyte_rep3 = fpkm2tpm(human_monocyte3)
data = data.frame(a$Geneid,human_monocyte_rep1,human_monocyte_rep2,human_monocyte_rep3)
write.table(data,'human_moncyte_TPM2.xls',sep = '\t',row.na
###############################
library(pheatmap)
a = read.csv('Y5_library.TPM.xls',sep = '\t',header = T)
a = a[2:4]
M=cor(log(a+1))
pheatmap(M,display_numbers = TRUE)
###############################
library(factoextra)
human_surface = read.csv('Human.blood.cell.surfece.all.xls',
                         header = T,row.names = 1,sep ='\t')
human = scale(human_surface)
human = human_surface
human = human[order(-human$H1),]
human = human[1:10000,]
human = t(human)
human = scale(human)
human = t(human)
dd <- dist(human, method = "manhattan") 
hc <- hclust(dd, method = "ward.D2") 
fviz_dend(hc, k = 4, 
          cex = 0.8, 
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, 
          rect_border = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          rect = TRUE, 
          repel = TRUE,
          #  k_colors = 'jco',
          type = "phylogenic",
          rect_fill = TRUE,
          phylo_layout = "layout_as_tree")
#########################
library(factoextra)
human_surface = read.csv('/home/xuchu/data/RNASeq/moncyte/mouse/combine/Mouse_surface_RNA.xls',
                         header = T,row.names = 1,sep ='\t')
human = scale(human_surface)
human = human_surface
human = human[order(-human$mouse.B.rep1),]
human = human[1:8000,]
human = t(human)
dd <- dist(human, method = "manhattan") 
hc <- hclust(dd, method = "ward.D2") 
fviz_dend(hc, k = 5, 
          cex = 0.8, 
          k_colors = c("#2E9FDF", "#00AFBB","cyan4", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, 
          rect_border = c("#2E9FDF", "#00AFBB","cyan4", "#E7B800", "#FC4E07"),
          rect = TRUE, 
          repel = TRUE,
          type ="phylogenic",
          rect_fill = TRUE
        )

######################

a = read.csv('overlap_matrix.xls',header = T,sep = '\t')
library(pheatmap)
data = a[4:9]
#data$Total = rowSums(data)
#data = data[order(-data$Total),]
#data = data[1:1000,]
data = data.frame(data$HEK293T.rep1,data$HEK293T.rep2,data$WGA.rep1,
                 data$WGA.rep2,data$eAMOUR.rep1,data$eAROUR.rep2)
M=cor(log(data+1))
pheatmap(M,display_numbers = TRUE,cluster_rows = F,cluster_cols = F)
col_dist = dist(M,method = 'euclidean')
hclust_1 <- hclust(col_dist)
manual_order = c("HEK293T.rep1", "HEK293T.rep2", 
   "WGA.rep1", "WGA.rep2", "eAMOUR.rep1",  "eAMOUR.rep2")
dend = reorder(as.dendrogram(hclust_1), 
    wts=order(match(manual_order, rownames(M)))
row_cluster = as.hclust(dend)
pheatmap(M, cluster_rows = row_cluster)
##############################################

human2 = human[which(rowSums(human)>150),]
human2 = log(human2+1)
pheatmap(human2,show_rownames = F)
#############################################
human_Hspc = rowMeans(data.frame(human$human.Hspc.rep1,human$human.Hspc.rep2,human$human.Hspc.rep3))
human_B = rowMeans(data.frame(human$human.B.rep1,human$human.B.rep2,human$human.B.rep3))
human_T = rowMeans(data.frame(human$human.T.rep1,human$human.T.rep2,human$human.T.rep3))
human_NK = rowMeans(data.frame(human$human.NK.rep1,human$human.NK.rep2,human$human.NK.rep3))
human_Mono = rowMeans(data.frame(human$human.Mono.rep1,human$human.Mono.rep2,human$human.Mono.rep3))

 mouse_Hspc = rowMeans(data.frame(human$mouse.Hspc.rep1,human$mouse.Hspc.rep2,human$mouse.Hspc.rep3))
mouse_B = rowMeans(data.frame(human$mouse.B.rep1,human$mouse.B.rep2,human$mouse.B.rep3))
mouse_T = rowMeans(data.frame(human$mouse.T.rep1,human$mouse.T.rep2,human$mouse.T.rep3))
mouse_NK = rowMeans(data.frame(human$mouse.NK.rep1,human$mouse.NK.rep2,human$mouse.NK.rep3))
mouse_Mono = rowMeans(data.frame(human$mouse.Mono.rep1,human$mouse.Mono.rep2,human$mouse.Mono.rep3))

human2 = data.frame(human_Hspc,human_B,human_T,human_NK,human_Mono,mouse_Hspc,mouse_B,mouse_T,mouse_NK,mouse_Mono)
human2 = human2[which(rowSums(human2)>5000),]
human2 = log(human2+1)
pheatmap(human2,show_rownames = T)
###############################################
human2 = human[which(rowSums(human)>5000),]
da = read.csv('Human_surface_RNA.xls',sep = '\t',header = T)
data = da[8:22]
rownames(data) = da$gene_name

####################
a =read.csv('significant.xls',header = F ,sep = '\t')
data = log(a[2:15]+1)
pheatmap(data,cluster_cols = F,cluster_rows = F)
#####################
data1 = read.csv('mouse_specfic.xls',header = F,sep = '\t',row.names = 1)
data = data1[1:5]
dd = t(data)
as = t(scale(dd))
pdf('celltype_specfic_mouse_tRNA.pdf',width=6,height=9)
annotation_row = data.frame(celltype =data1$V8)
annotation_col = data.frame(
  CellType = c('B cell','T cell','Hspc','NK','Monocyte')
)
rownames(annotation_col) <- colnames(as)
rownames(annotation_row) <- rownames(as)
pheatmap(as,cluster_rows = F,cluster_cols = F,
        labels_col = c('B cell','T cell','Hspc','NK','Monocyte'),
        annotation_col= annotation_col,annotation_row = annotation_row,show_rownames = F,
        color = colorRampPalette(colors = c("#51B187","white","#D96558"))(100)
)
dev.off()
####################
data = read.csv('matrix3',header = F,sep = '\t',row.names = 1)
data = data[1:15]
dd = t(data)
as = t(scale(dd))
pdf('celltype_specfic_surfaceRNA3.pdf',width=6,height=9)
pheatmap(as,cluster_rows = F,cluster_cols = F,
         labels_col = c('B cell rep1','B cell rep2','B cell rep3',
                        'T cell rep1','T cell rep2','T cell rep3',
                        'Hspc rep1','Hspc rep2','Hspc rep3',
                        'NK rep1','NK rep2','NK rep3',
                        'Monocyte rep1','Monocyte rep2','Monocyte rep3'),show_rownames = F)
dev.off()
######################################
a = read.csv('both_high_express.xls',header = T,sep = '\t')
human_B = rowMeans(data.frame(a$Human.Bcell.Rep1,a$Human.Bcell.Rep2,a$Human.Bcell.Rep3))
human_T = rowMeans(data.frame(a$Human.Tcell.Rep1,a$Human.Tcell.Rep2,a$Human.Tcell.Rep3))
human_Hspc = rowMeans(data.frame(a$Human.HSPC.Rep1,a$Human.HSPC.Rep2,a$Human.HSPC.Rep3))
human_NK = rowMeans(data.frame(a$Human.NK.Rep1,a$Human.NK.Rep2,a$Human.NK.Rep3))
human_Mono = rowMeans(data.frame(a$h.m.1,a$h.m.2,a$h.m.3))
mouse_B = rowMeans(data.frame(a$Mouse.Bcell.Rep1,a$Mouse.Bcell.Rep2,a$Mouse.Bcell.Rep3))
mouse_T = rowMeans(data.frame(a$Mouse.Tcell.Rep1,a$Mouse.Tcell.Rep2,a$Mouse.Tcell.Rep3))
mouse_Hspc = rowMeans(data.frame(a$Mouse.HSPC.Rep1,a$Mouse.HSPC.Rep2,a$Mouse.HSPC.Rep3))
mouse_NK = rowMeans(data.frame(a$Mouse.NK.Rep1,a$Mouse.NK.Rep2,a$Mouse.NK.Rep3))
mouse_Mono = rowMeans(data.frame(a$m.m.1,a$m.m.2,a$m.m.3))
name = paste(a$name,a$mouse.gene, sep = "/")
data = data.frame(a$name,human_B,human_T,human_Hspc,human_NK,human_Mono,
                  mouse_B,mouse_T,mouse_Hspc,mouse_NK,mouse_Mono)
#####################################
a = read.csv('overlap_high_GO.xls',sep = '\t',header = T)
Term = a$Description
q  = -a$Log10.q.
celltype = a$celltype
data = data.frame(Term,q,celltype)
data$Term = factor(Term,levels = rev(unique(Term)))
p = ggplot(data,aes(x=Term,y=q
  ,fill = celltype))+geom_bar(stat = 'identity')+theme_bw()+themes_pop()+coord_flip()+
  theme_classic(base_size = 14)+scale_fill_d3()+ylab('-log(qvalue)')+
  ggtitle('GO of High expressed in all celltype surface RNA')+xlab('GO Term')
p
ggsave('GO_of_High_expressed_in_all_celltype_surface_RNA.pdf',height = 2,width = 16)

######################entropy matrix
library(pheatmap)
a = read.csv('matrix2',header = F,sep = '\t')
data = t(scale(t(a[2:6])))
#pdf('celltype.pdf',width = 6,height = 9)
pheatmap(data,cluster_rows = F,cluster_cols = F,labels_col= c("B cell",
        "T cell","Hspc","NK","Monocyte"))
#dev.off()
#######################
a = read.csv('expression_matrix_celltype_specfic.xls',sep = '\t',head = F,row.names  = 1)
data = t(a)
exp = data
#######################
H2Aa = c(227.992,191.163,148.424,0.124172,0.147609,0.188432,0.0927966,0.244919,
          0.101991,28.4714,31.5784,33.5678,9.07326,3.03232,8.07538)
Grip2=c(242.084,379.082,408.427,49.2861,54.0561,42.9435,73.7005,107.943,87.9362,
        114.954,165.345,163.198,0.816918,0.851926,0.754341)
Snord11=c(38.0808,116.863,142.29,1000.93,707.741,631.294,31.1601,42.6529,40.0898,
          163.602,133.761,121.286,399.79,467.915,524.086)
Snord23=c(103.599,174.905,184.245,1210.63,1095.52,870.396,22.3953,41.8414,30.0359,
          251.892,270.58,162.884,40.8239,40.6642,37.5552)
Gm55936=c(104997,47251.6,35420.1,32394.8,46133.8,36765.8,277626,247783,264969,
          85016.6,101287,103083,40105.2,46733.2,58046.1)
Gm55286=c(24.2124,3.48334,2.64756,8.65107,18.5904,9.98871,86.5464,106.146,54.0974,
          59.0586,24.7012,35.5329,0,29.3046,27.0641)
Mir139=c(4.8606,3.29083,12.3398,100.196,48.6284,48.3691,1.56761,0,0,
         131.817,145.338,91.7238,17.723,0,0)
mtTl2=c(371.644,144.976,127.508,377.748,372.888,138.37,25.1607,73.969,65.6813,
         849.076,636.285,689.633,390.514,482.342,244.287)
Rny3=c(355.072,356.734,275.416,994.57,2202.54,902.498,289.22,492.57,392.183,
       1163.88,1731.8,1637.43,2515.46,2163.64,3103.62)
Mir8114=c(28.1562,60.209,67.4718,29.4721,34.8331,54.0407,15.5409,20.1804,21.4735,
          36.4057,39.5766,76.2391,382.093,206.84,230.753)
group = factor(c('B cell','T cell','HSPC','NK','Monocyte'),
               levels = c('B cell','T cell','HSPC','NK','Monocyte'))
anno = c('B cell','B cell','B cell','T cell','T cell','T cell',
         'HSPC','HSPC','HSPC','NK','NK','NK','Monocyte','Monocyte','Monocyte')
data = data.frame(H2Aa,Grip2,Snord11,Snord23,Gm55936,Gm55286,Mir139,mtTl2,Rny3,Mir8114,anno)
data$celltype = factor(anno,levels = group)
p1 = ggplot(data,aes(x = celltype,y = H2Aa,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('H2-Aa')+ylab('TPM')+scale_fill_d3("category10")
p2 = ggplot(data,aes(x = celltype,y = Grip2,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('Grip2')+ylab('TPM')+scale_fill_d3("category10")
p3 = ggplot(data,aes(x = celltype,y = Snord11,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('Snord11')+ylab('TPM')+scale_fill_d3("category10")
p4 = ggplot(data,aes(x = celltype,y = Snord23,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('Snord23')+ylab('TPM')+scale_fill_d3("category10")
p5 = ggplot(data,aes(x = celltype,y = Gm55936,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('Gm55936')+ylab('TPM')+scale_fill_d3("category10")
p6 = ggplot(data,aes(x = celltype,y = Gm55286,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('Gm55286')+ylab('TPM')+scale_fill_d3("category10")
p7 = ggplot(data,aes(x = celltype,y = Mir139,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('Mir139')+ylab('TPM')+scale_fill_d3("category10")
p8 = ggplot(data,aes(x = celltype,y = mtTl2,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('mt-Tl2')+ylab('TPM')+scale_fill_d3("category10")
p9 = ggplot(data,aes(x = celltype,y = Rny3,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('Rny3')+ylab('TPM')+scale_fill_d3("category10")
p10 = ggplot(data,aes(x = celltype,y = Mir8114,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('Mir8114')+ylab('TPM')+scale_fill_d3("category10")

p = (p1+p2)/(p3+p4)/(p5+p6)/(p7+p8)/(p9+p10)

ggsave('mouse_boxplot.pdf',height  = 20,width  = 14)
##########################
RNVU115=c(112.077,150.656,145.215,33.6861,19.5475,17.5826,18.3833,19.3319,19.5709,
           30.0869,13.5109,19.2072,119.258,57.3231,48.3355)
MIR339=c(17.0739,30.2411,19.7921,1.15216,1.86974,4.05956,13.5283,12.9821,7.24559,
         1.46652,0.442518,1.06715,9.09527,8.40906,5.19141)
RNA5SP442=c(9.54947,6.84101,3.9857,30.705,106.989,8.31612,11.4377,15.6626,7.14775,
            8.88849,0.708192,10.5825,28.3247,19.6464,29.6683)
RNVU128=c(945.149,555.053,203.477,4797.52,591.053,426.867,589.477,
           69.2388,233.317,1003.54,47.0567,101.829,192.011,78.0354,91.6854)
SNORD96A=c(299.686,848.995,463.926,30.9358,29.0142,783.611,1428.27,1102.1,1149.25,
           24.8589,54.9943,209.072,510.257,292.917,309.739)
SNORD15A=c(89.4583,47.2038,79.6789,140.279,59.0352,68.3203,632.034,964.781,857.264,
           112.176,83.1404,110.655,404.982,353.267,195.748)
#HBA2=c(12.888,7.85946,3.64277,7.44676,13.5699,0.568534,5.05843,6.25655,7.27975,
#       642.442,3414.78,573.66,23.1793,19.2216,5.728)
H3C14=c(0.738112,0.296366,0.271479,0.792184,0.778789,5.48E-05,1.50251,1.40981,0.850695,
        27.8341,21.6206,14.8433,13.374,12.4103,9.4492)

#NKG7=c(0.0505772,0.0547468,0.00104331,1.9823,0.702706,0.315329,0.277036,0.224131,0.0985902,
#       58.733,21.8002,23.556,8.36999,7.27282,11.1852)
IFITM1=c(0.138168,0.139093,0.389518,12.0487,6.30369,1.83195,0.611464,0.508495,0.944583,
         36.6079,17.5053,30.5318,10.3118,2.79301,7.55181)
Y_RNA=c(146.897,80.8119,180.453,217.691,27.mouse_celltype_avg_tRNA.xls2195,82.9347,99.6627,133.592,136.304,
        169.29,63.748,283.817,1032.68,734.212,753.618)
S100A8=c(1.33625,2.73409,3.95935,8.19662,0.522068,0.314217,0.835034,8.07262,
         10.1448,3.82544,4.60194,7.95546,182.767,408.067,321.749)

data = data.frame(RNVU115,MIR339,Snord11,Snord23,Gm55936,Gm55286,Mir139,mtTl2,Rny3,Mir8114,anno)
data$celltype = factor(anno,levels = group)
p1 = ggplot(data,aes(x = celltype,y = RNVU115,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('RNVU1-15')+ylab('TPM')+scale_fill_d3("category10")
p2 = ggplot(data,aes(x = celltype,y = MIR339,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('MIR339')+ylab('TPM')+scale_fill_d3("category10")
p3 = ggplot(data,aes(x = celltype,y = RNA5SP442,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('RNA5SP442')+ylab('TPM')+scale_fill_d3("category10")
p4 = ggplot(data,aes(x = celltype,y = RNVU128,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('RNVU1-28')+ylab('TPM')+scale_fill_d3("category10")
p5 = ggplot(data,aes(x = celltype,y = SNORD96A,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('SNORD96A')+ylab('TPM')+scale_fill_d3("category10")
p6 = ggplot(data,aes(x = celltype,y = SNORD15A,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('SNORD15A')+ylab('TPM')+scale_fill_d3("category10")
p7 = ggplot(data,aes(x = celltype,y = HBA2,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('HBA2')+ylab('TPM')+scale_fill_d3("category10")
p8 = ggplot(data,aes(x = celltype,y = IFITM1,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('IFITM1')+ylab('TPM')+scale_fill_d3("category10")
p9 = ggplot(data,aes(x = celltype,y = Y_RNA,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('Y_RNA')+ylab('TPM')+scale_fill_d3("category10")
p10 = ggplot(data,aes(x = celltype,y = S100A8,fill = celltype))+
  geom_boxplot()+theme_bw()+themes_pop()+ggtitle('S100A8')+ylab('TPM')+scale_fill_d3("category10")

p = (p1+p2)/(p3+p4)/(p5+p6)/(p7+p8)/(p9+p10)

ggsave('human_boxplot.pdf',height  = 20,width  = 14)

##########################
setwd('/home/xuchu/data/RNASeq/moncyte/mouse/combine')
a = read.csv('mouse_surface_T_TPM.xls',header=T,sep = '\t')
a$Mean = rowMeans(data.frame(a[,2],a[,3],a[,4]))
a = a[order(a$Mean, decreasing = TRUE),]
a = a[1:150,]
write.table(table(a$gene.type),'T_Top150_RNAtype.xls',sep = '\t')

##########################
a = read.csv('mouse_Top_150_percentage.xls',header = F,sep='\t')
data = data.frame(a)
le = c("snRNA",
       "rRNA",
       "Mt_tRNA",
       "protein_coding",
       "snoRNA",
       "miRNA",
       "misc_RNA",
       "Mt_rRNA",
       "nonsense_mediated_decay",
       "other"
)
data$anno = factor(data$V1,levels = le)
p = ggplot(data,aes(x=data$V3,y=data$V2
                    ,fill = data$anno))+geom_bar(stat = 'identity',position=position_dodge()
                    )+theme_bw()+themes_pop()+scale_fill_discrete(name = "hahah")
p = p+ylab('Proportion')+xlab('Sample name')+scale_fill_d3("category20")+
  scale_y_continuous(labels = scales::percent_format(scale = 100))+#,limits=c(0, 1))
  labs(fill = 'RNA type')+ggtitle('293T TPM > 50 gene proportion')
p
ggsave('HEK293T_rep1_rep2_proportion.pdf',width = 12,height = 8)

##########################
library(ComplexHeatmap)
library(pheatmap)
setwd('/home/xuchu/data/RNASeq/moncyte/mouse/combine/heatmap/uniquename/mouse_entrpy')
a = read.csv('mouse_entrpy_matrix.xls',sep='\t',header = F)
data = a[2:6]
data = t(scale(t(data)))
Heatmap(data,cluster_rows = FALSE,cluster_columns = FALSE,cluster_row_slices = TRUE,
        row_km = 2)
############################
library('eulerr')
library(RColorBrewer)
display.brewer.all()
########get colour
detail <- brewer.pal(n = 9, name = "Set1")
detail

completely_contained <- euler(c("AMOUR detect" = 374, "Ac4ManNAz-enriched RNAs" = 22,
                                "AMOUR detect&Ac4ManNAz-enriched RNAs" = 59))
pdf('veen.pdf',width = 8,height = 6)
plot(completely_contained,
     quantities = list(type = c("counts")),
     labels = list(col = c("black", "black", "black")),
     edges = list(col = "white", lex = 1),
     fills = c('#1f77b4','#d62728','#ff7f0e'))
#ggsave('veen.pdf',width = 8,height = 6)
dev.off()
################################
setwd('/home/xuchu/data/jiangx_10X/result/annotation_bySeurat_PBMC/GO/MonocyteGO')
a = read.csv('GO_matrix.xls',header = F,sep = '\t')
p <- ggplot(a,aes(x=a$V5,y=a$V1,colour=-1*log10(a$V4),size=a$V3))+
  geom_point()+
  scale_size(range=c(2, 8))+
  scale_colour_gradient(low = "blue",high = "red")+
  theme_bw()+
  ylab("GO Terms")+
  xlab("Y RNA")+
  labs(color=expression(-log[10](adj-PValue)))+
  labs(size=expression(GeneRatio))+
  theme(legend.title=element_text(size=14), legend.text = element_text(size=14))+
  theme(axis.title.y = element_text(margin = margin(r = 50)),axis.title.x = element_text(margin = margin(t = 20)))+
  theme(axis.text.x = element_text(face ="bold",color="black",angle=0,vjust=1))
p
plot
#保存图片
ggsave(plot,filename = "KEGG.pdf",width = 10,height = 6,dpi=300)
ggsave(plot,filename = "KEGG.png",width = 10,height = 6,dpi=300)
################################################
setwd('/home/xuchu/data/RNASeq/moncyte/mouse/combine/heatmap/uniquename/mouse_entrpy')
a = read.csv('mouse_entrpy_matrix.xls',sep='\t',header = F)
data = a[2:6]
data = t(scale(t(data)))
pheatmap(data, cluster_rows = F,cluster_cols = F)
##############################################
setwd('/home/xuchu/data/RNASeq/surfaceRNA/293T/comparecell')
data = read.csv('293T_detect3.xls',sep = '\t',header = F)
p = ggplot(data,aes(x=log(data$V3), y=log(data$V4),label=data$V5,col = data$V2#,size = data$V5
))+geom_point()+theme_bw()+themes_pop()+geom_text_repel()+ylim(0,10)+xlim(0,10)+
  theme(axis.text.x= element_text(size=15,colour = "black",angle = 0))+ylab('293T rep1 Log(TPM)')+xlab('293T rep2Log(TPM)')
p = p+scale_color_d3()+geom_hline(yintercept = log(10), linetype = "dashed")+ 
  geom_vline(xintercept = log(10), linetype = "dashed")+
  scale_color_manual(values = c('#1f77b4','#d62728','#ff7f0e','#2ca02c'))
p
ggsave('overlap_withcell4.pdf',width = 12,height = 10)

#############################################
library(pheatmap)
setwd('/home/xuchu/data/RNASeq/tRNA')
a = read.csv('human_celltype_avg_tRNA.xls',header = T,sep = '\t')
#colnames(a) = c('name','human B','human T','human Hspc','human NK','human Mono',
#                'mouse B','mouse T','mouse Hspc','mouse Nk','mouse Mono')
a = a[2:6]
pdf('both_tRNA_expression.pdf',width = 8,height=10)
pheatmap(log(data),show_rownames = F,show_colnames = T,
         color = colorRampPalette(colors = c("#51B187","white","#D96558"))(100))
dev.off()
############################################
data = read.csv('tRNA_scattor.xls',sep = '\t',header = F)
p = ggplot(data,aes(x=log(data$V2), y=log(data$V3),label=data$V5,col = data$V4#,size = data$V5
))+geom_point()+theme_bw()+themes_pop()+geom_text_repel()+ylim(0,10)+xlim(0,10)+
  theme(axis.text.x= element_text(size=15,colour = "black",angle = 0))+ylab('293T rep1 Log(TPM)')+xlab('293T rep2Log(TPM)')
p = p+scale_color_d3()+geom_hline(yintercept = log(10), linetype = "dashed")+ 
  geom_vline(xintercept = log(10), linetype = "dashed")+
  scale_color_manual(values = c('#1f77b4','#ff7f0e','#d62728','#2ca02c'))
p
ggsave('overlap_withcell_tRNA.pdf',width = 12,height = 10)
############################################
setwd('/home/xuchu/data/RNASeq/tRNA')
a = read.csv('mouse_celltype_avg_tRNA.xls',sep = '\t',row.names = 1,header =F)
entrpy = entropySpecificity(a)
a$entrpy = entrpy
write.table(a,'tRNA_entrpy_mouse.xls',sep = '\t')
#############################
a = read.csv('Human_Mouse_homologene_matrix.xls',header = T,sep = '\t')
human = a[3:32]
#human = human[which(rowSums(human)>10),]
#human = scale(human)
#human = t(human)
human_HSPC = rowMeans(data.frame(human$human.Hspc.rep1,human$human.Hspc.rep2,human$human.Hspc.rep3))
human_B = rowMeans(data.frame(human$human.B.rep1,human$human.B.rep2,human$human.B.rep3))
human_T = rowMeans(data.frame(human$human.T.rep1,human$human.T.rep2,human$human.T.rep3))
human_NK = rowMeans(data.frame(human$human.NK.rep1,human$human.NK.rep2,human$human.NK.rep3))
human_Mono = rowMeans(data.frame(human$human.Mono.rep1,human$human.Mono.rep2,human$human.Mono.rep3))

mouse_HSPC = rowMeans(data.frame(human$mouse.Hspc.rep1,human$mouse.Hspc.rep2,human$mouse.Hspc.rep3))
mouse_B = rowMeans(data.frame(human$mouse.B.rep1,human$mouse.B.rep2,human$mouse.B.rep3))
mouse_T = rowMeans(data.frame(human$mouse.T.rep1,human$mouse.T.rep2,human$mouse.T.rep3))
mouse_NK = rowMeans(data.frame(human$mouse.NK.rep1,human$mouse.NK.rep2,human$mouse.NK.rep3))
mouse_Mono = rowMeans(data.frame(human$mouse.Mono.rep1,human$mouse.Mono.rep2,human$mouse.Mono.rep3))

human2 = data.frame(human_HSPC,human_B,human_T,human_NK,human_Mono,mouse_HSPC,mouse_B,mouse_T,mouse_NK,mouse_Mono)
human2 = human2[which(rowSums(human2)>450),]
pdf('all.pdf',width=8,height=7)
pheatmap(log(human2+1),show_rownames = F,
  color = colorRampPalette(colors = c("#51B187","white","#D96558"))(100))
dev.off()

##############################################
library('eulerr')
library(RColorBrewer)
display.brewer.all()
########get colour
detail <- brewer.pal(n = 9, name = "Set1")
detail

completely_contained <- euler(c("AMOUR" = 9, "WGA" = 10, 'RT-AMP' = 2,
                                "AMOUR&RT-AMP" = 0,"WGA&RT-AMP" = 3,
                                "AMOUR&WGA" = 137,"AMOUR&WGA&RT-AMP"=233))
pdf('veen.pdf',width = 10,height = 10)
plot(completely_contained,
     quantities = list(type = c("counts")),
     labels = list(col = c("black", "black", "black")),
     edges = list(col = "white", lex = 1),
     fills = c('#87B5B2','#282A62','#2A5522','#A9C4E6','#E07E35','#F2CCA0','#B43970'))
#ggsave('veen.pdf',width = 8,height = 6)
dev.off()
###############
data = data.frame(group = c('AMOUR','WGA','RT-AMP'),tRNA = c(391,366,238))
data$Group = factor(data$group,levels = c('AMOUR','WGA','RT-AMP'))
p = ggplot(data,aes(x=Group,y=tRNA
                    ,fill = Group))+geom_bar(stat = 'identity')+theme_bw()+themes_pop()
p = p+scale_fill_d3()+ylab('tRNA Number')+xlab('Group')
p
ggsave('tRNA_number.pdf',width = 5,height = 8)
###############
f = read.table('mono_top.xls',header = T,sep ='\t')
data = f
data$name = factor(data$mono_mean,levels = data$mono_mean)
data
p = ggplot(data,aes(x=data$name,y=data$X
                    ,fill = name))+geom_bar(stat = 'identity')+theme_bw()+themes_pop()
p = p+scale_fill_manual(values = pal)+ylab('TPM')+xlab('Gene name')
p=p+theme(axis.text.x= element_text(size=15,colour = "black",angle = 90,vjust = 0.5))
p
ggsave('human_monocyte_TPM.pdf',width =8,height = 5)
