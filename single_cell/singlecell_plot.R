library(ggthemes)
library(ggrepel)
library(ggsci)
library(RColorBrewer)
library(wesanderson)
library(tidyverse)
library(scales)

setwd('/home/xuchu/data/jiangx_10X/result/annotation_bySeurat_PBMC/statcellnumber')
pal = wes_palette("Darjeeling1",8, type = "continuous")
a = read.csv('l2.type.percent.xls',sep = '\t',header = F)
data = data.frame(a)
data$CellType = data$V2
data$CellType[which(data$V5 < 0.01)] = 'Other'
data$CellType[which(data$V2 == 'HSPC')] = 'HSPC'
data$CellType[which(data$V2 == 'NK')] = 'NK'
data$CellType[which(data$V2 == 'B intermediate')] = 'Other'
data$CellType[which(data$V2 == 'Eryth')] = 'Other'
data$CellType = factor(data$CellType,levels = c('CD14 Mono',
                                          'CD16 Mono',
                                          'CD4 Naive',
                                          'CD4 TCM',
                                          'B naive',
                                          'CD8 Naive',
                                        #  'CD8 T',
                                          'cDC2',
                                          'NK',
                                          'HSPC',
                                          'Platelet',
                                          'Other'
                                          ))
p = ggplot(data,aes(x=data$V4,y=data$V5
                    ,fill = data$CellType))+geom_bar(stat = 'identity'#,position=position_dodge()
                    )+theme_bw()+themes_pop()+scale_fill_discrete(name = "hahah")
p = p+ylab('Proportion')+xlab('Sample name')+scale_fill_d3("category20")+
  scale_y_continuous(labels = scales::percent_format(scale = 100))#,limits=c(0, 1))
p
###############################
setwd('/home/xuchu/data/jiangx_10X/result/annotation_bySeurat_PBMC/statcellnumber')
pal = wes_palette("Darjeeling1",8, type = "continuous")
a = read.csv('l3.type.percent.xls',sep = '\t',header = F)
data = data.frame(a)
data$CellType = data$V2
data$CellType[which(data$V5 < 0.005)] = 'Other'
data$CellType[which(data$V2 == 'HSPC')] = 'HSPC'
data$CellType[which(data$V2 == 'NK')] = 'NK'
data$CellType[which(data$V2 == 'B intermediate')] = 'Other'
data$CellType[which(data$V2 == 'Eryth')] = 'Other'
data$CellType = factor(data$CellType,levels = c('CD14 Mono',
                                                'CD16 Mono',
                                                'CD4 Naive',
                                                'CD4 TCM',
                                                'B naive',
                                                'CD8 Naive',
                                                #  'CD8 T',
                                                'cDC2',
                                                'NK',
                                                'HSPC',
                                                'Platelet',
                                                'Other'
))
p = ggplot(data,aes(x=data$V4,y=data$V5
                    ,fill = data$CellType))+geom_bar(stat = 'identity'#,position=position_dodge()
                    )+theme_bw()+themes_pop()+scale_fill_discrete(name = "hahah")
p = p+ylab('Proportion')+xlab('Sample name')+scale_fill_d3("category20")+
  scale_y_continuous(labels = scales::percent_format(scale = 100))#,limits=c(0, 1))
p
#######################################3
setwd('/home/xuchu/data/jiangx_10X/result/annotation_bySeurat_PBMC/statcellnumber')
pal = wes_palette("Darjeeling1",8, type = "continuous")
a = read.csv('all.xls',sep = '\t',header = F)
data = data.frame(a)
data$CellType = data$V2
#data$CellType[which(data$V5 < 0.01)] = 'other'

data$CellType = factor(data$CellType,levels = c('Mono',
                                                'CD4 T',
                                                'B',
                                                'CD8 T',
                                                'DC',
                                                'NK',
                                                'other',
                                                'other T'
))
data$Group = factor(data$V4,levels = c('Y1','Y3','Y4','Y5','UCB'))
p = ggplot(data,aes(x=data$Group,y=data$V5
                    ,fill = data$CellType))+geom_bar(stat = 'identity'#,position=position_dodge()
                    )+theme_bw()+themes_pop()+scale_fill_discrete(name = "hahah")
p = p+ylab('Proportion')+xlab('Sample name')+scale_fill_d3("category20")+
  scale_y_continuous(labels = scales::percent_format(scale = 100))+labs(fill = 'Cell type')#,limits=c(0, 1))
p
