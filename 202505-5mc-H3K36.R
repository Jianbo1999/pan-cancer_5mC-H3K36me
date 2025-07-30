
#~/R/my_projects/202505-5mc-H3K36
if(T){
  rm(list = ls())
  setwd("~/R/my_projects/202505-5mc-H3K36")
  folder_path <- "./out"
  # 检查文件夹是否存在
  if (!dir.exists(folder_path)) {
    # 如果文件夹不存在，则创建它
    dir.create(folder_path)
    print(paste("Folder created at", folder_path))
  } else {
    print(paste("Folder already exists at", folder_path))
  }
  #c("#4766c2","#96dd73","#f9bc34","#e84545", "#51b1d6",'#59ab85','#fc7f4a')
  library(ggplot2);library(ggpubr)
  library(dplyr);library(ggplot2);library(viridis);library(pheatmap)
  colors <-c( '#D6E7A3', '#57C3F3', '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 

  #load("~/R/my_projects/202505-5mc-H3K36/202505-5mc-H3K36.RData")
}
#genes_5mc <-c("DNMT3A", "DNMT3B", "DNMT1")
#genes_H3K36me<-c("NSD1","NSD2","NSD3","ASH1L","SETMAR","SETD2")
genes_5mc <-c("DNMT3A", "DNMT3B", "DNMT1")
genes_H3K36me<-c("NSD1","WHSC1","WHSC1L1","ASH1L","SETMAR","SETD2")
#NSD2= WHSC1; NSD3 = "WHSC1L1";DNMT1= DNMT
genes <-c(genes_5mc,genes_H3K36me)#9-

#0. load data----
##setwd("D:\\研究生科研\\生物信息学\\phase\\data\\UCSC\\TCGA Pan-Cancer (PANCAN)22")
#load("~/R/ucsc_pancancer/data/panc_22_refined_fpkm_clinical.RData") #no
##load("~/R/ucsc_pancancer/data/UCSC.RData") #2.2G
exp<-panc_22#all samples=T+N #FPKM? TPM ?
# Normal  Tumor 
# 719  10295 
boxplot(panc_22[,1:10])
boxplot(log2(panc_22[,1:10]+1))#0-20
boxplot(scale(panc_22[,1:10]))#0-80
exp<-log2(exp+1)
boxplot(exp[,1:10])
###
genes_exp <-exp[rownames(exp) %in% genes,]#9-11069
#NSD2 NSD3="WHSC1L1" DNMT1=DNMT no match; change other name


for ( i in c("DNMT","Dnmt","dnmt","MCMT","ADCADN","HSN1E","HSN1E") ){
  print(rownames(exp)[grep(i,rownames(exp))])
} #no DNMT1/DNMT in panc_22_tumor; but all in panc_22!

#merge in meta
genes_exp <-as.data.frame(t(genes_exp))
genes_exp$sample <-rownames(genes_exp)

genes_exp_clinical <- merge(genes_exp,panc_22_clinical,by="sample")#11014-44

genes_exp_clinical_tumor <- genes_exp_clinical[genes_exp_clinical$type=="Tumor",]#10295
#1.0 expression----
##1.1 T vs N ?
#calculate  mean exp
# 使用dplyr计算平均值
library(dplyr)

result <-data.frame()
for (i in genes ){
  print(i)
  
  df <- genes_exp_clinical_tumor[,c(i,"cancer")] %>%  #tumor sample
    group_by(cancer) %>% summarise_each(funs(mean) )
    #summarise(mean_value = mean(genes_exp_clinical_tumor[,i], na.rm = TRUE))  # na.rm = TRUE表示忽略NA值
  #summarise(df,avg = mean(df[,i], na.rm = TRUE))  # na.rm = TRUE表示忽略NA值
  colnames(df)[2] <-"value"
  df$gene <-i
  result <-rbind(result,df)
} #

##dot plot
library(ggplot2);library(viridis)
ggplot()+
  #x和y代表位置坐标，z代表气泡的大小，group设置气泡的颜色类别，alpha设置透明度
  geom_point(data=result,aes(x=gene,y=cancer,color=cancer,#cancer value
                             size=value),alpha=0.8)+ 
  #设置气泡颜色
  scale_color_manual(values = colors )+
  #设置主题
  theme_bw()+
  #设置气泡大小范围
  #scale_size_continuous(range = c(1,5))+
  theme(#legend.position = 'none', #删除图例
        panel.grid.minor = element_blank(), #删除次grid
        panel.grid.major = element_line(linewidth = 0.4), #设置主grid
        axis.text = element_text(size=14), #坐标轴text设置
        axis.text.x  = element_text(angle = 90,hjust = 1),
        axis.title = element_text(size=16))+ #坐标轴lab设置
  labs(x='',y='')+ #x和y轴lab设置
  coord_flip()+ guides(color= "none")+labs(size="Exp")
ggsave(file='./out/dotpplot_exp.pdf',width = 8,height = 3,family="serif")
#
ggplot()+
  #x和y代表位置坐标，z代表气泡的大小，group设置气泡的颜色类别，alpha设置透明度
  geom_point(data=result,aes(x=gene,y=cancer,color=value,#cancer value
                             size=value),alpha=0.8)+ 
  #scale_color_manual(values = colors )+
  scale_color_viridis()+
  theme_bw()+
  theme(#legend.position = 'none', #删除图例
    panel.grid.minor = element_blank(), #删除次grid
    panel.grid.major = element_line(linewidth = 0.4), #设置主grid
    axis.text = element_text(size=14), #坐标轴text设置
    axis.text.x  = element_text(angle = 90,hjust = 1),
    axis.title = element_text(size=16))+ #坐标轴lab设置
  labs(x='',y='')+ #x和y轴lab设置
  coord_flip()#+ labs(size="Exp")#guides(color= "none")+
#ggsave(file='./out/dotpplot_exp-1.pdf',width = 8,height = 4,family="serif")

##热图
ggplot(data=result,aes(x=gene,y=cancer,color=value),alpha=0.8 ) +
  geom_tile(aes(fill = value), colour = "grey70") +
  # scale_fill_gradient(
  #   low = "white",
  #   high = list_temp$color,
  #   limits = c(0, 100)
  # ) +
  #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  scale_fill_gradient(low="skyblue",high = "pink")+
  theme_minimal() +
  #coord_fixed(ratio = 1) +
  labs(x = NULL, y = NULL, fill = "Exp") +
  theme(
    text = element_text(size = 12, colour = "black"),panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(file='./out/pheatmap_exp.pdf',width = 4,height = 5,family="serif")

##show number
library(reshape2)
result1 <- dcast(result,cancer~gene, value.var = "value")
#result1[2:ncol(result1)]<-lapply(result1[2:ncol(result1)],as.numeric())
result1<-as.data.frame(t(result1))
colnames(result1)<-result1[1,];result1<-result1[-1,]

result1[1:ncol(result1)]<-as.data.frame(lapply(result1[1:ncol(result1)],as.numeric))
#注释基因
result1$Type <-ifelse(rownames(result1) %in% genes_5mc,"5mC","H3K36me")

library(pheatmap)
#pheatmap(result1)
p1 <-pheatmap(result1[1:33],cluster_rows  = F, show_colnames = T,
              border_color = "white",#cellwidth=10,cellheight = 10,#gaps_col= 38, #table(exp_marker1$response)
              color=colorRampPalette( c("navy", "white", "firebrick3"))(10)#,
              #annotation_row= subset(result1, select = Type),annotation_names_row = F
)
p1
#保存热图
save_pheatmap_pdf <- function(x, filename, width, height) {
  library(grid)
  x$gtable$grobs[[1]]$gp <- gpar(lwd = 0.05)#修改聚类树线条宽度
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height,family = "serif")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}#保存热图的函数
save_pheatmap_pdf(p1, "./out/pheatmap_exp_2.pdf",6,3)

##
ann_colors = list(
  Type=c('5mC' ="#f9bc34",H3K36me="#96dd73")
)

p1 <-pheatmap(result1[1:33],cluster_rows  = F, show_colnames = T,
              border_color = "white",#cellwidth=10,cellheight = 10,#gaps_col= 38, #table(exp_marker1$response)
              color=colorRampPalette( c("navy", "white", "firebrick3"))(10),
              annotation_row= subset(result1, select = Type),annotation_names_row = F,
              annotation_colors = ann_colors
)
p1
save_pheatmap_pdf(p1, "./out/pheatmap_exp_3.pdf",6,3.2)
#+ number
p1 <-pheatmap(result1[1:33],cluster_rows  = F, show_colnames = T,
              border_color = "white",#cellwidth=10,cellheight = 10,#gaps_col= 38, #table(exp_marker1$response)
              color=colorRampPalette( c("skyblue", "white", "firebrick3"))(10),
              annotation_row= subset(result1, select = Type),annotation_names_row = F,
              annotation_colors = ann_colors,
              display_numbers=TRUE, number_format='%.1f'
)
p1
save_pheatmap_pdf(p1, "./out/pheatmap_exp_4.pdf",9,4)

##correlation基因相关性
library(ggcorrplot)
cormtcars <- round(cor(genes_exp_clinical_tumor1[2:10]), 3) #round()函数自定义小数点后位数
ggcorrplot(cormtcars,lab=T)
ggsave(file='./out/ggcorrplot.pdf',width = 5,height = 5,family="serif")

#默认是方形，修改为圆形显示，并标上相关系数
ggcorrplot(cormtcars,method = "circle",lab=T)
ggsave(file='./out/ggcorrplot-1.pdf',width = 5,height = 5,family="serif")

#使用ggcorrplot包的cor_pmat()函数计算p值：
pmtcars <- cor_pmat(genes_exp_clinical_tumor1[2:10])
ggcorrplot(cormtcars,hc.order = T,method = "circle",  #分等级聚类重排矩阵
           ggtheme = ggplot2::theme_void(base_size = 15), #主题修改
           colors = c("CornflowerBlue","white","Salmon"), #自定义颜色，看自己喜欢，或是参考好看的文献Figure用法。
           lab = T,lab_size = 4,    #相关系数文本字体大小
           tl.cex = 15,             #坐标轴字体大小
           p.mat = pmtcars,         #添加显著性信息
           sig.level = 0.01,        #显著性水平
           insig = "blank",pch = 4,                 #不够显著的色块进行标记，pch表示选择不同的标记方法，可以尝试其他数字表示什么标记方法
           pch.cex = 10)            #不显著标记的大小，使用insig = "blank"将不显著的空白处理
ggsave(file='./out/ggcorrplot-p.pdf',width = 5,height = 5,family="serif")

ggcorrplot(cormtcars,hc.order = T,method = "circle",  #分等级聚类重排矩阵
           ggtheme = ggplot2::theme_void(base_size = 15), #主题修改
           colors = c("CornflowerBlue","white","Salmon"), #自定义颜色，看自己喜欢，或是参考好看的文献Figure用法。
           lab = T,lab_size = 4,    #相关系数文本字体大小
           tl.cex = 15,             #坐标轴字体大小
           p.mat = pmtcars,         #添加显著性信息
           sig.level = 0.01,        #显著性水平
           pch = 4,                 #不够显著的色块进行标记，pch表示选择不同的标记方法，可以尝试其他数字表示什么标记方法
           pch.cex = 4)  
ggsave(file='./out/ggcorrplot-p-1.pdf',width = 5,height = 5,family="serif")

##显示一半
##corrplot
library(corrplot)
mycol2 <- colorRampPalette(c("skyblue","white", "pink"), alpha = TRUE)

pdf(file='./out/corrplot-p.pdf',width = 4,height = 4,family="serif")
corrplot(cormtcars, method = c('square'), type = c('lower'), 
         col = mycol2(100),
         outline = 'grey', #是否为图形添加外轮廓线，默认FLASE，可直接TRUE或指定颜色向量
         order = c('AOE'), #排序/聚类方式选择："original", "AOE", "FPC", "hclust", "alphabet"
         diag = FALSE, #是否展示对角线结果，默认TRUE
         tl.cex = 1.2, #文本标签大小
         tl.col = 'black', #文本标签颜色
         addgrid.col= 'grey' #格子轮廓颜色
)
dev.off()
#ggsave(file='./out/corrplot-p.pdf',width = 4,height = 4,family="serif")
# # 叠加相关性
# corrplot(cormtcars, method = c('ellipse'), 
#          type = c('lower'), 
#          col = mycol(100),
#          outline = 'grey',
#          order = c('AOE'),
#          diag = TRUE,
#          tl.cex = 1.2, 
#          tl.col = 'black',
#          addCoef.col = 'black', #在现有样式中添加相关性系数数字，并指定颜色
#          number.cex = 0.8, #相关性系数数字标签大小
# )
pdf(file='./out/corrplot-p-1.pdf',width = 4,height = 4,family="serif")
corrplot(cormtcars, method = c('square'), type = c('lower'), 
         col = mycol2(100),
         outline = 'grey', #是否为图形添加外轮廓线，默认FLASE，可直接TRUE或指定颜色向量
         order = c('AOE'), #排序/聚类方式选择："original", "AOE", "FPC", "hclust", "alphabet"
         diag = FALSE, #是否展示对角线结果，默认TRUE
         tl.cex = 1.2, #文本标签大小
         tl.col = 'black', #文本标签颜色
         addgrid.col= 'grey' #格子轮廓颜色
)
#下三角图添加不显著叉号：
corrplot(cormtcars, add = TRUE,
         #method = c('number'), 
         type = c('lower'),
         col = mycol(100),
         order = c('AOE'), 
         diag = FALSE, 
         number.cex = 0.9,#tl.cex=5,
         tl.pos = 'n', 
         cl.pos = 'n',
         p.mat = pmtcars,
         insig = "pch"
)
#ggsave(file='./out/corrplot-p-1.pdf',width = 4,height = 4,family="serif")
dev.off()

#2.0 correlation----
##2.1 score ----
#exp
##load("~/R/ucsc_pancancer/data/UCSC.RData") #2.2G
exp<-panc_22#

#geneset
gene_set<-subset(result1, select = Type)
gene_set$gene <-rownames(gene_set)
list<- split(as.matrix(gene_set)[,2], gene_set[,1])

library(GSVA)
#gsva_matrix<- gsva(as.matrix(exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

##First we should build a parameter object for the desired methodology.R
gsvaPar <- ssgseaParam(exprData = as.matrix(exp), 
                       geneSets = list,
                       normalize = TRUE)
##Second, we call the gsva() function with the parameter object as first argument. 
gsva_matrix <- gsva(gsvaPar, verbose = FALSE)
gsva_matrix<-as.data.frame(gsva_matrix);gsva_matrix<-as.data.frame(t(gsva_matrix ))
gsva_matrix$sample <-rownames(gsva_matrix)
#merge in meta
genes_exp_clinical_tumor1 <-merge(genes_exp_clinical_tumor,gsva_matrix,by="sample")
#"5mC"     "H3K36me"

##2.2 score exp ----
#平均值
library(dplyr)
result2 <-data.frame()
for (i in names(list)){
  print(i)
  
  df <- genes_exp_clinical_tumor1[,c(i,"cancer")] %>%  #tumor sample
    group_by(cancer) %>% summarise_each(funs(mean) )
  #summarise(mean_value = mean(genes_exp_clinical_tumor[,i], na.rm = TRUE))  # na.rm = TRUE表示忽略NA值
  #summarise(df,avg = mean(df[,i], na.rm = TRUE))  # na.rm = TRUE表示忽略NA值
  colnames(df)[2] <-"value"
  df$gene <-i
  result2 <-rbind(result2,df)
} #

#plot
ggplot()+
  #x和y代表位置坐标，z代表气泡的大小，group设置气泡的颜色类别，alpha设置透明度
  geom_point(data=result2,aes(x=gene,y=cancer,color=value,#cancer value
                             size=value),alpha=0.8)+ 
  #scale_color_manual(values = colors )+
  scale_color_viridis()+
  theme_bw()+
  theme(#legend.position = 'none', #删除图例
    panel.grid.minor = element_blank(), #删除次grid
    panel.grid.major = element_line(linewidth = 0.4), #设置主grid
    axis.text = element_text(size=14), #坐标轴text设置
    axis.text.x  = element_text(angle = 90,hjust = 1),
    axis.title = element_text(size=16))+ #坐标轴lab设置
  labs(x='',y='')#+ #x和y轴lab设置
  #coord_flip()#+ labs(size="Exp")#guides(color= "none")+
ggsave(file='./out/dotpplot_score-1.pdf',width = 2.5,height = 6.2,family="serif")

##热图
library(reshape2)
result3 <- dcast(result2,cancer~gene, value.var = "value")
#result1[2:ncol(result1)]<-lapply(result1[2:ncol(result1)],as.numeric())
result3<-as.data.frame(t(result3))
colnames(result3)<-result3[1,];result3<-result3[-1,]
result3[1:ncol(result3)]<-as.data.frame(lapply(result3[1:ncol(result3)],as.numeric))

p1 <-pheatmap(result3,cluster_rows  = F, show_colnames = T,
              border_color = "white",#cellwidth=10,cellheight = 10,#gaps_col= 38, #table(exp_marker1$response)
              gaps_row = 1,
              color=colorRampPalette( c("skyblue", "white", "firebrick3"))(10)#,
              #annotation_row= subset(result1, select = Type),annotation_names_row = F#,
              #annotation_colors = ann_colors,
              #display_numbers=TRUE, number_format='%.1f'
)
p1
save_pheatmap_pdf(p1, "./out/pheatmap_score.pdf",8,2)
#show number
p1 <-pheatmap(result3,cluster_rows  = F, show_colnames = T,#cluster_cols = F,
              border_color = "white",#cellwidth=10,cellheight = 10,#gaps_col= 38, #table(exp_marker1$response)
              gaps_row = 1,
              color=colorRampPalette( c("skyblue", "white", "firebrick3"))(10),
              #annotation_row= subset(result1, select = Type),annotation_names_row = F#,
              #annotation_colors = ann_colors,
              display_numbers=TRUE, number_format='%.1f'
)
p1
save_pheatmap_pdf(p1, "./out/pheatmap_score-2.pdf",8,2)

##2.3 score cor----
#Custom.color <-c("skyblue","pink")
df_cor <-as.data.frame(t(result3));df_cor$cancer <-rownames(df_cor)

library(ggplot2);library(ggpubr);library(ggsci)

ggplot(df_cor, aes(x=`5mC`, y=H3K36me)) + 
  #xlim(-20,15) + ylim(-15,10) +
  labs(x = "5mC score", y = "H3K36me score") +
  geom_point(aes(color=cancer),size = 3) +
  geom_smooth(method ='lm', size=0.5,color="skyblue") +
  stat_cor(method = "spearman",size = 6,color="red") +
  scale_colour_manual(values = colors ) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave(file='./out/cor_score.pdf',width = 6.5,height = 4.5,family="serif")

#no legend + text
ggplot(df_cor, aes(x=`5mC`, y=H3K36me)) + 
  #xlim(-20,15) + ylim(-15,10) +
  labs(x = "5mC score", y = "H3K36me score") +
  geom_point(color="navy",size = 3) + #aes(color=cancer)
  geom_smooth(method ='lm', size=0.5,color="skyblue") +
  ggrepel::geom_text_repel(aes(label=cancer))+ #geom_label_repel
  stat_cor(method = "spearman",size = 6,color="red") +
  scale_colour_manual(values = colors ) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave(file='./out/cor_score-1.pdf',width = 5,height = 4.5,family="serif")

##cor in cancer  各癌症类型中样本相关性
#genes_exp_clinical_tumor1
for (i in unique(genes_exp_clinical_tumor1$cancer) ){
  print(i)
  df <-genes_exp_clinical_tumor1[genes_exp_clinical_tumor1$cancer == i,]
  ggplot(df, aes(x=`5mC`, y=H3K36me)) + 
    #xlim(-20,15) + ylim(-15,10) +
    labs(x = "5mC score", y = "H3K36me score",title = i) +
    geom_point(color="navy",size = 1.5,alpha=0.5) + #aes(color=cancer)
    geom_smooth(method ='lm', size=0.5,color="skyblue") +
    #ggrepel::geom_text_repel(aes(label=cancer))+ #geom_label_repel
    stat_cor(method = "spearman",size = 6,color="red") +
    scale_colour_manual(values = colors ) +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 14), 
          axis.text.y = element_text(size = 14), 
          axis.title.x = element_text(size = 16), 
          axis.title.y = element_text(size = 16),panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 12))
  ggsave(file=paste0('./out/cor_score-',i,'.pdf'),width = 3.3,height = 3.3,family="serif")
  
}



#3.0 prognosis----

###3.1 COX Univariate Cox hazard analysis in all cancer---
library(survival)
data2<-genes_exp_clinical_tumor1
colnames(data2)[c(36,35)]<-c("futime","fustat")

input_gene= genes#c(genes,"H3K36me","5mC" )

uni_cox_fun <-function(data2,input_gene){
  outTab=data.frame()
  sigGenes=c("futime","fustat")
  for(i in input_gene ){ #trainset[,3:ncol(trainset)]
    cox <- coxph(Surv(futime, fustat) ~ data2[,i], data = data2)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    if(T){ #coxP<0.2
      sigGenes=c(sigGenes,i)
      outTab=rbind(outTab,
                   cbind(id=i,
                         HR= coxSummary$conf.int[,"exp(coef)"] ,
                         HR.95L=coxSummary$conf.int[,"lower .95"],
                         HR.95H=coxSummary$conf.int[,"upper .95"],
                         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
      )
    }
  }
  colnames(outTab)#;str(outTab)
  outTab[2:5]<-apply(outTab[2:5],2,as.numeric)
  outTab[2:4]<-round(outTab[2:4],2)
  outTab$Sign. <-ifelse(outTab$pvalue<0.001,"***",ifelse(outTab$pvalue<0.01,"**",ifelse(outTab$pvalue<0.05,"*","ns")))
  outTab$pvalue<-round(outTab$pvalue,3)
  # outTab<-dplyr:: arrange(outTab,desc(HR) )#HR排序 ?
  # outTab$id <-as.factor(outTab$id)
  outTab$color <-ifelse(outTab$HR < 1,"skyblue","pink")
  
  ##plot
  p <-ggplot(outTab, aes(HR, reorder(id,-HR) ))+
    geom_errorbarh(aes(xmax =HR.95H, xmin = HR.95L,color="black"),size= 0.5,height = 0.1) + #color=id
    #geom_point(aes(color=id),size=2.5)+
    geom_point(aes(size=HR,color=color)  )+ #,color="skyblue" id
    #scale_x_continuous(limits= c(-2.5, 2.5))+
    geom_vline(aes(xintercept = 1),color="gray",linetype="dashed", size = 0.5) +
    geom_text(size=4, aes(x = HR+0.1, label = Sign.),color="navy")+#,hjust = 1
    xlab('HR')+ ylab(' ')+
    theme_bw()+
    scale_color_manual(values = c("black","pink","skyblue") )+ #colors
    theme(axis.text.x = element_text(size = 14, color = "black"),axis.text.y = element_text(size = 14, color = "black"),
          panel.grid.minor = element_blank(),panel.grid.major  = element_blank(),
    )+
    theme(title=element_text(size=14),legend.position = "none")
  print(p)
}

uni_cox_fun(data2,genes)
#ggsave("./out/dot_uniCOX_HR_genes.pdf",width = 3.5,height = 3,family="serif")
uni_cox_fun(data2,c("H3K36me","5mC") )
#ggsave("./out/dot_uniCOX_HR_score.pdf",width = 3.5,height = 1.5,family="serif")

##genes km----
#all genes in all cancers
library(survminer)
for (j in c(genes,"H3K36me","5mC") ){
  print(j)
  data2["riskscore"] <-data2[j] 
  #
  dat <-data2
  cut_value <- surv_cutpoint(dat, #数据集
                             time = "futime", #生存状态
                             event = "fustat", #生存时间
                             variables = c("riskscore") )$cutpoint[,1] #best cutoff
  dat$Risk <- ifelse(dat$riskscore> cut_value,# median(dat$riskscore),#cut_value,
                     "High","Low") # 均值？
  dat$futime<-dat$futime/365
  
  fit <- survfit(Surv(futime, fustat) ~ Risk, data =dat)
  #pdf(paste0("./out/km_",i,"_",j,".pdf"),width=3,height=3.5,family = "serif")
  print(ggsurvplot_list(fit,data = dat,pval = T,pval.method = TRUE,##是否添加P值
                        conf.int = F,### 是否添加置信区间
                        legend.title = j,#title = i,#"Risk", # 设置图例标题
                        legend.labs = c("High", "Low"), #legend = c(0.8, 0.2),# 指定图例分组标签
                        risk.table = F, # 是否添加风险表
                        risk.table.col = "strata", 
                        censor=F,size=0.5,
                        surv.scale="percent",
                        ###linetype = "strata",
                        #surv.median.line = "hv", # 是否添加中位生存线
                        risk.table.y.text.col = F,risk.table.y.text = FALSE,
                        ggtheme = theme_bw(base_family = "serif")#+theme(legend.text = element_text(colour = c("red", "blue")))
                        +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
                        +theme(plot.title = element_text(hjust = 0.5,size = 16),axis.title.y.left = element_text(size = 16,vjust = 1),axis.title.x.bottom = element_text(size = 16,vjust = 0))
                        +theme(axis.text.x.bottom = element_text(size = 12,vjust = -0.8,colour = "black")) #,face = "bold"
                        +theme(axis.text.y.left = element_text(size = 12,vjust = 0.5,hjust = 0.5,angle = 90,colour = "black"))
                        +theme(legend.title = element_text(family = "Times",colour = "black",size = 14))
                        +theme(legend.text = element_text(family = "Times",colour = "black",size =12)), # Change ggplot2 theme
                        palette = c("#bc1e5d", "#0176bd"),#"lacent",#c("#bc1e5d", "#0176bd"),##如果数据分为两组，就写两种颜色，三组就写三种颜色，具体的颜色搭配参考上一期发布的R语言颜色对比表
                        #xlim=c(10,5000),##x轴的阈值，根据随访数据的天数进行更改
                        xlab = "Years")##随访的时间时天，就写Days，月份就写Months
        
        
  )
  #dev.off()
  ggsave(paste0("./out/km_all_cancer_",j,".pdf"),width=3,height=3.5,family = "serif")
  
}

#all genes in single cancer 
for (j in c(genes,"H3K36me","5mC") ){
  print(j)
  data2["riskscore"] <-data2[j] 
  #
  for (i in unique(data2$cancer)){
    print(i)
    dat <-data2[data2$cancer==i,]
    cut_value <- surv_cutpoint(dat, #数据集
                               time = "futime", #生存状态
                               event = "fustat", #生存时间
                               variables = c("riskscore") )$cutpoint[,1] #best cutoff
    dat$Risk <- ifelse(dat$riskscore> cut_value,# median(dat$riskscore),#cut_value,
                       "High","Low") # 均值？
    dat$futime<-dat$futime/365
    fit <- survfit(Surv(futime, fustat) ~ Risk, data =dat)
    #pdf(paste0("./out/km_",i,"_",j,".pdf"),width=3,height=3.5,family = "serif")
    print(ggsurvplot_list(fit,data = dat,pval = T,pval.method = TRUE,##是否添加P值
                          conf.int = F,### 是否添加置信区间
                          legend.title = j,#title = i,#"Risk", # 设置图例标题
                          legend.labs = c("High", "Low"), #legend = c(0.8, 0.2),# 指定图例分组标签
                          risk.table = F, # 是否添加风险表
                          risk.table.col = "strata", 
                          censor=F,size=0.5,
                          surv.scale="percent",
                          ###linetype = "strata",
                          #surv.median.line = "hv", # 是否添加中位生存线
                          risk.table.y.text.col = F,risk.table.y.text = FALSE,
                          ggtheme = theme_bw(base_family = "serif")#+theme(legend.text = element_text(colour = c("red", "blue")))
                          +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
                          +theme(plot.title = element_text(hjust = 0.5,size = 16),axis.title.y.left = element_text(size = 16,vjust = 1),axis.title.x.bottom = element_text(size = 16,vjust = 0))
                          +theme(axis.text.x.bottom = element_text(size = 12,vjust = -0.8,colour = "black")) #,face = "bold"
                          +theme(axis.text.y.left = element_text(size = 12,vjust = 0.5,hjust = 0.5,angle = 90,colour = "black"))
                          +theme(legend.title = element_text(family = "Times",colour = "black",size = 14))
                          +theme(legend.text = element_text(family = "Times",colour = "black",size =12)), # Change ggplot2 theme
                          palette = c("#bc1e5d", "#0176bd"),#"lacent",#c("#bc1e5d", "#0176bd"),##如果数据分为两组，就写两种颜色，三组就写三种颜色，具体的颜色搭配参考上一期发布的R语言颜色对比表
                          #xlim=c(10,5000),##x轴的阈值，根据随访数据的天数进行更改
                          xlab = "Years")##随访的时间时天，就写Days，月份就写Months
          
          
    )
    #dev.off()
    ggsave(paste0("./out/km_all_cancer_",i,"_",j,".pdf"),width=3,height=3.5,family = "serif")
    
  }
}

###3.2 bacth COX Univariate Cox hazard analysis in all single cancer---
#type<- c()
#out_tab<-data.frame()#(HR="",L95="",R95="",pval="",cancer="")
library(survminer)
for (j in c("H3K36me","5mC") ){
  print(j)
  data2["riskscore"] <-data2[j] 
  ##
  out_tab<-data.frame()
  for (i in  unique(data2$cancer) ) {
    print(i)
    dat <- data2[data2$cancer==i,]
    #dat <- na.omit(dat)
    #plot
    cut_value <- surv_cutpoint(dat, #数据集
                               time = "futime", #生存状态
                               event = "fustat", #生存时间
                               variables = c("riskscore") )$cutpoint[,1] #best cutoff
    dat$Risk <- ifelse(dat$riskscore> cut_value,# median(dat$riskscore),#cut_value,
                       "High","Low") # 均值？
    dat$futime<-dat$futime/365
    
    fit <- survfit(Surv(futime, fustat) ~ Risk, data =dat)
    pval <- survdiff(Surv(futime, fustat) ~ Risk, rho = 0, data =dat)
    pval <- pval$pvalue
    hr <- summary(coxph(Surv(futime, fustat) ~ Risk, data =dat))
    HR <-hr$conf.int[,"exp(coef)"]
    L95 <-hr$conf.int[,"lower .95"]
    R95 <-hr$conf.int[,"upper .95"]
    out_tab<-rbind(out_tab,
                   data.frame(HR=HR,L95=L95,R95=R95,pval=pval,cancer=i)
    )
    #plot km start
    pdf(paste0("./out/km_",i,"_",j,".pdf"),width=3,height=3.5,family = "serif")
    print(ggsurvplot_list(fit,data = dat,pval = T,pval.method = TRUE,##是否添加P值
                          conf.int = F,### 是否添加置信区间
                          legend.title = j,#title = i,#"Risk", # 设置图例标题
                          legend.labs = c("High", "Low"), #legend = c(0.8, 0.2),# 指定图例分组标签
                          risk.table = F, # 是否添加风险表
                          risk.table.col = "strata", 
                          censor=F,size=0.5,
                          surv.scale="percent",
                          ###linetype = "strata",
                          #surv.median.line = "hv", # 是否添加中位生存线
                          risk.table.y.text.col = F,risk.table.y.text = FALSE,
                          ggtheme = theme_bw(base_family = "serif")#+theme(legend.text = element_text(colour = c("red", "blue")))
                          +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
                          +theme(plot.title = element_text(hjust = 0.5,size = 16),axis.title.y.left = element_text(size = 16,vjust = 1),axis.title.x.bottom = element_text(size = 16,vjust = 0))
                          +theme(axis.text.x.bottom = element_text(size = 12,vjust = -0.8,colour = "black")) #,face = "bold"
                          +theme(axis.text.y.left = element_text(size = 12,vjust = 0.5,hjust = 0.5,angle = 90,colour = "black"))
                          +theme(legend.title = element_text(family = "Times",colour = "black",size = 14))
                          +theme(legend.text = element_text(family = "Times",colour = "black",size =12)), # Change ggplot2 theme
                          palette = c("#bc1e5d", "#0176bd"),#"lacent",#c("#bc1e5d", "#0176bd"),##如果数据分为两组，就写两种颜色，三组就写三种颜色，具体的颜色搭配参考上一期发布的R语言颜色对比表
                          #xlim=c(10,5000),##x轴的阈值，根据随访数据的天数进行更改
                          xlab = "Years")##随访的时间时天，就写Days，月份就写Months
          
          
    )
    dev.off()
    #plot km end 
  }
  if(T){
    #plot HR dot
    #plot https://github.com/yuyang3/pan-B/blob/main/Figure5_survival.R
    plot_df <-out_tab
    plot_df$R95[28]<-plot_df$L95[28]
    plot_df[1:3] <-log2(plot_df[1:3]+1)
    
    plot_df$color[plot_df$pval>= 0.05 & plot_df$HR < 1] <- "Better survival (P > 0.05)"
    plot_df$color[plot_df$pval< 0.05 & plot_df$HR < 1] <- "Better survival"
    plot_df$color[plot_df$pval >= 0.05 & plot_df$HR > 1] <- "Worse survival (P > 0.05)"
    plot_df$color[plot_df$pval < 0.05 & plot_df$HR > 1] <- "Worse survival"
    
    plot_df$color <- factor(plot_df$color, levels = c(
      "Better survival", "Better survival (P > 0.05)",
      "Worse survival (P > 0.05)", "Worse survival"
    ))
    mycolor <- c(
      "Better survival" = "#0F7B9F",
      "Better survival (P > 0.05)" = "#E0F3F8FF",
      "Worse survival (P > 0.05)" = "#FDDBC7FF",
      "Worse survival" = "#C3423F"
    )
    plot_df$significance <- ""
    plot_df$significance[plot_df$pval <= 0.05] <- "*"
    plot_df$significance[plot_df$pval <= 0.01] <- "**"
    plot_df$significance[plot_df$pval <= 0.001] <- "***"
    #plot_df$significance[plot_df$pval <= 0.0001] <- "#"
    
    #plot
    ggplot(data = plot_df, aes(x = reorder(cancer,-HR), y = HR, ymin = L95, ymax = R95, color = color)) +
      geom_pointrange() +
      geom_hline(yintercept = log2(1+1), colour = "grey40", linetype = "dashed", size = 0.2) + # add a dotted line at x=1 after flip
      geom_text(aes(x = cancer, y = max(R95), label = significance), show.legend = FALSE) +
      scale_color_manual(values = mycolor, name = "Survival association") +
      xlab("") +
      ylab("Hazard ratio (95% CI)") +
      theme_classic()+#cowplot::theme_cowplot() +
      theme(
        #text = element_text(size = 7, family = "ArialMT"),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
      ) #+ggtitle(paste0(j," score"))
     #ggsave(paste0("./out/HR_cancers_",j,"_score.pdf"),width = 9,height = 3,family="serif")
    #
  }
  ##
}

###3.3 m5c + H3K36me: double-group survival----
#A.batch group----
colors_4<-c("#0072B5","#20854E","#F37C95","#E18727")
dat <-data2
##all samples
for (i in c("H3K36me","5mC") ){
  print(i)
  group_name <- paste0(i,"_group")
   
  cut_value <- surv_cutpoint(dat, #数据集
                             time = "futime", #生存状态
                             event = "fustat", #生存时间
                             variables = i )$cutpoint[,1]  #c("riskscore") #best cut off
  dat[,group_name] <- ifelse(dat[,i]>cut_value ,#median(data2[,i]),
                               paste0(i,"_High"),paste0(i,"_Low") )
  dat$Group <- paste0(dat$`5mC_group`," + ",dat$H3K36me_group)
  ##plot all samples km
  fit <- survfit(Surv(futime/365, fustat) ~ Group, data =dat)
  pdf(paste0("./out/km_tow_group_",i,".pdf"),width=6,height=3.5,family = "serif")
  print(ggsurvplot_list(fit,data = dat,pval = T,pval.method = TRUE,##是否添加P值
                        conf.int = F,### 是否添加置信区间
                        legend.title = "",#j,#title = i,#"Risk", # 设置图例标题
                        legend.labs = unique(dat$Group),#c("High", "Low"), 
                        legend ="right",#legend = c(0.5, 0.8),# 指定图例分组标签
                        risk.table = F, # 是否添加风险表
                        risk.table.col = "strata", 
                        censor=F,size=0.5,
                        surv.scale="percent",
                        ###linetype = "strata",
                        #surv.median.line = "hv", # 是否添加中位生存线
                        risk.table.y.text.col = F,risk.table.y.text = FALSE,
                        ggtheme = theme_bw(base_family = "serif")#+theme(legend.text = element_text(colour = c("red", "blue")))
                        +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
                        +theme(plot.title = element_text(hjust = 0.5,size = 16),axis.title.y.left = element_text(size = 16,vjust = 1),axis.title.x.bottom = element_text(size = 16,vjust = 0))
                        +theme(axis.text.x.bottom = element_text(size = 12,vjust = -0.8,colour = "black")) #,face = "bold"
                        +theme(axis.text.y.left = element_text(size = 12,vjust = 0.5,hjust = 0.5,angle = 90,colour = "black"))
                        +theme(legend.title = element_text(family = "Times",colour = "black",size = 14))
                        +theme(legend.text = element_text(family = "Times",colour = "black",size =12)), # Change ggplot2 theme
                        palette = colors_4,#c("#bc1e5d", "#0176bd"),#"lacent",#c("#bc1e5d", "#0176bd"),##如果数据分为两组，就写两种颜色，三组就写三种颜色，具体的颜色搭配参考上一期发布的R语言颜色对比表
                        #xlim=c(10,5000),##x轴的阈值，根据随访数据的天数进行更改
                        xlab = "Years")##随访的时间时天，就写Days，月份就写Months
        
        
  )
  dev.off()
  
}

##single cancer
for (j in unique(data2$cancer)){
  print(j)
  dat <-data2[data2$cancer==j,]
  ##
  for (i in c("H3K36me","5mC") ){
    print(i)
    group_name <- paste0(i,"_group")
    
    cut_value <- surv_cutpoint(dat, #数据集
                               time = "futime", #生存状态
                               event = "fustat", #生存时间
                               variables = i )$cutpoint[,1]  #c("riskscore") #best cut off
    dat[,group_name] <- ifelse(dat[,i]>cut_value ,#median(data2[,i]),
                               paste0(i,"_High"),paste0(i,"_Low") )
    dat$Group <- paste0(dat$`5mC_group`," + ",dat$H3K36me_group)
    print(table(dat$Group))
    
  }
  
  ##plot all samples km
  fit <- survfit(Surv(futime/365, fustat) ~ Group, data =dat)
  pdf(paste0("./out/km_tow_group_",j,"_.pdf"),width=6,height=3.5,family = "serif")
  print(ggsurvplot_list(fit,data = dat,pval = T,pval.method = TRUE,##是否添加P值
                        conf.int = F,### 是否添加置信区间
                        legend.title = j,#title = i,#"Risk", # 设置图例标题
                        legend.labs = unique(dat$Group),#c("High", "Low"), 
                        legend ="right",#legend = c(0.5, 0.8),# 指定图例分组标签
                        risk.table = F, # 是否添加风险表
                        risk.table.col = "strata", 
                        censor=F,size=0.5,
                        surv.scale="percent",
                        ###linetype = "strata",
                        #surv.median.line = "hv", # 是否添加中位生存线
                        risk.table.y.text.col = F,risk.table.y.text = FALSE,
                        ggtheme = theme_bw(base_family = "serif")#+theme(legend.text = element_text(colour = c("red", "blue")))
                        +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
                        +theme(plot.title = element_text(hjust = 0.5,size = 16),axis.title.y.left = element_text(size = 16,vjust = 1),axis.title.x.bottom = element_text(size = 16,vjust = 0))
                        +theme(axis.text.x.bottom = element_text(size = 12,vjust = -0.8,colour = "black")) #,face = "bold"
                        +theme(axis.text.y.left = element_text(size = 12,vjust = 0.5,hjust = 0.5,angle = 90,colour = "black"))
                        +theme(legend.title = element_text(family = "Times",colour = "black",size = 14))
                        +theme(legend.text = element_text(family = "Times",colour = "black",size =12)), # Change ggplot2 theme
                        palette = colors_4,#c("#bc1e5d", "#0176bd"),#"lacent",#c("#bc1e5d", "#0176bd"),##如果数据分为两组，就写两种颜色，三组就写三种颜色，具体的颜色搭配参考上一期发布的R语言颜色对比表
                        #xlim=c(10,5000),##x轴的阈值，根据随访数据的天数进行更改
                        xlab = "Years")##随访的时间时天，就写Days，月份就写Months
        
        
  )
  dev.off()
}

#10295 tumor
# ACC BLCA BRCA CESC CHOL COAD DLBC ESCA  GBM HNSC KICH KIRC KIRP LAML  LGG LIHC LUAD 
# 79  408 1102  306   36  451   48  185  166  522   66  534  291  173  529  373  517 
# LUSC MESO   OV PAAD PCPG PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC  UCS  UVM 
# 501   87  308  179  184  498  160  263  472  415  139  513  120  533   57   80 

##normal 
table(panc_22_clinical$type)
# Normal  Tumor 
# 1413  11178 
count_data <-panc_22_clinical[panc_22_clinical$type=="Normal",] %>%
   group_by(cancer) %>% #group_by(type) %>%
  summarise(count = n())
print(count_data,n=nrow(count_data) )
# cancer count
# <chr>  <int>
#   1 BLCA      23
# 2 BRCA     132
# 3 CESC       3
# 4 CHOL       9
# 5 COAD      85
# 6 ESCA      18
# 7 GBM        1
# 8 HNSC      74
# 9 KICH      25
# 10 KIRC     407
# 11 KIRP      60
# 12 LIHC      59
# 13 LUAD     120
# 14 LUSC     119
# 15 OV         4
# 16 PAAD      10
# 17 PCPG       3
# 18 PRAD      67
# 19 READ      16
# 20 SARC       6
# 21 SKCM       2
# 22 STAD      68
# 23 THCA      65
# 24 THYM       2
# 25 UCEC      35