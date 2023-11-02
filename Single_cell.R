############ Star analysis ##############
library(Seurat)
library(dplyr)
library(patchwork)
library(SeuratDisk)
library(DoubletFinder)
library(pals)
library(cowplot)
library(ggplot2)
library(SeuratObject)
set.seed(9999)
setwd('F:/analysis/Cecum_10X+C4')

#color scheme
library(RColorBrewer)
zzcolors <- c("#8DD3C7", "#F4F4B9", "#CFCCCF", "#D1A7B9", "#F4867C", "#C0979F", "#86B1CD", "#CEB28B", 
              "#EDBC63", "#C2D567", "#F8CDDE", "#E9D3DE", "#D5CFD6", "#C59CC5","#CFECBB", "#CDD796",
              "#E28BC3", "#D2A29F", "#BABF77", "#AAD852", "#CBD844", "#ECD836", "#FAD53E", "#F1CD64")
zz2colors <- c("#66C2A5", "#9DAE8C", "#D49A73", "#F08F6D", "#C79693", "#9E9DBA", "#9F9BC9", "#C193C6")
webcolors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175')














############ color scheme####################
#color <- readRDS('F:/analysis/Cecum_10X+C4/intestine_color.rds')
library(hash)
#Level1 Cell Color Scheme
h <- hash()
h = hash(keys=c('Epithelial','T/NK cell','Fibroblast','B cell','Myeloid','SMC','LEC','Endothelial','Terminal glial cell','Mesothelial cell'),values=c('#FF4A46','#FF8C00','#FFFF00','#06D6A0','#1CE6FF','#FFE4E1','#9395E7','#B1CE46','#F8AC8C','#FF90C9'))
saveRDS(h,'F:/analysis/Cecum_10X+C4/color/Cecum_color.rds')
#Level2 Cell Color Scheme
h <- hash()
h = hash(keys=c('Enterocyte','Goblet','Tuft','Enteroendocrine','NKT','Naive Cd4 T',
                'Cd4 T','NK','Fibroblast(Gsn high)','Fibroblast(Adamdec1 high)','B cell',
                'Iga_Plasma','B_Cycling','Macrophage','Mast cell',
                'SMC','LEC','Endothelial','Terminal glial cell','Mesothelial cell'),
         values=c('#FF8C00','#FFFF00','#0000A6','#9395E7','#B1CE46','#F8AC8C',
                '#63E398','#2F7FC1','#FF4A46','#1CE6FF','#C4A5DE',
                '#06D6A0','#0F4C5C','#F1C395','#BC6C25',
                '#808000','#71C33A','#F5F5DC','#927D85','#A30059'))
saveRDS(h,'F:/analysis/Cecum_10X+C4/color/Cecum_color_Level2.rds')
#Immune cell color scheme
h <- hash()
h = hash(keys=c('NKT','Naive Cd4 T','Cd4 T','NK','B cell','Iga_Plasma','B_Cycling','Macrophage','Mast cell'),
         values=c('#B1CE46','#F8AC8C','#63E398','#2F7FC1','#C4A5DE','#06D6A0',
                  '#0F4C5C','#F1C395','#BC6C25'))
#Epithelial cell subclass color scheme
Ep <- hash()
Ep = hash(keys=c('EC(Saa1 high)','EC(Selenbp1 high)','EC(Hmgb2 high)','Goblet','mLTo','EC(Rpl3 high)','Tuft','EEC','EC(Dying)','Ep_Fibroblast','EC(Reg3b high)'),values=c('#FF4A46','#FF8C00','#FFFF00','#06D6A0','#1CE6FF','#E14AEC','#FF90C9','#FFE4E1','#9395E7','#B1CE46','#F8AC8C'))
saveRDS(Ep,'Ep_color.rds')
#T/NK cell subclass color scheme
T <- hash()
T = hash(keys=c('Immature NKT','Naive Cd4 T','T cell','Treg','T_Dying','Mature NKT','Th17','Mature B'),values=c('#FF4A46','#FF8C00','#FFFF00','#06D6A0','#1CE6FF','#E14AEC','#FF90C9','#FFE4E1'))
saveRDS(T,'T_color.rds')
#Fibroblast color scheme
h <- hash()
h = hash(keys=c('Fibroblast(Gsn high)','Fibroblast(Adamdec1 high)','Mesothelial cell'),values=c('#FF4A46','#1CE6FF','#A30059'))
############ data integration ##########
### START ANALYSING
library(Seurat)
library(dplyr)
library(patchwork)
library(SeuratDisk)
library(DoubletFinder)
library(pals)
library(cowplot)
library(ggplot2)
set.seed(9999)
setwd('D:/analysis')

Find_doublet <- function(data){
  sweep.res.list <- paramSweep_v3(data, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  nExp_poi <- round(0.05*ncol(data))
  p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  data <- doubletFinder_v3(data, PCs = 1:10, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  data
}

bm_GF.data <- Read10X(data.dir = 'D:/analysis/GF_HJ/')
bm_SPF.data <- Read10X(data.dir = 'D:/analysis/SPF_HJ/')
bm_GFM.data <- Read10X(data.dir = 'D:/analysis/GF_M/')
bm_SPFM.data <- Read10X(data.dir = 'D:/analysis/SPF_M/')

bm_GF <- CreateSeuratObject(counts = bm_GF.data, project = "GF", min.cells = 3, min.features = 100)
bm_SPF <- CreateSeuratObject(counts = bm_SPF.data, project = "SPF", min.cells = 3, min.features = 100)
bm_GFM <- CreateSeuratObject(counts = bm_GFM.data, project = "GFM", min.cells = 3, min.features = 100)
bm_SPFM <- CreateSeuratObject(counts = bm_SPFM.data, project = "SPFM", min.cells = 3, min.features = 100)
#
bm_GF$Label <- "GF"
bm_GF <- subset(bm_GF, subset = nFeature_RNA > 100)
bm_GF <- NormalizeData(bm_GF, verbose = FALSE)
bm_GF <- FindVariableFeatures(bm_GF, selection.method = "vst", nfeatures = 2000)
bm_GF <- ScaleData(bm_GF)
bm_GF <- RunPCA(bm_GF)
bm_GF <- RunUMAP(bm_GF, dims = 1:10)
bm_GF <- Find_doublet(bm_GF)
bm_GF <- subset(bm_GF,subset=doublet_info=="Singlet")

bm_SPF$Label <- "SPF"
bm_SPF <- subset(bm_SPF, subset = nFeature_RNA > 100)
bm_SPF <- NormalizeData(bm_SPF, verbose = FALSE)
bm_SPF <- FindVariableFeatures(bm_SPF, selection.method = "vst", nfeatures = 2000)
bm_SPF <- ScaleData(bm_SPF)
bm_SPF <- RunPCA(bm_SPF)
bm_SPF <- RunUMAP(bm_SPF, dims = 1:10)
bm_SPF <- Find_doublet(bm_SPF)
bm_SPF <- subset(bm_SPF,subset=doublet_info=="Singlet")

bm_GFM$Label <- "GFM"
bm_GFM <- subset(bm_GFM, subset = nFeature_RNA > 100)
bm_GFM <- NormalizeData(bm_GFM, verbose = FALSE)
bm_GFM <- FindVariableFeatures(bm_GFM, selection.method = "vst", nfeatures = 2000)
bm_GFM <- ScaleData(bm_GFM)
bm_GFM <- RunPCA(bm_GFM)
bm_GFM <- RunUMAP(bm_GFM, dims = 1:10)
bm_GFM <- Find_doublet(bm_GFM)
bm_GFM <- subset(bm_GFM,subset=doublet_info=="Singlet")

bm_SPFM$Label <- "SPFM"
bm_SPFM <- subset(bm_SPFM, subset = nFeature_RNA > 100)
bm_SPFM <- NormalizeData(bm_SPFM, verbose = FALSE)
bm_SPFM <- FindVariableFeatures(bm_SPFM, selection.method = "vst", nfeatures = 2000)
bm_SPFM <- ScaleData(bm_SPFM)
bm_SPFM <- RunPCA(bm_SPFM)
bm_SPFM <- RunUMAP(bm_SPFM, dims = 1:10)
bm_SPFM <- Find_doublet(bm_SPFM)
bm_SPFM <- subset(bm_SPFM,subset=doublet_info=="Singlet")
saveRDS(bm_GF,'bm_GF_normalize.rds')
saveRDS(bm_SPF,'bm_SPF_normalize.rds')
saveRDS(bm_GFM,'bm_GFM_normalize.rds')
saveRDS(bm_SPFM,'bm_SPFM_normalize.rds')
##
all_features <- lapply(list(bm_GF, bm_SPF,bm_GFM,bm_SPFM), row.names) %>% Reduce(intersect, .)
immune.anchors <- FindIntegrationAnchors(object.list = list(bm_GF, bm_SPF,bm_GFM,bm_SPFM), dims = 1:30 ,anchor.features = 3000)
bm <- IntegrateData(anchorset = immune.anchors, dims = 1:30, features.to.integrate = all_features)
saveRDS(bm,'bm_dataintegration_V3.rds')
############ filter ##########
bm <- readRDS('bm_dataintegration_V3.rds')
DefaultAssay(bm) <- "RNA"
#
bm[["percent.mt"]] <- PercentageFeatureSet(bm, pattern = "^mt-")
png("VlnPlot.png",width = 3000,height = 1000)
VlnPlot(bm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
#
bm_SPFM<-readRDS('bm_SPFM_normalize.rds')
bm_SPFM[["percent.mt"]] <- PercentageFeatureSet(bm_SPFM, pattern = "^mt-")
png("VlnPlot_SPFM.png",width = 3000,height = 1000)
VlnPlot(bm_SPFM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 3)
dev.off()
QC.data=data.frame(Barcode=rownames(bm_SPFM[['nCount_RNA']]),nCount_RNA=bm_SPFM[['nCount_RNA']],nFeature_RNA=bm_SPFM[['nFeature_RNA']],percent.mt=bm_SPFM[['percent.mt']],percent.ribo=bm_SPFM[['percent.ribo']])
write.table(QC.data,file="D:/analysis/qc.xls",sep="\t",quote=F,row.names=F)
#
mean(bm$nFeature_RNA)
median(bm$nFeature_RNA)
var(bm$nFeature_RNA)
sd(bm$nFeature_RNA)
#
mean(bm$nCount_RNA)
median(bm$nCount_RNA)
var(bm$nCount_RNA)
sd(bm$nCount_RNA)
#
plot1 <- FeatureScatter(bm, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png("Scatter.png",width = 2000,height = 1000)
plot1 + plot2
dev.off()
#
bm <- subset(bm, subset = nFeature_RNA > 200 & percent.mt < 50 )
saveRDS(bm,'bm_afterfilter.rds')
############ clustering ############
rds='bm_afterfilter.rds'
i=0.5
outpathway='./00_all_umap/'

bm<-readRDS(rds)
DefaultAssay(bm) <- "integrated"
bm <- ScaleData(bm, verbose = FALSE)
bm <- RunPCA(bm, npcs = 30, verbose = FALSE)
bm <- FindNeighbors(bm, reduction = "pca", dims = 1:30)
bm <- FindClusters(bm, resolution = i)
bm <- RunUMAP(bm, umap.method = 'uwot', reduction = "pca", dims = 1:30, metric = "correlation" )
bm <- RunTSNE(bm, reduction = "pca", dims = 1:30)

saveRDS(bm,paste0(outpathway,'all_afterumap_',i,'.rds'))

pdf(paste0(outpathway,"umap_",i,".pdf"),width = 10,height = 8)
DimPlot(bm, reduction = "umap", pt.size = 1.2, label.size = 6, label = TRUE)
dev.off()

pdf(paste0(outpathway,"tsne_",i,".pdf"),width = 10,height = 8)
DimPlot(bm, reduction = "tsne", pt.size = 1.2, label.size = 6, label = TRUE)
dev.off()

DefaultAssay(bm) <- "RNA"
library(future)
plan(strategy = "multicore", workers = 4)
options(future.globals.maxSize = 8000 * 1024^2)
inputmarkers <- FindAllMarkers(bm)
write.table(inputmarkers,paste0(outpathway,'marker_resolution_',i,'.xls'),sep="\t",quote = FALSE)
############ annotation #############
#Level1 annotation
bm <- readRDS('F:/analysis/Cecum_10X+C4/data/cecum_IntegrateP.rds')
DefaultAssay(bm) <- "RNA"
bm <- RenameIdents(bm,  '0' = "B cell", '1' = "Epithelial", '2' = "Fibroblast", `3` = "T/NK cell", `4` = "Epithelial", `5` = "T/NK cell", `6` = "Epithelial", `7` = "Epithelial",
                   '8' = "Myeloid", '9' = "Neuronal", '10' = "T/NK cell", `11` = "T/NK cell", `12` = "Epithelial", `13` = "unknown", `14` = "Fibroblast", `15` = "SMC",
                   '16' = "B cell", '17' = "Epithelial", '18' = "B cell", `19` = "LEC", `20` = "Endothelial", `21` = "Epithelial", `22` = "Terminal glial cell", `23` = "Myeloid",
                   '24' = "Mesothelial cell", '25' = "Epithelial", '26' = "T/NK cell")
bm@meta.data$celltypes<-bm@active.ident
saveRDS(bm,'cecum_IntegrateA.rds')
#Level2 annotation
bm <- readRDS('F:/analysis/Cecum_10X+C4/data/cecum_IntegrateA2.rds')
DefaultAssay(bm) <- "RNA"
bm <- RenameIdents(bm,  '0' = "B cell", '1' = "Enterocyte", '2' = "Fibroblast(Gsn high)", `3` = "NKT", `4` = "Enterocyte", `5` = "Naive Cd4 T", `6` = "Goblet", `7` = "Enterocyte",
                     '8' = "Macrophage", '9' = "Neuronal", '10' = "Cd4 T", `11` = "NK", `12` = "Enterocyte", `13` = "unknown", `14` = "Fibroblast(Adamdec1 high)", `15` = "SMC",
                     '16' = "Iga_Plasma", '17' = "Enterocyte", '18' = "B_Cycling", `19` = "LEC", `20` = "Endothelial", `21` = "Tuft", `22` = "Terminal glial cell", `23` = "Mast cell",
                     '24' = "Mesothelial cell", '25' = "Enteroendocrine", '26' = "NKT")
bm@meta.data$celltypes<-bm@active.ident
bm_subset<-subset(bm,idents=c('B cell','Naive Cd4 T','Fibroblast(Gsn high)','Macrophage','Goblet','Enterocyte','NK','Enteroendocrine','NKT','Cd4 T','Tuft','Fibroblast(Adamdec1 high)','SMC','Mast cell','Iga_Plasma','B_Cycling'))
bm_subset@meta.data$celltypes<-bm_subset@active.ident
saveRDS(bm_subset,'F:/analysis/Cecum_10X+C4/data/cecum_IntegrateA2_scmeta.rds')
#T/NK Subclass Annotations
bm_T <- readRDS('F:/analysis/Cecum_10X+C4/data/cecum_T_NK_cell_Integrate.rds')
DefaultAssay(bm_T) <- "RNA"
bm_T <- RenameIdents(bm_T,  '0' = "Immature NKT", '1' = "Naive Cd4 T", '2' = "T cell", `3` = "Treg", `4` = "T_Dying", `5` = "Mature NKT", `6` = "Th17", `7` = "Mature B")
bm_T@meta.data$Level1<-bm_T@meta.data$celltypes
bm_T@meta.data$celltypes<-bm_T@active.ident
saveRDS(bm_T,'cecum_T_NK_cell_IntegrateA.rds')
#EP Subclass Annotations
bm_Ep <- readRDS('F:/analysis/Cecum_10X+C4/data/cecum_Epithelial_Integrate.rds')
DefaultAssay(bm_Ep) <- "RNA"
bm_Ep <- RenameIdents(bm_Ep,  '0' = "EC(Saa1 high)", '1' = "EC(Selenbp1 high)", '2' = "EC(Hmgb2 high)", `3` = "Goblet", `4` = "mLTo", `5` = "EC(Rpl3 high)", `6` = "Tuft", `7` = "EEC",`8` = "EC(Dying)",`9` = "Ep_Fibroblast",`10` = "EC(Reg3b high)")
bm_Ep@meta.data$Level1<-bm_Ep@meta.data$celltypes
bm_Ep@meta.data$celltypes<-bm_Ep@active.ident
saveRDS(bm_Ep,'cecum_Epithelial_IntegrateA.rds')

############ Cell ratio_ new #############
suppressMessages(library(ggrepel))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(hash))
library(rstatix)
library(ggpubr)
setwd('F:/analysis/Cecum_10X+C4/Ep')
# ReadRDS #
RDS <- readRDS('F:/analysis/Cecum_10X+C4/data/cecum_Epithelial_IntegrateA.rds')
# celltypes #
celltypes <- read.csv('F:/analysis/Cecum_10X+C4/bin/Ep_celltypes_plotindex.txt',sep = '\t')
len <- length(unique(celltypes$annotation))
# color #
color <- hash()
color<-h
color<-Ep
### Function ###
Do_stat <- function(object){
  stat <- table(object@meta.data$celltypes,object@meta.data$mice)
  stat <- as.data.frame(stat)
  colnames(stat) <- c("celltypes","mice","counts")
  stat3 <- stat %>% group_by(mice) %>% summarise(sum = sum(counts))
  stat4 <- left_join(stat,stat3)  %>% mutate(diff = sum - counts) %>% 
    mutate(percentage = counts/sum * 100)
  return(stat4)
}
# output a signif data
Do_FisterTest <-  function(data,celltypes){
  celltype <- c()
  pvalue <-c()
  y <- c()
  for(i in celltypes){
    data2 <- data %>% filter(celltypes == i)
    rownames(data2) <- data2$mice
    yp <- max(data2$percentage) + 2
    data2 <- data2 %>% select(counts,diff)
    data2 <- t(as.matrix(data2))
    p <- fisher.test(data2)
    celltype <- c(celltype,i)
    pvalue <- c(pvalue,p$p.value)
    y <- c(y,yp)
  }
  print(y) 
  print(celltype)
  print(pvalue)
  dat <- data.frame(celltypes = celltype , pvalue = pvalue,ypostion = y)
  return(dat)
}
# Do a signif to signif.data
Do_signif <- function(plot.data,signif.data){
  signif.data <- signif.data %>% add_significance()
  GF <- plot.data %>% filter(mice == "GF") %>% select(celltypes,mice)
  colnames(GF) <- c("celltypes","group1")
  SPF <- plot.data %>% filter(mice == "SPF") %>% select(celltypes,mice)
  colnames(SPF) <- c("celltypes","group2")
  group <- left_join(GF,SPF)
  signif.data <- left_join(group,signif.data)
  signif.data$xmin = 1
  signif.data$xmax = 2
  signif.data$celltypes <- factor(signif.data$celltypes,levels = unique(plot.data$celltypes))
  return(signif.data)
}
# order of bar position 
Do_OrderPlotData <- function(plot.data,signif.data){
  order_index <- order(signif.data$ypostion)
  signif.data <- signif.data[order_index,]
  plot.data$celltypes <- factor(plot.data$celltypes, 
                                levels= signif.data$celltypes)
  return(plot.data)
}
# order of Color Index
Do_OrderColorIndex <- function(celltypes,hash){
  celltypes <- celltypes %>% mutate(cell_type = paste(clusters,annotation,sep = " "))
  for(i in 1:length(celltypes$annotation)){
    celltypes$color[i] <- hash[[celltypes$annotation[i]]]
  }
  celltypes <- celltypes %>% select(annotation,color)
  Index <- celltypes %>% duplicated()
  celltypes <- celltypes[!Index,]
  IndexColor <- celltypes$color
  names(IndexColor) <-celltypes$annotation
  return(IndexColor)
}
# bar color 
Do_Plot3 <- function(data,signif,Colorplate,font.family = "Times",axis.text.fsize=20,axis.title.fsize=25,
                     strip.fsize =20,strip.text.angle = 45){
  p <- ggplot(data,aes(x =  mice,y = percentage,fill = celltypes)) + 
    geom_bar(stat = "identity",position = "dodge",width = 0.7,
             colour = "black") + 
    facet_grid(. ~ celltypes) + theme_classic() + 
    ylab("The percentage of celltypes") + xlab("") + 
    theme(plot.margin = unit(c(1,1,1,1.5),"cm")) +
    theme(legend.position = 'none') + 
    scale_fill_manual(values = Colorplate) +
    theme( axis.text.x = element_text(vjust = -2))+
    theme(axis.text = element_text(size = axis.text.fsize),
          axis.title  = element_text(size = axis.title.fsize) , 
          axis.title.y  = element_text(vjust =  7) ,
          axis.line = element_line(linetype = 1,color= "black",size = 1),
          axis.ticks  = element_line(color = "black",size = 1,lineend = 2)) +
    theme(axis.text.y.left = element_text(hjust = 1, size=axis.text.fsize,
                                          margin = margin(0,18,0,0)))
  p2 <- p + stat_pvalue_manual(data = signif,label = "pvalue.signif",
                               y.position = "ypostion",size = 6) + 
    theme(strip.background = element_blank(),
          strip.text =  element_text(size = strip.fsize,
                                     angle = strip.text.angle)) + 
    theme(text = element_text(family = font.family))
  return(p2)
}
Do_CelltypePercentage <- function(Seurat,celltypes,color,font.family,axis.text.fsize,axis.title.fsize,strip.text.angle,strip.fsize){
  stat <- Do_stat(Seurat)
  celltype <- unique(Seurat@meta.data$celltypes)
  test <- Do_FisterTest(data =  stat,celltypes = celltype)
  signif <- Do_signif(stat,test)
  stat  <- Do_OrderPlotData(stat,signif)
  Color <- Do_OrderColorIndex(celltypes = celltypes,hash = color)
  p <- Do_Plot3(data = stat, signif = signif,Colorplate =  Color ,
                font.family = font.family , 
                axis.text.fsize = axis.text.fsize,
                axis.title.fsize = axis.title.fsize,
                strip.text.angle = strip.text.angle,
                strip.fsize = strip.fsize )
  return(p)
}

pdf("Ep_CelltypesPercentage1.pdf",w = 20 ,h = 5)
Do_CelltypePercentage(Seurat = RDS, celltypes = celltypes, color = color , font.family = "sans" , axis.text.fsize = 25 , axis.title.fsize = 25, strip.text.angle = 0 , strip.fsize = 13)
dev.off()

############ Annotated figure #############
suppressMessages(library(ggrepel))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(hash))
library(rstatix)
library(ggpubr)
setwd('F:/analysis/Cecum_10X+C4/Level2')

# ReadRDS #
RDS <- readRDS('F:/analysis/Cecum_10X+C4/data/cecum_IntegrateA2_New.rds')
# celltypes #
celltypes <- read.csv('F:/analysis/Cecum_10X+C4/bin/Level2_celltypes_plotindex_New.txt',sep = '\t')
len <- length(unique(celltypes$annotation))
print(len)
# color #
color <- hash()
color <- h
### Function ###
GetUmapEmbeddings <- function(Seurat){ 
  umap = Seurat@reductions$umap@cell.embeddings %>%  
    as.data.frame() %>% cbind(cell_type = Seurat@meta.data$celltypes) %>% 
    cbind(cluster = Seurat@meta.data$seurat_clusters)
  return(umap)
}
Do_hashTable <- function(celltypes){
  celltypes <- celltypes %>% select(plotindex,annotation)
  Index <- celltypes %>% duplicated()
  celltypes <- celltypes[!Index,] 
  celltypes <- celltypes %>% mutate(factor = paste(plotindex,annotation,sep = " ")) %>% arrange(plotindex)
  return(celltypes)
}
Do_UmapPlotData <- function(umap,hashtable){
  umap$cell_type <- as.vector(umap$cell_type)
  factor_hash.list <- hash(hashtable$annotation,hashtable$factor)
  print(factor_hash.list)
  num_hash.list <- hash(hashtable$annotation,hashtable$plotindex)
  print(num_hash.list)
  for(i in 1:length(umap$cell_type)){
    umap$annot[i] <- factor_hash.list[[umap$cell_type[i]]]
  }
  for(i in 1:length(umap$cell_type)){
    umap$plotindex[i] <- num_hash.list[[umap$cell_type[i]]]
  }
  umap$annot <- factor(umap$annot,level = hashtable$factor)
  return(umap)
}
Do_ColorIndex <- function(celltype,color){
  for(i in 1:length(celltype$annotation)){
    celltype$color[i] <- color[[celltype$annotation[i]]]
  }
  A <- celltype$color
  names(A) <- celltype$factor
  return(A)
}
GetCelltypesCenter <- function(umap){
  celltypesCenter <- umap %>% group_by(plotindex) %>%summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  return(celltypesCenter)
}
DoUmapPlot2 <- function(umap,CelltypesCenter,Colorplate,legend.fsize = 20,pindex.fsize = 7,font.family = "Times",background = "white"){
  p1 <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = annot)) +  geom_point(size = 0.2 , alpha =1 )  + 
    scale_color_manual(values = Colorplate) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          axis.title = element_blank(),  
          axis.text = element_blank(), 
          axis.ticks = element_blank())           
  if (background == "white"){
    p2 <- p1 + geom_segment(aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2) ,
                                xend = min(umap$UMAP_1) +3, yend = min(umap$UMAP_2) ),
                            colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+
      geom_segment(aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,
                       xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 3),
                   colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
      annotate("text", x = min(umap$UMAP_1) +1.5, y = min(umap$UMAP_2) -1, label = "UMAP_1",
               color="black",size = 5, fontface="bold" ) +
      annotate("text", x = min(umap$UMAP_1) -1, y = min(umap$UMAP_2) + 1.5, label = "UMAP_2",
               color="black",size = 5, fontface="bold" ,angle=90) +
      geom_text_repel(aes(label=plotindex), fontface="bold",data = CelltypesCenter,
                      point.padding=unit(2, "lines"),color = "black",size = pindex.fsize) + theme(text = element_text(family = font.family))
    p3 <- p2 +  theme(panel.background = element_rect(fill = 'white'), plot.background = element_rect(fill = 'white'),  
                      legend.title = element_blank(),
                      legend.key=element_rect(fill='white'),
                      legend.text = element_text(size=legend.fsize),
                      legend.key.size=unit(1,'cm')) + guides(color = guide_legend(override.aes = list(size=5)))
  } else if (background == "black"){
    p2 <- p1 + geom_segment(aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2) ,
                                xend = min(umap$UMAP_1) +3, yend = min(umap$UMAP_2) ),
                            colour = "white", size=1,arrow = arrow(length = unit(0.3,"cm")))+
      geom_segment(aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,
                       xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 3),
                   colour = "white", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
      annotate("text", x = min(umap$UMAP_1) +1.5, y = min(umap$UMAP_2) -1, label = "UMAP_1",
               color="white",size = 5, fontface="bold" ) +
      annotate("text", x = min(umap$UMAP_1) -1, y = min(umap$UMAP_2) + 1.5, label = "UMAP_2",
               color="white",size = 5, fontface="bold" ,angle=90) +
      geom_text_repel(aes(label=plotindex), fontface="bold",data = CelltypesCenter,
                      point.padding=unit(2, "lines"),color = "white",size = pindex.fsize) + theme(text = element_text(family = font.family))
    p3 <- p2 + theme(panel.background  = element_rect(fill = "black"), plot.background = element_rect(fill = 'black'),
                     legend.background = element_rect(fill = "black"),
                     legend.title  = element_blank(),
                     legend.key=element_rect(fill='black'),
                     legend.text = element_text(size=legend.fsize,color = "white"),
                     legend.key.size=unit(1,'cm')) + guides(color = guide_legend(override.aes = list(size=5)))
  } else {
    stop("send a message to Lwm with Wechat!")
  }
  return(p3)
}
Do_CelltypesUmap <- function(Seurat,celltypes,color,legend.fsize,pindex.fsize,font.family,background){
  umap <- GetUmapEmbeddings(Seurat = Seurat)
  celltype <- Do_hashTable(celltypes = celltypes)
  umap <- Do_UmapPlotData(umap = umap, hashtable = celltype)
  ColorIndex <- Do_ColorIndex(celltype,color = color)
  CelltypesCenter <- GetCelltypesCenter(umap)
  p <- DoUmapPlot2(umap,CelltypesCenter,ColorIndex,legend.fsize,pindex.fsize,font.family,background)
  return(p)
}

# plot #
pdf("Level2_CelltypesUmap_black_New.pdf",w = 13,h = 8)
Do_CelltypesUmap(Seurat = RDS , celltypes = celltypes , color = color , legend.fsize = 20 , pindex.fsize = 7 ,
                 font.family = "Times",background = "black")
dev.off()


############ Marker genes for each cluster################
gene <- read.table("F:/analysis/Cecum_10X+C4/gene_Marker.txt")
target_gene <- as.character(gene[,1])
png("FeaturePlot_LEC.png",width = 800,height = 900)
FeaturePlot(bm, features = target_gene, min.cutoff = "q9")
dev.off()

pdf("FeaturePlot_EP.pdf",width = 8,height = 9)
FeaturePlot(bm, features = target_gene, min.cutoff = "q9")
dev.off()
############ Gene express ################
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(ggthemes)
suppressMessages(library(ggsci))
setwd("F:/analysis/Cecum_10X+C4")
bm <- readRDS('F:/analysis/Cecum_10X+C4/data/cecum_T_NK_cell_IntegrateA.rds')
bm_NK <- subset(bm,idents=c('Resting NKT','Activated NKT'))
bm_NK@meta.data$celltypes<-bm_NK@active.ident
p = VlnPlot(bm_NK, features = 'Gzmk',
            pt.size = 0,
            group.by = 'celltypes',
            split.by = 'mice',
            ncol = 5
)
p


celltype='Epithelial'
Features=c('Gzma','Gzmb')
for(i in Features){
  options(repr.plot.width = 8,repr.plot.height = 5)
  DefaultAssay(bm_Ep)
  str(bm_Ep)
  Idents(bm_Ep)<-"Label"
  p = VlnPlot(bm_Ep, features = i,
              pt.size = 0,
              group.by = NULL,
              ncol = 5
  ) & geom_boxplot(width=.1,col="black",fill="white")& stat_compare_means(label.y =6.5)&ylim(0,7)&labs(title= i,x= celltype,y="Expression Level")
  #p
  pdf(paste0(i,'_',celltype,'.pdf'),width = 11,height = 5)
  print(p)
  dev.off()
}

############ EP Difference Analysis###############
setwd('F:/analysis/Cecum_10X+C4/Ep/MAST_0.25')
bm_Ep <- readRDS('F:/analysis/Cecum_10X+C4/data/cecum_Epithelial_IntegrateA.rds')

library(Seurat)
library(ggpubr)
library(ggthemes)

x1=unique(as.character(bm_Ep$celltypes))
meta
for(i in unique(x1)){
  print(i)
  markers <- FindMarkers(bm_Ep, test.use = "MAST" ,ident.1 = "GF",ident.2 = "SPF" , group.by = 'mice', subset.ident = i,logfc.threshold = 0.15)
  
  markers$gene <- row.names(markers)
  markers$logP <- -log10(markers$p_val_adj)
  
  markers$group="not-significant"
  markers$group[which( (markers$p_val_adj < 0.05) & (markers$avg_log2FC > 0.25 ) )] = "up-regulated"
  markers$group[which( (markers$p_val_adj < 0.05) & (markers$avg_log2FC < -0.25 ) )] = "down-regulated"
  
  markers$Label = ""
  markers <- markers[order(markers$p_val_adj),]
  up.genes <-head(markers$gene[which(markers$group == "up-regulated")],30)
  down.genes <-head(markers$gene[which(markers$group == "down-regulated")],30)
  markers.top30<-c(as.character(up.genes), as.character((down.genes)))
  markers$Label[match(markers.top30, markers$gene)] <- markers.top30
  
  write.table(markers,paste0("MAST_",i,"_GF_SPF.txt"),sep='\t',quote=FALSE)
  
  pdf(paste0("MAST_",i,"_GF_SPF.pdf"),width = 9,height = 6)
  p<-ggscatter(markers,x="avg_log2FC",y="logP",color = "group",palette = c("#2f5688","#BBBBBB","#CC0000"), 
               size=0.5 ,label = markers$Label, font.label = 8 ,repel = F, xlab = "avg_log2FoldChange", 
               ylab = "-log10(Adjusted P-value)",title = i) +
    theme_base()+ scale_x_continuous(limits = c(-2.5, 2.5)) + scale_y_continuous(limits = c(-2,max(markers$logP)+20))+
    geom_hline(yintercept = 1.30, linetype="dashed")+
    geom_vline(xintercept = c(-0.58,0.58), linetype="dashed") 
  print(p)
  dev.off()
}
############ TNK difference analysis ##############
bm_T <- readRDS('F:/analysis/Cecum_10X+C4/data/cecum_T_NK_cell_IntegrateA.rds')
setwd('F:/analysis/Cecum_10X+C4/TNK')
library(Seurat)
library(ggpubr)
library(ggthemes)

x1=unique(as.character(bm_T$celltypes))
meta<-t(table(bm_T@meta.data$mice,bm_T@active.ident))
meta
for(i in unique(x1)){
  print(i)
  markers <- FindMarkers(bm_T, test.use = "MAST" ,ident.1 = "GF",ident.2 = "SPF" , group.by = 'mice', subset.ident = i,logfc.threshold = 0.15)
  
  markers$gene <- row.names(markers)
  markers$logP <- -log10(markers$p_val_adj)
  
  markers$group="not-significant"
  markers$group[which( (markers$p_val_adj < 0.05) & (markers$avg_log2FC > 0.58 ) )] = "up-regulated"
  markers$group[which( (markers$p_val_adj < 0.05) & (markers$avg_log2FC < -0.58 ) )] = "down-regulated"
  
  markers$Label = ""
  markers <- markers[order(markers$p_val_adj),]
  up.genes <-head(markers$gene[which(markers$group == "up-regulated")],30)
  down.genes <-head(markers$gene[which(markers$group == "down-regulated")],30)
  markers.top30<-c(as.character(up.genes), as.character((down.genes)))
  markers$Label[match(markers.top30, markers$gene)] <- markers.top30
  
  write.table(markers,paste0("MAST_",i,"_GF_SPF.txt"),sep='\t',quote=FALSE)
  
  pdf(paste0("MAST_",i,"_GF_SPF.pdf"),width = 9,height = 6)
  p<-ggscatter(markers,x="avg_log2FC",y="logP",color = "group",palette = c("#2f5688","#BBBBBB","#CC0000"), 
               size=0.5 ,label = markers$Label, font.label = 8 ,repel = F, xlab = "avg_log2FoldChange", 
               ylab = "-log10(Adjusted P-value)",title = i) +
    theme_base()+ scale_x_continuous(limits = c(-2.5, 2.5)) + scale_y_continuous(limits = c(-2,max(markers$logP)+20))+
    geom_hline(yintercept = 1.30, linetype="dashed")+
    geom_vline(xintercept = c(-0.58,0.58), linetype="dashed") 
  print(p)
  dev.off()
}


############ Level 1&2 Difference Analysis ################
library(Seurat)
library(ggpubr)
library(ggthemes)
library(SeuratObject)
bm <- readRDS('F:/analysis/Cecum_10X+C4/data/cecum_IntegrateA2.rds')
DefaultAssay(bm) <- "RNA"
#Level2
bm <- RenameIdents(bm,  '0' = "B cell", '1' = "Enterocyte", '2' = "Fibroblast(Gsn high)", `3` = "NKT", `4` = "Enterocyte", `5` = "Naive Cd4 T", `6` = "Goblet", `7` = "Enterocyte", `8` = "Macrophage", `9` = "Neuronal",
                   `10` = "Cd4 T", `11` = "NK", `12` = "Enterocyte",`13`= "unknown", `14` = "Fibroblast(Adamdec1 high)",`15` = "SMC",`16` = "Iga_Plasma",`17` = "Enterocyte",`18` = "B_Cycling",
                   `19` = "LEC",`20` = "Endothelial",`21` = "Tuft",`22` = "Terminal glial cell",`23` = "Mast cell",`24` = "Mesothelial cell",`25` = "Enteroendocrine",`26` = "NKT")
bm@meta.data$level2 <- bm@active.ident
#Level1
bm <- RenameIdents(bm,  '0' = "B cell", '1' = "Epithelial", '2' = "Fibroblast", `3` = "T/NK cell", `4` = "Epithelial", `5` = "T/NK cell", `6` = "Epithelial", `7` = "Epithelial", `8` = "Myeloid", `9` = "Neuronal",
                   `10` = "T/NK cell", `11` = "T/NK cell", `12` = "Epithelial",`13`= "unknown", `14` = "Fibroblast",`15` = "SMC",`16` = "B cell",`17` = "Epithelial",`18` = "B cell",
                   `19` = "LEC",`20` = "Endothelial",`21` = "Epithelial",`22` = "Terminal glial cell",`23` = "Myeloid",`24` = "Mesothelial cell",`25` = "Epithelial",`26` = "T/NK cell")
bm@meta.data$level1 <- bm@active.ident

x1=unique(as.character(bm$celltypes))
for(i in unique(x1)){
  print(i)
  markers <- FindMarkers(bm, test.use = "MAST" ,ident.1 = "GF",ident.2 = "SPF" , group.by = 'mice', subset.ident = i,logfc.threshold = 0.15)
  
  markers$gene <- row.names(markers)
  markers$logP <- -log10(markers$p_val_adj)
  
  markers$group="not-significant"
  markers$group[which( (markers$p_val_adj < 0.05) & (markers$avg_log2FC > 0.58 ) )] = "up-regulated"
  markers$group[which( (markers$p_val_adj < 0.05) & (markers$avg_log2FC < -0.58 ) )] = "down-regulated"
  
  markers$Label = ""
  markers <- markers[order(markers$p_val_adj),]
  up.genes <-head(markers$gene[which(markers$group == "up-regulated")],30)
  down.genes <-head(markers$gene[which(markers$group == "down-regulated")],30)
  markers.top30<-c(as.character(up.genes), as.character((down.genes)))
  markers$Label[match(markers.top30, markers$gene)] <- markers.top30
  
  write.table(markers,paste0("MAST_",i,"_GF_SPF.txt"),sep='\t',quote=FALSE)
  
  pdf(paste0("MAST_",i,"_GF_SPF.pdf"),width = 9,height = 6)
  p<-ggscatter(markers,x="avg_log2FC",y="logP",color = "group",palette = c("#2f5688","#BBBBBB","#CC0000"), 
               size=0.5 ,label = markers$Label, font.label = 8 ,repel = F, xlab = "avg_log2FoldChange", 
               ylab = "-log10(Adjusted P-value)",title = i) +
    theme_base()+ scale_x_continuous(limits = c(-2.5, 2.5)) + scale_y_continuous(limits = c(-2,max(markers$logP)+20))+
    geom_hline(yintercept = 1.30, linetype="dashed")+
    geom_vline(xintercept = c(-0.58,0.58), linetype="dashed") 
  print(p)
  dev.off()
}



############ Differential gene map############
library(tidyverse)
library(ggrepel)
setwd('F:/analysis/Cecum_10X+C4/Ep/Mast')
x1 = c('EC(Saa1 high)','EC(Selenbp1 high)','EC(Hmgb2 high)','Goblet','mLTo','EC(Rpl3 high)','Tuft','EEC','EC(Dying)','Ep_Fibroblast','EC(Reg3b high)')

alltype = list()
for(i in x1){
  print(i)
  sce.markers <- read.delim(file = paste0("MAST_",i,"_GF_SPF.txt"),sep = '\t')
  sce.markers$celltype <- i
  #rbind(df,sce.markers)
  alltype[[i]] = sce.markers
}
df <- do.call(rbind, alltype)
unique(df$celltype)

df$label <- ifelse(df$p_val_adj<0.05,"adjusted P-val<0.05","adjusted P-val>=0.05")
unique(df$label)

top10sig0 <- filter(df,celltype =="EC(Saa1 high)") %>% distinct(gene,.keep_all = T) %>% top_n(6,abs(avg_log2FC))
top10sig1 <- filter(df,celltype =="EC(Selenbp1 high)") %>% distinct(gene,.keep_all = T) %>% top_n(6,abs(avg_log2FC))
top10sig2 <- filter(df,celltype =="EC(Hmgb2 high)") %>% distinct(gene,.keep_all = T) %>% top_n(6,abs(avg_log2FC))
top10sig3 <- filter(df,celltype =="Goblet") %>% distinct(gene,.keep_all = T) %>% top_n(6,abs(avg_log2FC))
top10sig4 <- filter(df,celltype =="mLTo") %>% distinct(gene,.keep_all = T) %>% top_n(6,abs(avg_log2FC))
top10sig5 <- filter(df,celltype =="EC(Rpl3 high)") %>% distinct(gene,.keep_all = T) %>% top_n(6,abs(avg_log2FC))
top10sig6 <- filter(df,celltype =="Tuft") %>% distinct(gene,.keep_all = T) %>% top_n(6,abs(avg_log2FC))
top10sig7 <- filter(df,celltype =="EEC") %>% distinct(gene,.keep_all = T) %>% top_n(6,abs(avg_log2FC))
top10sig8 <- filter(df,celltype =="EC(Dying)") %>% distinct(gene,.keep_all = T) %>% top_n(6,abs(avg_log2FC))
top10sig9 <- filter(df,celltype =="Ep_Fibroblast") %>% distinct(gene,.keep_all = T) %>% top_n(6,abs(avg_log2FC))
top10sig10 <- filter(df,celltype =="EC(Reg3b high)") %>% distinct(gene,.keep_all = T) %>% top_n(6,abs(avg_log2FC))
top10sig <- rbind(top10sig0,top10sig1,top10sig2,top10sig3,top10sig4,top10sig5,top10sig6,top10sig7,top10sig8,top10sig9,top10sig10) 

df$size <- case_when(!(df$gene %in% top10sig$gene)~ 1,
                     df$gene %in% top10sig$gene ~ 2)


dt <- filter(df,size==1)
unique(dt$celltype)

options(repr.plot.width = 15,repr.plot.height =10)
p <- ggplot()+
  geom_jitter(data = df,aes(x = celltype, y = avg_log2FC, color = label),size = 0.4,width = 0.4)
p


options(repr.plot.width = 15,repr.plot.height =10)
p <- ggplot()+
  geom_jitter(data = dt,aes(x = celltype, y = avg_log2FC, color = label),size = 0.8,width =0.4)+
  geom_jitter(data = top10sig,aes(x = celltype, y = avg_log2FC, color = label),size = 1.5,width =0.4)
p

rmax <- vector()
rmin <- vector()
for(i in x1){
  print(i)
  A <- filter(df,celltype == i) %>% distinct(gene,.keep_all = T) %>% top_n(10,abs(avg_log2FC))
  rmax = c(rmax,max(A$avg_log2FC)+0.1)
  rmin = c(rmin,min(A$avg_log2FC)-0.1)
}

options(repr.plot.width = 15,repr.plot.height =10)
dfbar<-data.frame(x=c('EC(Saa1 high)','EC(Selenbp1 high)','EC(Hmgb2 high)','Goblet','mLTo','EC(Rpl3 high)','Tuft','EEC','EC(Dying)','Ep_Fibroblast','EC(Reg3b high)')
                  ,y=c(rmax[1],rmax[2],rmax[3],rmax[4],0.8,rmax[6],
                       rmax[7],rmax[8],rmax[9],rmax[10],rmax[11]))
dfbar1<-data.frame(x=c('EC(Saa1 high)','EC(Selenbp1 high)','EC(Hmgb2 high)','Goblet','mLTo','EC(Rpl3 high)','Tuft','EEC','EC(Dying)','Ep_Fibroblast','EC(Reg3b high)'),
                   y=c(rmin[1],rmin[2],rmin[3],rmin[4],rmin[5],rmin[6],
                       rmin[7],rmin[8],rmin[9],rmin[10],rmin[11]))


p1 <- ggplot()+
  geom_col(data = dfbar,mapping = aes(x = x,y = y),fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,mapping = aes(x = x,y = y),fill = "#dcdcdc",alpha = 0.6)
p1


options(repr.plot.width = 15,repr.plot.height =10)
p2 <- ggplot()+
  geom_col(data = dfbar,mapping = aes(x = x,y = y),fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,mapping = aes(x = x,y = y),fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = dt,aes(x = celltype, y = avg_log2FC, color = label),size = 0.8,width =0.4)+
  geom_jitter(data = top10sig,aes(x = celltype, y = avg_log2FC, color = label),size = 1.5,width =0.4)
p2


dfcol<-data.frame(x=c('EC(Saa1 high)','EC(Selenbp1 high)','EC(Hmgb2 high)','Goblet','mLTo','EC(Rpl3 high)','Tuft','EEC','EC(Dying)','Ep_Fibroblast','EC(Reg3b high)'),y=0,label=c('EC(Saa1 high)','EC(Selenbp1 high)','EC(Hmgb2 high)','Goblet','mLTo','EC(Rpl3 high)','Tuft','EEC','EC(Dying)','Ep_Fibroblast','EC(Reg3b high)'))
mycol <- c('#FF4A46','#FF8C00','#FFFF00','#06D6A0','#1CE6FF','#E14AEC','#FF90C9','#FFE4E1','#9395E7','#B1CE46','#F8AC8C')
p3 <- p2 + geom_tile(data = dfcol,aes(x=x,y=y),height=0.3,color = "black",fill = mycol,alpha = 0.6,show.legend = F)
p3

p4 <- p3+ 
  geom_text_repel(data=top10sig, aes(x=celltype,y=avg_log2FC,label=gene),size=4, min.segment.length = 0,
                  fontface="italic", segment.color="grey30",position= position_jitter(width=0.1, seed=1),box.padding=0,point.padding = 0,segment.size=0.2)
p4

options(repr.plot.width = 15,repr.plot.height =10)
p5 <- p4+
  labs(x="Cluster",y="average logFC")+geom_text(data=dfcol,
                                                aes(x=x,y=y,label=label),size =3,color ="black")
p5

p6 <- p5+
  theme_minimal()+theme(axis.title = element_text(size = 13,color = "black",face = "bold"),
                        axis.line.y = element_line(color = "black",size = 1.2),axis.line.x = element_blank(),axis.text.x = element_blank(),panel.grid = element_blank(),
                        legend.position = "top",legend.direction = "vertical",legend.justification = c(1,0),legend.text = element_text(size = 15)
  )
p6

pdf("DEG_Ep_color_top6.pdf",width = 12,height = 8)
print(p6)
dev.off()
############ Marker gene dotplot #############
setwd('F:/analysis/Cecum_10X+C4/Fibroblast')
bm <- readRDS('F:/analysis/Cecum_10X+C4/data/cecum_IntegrateA2.rds')
DefaultAssay(bm)<-"RNA"
A<-as.character(unique(bm$celltypes))
A<-A[-11]
A<-A[-21]
bm<-subset(bm,idents=A)
bm$celltypes <- bm@active.ident

markers.to.plot <-c('Cd160','Xcl1',
                    'Igfbp4','Sell',
                    'Ncr1','Klrk1',
                    'Ctla4','Icos',
                    'mt-Co2','mt-Nd4',
                    'Gzma','Gzmb',
                    'Ncoa7','Sdc4',
                    'Cd79a','Ebf1')
bm@meta.data <- bm@meta.data %>% mutate(haha = paste(mice,celltypes,sep = "_"))
unique(bm$haha)
bm$haha <- factor(bm$haha, levels=c("GF_B cell","SPF_B cell","GF_Naive Cd4 T","SPF_Naive Cd4 T","GF_Fibroblast(Gsn high)","SPF_Fibroblast(Gsn high)","GF_Macrophage","SPF_Macrophage","GF_Goblet","SPF_Goblet",
                                    "GF_Enterocyte","SPF_Enterocyte",'GF_NK','SPF_NK','GF_Enteroendocrine','SPF_Enteroendocrine','GF_NKT','SPF_NKT','GF_Cd4 T','SPF_Cd4 T','GF_Tuft','SPF_Tuft','GF_Fibroblast(Adamdec1 high)','SPF_Fibroblast(Adamdec1 high)',
                                    'GF_SMC','SPF_SMC','GF_Mast cell','SPF_Mast cell','GF_Iga_Plasma','SPF_Iga_Plasma','GF_Endothelial','SPF_Endothelial','GF_B_Cycling','SPF_B_Cycling','GF_LEC','SPF_LEC','GF_Terminal glial cell','SPF_Terminal glial cell',
                                    'GF_Mesothelial cell','SPF_Mesothelial cell'), ordered=TRUE)
p1<-DotPlot(bm, features = c('Tgfb1','Tgfb2','Tgfb3'), group.by = 'haha') + RotatedAxis()
pdf(file="DotPlot.pdf",width = 6,height = 10)
print(p1)
dev.off()

ggplot(aes(x=celltypes,y=gene,size=ratio,color=exp))+geom_point()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  scale_size_continuous(limits = c(0,1))+theme_bw()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )



############ UMAP of DEG number ##########
library(plyr)
bm <- readRDS('F:/analysis/Cecum_10X+C4/data/cecum_IntegrateA_New.rds')
setwd('F:/analysis/Cecum_10X+C4/Level1/MAST')
bm <- RenameIdents(bm,  'T/NK cell' = "TNK")
bm$cell.ident.subset<-bm@active.ident
bm$cell.ident.main<-bm@active.ident
bm@meta.data$celltypes<-bm@active.ident
bm$cell.ident.subset <- factor(bm$cell.ident.subset, levels=c("B cell","TNK","Fibroblast","Myeloid","Epithelial","SMC","Endothelial","LEC","Terminal glial cell",
                                                              "Mesothelial cell"), ordered=TRUE)
genelist = data.frame(celltype = c("B cell","TNK","Fibroblast","Myeloid","Epithelial","SMC","Endothelial","LEC","Terminal glial cell",
                                   "Mesothelial cell"),number = rep(0,10))
for (i in dir(pattern='.txt$')){
  DEG <- read.delim(file = i,sep = '\t')
  i <- gsub('_GF_SPF.txt','',gsub('MAST_','',i))
  genelist[which(genelist$celltype == i),]$number <- length(DEG[which(DEG$group == 'down-regulated' | DEG$group == 'up-regulated'),]$gene)
}
bm$nDEG <- bm$cell.ident.subset
bm$nDEG <- mapvalues(bm$nDEG, from = as.character(genelist$celltype), to = genelist$number)
bm$nDEG <- as.integer(as.character(bm$nDEG))
max(bm$nDEG)
pE<-FeaturePlot(bm,features = 'nDEG',raster = F)
pE
ggsave('umap_DEG_New.pdf',plot = pE,width = 5,height = 4)
############ GO enrichment ##############
library(Seurat)
library(dplyr)
setwd('F:/analysis/Cecum_10X+C4/Ep/MAST_0.25')
bm <- readRDS('F:/analysis/Cecum_10X+C4/data/cecum_Epithelial_IntegrateA.rds')

EnrichGeneGo <- function(gene_set){
  Go_anno <- enrichGO(gene_set, OrgDb = go_organism, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = 'ENTREZID')
  Go_anno <- as.data.frame(Go_anno)
  return(Go_anno)
}

MAST_name=unique(as.character(bm$celltypes))

for(i in MAST_name){
  print(i)
  sce.markers <- read.delim(file = paste0("MAST_",i,"_GF_SPF.txt"),sep = '\t')
  colnames(sce.markers)[6] <- "Gene"
  sce.markers <- sce.markers %>% tibble::rownames_to_column('gene')
  head(sce.markers)

  library(ggpubr)
  library(clusterProfiler)
  options(connectionObserver = NULL) #加载org.Hs.eg.db失败时的解决方法
  options(stringsAsFactors = F)
  suppressMessages(library('org.Mm.eg.db'))
  keytypes(org.Mm.eg.db)
  ##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
  ##  [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
  ## [11] "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"        
  ## [16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
  ## [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"       "UNIGENE"     
  ## [26] "UNIPROT"
  
  go_organism <- "org.Mm.eg.db"      # org.Hs.eg.db/org.Mm.eg.db
  kegg_organism <- 'mmu'             # hsa/mmu
  
  ids = bitr(sce.markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=go_organism)
  sce.markers = merge(sce.markers, ids, by.x='gene', by.y='SYMBOL')
  sce.markers <- sce.markers[which((sce.markers$p_val_adj < 0.05) & ((sce.markers$avg_log2FC > 0.25 ) | (sce.markers$avg_log2FC < -0.25 ))),]
  sce.markers$group <- factor(ifelse(sce.markers$avg_log2FC < 0, -1, 1), levels = c(-1, 1))
  gcSample = split(sce.markers$ENTREZID, sce.markers$group)
  up_genes <- subset(sce.markers, group==1)$ENTREZID
  down_genes <- subset(sce.markers, group==-1)$ENTREZID
  
  up_GO <- EnrichGeneGo(gene_set=up_genes) %>% dplyr::mutate(group=rep(1, n()))
  down_GO <- EnrichGeneGo(gene_set=down_genes) %>% dplyr::mutate(group=rep(-1, n()))
  dat <- rbind(up_GO, down_GO)
  
  sample_names <- c(paste0('GF'," up-regulated"), paste0('SPF'," up-regulated"))
  dat <- dat %>%
    dplyr::mutate(group_type = factor(ifelse(group == 1, sample_names[1], sample_names[2]), levels = sample_names)) %>%
    dplyr::mutate(Gene_Number = Count * group)
  
  dat <- dat %>% dplyr::group_by(group) %>% dplyr::do(head(., n = 10))
  dat<-dat[c(1,7,13,14,17,20),]
  
  library(ggplot2)
  library(cowplot)
  th <- theme(axis.text.y = element_blank(), 
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              panel.background = element_rect(fill = 'white'), 
              panel.border = element_rect(color = 'black', fill = NA),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
              panel.grid.minor.x = element_blank())
  
  get_wraper <- function(width) {
    function(x) {
      lapply(strwrap(x, width = width, simplify = FALSE), paste, collapse="\n")
    }
  }
  
  pdf(paste0("GO_",i,"_1.pdf"),width = 9,height = 6)
  q <- ggplot(dat, aes(x = reorder(Description, Gene_Number), y = Gene_Number)) +
    geom_bar(aes(fill = group_type), stat = "identity") + 
    labs(x = '', y = 'Gene Number', fill = i) +
    scale_x_discrete(labels= get_wraper(100))+
    theme(axis.text = element_text(size=3),)+
    coord_flip() + theme(legend.position = 'bottom') + guides(fill = guide_legend(ncol = 1))
  print(q)
  dev.off()
  
  pdf(paste0("GO_",i,"_2.pdf"),width = 9,height = 6)
  p <- ggplot(dat, aes(x = reorder(Description, Gene_Number), y = Gene_Number)) +
    geom_bar(aes(fill = group_type), stat = "identity") + 
    labs(y = 'Gene Number', fill = i) +
    coord_flip() + 
    geom_text(aes(y = 0, label = Description, hjust = as.numeric(Gene_Number > 0))) +  # label text based on value
    th + theme(legend.position = 'bottom') + guides(fill = guide_legend(ncol = 2))
  print(p)
  dev.off()
  
  nbreaks = 10
  minimum <- floor(min(dat$Gene_Number)/nbreaks)*nbreaks
  maximum <- ceiling(max(dat$Gene_Number)/nbreaks)*nbreaks
  pdf(paste0("GO_",i,"_3.pdf"),width = 9,height = 6)
  r <- ggplot(dat, aes(x = reorder(Description, Gene_Number), y = Gene_Number)) +
    geom_bar(aes(fill = Gene_Number), stat = "identity") + 
    labs(y = 'Gene Number', fill = i) +
    scale_x_discrete(labels= get_wraper(35))+
    coord_flip() +  
    geom_text(aes(y = 0, label = Description, hjust = as.numeric(Gene_Number > 0))) +  # label text based on value
    th + 
    scale_y_continuous(limits = c(minimum, maximum), breaks = seq(minimum, maximum, nbreaks)) +
    scale_fill_gradient2(low = 'darkblue', high = 'red', mid = 'white',
                         limits = c(minimum, maximum), breaks = seq(minimum, maximum, nbreaks))
  print(r)
  dev.off()
}

############ CCI-Level1 #############
setwd('F:/analysis/Cecum_10X+C4/CCI_Level1_New')
bm <- readRDS('F:/analysis/Cecum_10X+C4/data/cecum_IntegrateA_New.rds')
DefaultAssay(bm)<-'RNA'
bm_GF<-subset(bm,subset = mice == "GF")
bm_SPF<-subset(bm,subset = mice == "SPF")
saveRDS(bm_GF,'Level1_GF.rds')
saveRDS(bm_SPF,'Level1_SPF.rds')

write.table(data.frame(bm_GF@active.ident),'GF.bm.celltype.txt',
            row.names=TRUE,col.names=TRUE,sep="\t",quote = FALSE)
write.table(data.frame(bm_SPF@active.ident),'SPF.bm.celltype.txt',
            row.names=TRUE,col.names=TRUE,sep="\t",quote = FALSE)
bm_GF <- readRDS('F:/analysis/Cecum_10X+C4/CCI/Level1_GF.rds')
write.table(GetAssayData(object = bm_GF, slot = "data"),'GF.norm.matrix',
            row.names=TRUE,col.names=TRUE,sep="\t",quote = FALSE)
bm_SPF <- readRDS('F:/analysis/Cecum_10X+C4/CCI/Level1_SPF.rds')
write.table(GetAssayData(object = bm_SPF, slot = "data"),'SPF.norm.matrix',
            row.names=TRUE,col.names=TRUE,sep="\t",quote = FALSE)
#CellChat
library(CellChat)
library(patchwork)

DefaultAssay(bm)<-'RNA'

# Prepare GF data for CelChat analysis
data.input = GetAssayData(object = bm,slot = "data") # normalized data matrix
meta = bm@meta.data # a dataframe with rownames containing cell mata data
cell.use = rownames(meta)[meta$mice == "GF"] # extract the cell names from GF data
data.input = data.input[, cell.use]
meta$labels=bm@active.ident
meta = meta[cell.use, ]
unique(meta$labels) # check the cell labels
cellchat_GF <- createCellChat(object = data.input, meta = meta, group.by = "labels")

cellchat_GF@DB <- CellChatDB.mouse
cellchat_GF <- subsetData(cellchat_GF) # This step is necessary even if using the whole database
cellchat_GF <- identifyOverExpressedGenes(cellchat_GF)
cellchat_GF <- identifyOverExpressedInteractions(cellchat_GF)
cellchat_GF <- projectData(cellchat_GF, PPI.mouse)
cellchat_GF <- computeCommunProb(cellchat_GF)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_GF <- filterCommunication(cellchat_GF, min.cells = 10)
cellchat_GF <- computeCommunProbPathway(cellchat_GF)
cellchat_GF <- aggregateNet(cellchat_GF)
cellchat_GF <- netAnalysis_computeCentrality(cellchat_GF)
saveRDS(cellchat_GF,'cellchat_GF.rds')

# Prepare SPF data for CelChat analysis
data.input = GetAssayData(object = bm, slot = "data") # normalized data matrix
meta = bm@meta.data # a dataframe with rownames containing cell mata data
cell.use = rownames(meta)[meta$mice == "SPF"] # extract the cell names from SPF data
data.input = data.input[, cell.use]
meta$labels=bm@active.ident
meta = meta[cell.use, ]
#meta = data.frame(labels = meta$labels[cell.use], row.names = colnames(data.input)) # manually create a dataframe consisting of the cell labels
unique(meta$labels) # check the cell labels
cellchat_SPF <- createCellChat(object = data.input, meta = meta, group.by = "labels")


#setwd('D:/analysis/Cecum/CCI/cellchat/SPF')
cellchat_SPF@DB <- CellChatDB.mouse
cellchat_SPF <- subsetData(cellchat_SPF) # This step is necessary even if using the whole database
#future::plan("multiprocess", workers = 4) # do parallel
cellchat_SPF <- identifyOverExpressedGenes(cellchat_SPF)
cellchat_SPF <- identifyOverExpressedInteractions(cellchat_SPF)
cellchat_SPF <- projectData(cellchat_SPF, PPI.mouse)
cellchat_SPF <- computeCommunProb(cellchat_SPF)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_SPF <- filterCommunication(cellchat_SPF, min.cells = 10)
cellchat_SPF <- computeCommunProbPathway(cellchat_SPF)
cellchat_SPF <- aggregateNet(cellchat_SPF)
cellchat_SPF <- netAnalysis_computeCentrality(cellchat_SPF)
saveRDS(cellchat_SPF,'cellchat_SPF.rds')

##################################### start analysis
library(CellChat)
library(patchwork)
setwd('F:/analysis/Cecum_10X+C4/CCI_Level1_New')
cellchat_GF <- readRDS('F:/analysis/Cecum_10X+C4/CCI_Level1_New/cellchat_GF.rds')
cellchat_SPF <- readRDS('F:/analysis/Cecum_10X+C4/CCI_Level1_New/cellchat_SPF.rds')

########################### CellChatDB
CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)
########################### CellChat GF

groupSize <- as.numeric(table(cellchat_GF@idents))
pdf("netVisual_circle_GF.pdf",width = 18,height = 9)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_GF@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_GF@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat <- cellchat_GF@net$weight
pdf("netVisual_circle_sparate_GF.pdf",width = 21,height = 14)
par(mfrow = c(4,6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

pathways.show <- c("MPZ") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,14) # a numeric vector. 
#netVisual_aggregate(cellchat_GF, signaling = pathways.show,  vertex.receiver = vertex.receiver)
pdf(paste0("pathway_hierarchy_GF_",pathways.show,".pdf"),width = 16,height = 10)
netVisual_aggregate(cellchat_GF, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout ="hierarchy")
dev.off()
# Circle plot
pdf(paste0("pathway_circle_GF_",pathways.show,".pdf"),width = 8,height = 10)
netVisual_aggregate(cellchat_GF, signaling = pathways.show, layout = "circle")
dev.off()
# Chord diagram
pdf(paste0("pathway_chord_GF_",pathways.show,".pdf"),width = 8,height = 8)
netVisual_aggregate(cellchat_GF, signaling = pathways.show, layout = "chord")
dev.off()
# Heatmap
pdf(paste0("pathway_Reds_GF_",pathways.show,".pdf"),width = 8,height = 5.5)
netVisual_heatmap(cellchat_GF, signaling = pathways.show, color.heatmap = "Reds")
dev.off()

####################### CellChat SPF

groupSize <- as.numeric(table(cellchat_SPF@idents))
pdf("netVisual_circle_SPF.pdf",width = 18,height = 9)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_SPF@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_SPF@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat <- cellchat_SPF@net$weight
pdf("netVisual_circle_sparate_SPF.pdf",width = 21,height = 14)
par(mfrow = c(4,6), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

pathways.show <- c("MPZ") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(from=1,to=12,by=1)# a numeric vector. 
pdf(paste0("pathway_hierarchy_SPF_",pathways.show,".pdf"),width = 16,height = 10)
netVisual_aggregate(cellchat_SPF, signaling = pathways.show,vertex.receiver = vertex.receiver,layout ="hierarchy")
dev.off()
# Circle plot
pdf(paste0("pathway_circle_SPF_",pathways.show,".pdf"),width = 8,height = 10)
netVisual_aggregate(cellchat_SPF, signaling = pathways.show, layout = "circle")
dev.off()
# Chord diagram
pdf(paste0("pathway_chord_SPF_",pathways.show,".pdf"),width = 8,height = 8)
netVisual_aggregate(cellchat_SPF, signaling = pathways.show, layout = "chord")
dev.off()
# Heatmap
pdf(paste0("pathway_Reds_SPF_",pathways.show,".pdf"),width = 8,height = 5.5)
#netVisual_heatmap(cellchat_SPF, signaling = pathways.show, color.heatmap = "Reds",slot.name = "data")
netVisual_heatmap(cellchat_SPF, signaling = pathways.show, color.heatmap = "Reds")
dev.off()

######################### CellChat compare
object.list <- list(GF = cellchat_GF, SPF = cellchat_SPF)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(cellchat,'cellchat_compare.rds')

###### Part 0: Load CellChat object of each dataset and then merge together
cellchat_GF <- readRDS('F:/analysis/Cecum_10X+C4/CCI_Level1_New/cellchat_GF.rds')
cellchat_SPF <- readRDS('F:/analysis/Cecum_10X+C4/CCI_Level1_New/cellchat_SPF.rds')
cellchat <- readRDS('F:/analysis/Cecum_10X+C4/CCI_Level1_New/cellchat_compare.rds')
CellChatDB <- CellChatDB.mouse 

###### Part I: Predict general principles of cell-cell communication
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf("CompareInteraction_bar.pdf",width = 5,height = 3)
gg1+gg2
dev.off()


pdf("CompareInteraction_number&strength.pdf",width = 14,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
pdf("CompareInteraction_heatmap.pdf",width = 11,height = 6)
gg1 + gg2
dev.off()

pdf("CompareInteraction_separate.pdf",width = 11,height = 6)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
pdf("CompareInteraction_Scatter.pdf",width = 11,height = 6)
patchwork::wrap_plots(plots = gg)
dev.off()

pdf("CompareInteraction_Scatter_Fib_Ep.pdf",width = 11,height = 6)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Fibroblast", signaling.exclude = "MPZ")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Epithelial", signaling.exclude = c("MPZ"))
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()

############## Part II: Identify the conserved and context-specific signaling pathways识别保守的和上下文特定的信号通路
#functional
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
# Visualization in 2D-space
pdf("Identify_Pathway_functional_similarity.pdf",width = 8,height = 7)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 2)
dev.off()
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)

#structural
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
pdf("Identify_Pathway_structural_similarity.pdf",width = 8,height = 7)
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
dev.off()

pdf("Identify_Pathway_structural_ZoomIn.pdf",width = 8,height = 7)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()

pdf("RankSimilarity.pdf",width = 5,height = 7)
rankSimilarity(cellchat, type = "functional")
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
pdf("RankSimilarity_compare.pdf",width = 9,height = 12)
gg1 + gg2
dev.off()

library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 7, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 7, height = 20)
pdf("RankSimilarity_Heatmap.pdf",width = 8,height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

#####Identify and visualize the outgoing and incoming communication patterns of target cells—SPF
cellchat_SPF <- readRDS('D:/analysis/Cecum/CCI/cellchat/SPF/cellchat_SPF.rds')
#outgoing
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
setwd('F:/analysis/Cecum_10X+C4/CCI')
p<-selectK(cellchat_SPF, pattern = "outgoing")
pdf("selectK_outgoing_SPF.pdf",width = 8,height = 6)
p
dev.off()

nPatterns = 4
cellchat_SPF <- identifyCommunicationPatterns(cellchat_SPF, pattern = "outgoing",
                                              k = nPatterns,width = 16,height = 20,font.size = 8)

# river plot
p1 <- netAnalysis_river(cellchat_SPF, pattern = "outgoing")
pdf("river plot_outgoing_SPF.pdf",width = 8,height = 14)
p1
dev.off()
# dot plot
p2 <- netAnalysis_dot(cellchat_SPF, pattern = "outgoing")
pdf("dot plot_outgoing_SPF.pdf",width = 16,height = 8)
p2
dev.off()

######incoming
p<-selectK(cellchat_SPF, pattern = "incoming")
pdf("selectK_incoming_SPF.pdf",width = 8,height = 6)
p
dev.off()

nPatterns = 4
cellchat_SPF <- identifyCommunicationPatterns(cellchat_SPF, pattern = "incoming",
                                              k = nPatterns,width = 16,height = 20,font.size = 8)
# river plot
p1 <- netAnalysis_river(cellchat_SPF, pattern = "incoming")
pdf("river plot_incoming_SPF.pdf",width = 8,height = 14)
p1
dev.off()
# dot plot
p2 <- netAnalysis_dot(cellchat_SPF, pattern = "incoming")
pdf("dot plot_incoming_SPF.pdf",width = 16,height = 8)
p2
dev.off()

#####Identify and visualize the outgoing and incoming communication patterns of target cells—GF
cellchat_GF <- readRDS('D:/analysis/Cecum/CCI/cellchat/GF/cellchat_GF.rds')
#outgoing
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)
setwd('F:/analysis/Cecum_10X+C4/CCI')
p<-selectK(cellchat_GF, pattern = "outgoing")
pdf("selectK_outgoing_GF.pdf",width = 8,height = 6)
p
dev.off()

nPatterns = 3
cellchat_GF <- identifyCommunicationPatterns(cellchat_GF, pattern = "outgoing",
                                             k = nPatterns,width = 16,height = 20,font.size = 8)
# river plot
p1 <- netAnalysis_river(cellchat_GF, pattern = "outgoing")
pdf("river plot_outgoing_GF.pdf",width = 8,height = 14)
p1
dev.off()
# dot plot
p2 <- netAnalysis_dot(cellchat_GF, pattern = "outgoing")
pdf("dot plot_outgoing_GF.pdf",width = 16,height = 8)
p2
dev.off()

######incoming
p<-selectK(cellchat_GF, pattern = "incoming")
pdf("selectK_incoming_GF.pdf",width = 8,height = 6)
p
dev.off()

nPatterns = 4
cellchat_GF <- identifyCommunicationPatterns(cellchat_GF, pattern = "incoming",
                                             k = nPatterns,width = 16,height = 20,font.size = 8)
# river plot
p1 <- netAnalysis_river(cellchat_GF, pattern = "incoming")
pdf("river plot_incoming_GF.pdf",width = 8,height = 14)
p1
dev.off()
# dot plot
p2 <- netAnalysis_dot(cellchat_GF, pattern = "incoming")
pdf("dot plot_incoming_GF.pdf",width = 16,height = 8)
p2
dev.off()

############ all_cluster_GO_enrichment_dotplot ##############
setwd('F:/analysis/Cecum_10X+C4/Ep/each cluster GO')
bm <- readRDS('F:/analysis/Cecum_10X+C4/data/cecum_IntegrateA2.rds')
bm<- subset(bm,idents=c('Fibroblast(Gsn high)','Fibroblast(Adamdec1 high)','Mesothelial cell'))
markers <- FindAllMarkers(bm)
#markers <- read.delim('bm.markers.txt', header = TRUE, sep = '\t')
write.table(markers,'bm.level2.Fib.markers.txt',sep='\t',quote=FALSE)
myData <- read.table('bm.Ep.markers.txt', header = T, row.names = 1, sep = "\t")
myData.sort <- myData[myData$p_val_adj < 0.05 &  myData$avg_log2FC > 0.5,]
ref_db <- "org.Mm.eg.db"

clusters <- c('EC(Saa1 high)','EC(Selenbp1 high)','EC(Hmgb2 high)','Goblet','mLTo','Ep_Fibroblast')
#unique(myData$cluster)
get_go <- function(genelist){
  suppressMessages(library(clusterProfiler))
  suppressMessages(library(AnnotationDbi))
  suppressMessages(library(org.Hs.eg.db))
  suppressMessages(library(org.Mm.eg.db))
  suppressMessages(library(DOSE))
  suppressMessages(library(topGO))
  
  eg <- suppressMessages(bitr(genelist, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb=ref_db))
  genelist <- eg$ENTREZID
  genelist  <- genelist[!duplicated(genelist)]
  go.BP <- enrichGO(genelist, OrgDb = ref_db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = 'ENTREZID')
  return(list(go.BP, genelist))
}

gc <- list()

cl_name <- c()
for (i in clusters) {
  cluster.geneList <- myData.sort$gene[myData.sort$cluster == i]
  try(go <- suppressWarnings(get_go(cluster.geneList)),silent = TRUE)
  if ('NULL' %in% go) {
    next
  } else{ 
    go_enrichment <- suppressMessages(go[[1]]@result)
    cl_name <- c(cl_name, i)
    gc <- c(gc,list(go[[2]]))
  }
}
names(gc) <- cl_name
#define a function to adjust the y label of dotplot
get_wraper <- function(width) {
  function(x) {
    lapply(strwrap(x, width = width, simplify = FALSE), paste, collapse="\n")
  }
}
#
x <- compareCluster(gc,fun='enrichGO',OrgDb=ref_db,ont="BP")
p <- dotplot(x)
p <- p+ aes(color=-log10(p.adjust))+
  #scale_color_continuous(low="darkred", high="red", guide=guide_colorbar(reverse=TRUE))+
  scale_color_continuous(low="red", high="green", guide=guide_colorbar(reverse=TRUE))+
  theme(axis.text.x = element_text(size = 12,angle = 30,hjust = 1),
        axis.text.y = element_text(size = 12)) + 
  labs(color=expression(p.adjust),size="Gene Ratio",x="",y="")+
  scale_size(range=c(1, 8))+
  scale_y_discrete(labels= get_wraper(60))+
  scale_x_discrete(labels= get_wraper(50))
  

# output the dotplot
pdf(file="EC3_GO_enrichment_dotplot.pdf",height = 8, width = 8)
print(p)
while (!is.null(dev.list()))  dev.off()
############ scMetabolism ##########
.libPaths(c("/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R411/lib/R/library", .libPaths()))
library(homologene)
library(scMetabolism)
library(ggplot2)
library(rsvd)

#' scMetabolism
#' scMetabolism
#' @param obj
#' @keywords scMetabolism
#' @examples
#' sc.metabolism.Seurat()
#' @export sc.metabolism.Seurat
sc.metabolism.Seurat <- function(obj, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG",threshold=0.01,projection_genes="threshold") {
  countexp<-obj@assays$RNA@counts
  countexp<-data.frame(as.matrix(countexp))
  signatures_KEGG_metab <- system.file("data", "mouse_KEGG_metabolism_nc.gmt", package = "scMetabolism")
  signatures_REACTOME_metab <- system.file("data", "mouse_REACTOME_metabolism.gmt", package = "scMetabolism")
  if (metabolism.type == "KEGG")  {gmtFile<-signatures_KEGG_metab; cat("Your choice is: KEGG\n")}
  if (metabolism.type == "REACTOME")  {gmtFile<-signatures_REACTOME_metab; cat("Your choice is: REACTOME\n")}
  #imputation
  if (imputation == F) {
    countexp2<-countexp
  }
  if (imputation == T) {
    cat("Start imputation...\n")
    #Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588
    #Github: https://github.com/KlugerLab/ALRA
    cat("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588 \n")
    result.completed <- alra(as.matrix(countexp))
    countexp2 <- result.completed[[3]]; row.names(countexp2) <- row.names(countexp)
  }
  #signature method
  cat("Start quantify the metabolism activity...\n")
  #VISION
  if (method == "VISION") {
    library(VISION)
    n.umi <- colSums(countexp2)
    scaled_counts <- t(t(countexp2) / n.umi) * median(n.umi)
    vis <- Vision(scaled_counts, signatures = gmtFile,threshold = threshold,projection_genes=c(projection_genes))
    options(mc.cores = ncores)
    vis <- analyze(vis)
    signature_exp<-data.frame(t(vis@SigScores))
  }
  #AUCell
  if (method == "AUCell") {
    library(AUCell)
    library(GSEABase)
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), nCores=ncores, plotStats=F) #rank
    geneSets <- getGmt(gmtFile) #signature read
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) #calc
    signature_exp <- data.frame(getAUC(cells_AUC))
  }
  #ssGSEA
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile) #signature read
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method=c("ssgsea"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)
  }
  #GSVA
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile) #signature read
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method=c("gsva"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)
  }
  #obj@assays$METABOLISM$score<-signature_exp
  signature_exp
}

#run
obj=readRDS("/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/chaitailiang/10X_C4/data/cecum_Epithelial_IntegrateA_scmeta.rds")
kegg=sc.metabolism.Seurat(obj, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG",threshold=0.01,projection_genes="threshold")
reactome=sc.metabolism.Seurat(obj, method = "VISION", imputation = F, ncores = 2, metabolism.type = "REACTOME",threshold=0.01,projection_genes="threshold")


outdir="/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/chaitailiang/10X_C4/Result/scMetabolism_Ep"
write.table(kegg,file=paste(outdir,"/kegg.xls",sep=""),row.names=rownames(kegg),col.names=colnames(kegg),sep="\t",quote=F)
write.table(reactome,file=paste(outdir,"/reactome.xls",sep=""),row.names=rownames(reactome),col.names=colnames(reactome),sep="\t",quote=F)


############ scMetabolism_plot ##########
library(Seurat)
library(dplyr)
library(patchwork)
library(SeuratDisk)
library(DoubletFinder)
library(pals)
library(cowplot)
library(ggplot2)


DotPlot.vs.scmetabolism <- function (obj, pathway, phenotype, norm = "y") {
  input.norm = norm
  input.pathway <- pathway
  input.parameter <- phenotype
  metadata <- obj@meta.data
  metabolism.matrix <- obj@assays$METABOLISM$score
  metadata[, input.parameter] <- as.character(metadata[, input.parameter])
  metabolism.matrix_sub <- t(metabolism.matrix[input.pathway, ])
  gg_table <- c()
  for (i in 1:length(input.pathway)) { #???
    gg_table <- rbind(gg_table, cbind(metadata[, input.parameter], input.pathway[i], metabolism.matrix_sub[, i]))
  }
  gg_table <- data.frame(gg_table)
  gg_table_median <- c()
  input.group.x <- unique(as.character(gg_table[, 1]))
  input.group.y <- unique(as.character(gg_table[, 2]))
  for (x in 1:length(input.group.x)) {
    for (y in 1:length(input.group.y)) {
      gg_table_sub <- subset(gg_table, gg_table[, 1] == input.group.x[x] & gg_table[, 2] == input.group.y[y])
      gg_table_median <- rbind(gg_table_median, cbind(input.group.x[x], input.group.y[y], 
                                                      median(as.numeric(as.character(gg_table_sub[, 3])))))
    }
  }
  gg_table_median <- data.frame(gg_table_median)
  gg_table_median[, 3] <- as.numeric(as.character(gg_table_median[, 
                                                                  3]))
  gg_table_median_norm <- c()
  input.group.x <- unique(as.character(gg_table[, 1]))
  input.group.y <- unique(as.character(gg_table[, 2]))
  range01 <- function(x) {
    (x - min(x))/(max(x) - min(x))
  }
  if (input.norm == "y") 
    for (y in 1:length(input.group.y)) {
      gg_table_median_sub <- subset(gg_table_median, gg_table_median[, 2] == input.group.y[y])
      norm_value <- range01(as.numeric(as.character(gg_table_median_sub[, 3])))
      gg_table_median_sub[, 3] <- norm_value
      gg_table_median_norm <- rbind(gg_table_median_norm, 
                                    gg_table_median_sub)
    }
  if (input.norm == "x") 
    for (x in 1:length(input.group.x)) {
      gg_table_median_sub <- subset(gg_table_median, gg_table_median[, 1] == input.group.x[x])
      norm_value <- range01(as.numeric(as.character(gg_table_median_sub[, 3])))
      gg_table_median_sub[, 3] <- norm_value
      gg_table_median_norm <- rbind(gg_table_median_norm, 
                                    gg_table_median_sub)
    }
  if (input.norm == "na") 
    gg_table_median_norm <- gg_table_median
  gg_table_median_norm <- data.frame(gg_table_median_norm)
  gg_table_median_norm[, 3] <- as.numeric(as.character(gg_table_median_norm[, 3]))
  library(wesanderson)
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  # order x axis
  gg_table_median_norm[, 1] <- factor(gg_table_median_norm[, 1],levels =unique(gg_table_median_norm[, 1]))
  ggplot(data = gg_table_median_norm, aes(x = gg_table_median_norm[, 1], y = gg_table_median_norm[, 2],
                                          color = gg_table_median_norm[, 3])) + 
    geom_point(data = gg_table_median_norm, aes(size = gg_table_median_norm[, 3])) + 
    ylab("Metabolic Pathway") + xlab("Celltypes") + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                       panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
    scale_color_gradientn(colours = pal) + labs(color = "Value", size = "Value") + NULL
}
#
DoMelt <- function(object){
  object@meta.data <- object@meta.data %>% arrange(celltypes) %>%
    mutate(scMeta = paste(mice,tissue,celltypes,sep = "_"))
  cellid <- rownames(object@meta.data)
  #object@assays$METABOLISM$score <- object@assays$METABOLISM$score[,cellid]
  return(object)
}
#haha <- colon_meta@assays$METABOLISM$score

# main #
setwd('F:/analysis/Cecum_10X+C4/scMeta')
colon_meta <- readRDS("F:/analysis/Cecum_10X+C4/scMeta/cecum_EP_kegg_scMeta.rds")
colon_meta<-RenameIdents(colon_meta,'Goblet'='Epithelial','EC(Hmgb2 high)'='Epithelial','EEC'='Epithelial','EC(Saa1 high)'='Epithelial',
                         'EC(Selenbp1 high)'='Epithelial','EC(Rpl3 high)'='Epithelial','mLTo'='Epithelial','Tuft'='Epithelial','EC(Dying)'='Epithelial','Ep_Fibroblast'='Epithelial','EC(Reg3b high)'='Epithelial')
colon_meta$celltypes<-colon_meta@active.ident
unique(colon_meta@meta.data$celltypes)
colon_meta <- DoMelt(colon_meta)
unique(colon_meta@meta.data$scMeta)
pathway <- rownames(colon_meta@assays$METABOLISM$score)
DotPlot.vs.scmetabolism(obj = colon_meta, pathway = pathway, phenotype = 'scMeta', norm = "y")
