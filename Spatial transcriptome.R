############ gem_to_rds ###########
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))
suppressMessages(library(SeuratDisk))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(Matrix))
suppressMessages(library(rjson))
suppressMessages(library(RColorBrewer))
setwd('F:/analysis/ST/gem to rds/RDS')

gem_to_seuratObject <- function(gem, prefix, binsize){
  #' group counts into bins
  data=fread(file=gem)
  colnames(data)=c("geneID","x","y","MIDCounts")
  
  
  if ('MIDCounts' %in% colnames(data)) {
    
    data <- data[, .(counts=sum(MIDCounts)), by = .(geneID, x, y)]
  } else {
    data <- data[, .(counts=sum(UMICount)), by = .(geneID, x, y)]
  }
  
  #' create sparse matrix from stereo
  data$cell <- paste0(prefix, ':', data$x, '_', data$y)
  data$geneIdx <- match(data$geneID, unique(data$geneID))
  data$cellIdx <- match(data$cell, unique(data$cell))
  
  mat <- sparseMatrix(i = data$geneIdx, j = data$cellIdx, x = data$counts, 
                      dimnames = list(unique(data$geneID), unique(data$cell)))
  
  cell_coords <- unique(data[, c('cell', 'x', 'y')])
  
  rownames(cell_coords) <- cell_coords$cell
  
  seurat_spatialObj <- CreateSeuratObject(counts = mat, project = 'Stereo', assay = 'Spatial', 
                                          names.delim = ':', meta.data = cell_coords)
  
  
  #' create pseudo image
  cell_coords$x <- cell_coords$x - min(cell_coords$x) + 1
  cell_coords$y <- cell_coords$y - min(cell_coords$y) + 1
  
  tissue_lowres_image <- matrix(1, max(cell_coords$y), max(cell_coords$x))
  
  tissue_positions_list <- data.frame(row.names = cell_coords$cell,
                                      tissue = 1,
                                      row = cell_coords$y, col = cell_coords$x,
                                      imagerow = cell_coords$y, imagecol = cell_coords$x)
  
  
  scalefactors_json <- toJSON(list(fiducial_diameter_fullres = binsize,
                                   tissue_hires_scalef = 1,
                                   tissue_lowres_scalef = 1))
  #' function to create image object
  generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE){
    if (filter.matrix) {
      tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
    }
    
    unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
    
    spot.radius <- unnormalized.radius / max(dim(x = image))
    
    return(new(Class = 'VisiumV1', 
               image = image, 
               scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
                                            fiducial = scale.factors$fiducial_diameter_fullres, 
                                            hires = scale.factors$tissue_hires_scalef, 
                                            lowres = scale.factors$tissue_lowres_scalef), 
               coordinates = tissue.positions, 
               spot.radius = spot.radius))
  }
  
  spatialObj <- generate_spatialObj(image = tissue_lowres_image, 
                                    scale.factors = fromJSON(scalefactors_json), 
                                    tissue.positions = tissue_positions_list)
 
  
  #' import image into seurat object
  spatialObj <- spatialObj[Cells(x = seurat_spatialObj)]
  DefaultAssay(spatialObj) <- 'Spatial'
  
  seurat_spatialObj[['slice1']] <- spatialObj
  rm("spatialObj")
  rm("data")
  rm("mat")
  
  return(seurat_spatialObj)
}

gem='F:/analysis/ST/gem to rds/gemChip/GF_CC_3.gem'
prefix='GF_CC_3'
binsize=50

seurat_spatialObj=gem_to_seuratObject(gem,prefix,binsize)
saveRDS(seurat_spatialObj,paste0(prefix,'.rds'))



############ BayesSpace ##########
library(BayesSpace)
library(SingleCellExperiment)
library(ggplot2)
library(Seurat)
library(Matrix)
setwd('F:/analysis/ST')
#
Stdata=readRDS('F:/analysis/ST/GFCUT_filt_norm_cluster_seuratObject.rds')
df=data.frame(X=Stdata@meta.data$x,Y=Stdata@meta.data$y)
df$spot_id=colnames(Stdata@assays$Spatial@counts)
df<-df[, c("spot_id","X","Y")]
write.table(df,'00.bin.info_GFCUT.tsv',quote = FALSE,col.names = FALSE,row.names = FALSE,sep='\t')
#
counts = as(Stdata@assays$Spatial@counts, "dgCMatrix")
colData = read.table('00.bin.info_GFCUT.tsv', stringsAsFactors=F, sep="\t",header = F)
names(colData) = c( "spot", "row", "col" )
rownames(colData) = colData$spot
genes = rownames(counts)
rowData = data.frame(genes)
sce <- SingleCellExperiment(assays=list(counts=counts),
                            rowData=rowData,
                            colData=colData)
set.seed(102)
sce <- spatialPreprocess(sce, platform="ST",
                         n.PCs=30, n.HVGs=3000, log.normalize=TRUE)
saveRDS(sce,'BayesSpa_GFCUT.rds')

#
clusternum=4 
nrep.num=10000 

sce <- readRDS('BayesSpa_GFCUT.rds')
sce <- spatialCluster(sce, q=clusternum, platform="ST", d=30,
                      init.method="mclust", model="t", gamma=2,
                      nrep=nrep.num, burn.in=200,
                      save.chain=TRUE)
out.df = colData(sce)
write.table(out.df, file=paste0('02.cluster_q',clusternum,'.info.txt'), quote=F, row.names=T, col.names=F, sep="\t")

#define a function for color palette
ColorPalette <- function(number) {
  print(number)
  if (number < 25) {
    colorScheme = c("#999999","#FF0099","#E69F00","#56B4E9","#009E73","#F0E442",
                    "#0072B2","#D55E00","#CC79A7","#990000","#9900cc","#66FF66",
                    "#663300","#0000FF","#CC0033","#FF0000","#000099","#660066",
                    "#333333","#FFCCCC","#993366","#33CC33","#000099","#CC9900")
  } else {
    suppressMessages(library(RColorBrewer))
    qual_col_pals <-
      brewer.pal.info[brewer.pal.info$category == 'qual', ]
    colorScheme <-
      unique(unlist(mapply(
        brewer.pal,
        qual_col_pals$maxcolors,
        rownames(qual_col_pals)
      )))
    return (colorScheme)
  }
}
#' define a function to visualize the spatial cluster
df <- read.table(file=paste0('02.cluster_q',clusternum,'.info.txt'), header = F, row.names = 1)
colnames(df) <- c('bin_id', 'x', 'y', 'B', 'cluster1', 'cluster2')
df$y <- -(df$y)  ## you can comment this for your need
df$cluster2 <- as.character(df$cluster2)
colorPalette <- ColorPalette(length(unique(df$cluster2)))

point_size=0.4
p0 <-ggplot(df, aes(x = x, y = y, color = cluster2)) +
  geom_point(shape = 19, size = point_size) +  ## you can modify the point size to avoid overlap
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    legend.position = 'right'
  ) +
  coord_fixed()  +
  scale_color_manual(values = colorPalette) +
  guides(colour = guide_legend(
    override.aes = list(size = 3),
    nrow = 15,
    title = 'Clusters'
  ))  +
  theme_void()

pdf(file=paste0('03.BayesSpace.spatial.cluster_',clusternum,'.pdf'), width = 8, height = 6)
p0
dev.off()

############ Spatial distribution of marker genes  ############
setwd('F:/analysis/ST/Marker')
Stdata=readRDS('F:/analysis/ST/GF_filt_norm_cluster_seuratObject.rds')
Stdata@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = Stdata@meta.data[, c('y', 'x')]
)
#
library(gridExtra)
target_gene <- c('Slc10a2','Slc51a','Slc51b')
title <- c('Slc10a2','Slc51a','Slc51b')
options(repr.plot.width =7.5,repr.plot.height = 5)#*(length(target_gene)%/%4))
for(i in 1:length(target_gene))
{
  plot <-paste("P",i,sep="")
  run  <- paste (plot," <- SpatialFeaturePlot(Stdata, features = '",target_gene[i],"', alpha = c(0.7, 1),ncol=5,stroke=NA,images='image',pt.size.factor= 1)+ labs(title='",title[i],"')+ coord_flip() +scale_y_reverse() + theme(legend.text=element_text(size=7),legend.title=element_text(size=7),legend.position = 'right',legend.key.width=unit(0.3,'cm'),legend.key.height=unit(0.3,'cm'),plot.title=element_text(size=7),legend.margin=margin(0,0,0,0),legend.box.margin=margin(-5,-5,0,-5))",sep="")
  eval(parse(text = run))
}
pdf("bile acide.pdf",width=7.5,height=5)
grid.arrange(P1,P2,P3,ncol=5,nrow=4)
dev.off()
############ Addmodulescore ###########
library(Seurat)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(ggpubr)
setwd('F:/analysis/ST/Addmodulescore')
GF=readRDS('F:/analysis/ST/data/GF_filt_norm_cluster_seuratObject.rds')
SPF=readRDS('F:/analysis/ST/data/SPF3_filt_norm_cluster_seuratObject.rds')

GF@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = GF@meta.data[, c('y', 'x')]
)
SPF@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = SPF@meta.data[, c('y', 'x')]
)
combined <- merge(GF, y = SPF, project = "combined")
combined <- NormalizeData(combined , verbose = FALSE)

#Addmodulescore

#Response to Unfolded Proteins
gene <- list(c("Atf3","Atf4","Bag3","Hsp90aa1","Hspa5","Hspa9","Hspb1","Ppp1r15a","Thbs1"))
dfscore_list <- AddModuleScore(combined, features=gene, nbin = 24, ctrl = 100, name = "Unfolded", seed = 1)
head(dfscore_list)
unique(dfscore_list$orig.ident)

p<-VlnPlot(dfscore_list,features = 'Unfolded1',group.by='orig.ident',pt.size = 0) & stat_compare_means(label.y =5)&ylim(-5,5) & geom_boxplot(width=.1,col="black",fill="white",outlier.shape=NA)

pdf('Unfolded_gene.pdf',width = 4,height = 5)
print(p)
dev.off()

