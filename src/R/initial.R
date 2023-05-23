library(qs)
library(ggplot2)
library(magrittr)
library(dplyr)
install.packages("viridis")
library(viridis)
library(Seurat)
library(readr)
require(cowplot)

beadplot <- function(sb.data){
  ggplot(sb.data,aes(x=x_um,y=y_um,col=umi)) +
    geom_point(size=1) +
    coord_fixed() +
    theme_classic() +
    scale_color_viridis(trans="log", option="B") +
    labs(x="x (um)", y="y (um)") +
    ggtitle("Spatial Barcode nUMIs (CB grouped)")
}

# v1 format data - spatialdata is in qs format
obj1 = '/Users/mraj/Desktop/work/data/slide_tag/v1/seurat.qs'
obj2 = '/Users/mraj/Desktop/work/data/slide_tag/v1/spatialdata.qs'
obj <- qread(obj1)
data <- qread(obj2)
# v1 pre plot data procesing
sb.data = data %>% group_by(sb) %>% summarize(unique.cb=n(),umi=sum(umi),reads=sum(reads),x_um=mean(x_um),y_um=mean(y_um)) %>% {.[order(.$umi),]}

# v2 format data - spatial data is in csv format
o2a = '/Users/mraj/Desktop/work/data/slide_tag/v2/1915seurat.qs'
o2b = '/Users/mraj/Desktop/work/data/slide_tag/v2/1915coords.csv'
obj <- qread(o2a)
data <- read_csv(o2b)
# modified version of preplot data processing for v2 data
sb.data = data %>% group_by(sb) %>% summarize(unique.cb=n(),umi=sum(umi),reads=sum(reads),x_um=mean(x),y_um=mean(y)) %>% {.[order(.$umi),]}


# v2 format data - spatial data is in csv format
o3a = '/Users/mraj/Desktop/work/data/slide_tag/v2/6214seurat.qs'
o3b = '/Users/mraj/Desktop/work/data/slide_tag/v2/6214coords.csv'
obj <- qread(o3a)
data <- read_csv(o3b)
# modified version of preplot data processing for v2 data
sb.data = data %>% group_by(sb) %>% summarize(unique.cb=n(),umi=sum(umi),reads=sum(reads),x_um=mean(x),y_um=mean(y)) %>% {.[order(.$umi),]}

beadplot(sb.data)
DimPlot(obj, reduction="spatial")


# Sample code for plotting with no margins
# https://stackoverflow.com/questions/31254533/when-using-ggplot-in-r-how-do-i-remove-margins-surrounding-the-plot-area
qplot(1:10, (1:10)^2, geom='point', color=I('red'), size=I(1)) + theme_nothing() + 
  scale_x_continuous(expand=c(0,0), limits=c(1,50)) +
  scale_y_continuous(expand=c(0,0), limits=c(1,50)) +
  labs(x = NULL, y = NULL)
impath = '/Users/mraj/Desktop/tmp.png' 
ggsave(impath, dpi=75)

# obj@reductions$spatial@cell.embeddings or obj[['spatial']]@cell.embeddings
# https://satijalab.org/seurat/archive/v3.0/dim_reduction_vignette.html
# obj$seurat_clusters
# plotting by color - http://www.sthda.com/english/wiki/qplot-quick-plot-with-ggplot2-r-software-and-data-visualization
# filtering clusters https://www.biostars.org/p/483016/
