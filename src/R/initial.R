library(qs)
library(ggplot2)
library(magrittr)
library(dplyr)
install.packages("viridis")
library(viridis)
library(Seurat)

obj1 = '/Users/mraj/Desktop/work/data/slide_tag/v1/seurat.qs'
obj2 = '/Users/mraj/Desktop/work/data/slide_tag/v1/spatialdata.qs'

obj <- qread(obj1)
data <- qread(obj2)

beadplot <- function(sb.data){
  ggplot(sb.data,aes(x=x_um,y=y_um,col=umi)) +
    geom_point(size=1) +
    coord_fixed() +
    theme_classic() +
    scale_color_viridis(trans="log", option="B") +
    labs(x="x (um)", y="y (um)") +
    ggtitle("Spatial Barcode nUMIs (CB grouped)")
}


sb.data = data %>% group_by(sb) %>% summarize(unique.cb=n(),umi=sum(umi),reads=sum(reads),x_um=mean(x_um),y_um=mean(y_um)) %>% {.[order(.$umi),]}
beadplot(sb.data)
