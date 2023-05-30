library(qs)
library(ggplot2)
library(magrittr)
library(dplyr)
# install.packages("viridis")
library(viridis)
library(Seurat)
library(readr)
require(cowplot)

inpdir <- snakemake@input[['raw_dir']]
opdir <- snakemake@output[['stag_imgs_dir']]

print("Creating output directory:")
print(opdir)
dir.create(file.path(opdir), recursive = TRUE)


# if coords csv file exists in input dir, create plot and save in output dir
data_folder_path <- file.path(inpdir, snakemake@config['dataname'])
files <- list.files(path=data_folder_path, pattern = "\\.csv$")
print(files[[1]])
print(length(files))

if (length(files)==1){
  coords_file <- file.path(data_folder_path, files[[1]])
  print(coords_file)

  # v2 format data - spatial data is in csv format
  o3b <- coords_file
  data <- read_csv(o3b)
  # modified version of preplot data processing for v2 data
  sb.data = data %>% group_by(sb) %>% summarize(unique.cb=n(),umi=sum(umi),reads=sum(reads),x_um=mean(x),y_um=mean(y)) %>% {.[order(.$umi),]}

  # for extents
  xlims <- c(min(sb.data$x_um), max(sb.data$x_um))
  ylims <- c(min(sb.data$y_um), max(sb.data$y_um))

  qplot(sb.data$x_um, sb.data$y_um, geom='point', color=I('red'), size=I(0.1)) + theme_nothing() + 
    scale_x_continuous(expand=c(0,0), limits=xlims) +
    scale_y_continuous(expand=c(0,0), limits=ylims) +
    labs(x = NULL, y = NULL)

  # impath = '/Users/mraj/Desktop/tmp.png' 
  impath <- file.path(opdir, 'bead_plot.png')
  ggsave(impath, dpi=75)


}


# if seurat.qs file exists in input dir, create plots and save in output dir

