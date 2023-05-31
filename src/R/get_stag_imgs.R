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
opdir_coords <- snakemake@output[['stag_coords_dir']]

print("Creating output directory:")
print(opdir)
dir.create(file.path(opdir), recursive = TRUE)

print("Creating output directory for coords:")
print(opdir_coords)
dir.create(file.path(opdir_coords), recursive = TRUE)

# if coords csv file exists in input dir and if seurat.qs file exists in input dir, create plots and save in output dir, also store coords in interim/stags_coords folder
data_folder_path <- file.path(inpdir, snakemake@config['dataname'])
files <- list.files(path=data_folder_path, pattern = "\\.csv$")
print(files[[1]])
print(length(files))
files_qs <- list.files(path=data_folder_path, pattern = "\\qs$")

if (length(files)==1 && length(files_qs)==1){
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

  # save coords in opdir_coords folder in csv file named bead_coords.csv and only x and y columns
  coords_file <- file.path(opdir_coords, 'bead_coords.csv')
  write.table(sb.data[,c('x_um', 'y_um', 'umi')], coords_file, row.names = FALSE, col.names = FALSE, sep = ",")


  print ("Creating plots for seurat object")

  seurat_file <- file.path(data_folder_path, files_qs[[1]])
  obj  <- qread(seurat_file)
  b <- table(obj$seurat_clusters)
  c <- sort(b, decreasing=TRUE)
  num_clusters <- min(c(length(names(c)), 10)) # select top 5 or less number of cells to generate imgs
  print(c('num_clusters',num_clusters))
  print(c)
  
  # loop over num_clusters
  for (i in 1:num_clusters){

    locidx <- strtoi(names(c)[i])
    inds <- obj$seurat_clusters==locidx
    print(c('locidx',locidx, 'num cells', sum(inds, na.rm=TRUE)))
    d <- obj@reductions$spatial@cell.embeddings
    z <- na.omit(as.data.frame(d[inds,]))

    # keep extents same as bead plot

    # plot using qplot
    qplot(z$s_1, z$s_2, geom='point', color=I('red'), size=I(1)) + theme_nothing() + 
      scale_x_continuous(expand=c(0,0), limits=xlims) +
      scale_y_continuous(expand=c(0,0), limits=ylims) +
      labs(x = NULL, y = NULL)

    # save plot in output dir
    impath <- file.path(opdir, paste0('stag_', locidx,'.png'))
    ggsave(impath, dpi=75)

    # save  coords in opdir_coords folder
    coords_file <- file.path(opdir_coords, paste0('stag_', locidx,'_coords.csv'))
    write.table(z[,c('s_1', 's_2')], coords_file, row.names = FALSE, col.names = FALSE, sep = ",")

  }   
  
}


