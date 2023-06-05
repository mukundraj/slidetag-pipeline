
library(qs)
library(ggplot2)
library(magrittr)
library(dplyr)
# install.packages("viridis")
library(viridis)
library(Seurat)
library(readr)
require(cowplot)
library(readbitmap)

inpdir <- snakemake@input[['raw_dir']]
opdir <- snakemake@output[['stag_imgs_dir']]
opdir_coords <- snakemake@output[['stag_coords_dir']]

print("Creating output directory:")
print(opdir)
dir.create(file.path(opdir), recursive = TRUE)

print("Creating output directory for coords:")
print(opdir_coords)
dir.create(file.path(opdir_coords), recursive = TRUE)

data_folder_path <- file.path(inpdir, snakemake@config['dataname'])
files_qs <- list.files(path=data_folder_path, pattern = "\\qs$")

if (length(files_qs)==1){

  seurat_file <- file.path(data_folder_path, files_qs[[1]])

  sobj <- qread(seurat_file)
  # sobj <- qread('/Users/mraj/Desktop/work/data/slide_tag/v3/sample01.qs')
  cells <- unique(sobj$type)
  cellembs <- as.data.frame(sobj@reductions$spatial@cell.embeddings)
  cellembs_rnames <- rownames(cellembs)
  typ <- character(length(cellembs_rnames))
  cellembs$typ <- typ
  for (rid in 1:nrow(cellembs)){
    cellembs$typ[rid] = sobj$type[cellembs_rnames[rid]]
  }
  cellembs_nona <- na.omit(cellembs)

  xlims <- c(min(cellembs_nona$s_1), max(cellembs_nona$s_1))
  ylims <- c(min(cellembs_nona$s_2), max(cellembs_nona$s_2))

  # print xlims and ylims
  print('xlims')
  print(xlims)
  print('ylims')
  print(ylims)

  qplot(cellembs_nona$s_1, cellembs_nona$s_2, geom='point', color=I('red'), size=I(0.1)) + theme_nothing() +
    scale_x_continuous(expand=c(0,0), limits=xlims) +
    scale_y_continuous(expand=c(0,0), limits=ylims) +
    labs(x = NULL, y = NULL)&coord_fixed()
  impath <- file.path(opdir, paste0('bead_plot.png'))
  print(impath)
  ggsave(impath, dpi=75)

  # get img dims of image
  img_dims <- dim(read.bitmap(impath))
  print('img_dims')
  print(img_dims)


  for (cell in cells){
    #print(sum(cellembs_nona$typ==cell))
    sub_cellembs <- cellembs_nona[cellembs_nona$typ==cell,]
    print(dim(sub_cellembs))

    qplot(sub_cellembs$s_1, sub_cellembs$s_2, geom='point', color=I('red'), size=I(0.1)) + theme_nothing() +
      scale_x_continuous(expand=c(0,0), limits=xlims) +
      scale_y_continuous(expand=c(0,0), limits=ylims) +
      labs(x = NULL, y = NULL)&coord_fixed()

    # save plot in output dir
    impath <- file.path(opdir, paste0('stag_', cell,'.png'))
    print(impath)
  ggsave(impath, dpi=75)

    # normalize coords using extents before writing and normalize
    sub_cellembs$s_1 <- (sub_cellembs$s_1 - xlims[1])/(xlims[2]-xlims[1]) * img_dims[1]
    sub_cellembs$s_2 <- (sub_cellembs$s_2 - ylims[1])/(ylims[2]-ylims[1]) * img_dims[2]

    coords_file <- file.path(opdir_coords, paste0('stag_', cell,'_coords.csv'))

    write.table(sub_cellembs[,c('s_1', 's_2')], coords_file, row.names = FALSE, col.names = FALSE, sep = ",")

  }


  ###




  # normalize coords using extents before writing and normalize
  cellembs_nona$s_1 <- (cellembs_nona$s_1 - xlims[1])/(xlims[2]-xlims[1]) * img_dims[1]
  cellembs_nona$s_2 <- (cellembs_nona$s_2 - ylims[1])/(ylims[2]-ylims[1]) * img_dims[2]


  # save coords in opdir_coords folder in csv file named bead_coords.csv and only x and y columns
  coords_file <- file.path(opdir_coords, 'bead_coords.csv')

  write.table(cellembs_nona[,c('s_1', 's_2')], coords_file, row.names = FALSE, col.names = FALSE, sep = ",")



}


