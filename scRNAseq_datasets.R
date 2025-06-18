# Reading in multiple GEO dataset files and creating a Seurat object
# https://github.com/crazyhottommy/compbio_tutorials/blob/main/scripts/04_create_seurat_object_from_GEO.Rmd


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/GSE116256")

# download data  
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116256
# 43 samples: download the data, create a Seurat object for each file and then merge them into a single Seurat objet

# bash code
# cd /Users/kristinvandenham/kmvanden/RStudio/GSE116256
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116256/suppl/GSE116256_RAW.tar
# tar xvf GSE116256_RAW.tar
# rm  GSE116256_RAW.tar

# load libraries
library(stringr)
library(dplyr)
library(ggplot2)
library(Seurat)
library(purrr)
library(readr)
library(harmony)
library(scCustomize)
library(SeuratDisk)
library(Matrix) # for sparse matrices


# function to read in the count files as matrices
read_counts <- function(file){
  x <- read_tsv(file) # read in the files
  x <- as.data.frame(x) # convert into a data.frame
  genes <- x$Gene # save Gene column
  x <- x[, -1] # remove first column (Gene)
  rownames(x) <- genes # make Genes the rownames 
  return(as.matrix(x)) # return a matrix
}

# list the names of the counts files
counts_files <- list.files(full.names = TRUE, pattern = "*dem.txt.gz")
# extract the last part of the file name for the sample names 
samples <- map_chr(counts_files, basename) 
# remove "dem.txt.gz" from the sample names
samples <- str_replace(samples, "(GSM[0-9]+_.+).dem.txt.gz", "\\1")
# add the sample names to count_files
names(counts_files) <- samples

# read in the counts files as matrices using the created function and purrr
counts <- purrr::map(counts_files, read_counts)
counts$`GSM3587923_AML1012-D0`


# function to read in the metadata files
read_meta <- function(file){
  y <- read_tsv(file) # read in the files
  y <- as.data.frame(y) # convert the files into data.frames
  cells <- y$Cell # save the Cell column
  y <- y[,-1] # remove the first column (Cell)
  rownames(y) <- cells # add Cell as rownames
  return(y)
}

# list the names of the metadata files
meta_files <- list.files(full.names = TRUE, pattern = "*anno.txt.gz")
# extract the last part of the file name for the sample names 
meta_names <- map_chr(meta_files, basename)
# remove "anno.txt.gz" from the sample names
meta_names <- str_replace(meta_names, "(GSM[0-9]+_.+).anno.txt.gz", "\\1")

# check that the counts files and metadata files are in the same order 
sum(str_replace(samples, "GSM[0-9]+_(.+)", "\\1") == str_replace(meta_names, "GSM[0-9]+_(.+)", "\\1")) # 43

# add the meta names to meta_files
names(meta_files) <- meta_names

# read in the metadata files using the created function and purrr
meta <- purrr::map(meta_files, read_meta)
meta$`GSM3587924_AML1012-D0`


### create the Seurat object
# .x refers to every element in the counts list and .y refers to every element in the meta list
objs <- purrr::map2(counts, meta,  
                   ~CreateSeuratObject(counts = as(.x, "sparseMatrix"), 
                                       meta.data = .y))
objs # 43 Seurat objects

# merge the 43 Seurat objects into a single Seurat object 
merged_seurat <- purrr::reduce(objs, function(x,y) {merge(x,y)})
merged_seurat


### remove the large files that were created to free up memory
rm(counts)
rm(objs)
rm(meta)
gc()


### follow the standard Seurat preprocessing workflow
merged_seurat <- merged_seurat %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures( selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony(group.by.vars = "orig.ident", dims.use = 1:30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.6)


### visualization
DimPlot_scCustom(seurat_object = merged_seurat)
head(merged_seurat@meta.data)


sessionInfo()
# R version 4.5.0 (2025-04-11)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.5
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/Edmonton
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] future_1.58.0         SeuratDisk_0.0.0.9021 scCustomize_3.0.1    
# [4] harmony_1.2.3         Rcpp_1.0.14           readr_2.1.5          
# [7] purrr_1.0.4           Seurat_5.3.0          SeuratObject_5.1.0   
# [10] sp_2.2-0              ggplot2_3.5.2         dplyr_1.1.4          
# [13] stringr_1.5.1         here_1.0.1           
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     rstudioapi_0.17.1      jsonlite_2.0.0        
# [4] shape_1.4.6.1          magrittr_2.0.3         spatstat.utils_3.1-4  
# [7] ggbeeswarm_0.7.2       farver_2.1.2           GlobalOptions_0.1.2   
# [10] vctrs_0.6.5            ROCR_1.0-11            spatstat.explore_3.4-3
# [13] paletteer_1.6.0        janitor_2.2.1          htmltools_0.5.8.1     
# [16] forcats_1.0.0          sctransform_0.4.2      parallelly_1.45.0     
# [19] KernSmooth_2.23-26     htmlwidgets_1.6.4      ica_1.0-3             
# [22] plyr_1.8.9             lubridate_1.9.4        plotly_4.10.4         
# [25] zoo_1.8-14             igraph_2.1.4           mime_0.13             
# [28] lifecycle_1.0.4        iterators_1.0.14       pkgconfig_2.0.3       
# [31] Matrix_1.7-3           R6_2.6.1               fastmap_1.2.0         
# [34] snakecase_0.11.1       fitdistrplus_1.2-2     shiny_1.10.0          
# [37] digest_0.6.37          colorspace_2.1-1       rematch2_2.1.2        
# [40] patchwork_1.3.0        rprojroot_2.0.4        tensor_1.5            
# [43] prismatic_1.1.2        RSpectra_0.16-2        irlba_2.3.5.1         
# [46] labeling_0.4.3         progressr_0.15.1       timechange_0.3.0      
# [49] spatstat.sparse_3.1-0  httr_1.4.7             polyclip_1.10-7       
# [52] abind_1.4-8            compiler_4.5.0         bit64_4.6.0-1         
# [55] withr_3.0.2            doParallel_1.0.17      fastDummies_1.7.5     
# [58] MASS_7.3-65            tools_4.5.0            vipor_0.4.7           
# [61] lmtest_0.9-40          beeswarm_0.4.0         httpuv_1.6.16         
# [64] future.apply_1.20.0    goftest_1.2-3          glue_1.8.0            
# [67] nlme_3.1-168           promises_1.3.3         grid_4.5.0            
# [70] Rtsne_0.17             cluster_2.1.8.1        reshape2_1.4.4        
# [73] generics_0.1.4         hdf5r_1.3.12           gtable_0.3.6          
# [76] spatstat.data_3.1-6    tzdb_0.5.0             tidyr_1.3.1           
# [79] data.table_1.17.4      hms_1.1.3              spatstat.geom_3.4-1   
# [82] RcppAnnoy_0.0.22       ggrepel_0.9.6          RANN_2.6.2            
# [85] foreach_1.5.2          pillar_1.10.2          vroom_1.6.5           
# [88] ggprism_1.0.6          spam_2.11-1            RcppHNSW_0.6.0        
# [91] later_1.4.2            circlize_0.4.16        splines_4.5.0         
# [94] lattice_0.22-7         bit_4.6.0              survival_3.8-3        
# [97] deldir_2.0-4           tidyselect_1.2.1       miniUI_0.1.2          
# [100] pbapply_1.7-2          gridExtra_2.3          scattermore_1.2       
# [103] RhpcBLASctl_0.23-42    spacexr_2.2.1          matrixStats_1.5.0     
# [106] stringi_1.8.7          lazyeval_0.2.2         codetools_0.2-20      
# [109] tibble_3.3.0           cli_3.6.5              uwot_0.2.3            
# [112] xtable_1.8-4           reticulate_1.42.0      dichromat_2.0-0.1     
# [115] globals_0.18.0         spatstat.random_3.4-1  png_0.1-8             
# [118] ggrastr_1.0.2          spatstat.univar_3.1-3  parallel_4.5.0        
# [121] dotCall64_1.2          listenv_0.9.1          viridisLite_0.4.2     
# [124] scales_1.4.0           ggridges_0.5.6         crayon_1.5.3          
# [127] rlang_1.1.6            cowplot_1.1.3  

