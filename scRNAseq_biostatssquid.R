# Standard scRNAseq preprocessing workflow with Seurat
# https://biostatsquid.com/scrnaseq-preprocessing-workflow-seurat/


### clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# load libraries 
library(tidyverse)
library(Seurat)

### setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

### set seed
set.seed(42)

### load data
# dataset from 10x genomics: https://www.10xgenomics.com/datasets/40-k-mixture-of-nsclc-dt-cs-from-7-donors-3-ht-v-3-1-3-1-high-6-1-0
nsclc_sm <- Read10X_h5("40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")
# the object is a list of matrices containing different types of data:  gene expression matrix, antibody capture, and multiplexing capture
str(nsclc_sm)
cts <- nsclc_sm$`Gene Expression` # sparse matrix 
class(cts) # dgCMatrix


### create Seurat object (raw counts)
# features that are in at least 3 cells and cells that have a minimum of 200 features
nsclc_seu <- CreateSeuratObject(counts = cts, project = 'NSCLC', min.cells = 3, min.features = 200)
str(nsclc_seu)
class(nsclc_seu) # SeuratObject


### quality control
## percent mitochondrial features
nsclc_seu[['percent_mt']] <- PercentageFeatureSet(nsclc_seu, pattern = '^MT-')
head(nsclc_seu@meta.data)

# visualize quality metrics to determine thresholds
VlnPlot(nsclc_seu, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"), ncol = 3)
FeatureScatter(nsclc_seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')


### filtering based on quality control thresholds
nsclc_seu <- subset(nsclc_seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mt < 5)


### normalization
# default: global-scaling normalization method (LogNormalize)
# normalizes the feature expression measurement for each cell by the total expression
# then multiples this by a scale factor (10,000) and log-transforms the result
nsclc_seu <- NormalizeData(nsclc_seu)
str(nsclc_seu)


### identify highly-variable features (features that have a high cell-to-cell variation)
# cells with low cell-to-cell variation (i.e., housekeeping genes) are not very informative toward identifying cell subsets
nsclc_seu <- FindVariableFeatures(nsclc_seu, selection.method =  'vst', nfeatures = 2000)
# Identify the top 10 HVGs
top10 <- head(VariableFeatures(nsclc_seu), 10)
top10_plot <- VariableFeaturePlot(nsclc_seu)
LabelPoints(plot = top10_plot, points = top10, repel = TRUE)


# scaling
all_genes <- rownames(nsclc_seu)
nsclc_seu <- ScaleData(nsclc_seu, features = all_genes)
head(nsclc_seu@assays$RNA)


# dimensionality reduction (PCA)
# summarizes multidimensionality data into principle components (PCs)
# cells with similar expression profiles will cluster together
nsclc_seu <- RunPCA(nsclc_seu, features = VariableFeatures(nsclc_seu))
print(nsclc_seu[['pca']], dims = 1:5, nfeatures = 5)

# visualize the PCs
DimHeatmap(nsclc_seu, dims = 1, cells = 500, balanced = TRUE)
DimPlot(nsclc_seu, reduction = "pca") + NoLegend()

# determine number of PCs to retain to explain the variation in the dataset
ElbowPlot(nsclc_seu, ndims = 50)


### clustering
# identify subgroups in the data so that the cells belonging to the same subgroup (cluster) are similar and cells in different clusters are very different
nsclc_seu <- FindNeighbors(nsclc_seu, dims = 1:15) # compute the k-nearest neighbors of each cell using the number of PCs chosen
nsclc_seu <- FindClusters(nsclc_seu, resolution = c(0.1, 0.3, 0.5, 0.7, 1)) # graph clustering assigns each cell a number (the cluster they belong to)
# number of clusters determined by the resolution parameter (higher the resolution = more clusters)
head(nsclc_seu@meta.data)

DimPlot(nsclc_seu, group.by = 'RNA_snn_res.0.1', label = TRUE)
Idents(nsclc_seu) <- 'RNA_snn_res.0.1' # set identity of clusters


# dimensionality reduction (UMAP) to visualize data
nsclc_seu <- RunUMAP(nsclc_seu, dims = 1:15)
DimPlot(nsclc_seu, reduction = 'umap')


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
#   [1] future_1.58.0      Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0          
# [5] lubridate_1.9.4    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4       
# [9] purrr_1.0.4        readr_2.1.5        tidyr_1.3.1        tibble_3.3.0      
# [13] ggplot2_3.5.2      tidyverse_2.0.0   
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     rstudioapi_0.17.1      jsonlite_2.0.0        
# [4] magrittr_2.0.3         spatstat.utils_3.1-4   ggbeeswarm_0.7.2      
# [7] farver_2.1.2           vctrs_0.6.5            ROCR_1.0-11           
# [10] spatstat.explore_3.4-3 htmltools_0.5.8.1      sctransform_0.4.2     
# [13] parallelly_1.45.0      KernSmooth_2.23-26     htmlwidgets_1.6.4     
# [16] ica_1.0-3              plyr_1.8.9             plotly_4.10.4         
# [19] zoo_1.8-14             igraph_2.1.4           mime_0.13             
# [22] lifecycle_1.0.4        iterators_1.0.14       pkgconfig_2.0.3       
# [25] Matrix_1.7-3           R6_2.6.1               fastmap_1.2.0         
# [28] fitdistrplus_1.2-2     shiny_1.10.0           digest_0.6.37         
# [31] colorspace_2.1-1       patchwork_1.3.0        tensor_1.5            
# [34] RSpectra_0.16-2        irlba_2.3.5.1          labeling_0.4.3        
# [37] progressr_0.15.1       spatstat.sparse_3.1-0  timechange_0.3.0      
# [40] mgcv_1.9-3             httr_1.4.7             polyclip_1.10-7       
# [43] abind_1.4-8            compiler_4.5.0         bit64_4.6.0-1         
# [46] withr_3.0.2            doParallel_1.0.17      fastDummies_1.7.5     
# [49] MASS_7.3-65            tools_4.5.0            vipor_0.4.7           
# [52] lmtest_0.9-40          beeswarm_0.4.0         httpuv_1.6.16         
# [55] future.apply_1.20.0    goftest_1.2-3          glue_1.8.0            
# [58] nlme_3.1-168           promises_1.3.3         grid_4.5.0            
# [61] Rtsne_0.17             cluster_2.1.8.1        reshape2_1.4.4        
# [64] generics_0.1.4         hdf5r_1.3.12           gtable_0.3.6          
# [67] spatstat.data_3.1-6    tzdb_0.5.0             data.table_1.17.4     
# [70] hms_1.1.3              spatstat.geom_3.4-1    RcppAnnoy_0.0.22      
# [73] ggrepel_0.9.6          RANN_2.6.2             foreach_1.5.2         
# [76] pillar_1.10.2          spam_2.11-1            RcppHNSW_0.6.0        
# [79] later_1.4.2            splines_4.5.0          lattice_0.22-7        
# [82] survival_3.8-3         bit_4.6.0              deldir_2.0-4          
# [85] tidyselect_1.2.1       miniUI_0.1.2           pbapply_1.7-2         
# [88] gridExtra_2.3          scattermore_1.2        spacexr_2.2.1         
# [91] matrixStats_1.5.0      stringi_1.8.7          lazyeval_0.2.2        
# [94] codetools_0.2-20       cli_3.6.5              uwot_0.2.3            
# [97] xtable_1.8-4           reticulate_1.42.0      dichromat_2.0-0.1     
# [100] Rcpp_1.0.14            globals_0.18.0         spatstat.random_3.4-1 
# [103] png_0.1-8              ggrastr_1.0.2          spatstat.univar_3.1-3 
# [106] parallel_4.5.0         dotCall64_1.2          listenv_0.9.1         
# [109] viridisLite_0.4.2      scales_1.4.0           ggridges_0.5.6        
# [112] crayon_1.5.3           rlang_1.1.6            cowplot_1.1.3 

