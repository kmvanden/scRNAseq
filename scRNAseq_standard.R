# Walk-through of standard Seurat workflow to analyze single-cell RNA sequencing
# https://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_standard_workflow.R


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# load libraries
library(Seurat)
library(tidyverse)

# data: 20k Mixture of NSCLC DTCs from 7 donors, 3' v3.1 (raw output)
# data source: https://www.10xgenomics.com/datasets/20k-mixture-of-nsclc-dtcs-from-7-donors-3-v3-1-with-intronic-reads-3-1-standard

# load the NSCLC dataset
nsclc.sparse.m <- Read10X_h5(filename = "20k_NSCLC_DTC_3p_nextgem_intron_Multiplex_count_raw_feature_bc_matrix.h5")
# Genome matrix has multiple modalities, returning a list of matrices for this genome
str(nsclc.sparse.m)
cts <- nsclc.sparse.m$`Gene Expression` # use counts of gene expression for the tutorial
cts[1:10, 1:10]


# read the counts (raw, non-normalized data) into a Seurat object 
# keep all features that are expressed in at least 3 cells and keep all the cells that have at least 200 features
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)
str(nsclc.seurat.obj)
nsclc.seurat.obj
# An object of class Seurat 
# 32978 features across 71880 samples within 1 assay 
# Active assay: RNA (32978 features, 0 variable features)
# 1 layer present: counts


### quality control
head(nsclc.seurat.obj@meta.data)

# create column of percent MT reads
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
head(nsclc.seurat.obj@meta.data)

# look at distribution of: number of features/genes per cell, number of counts per cell, and percentage of mitochondrial reads 
# poor quality cells have low numbers of features and low counts and a high percentage of MT reads
# high numbers of features and/counts could indicate doublets/multiple cells
VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')


### filter based on quality control
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)


### normalize data
# normalization.method = "LogNormalize" | scale.factor = 10000
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
nsclc.seurat.obj@commands


### identify highly variable features
# select a subset of features that exhibit high cell to cell variation
# highlights the biological signal in the dataset
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj)

top10 <- head(VariableFeatures(nsclc.seurat.obj), 10) # top 10 most highly variable genes

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


### scaling 
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj)
str(nsclc.seurat.obj)


### linear dimensionality reduction (PCA)
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj)

print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5) # top 5 features (positive and negaive for the first 5 PCs)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE) # visualize PCA results

ElbowPlot(nsclc.seurat.obj, ndims = 50) # determine dimensionality of the data


### find neighbors and find clusters clustering
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:30)
# resolution defines the granularity of the clusters
# the lower the number, the fewer the clusters
# resolution = 0.8 is the default
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
head(nsclc.seurat.obj@meta.data)


table(nsclc.seurat.obj@meta.data$RNA_snn_res.0.1) # 8 clusters
table(nsclc.seurat.obj@meta.data$RNA_snn_res.1) # 28 clusters

DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.1", label = TRUE)
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.1", label = TRUE)

# setting identity of the cells to the clusters
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.7"
Idents(nsclc.seurat.obj)


## non-linear dimensionality reduction (UMAP)
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:30)
DimPlot(nsclc.seurat.obj, reduction = "umap")


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
#   [1] future_1.58.0      lubridate_1.9.4    forcats_1.0.0      stringr_1.5.1     
# [5] dplyr_1.1.4        purrr_1.0.4        readr_2.1.5        tidyr_1.3.1       
# [9] tibble_3.3.0       ggplot2_3.5.2      tidyverse_2.0.0    Seurat_5.3.0      
# [13] SeuratObject_5.1.0 sp_2.2-0          
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     rstudioapi_0.17.1      jsonlite_2.0.0         magrittr_2.0.3        
# [5] spatstat.utils_3.1-4   ggbeeswarm_0.7.2       farver_2.1.2           vctrs_0.6.5           
# [9] ROCR_1.0-11            spatstat.explore_3.4-3 htmltools_0.5.8.1      sctransform_0.4.2     
# [13] parallelly_1.45.0      KernSmooth_2.23-26     htmlwidgets_1.6.4      ica_1.0-3             
# [17] plyr_1.8.9             plotly_4.10.4          zoo_1.8-14             igraph_2.1.4          
# [21] mime_0.13              lifecycle_1.0.4        iterators_1.0.14       pkgconfig_2.0.3       
# [25] Matrix_1.7-3           R6_2.6.1               fastmap_1.2.0          fitdistrplus_1.2-2    
# [29] shiny_1.10.0           digest_0.6.37          colorspace_2.1-1       patchwork_1.3.0       
# [33] tensor_1.5             RSpectra_0.16-2        irlba_2.3.5.1          labeling_0.4.3        
# [37] progressr_0.15.1       spatstat.sparse_3.1-0  timechange_0.3.0       mgcv_1.9-3            
# [41] httr_1.4.7             polyclip_1.10-7        abind_1.4-8            compiler_4.5.0        
# [45] bit64_4.6.0-1          withr_3.0.2            doParallel_1.0.17      fastDummies_1.7.5     
# [49] MASS_7.3-65            tools_4.5.0            vipor_0.4.7            lmtest_0.9-40         
# [53] beeswarm_0.4.0         httpuv_1.6.16          future.apply_1.20.0    goftest_1.2-3         
# [57] glue_1.8.0             nlme_3.1-168           promises_1.3.3         grid_4.5.0            
# [61] Rtsne_0.17             cluster_2.1.8.1        reshape2_1.4.4         generics_0.1.4        
# [65] hdf5r_1.3.12           gtable_0.3.6           spatstat.data_3.1-6    tzdb_0.5.0            
# [69] data.table_1.17.4      hms_1.1.3              spatstat.geom_3.4-1    RcppAnnoy_0.0.22      
# [73] ggrepel_0.9.6          RANN_2.6.2             foreach_1.5.2          pillar_1.10.2         
# [77] spam_2.11-1            RcppHNSW_0.6.0         later_1.4.2            splines_4.5.0         
# [81] lattice_0.22-7         survival_3.8-3         bit_4.6.0              deldir_2.0-4          
# [85] tidyselect_1.2.1       miniUI_0.1.2           pbapply_1.7-2          gridExtra_2.3         
# [89] scattermore_1.2        spacexr_2.2.1          matrixStats_1.5.0      stringi_1.8.7         
# [93] lazyeval_0.2.2         codetools_0.2-20       cli_3.6.5              uwot_0.2.3            
# [97] xtable_1.8-4           reticulate_1.42.0      dichromat_2.0-0.1      Rcpp_1.0.14           
# [101] globals_0.18.0         spatstat.random_3.4-1  png_0.1-8              ggrastr_1.0.2         
# [105] spatstat.univar_3.1-3  parallel_4.5.0         dotCall64_1.2          listenv_0.9.1         
# [109] viridisLite_0.4.2      scales_1.4.0           ggridges_0.5.6         crayon_1.5.3          
# [113] rlang_1.1.6            cowplot_1.1.3   

