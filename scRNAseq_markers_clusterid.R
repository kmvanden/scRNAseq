# Find markers and cluster identification in scRNAseq data using Seurat
# https://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_CI_markers.R
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

# cluster identification: requires multiple iterations, modifying the cluster resolution and parameter thresholds in FindMarkers() function


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# set seed
set.seed(1234)

# load libraries
library(Seurat)
library(tidyverse)

# load the data (use ifb_harmony.rds file generated during scRNAseq_integration.R script
ifnb_harmony <- readRDS(file = "ifnb_harmony.rds")
str(ifnb_harmony)
head(ifnb_harmony@meta.data)

# visualize the data
clusters <- DimPlot(ifnb_harmony, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
condition <- DimPlot(ifnb_harmony, reduction = "umap", group.by = "stim")

condition|clusters


# # FindAllMarkers(): compare genes differentially expressed between different clusters -> when you have one condition
# # iteratively comapres each cluster to every other cluster
# DefaultAssay(ifnb_harmony) # needs to be "RNA" (will be changed by other integration methods, but not by Harmony)
# # DefaultAssay(ifnb_harmony) <- "RNA"
# diff_clusters <- FindAllMarkers(ifnb_harmony, logfc.threshold = 0.25, min.pct = 0.1,
#                                 only.pos = TRUE, test.use = "DESeq2", slot = "counts")
# # DESeq2 uses raw counts (default for FindAllMarkers is the data slot that stores the normalized data)


##########################################################
### IDENTIFY CELL TYPES IN EACH CLUSTER - MARKER GENES ###
##########################################################

# FindConservedMarkers(): identify the cell types present in each cluster (markers that are highly expressed in each cluster) --> when you have multiple conditions --> separates by condition
# example what is the identity of cells in cluster 3 (compare to all other clusters; can also use ident.2 = 11 to compare to cluster 11)
markers_cluster3 <- FindConservedMarkers(ifnb_harmony, ident.1 = 3, grouping.var = "stim")
head(markers_cluster3)
rownames(markers_cluster3)

# visualize the top features in markers_cluster3 (top genes the same in control an stimulated cells)
FeaturePlot(ifnb_harmony, features = c("FCGR3A"), min.cutoff = "q10")
FeaturePlot(ifnb_harmony, features = c("VMO1"), min.cutoff = "q10")
FeaturePlot(ifnb_harmony, features = c("MS4A7"), min.cutoff = "q10")
FeaturePlot(ifnb_harmony, features = c("CXCL16"), min.cutoff = "q10")
# cells that have <= tenth quantile expression of the genes will be in grey --> increases the contrast


# rename the Idents for cluster 3 to CD16 monocytes (FCGR3A is a marker for CD16 monocytes)
Idents(ifnb_harmony)
ifnb_harmony <- RenameIdents(ifnb_harmony, "3" = "CD16 Mono")
DimPlot(ifnb_harmony, reduction = "umap", label = TRUE) # cluster named CD16 Mono on graph

### to identify the other clusters, would need to repeat this process for each of the clusters


### cells already have annotations provided in the metadata --> seurat_annotations (data from published article)
head(ifnb_harmony@meta.data)

# set Idents to seurat_annotations
Idents(ifnb_harmony) <- ifnb_harmony@meta.data$seurat_annotations
Idents(ifnb_harmony)

DimPlot(ifnb_harmony, reduction = "umap", label = TRUE)


#######################################################################################
### GENES THAT ARE DIFFERNTIALLY EXPRESSED BETWEEN CONDITIONS FOR A GIVEN CELL TYPE ###
#######################################################################################

# create a column with both cell type and stimulation information
ifnb_harmony$celltype.cnd <- paste0(ifnb_harmony$seurat_annotations,'_', ifnb_harmony$stim)
head(ifnb_harmony@meta.data)
Idents(ifnb_harmony) <- ifnb_harmony$celltype.cnd

DimPlot(ifnb_harmony, reduction = 'umap', label = TRUE)

# look at the change in expression of genes in CD16 monocytes after IFNb stimulation (same cell type before and after stimulation)
b.interferon.response <- FindMarkers(ifnb_harmony, ident.1 = "CD16 Mono_STIM", ident.2 = "CD16 Mono_CTRL")
head(b.interferon.response)

# plot using the top two features from FindConservedMarkers() and the top two features from FindMarkers()

FeaturePlot(ifnb_harmony, features = c("FCGR3A", "VMO1", "IFIT1", "ISG15"), split.by = "stim", min.cutoff = "q10")
# FindConservedMarkers(): high in one cluster/cell type in both stim and control conditions
# FindMarkers(): high in CD16 monocytes that were stimulated (cell type that was chosen for FindMarkers() above)
    # possibly increased in multiple clusters/cell types

# can be used to find marker genes for one cluster compared to another cluster or all other clusters (if cells have been stimulated, the markers might be due to stimulation and not cell type markers)
    # comparison based on what is chosen for ident.1 and ident.2 in FindMarkers()


### can examine results using several different types of plots in addition to the feature plot
    # plot using the top two features from FindConservedMarkers() and the top two features from FindMarkers()

# FeaturePlot()
FeaturePlot(ifnb_harmony, features = c("FCGR3A", "VMO1", "IFIT1", "ISG15"), split.by = "stim", min.cutoff = "q10")
# VlnPlot()
vln <- VlnPlot(ifnb_harmony, features = c("FCGR3A", "VMO1", "IFIT1", "ISG15"), split.by = "stim")
(vln[[1]] | vln[[2]]) / (vln[[3]] | vln[[4]])
#DotPlot()
DotPlot(ifnb_harmony, features = c("FCGR3A", "VMO1", "IFIT1", "ISG15"), split.by = "stim")


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
#   [1] future_1.58.0             lubridate_1.9.4           forcats_1.0.0            
# [4] stringr_1.5.1             dplyr_1.1.4               purrr_1.0.4              
# [7] readr_2.1.5               tidyr_1.3.1               tibble_3.3.0             
# [10] ggplot2_3.5.2             tidyverse_2.0.0           stxBrain.SeuratData_0.1.2
# [13] ssHippo.SeuratData_3.1.4  pbmcsca.SeuratData_3.0.0  pbmcref.SeuratData_1.0.0 
# [16] pbmc3k.SeuratData_3.1.4   panc8.SeuratData_3.0.2    ifnb.SeuratData_3.1.0    
# [19] SeuratData_0.2.2.9002     Seurat_5.3.0              SeuratObject_5.1.0       
# [22] sp_2.2-0                  harmony_1.2.3             Rcpp_1.0.14              
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22            splines_4.5.0               later_1.4.2                
# [4] polyclip_1.10-7             fastDummies_1.7.5           lifecycle_1.0.4            
# [7] Rdpack_2.6.4                doParallel_1.0.17           globals_0.18.0             
# [10] lattice_0.22-7              MASS_7.3-65                 magrittr_2.0.3             
# [13] limma_3.64.1                plotly_4.10.4               plotrix_3.8-4              
# [16] qqconf_1.3.2                httpuv_1.6.16               sn_2.1.1                   
# [19] sctransform_0.4.2           spam_2.11-1                 spatstat.sparse_3.1-0      
# [22] reticulate_1.42.0           cowplot_1.1.3               pbapply_1.7-2              
# [25] RColorBrewer_1.1-3          multcomp_1.4-28             abind_1.4-8                
# [28] Rtsne_0.17                  GenomicRanges_1.60.0        presto_1.0.0               
# [31] BiocGenerics_0.54.0         TH.data_1.1-3               sandwich_3.1-1             
# [34] rappdirs_0.3.3              GenomeInfoDbData_1.2.14     IRanges_2.42.0             
# [37] S4Vectors_0.46.0            ggrepel_0.9.6               irlba_2.3.5.1              
# [40] listenv_0.9.1               spatstat.utils_3.1-4        TFisher_0.2.0              
# [43] goftest_1.2-3               RSpectra_0.16-2             spatstat.random_3.4-1      
# [46] fitdistrplus_1.2-2          parallelly_1.45.0           codetools_0.2-20           
# [49] DelayedArray_0.34.1         tidyselect_1.2.1            UCSC.utils_1.4.0           
# [52] farver_2.1.2                matrixStats_1.5.0           stats4_4.5.0               
# [55] spatstat.explore_3.4-3      mathjaxr_1.8-0              jsonlite_2.0.0             
# [58] multtest_2.64.0             progressr_0.15.1            ggridges_0.5.6             
# [61] survival_3.8-3              iterators_1.0.14            foreach_1.5.2              
# [64] tools_4.5.0                 ica_1.0-3                   glue_1.8.0                 
# [67] mnormt_2.1.1                gridExtra_2.3               SparseArray_1.8.0          
# [70] metap_1.12                  DESeq2_1.48.1               MatrixGenerics_1.20.0      
# [73] GenomeInfoDb_1.44.0         withr_3.0.2                 numDeriv_2016.8-1.1        
# [76] fastmap_1.2.0               digest_0.6.37               timechange_0.3.0           
# [79] R6_2.6.1                    mime_0.13                   colorspace_2.1-1           
# [82] scattermore_1.2             tensor_1.5                  dichromat_2.0-0.1          
# [85] spatstat.data_3.1-6         RhpcBLASctl_0.23-42         generics_0.1.4             
# [88] data.table_1.17.4           httr_1.4.7                  htmlwidgets_1.6.4          
# [91] S4Arrays_1.8.1              uwot_0.2.3                  pkgconfig_2.0.3            
# [94] gtable_0.3.6                lmtest_0.9-40               XVector_0.48.0             
# [97] htmltools_0.5.8.1           dotCall64_1.2               scales_1.4.0               
# [100] Biobase_2.68.0              png_0.1-8                   spatstat.univar_3.1-3      
# [103] rstudioapi_0.17.1           tzdb_0.5.0                  reshape2_1.4.4             
# [106] nlme_3.1-168                zoo_1.8-14                  spacexr_2.2.1              
# [109] KernSmooth_2.23-26          parallel_4.5.0              miniUI_0.1.2               
# [112] pillar_1.10.2               grid_4.5.0                  vctrs_0.6.5                
# [115] RANN_2.6.2                  promises_1.3.3              xtable_1.8-4               
# [118] cluster_2.1.8.1             mvtnorm_1.3-3               cli_3.6.5                  
# [121] locfit_1.5-9.12             compiler_4.5.0              rlang_1.1.6                
# [124] crayon_1.5.3                future.apply_1.20.0         mutoss_0.1-13              
# [127] labeling_0.4.3              plyr_1.8.9                  stringi_1.8.7              
# [130] viridisLite_0.4.2           deldir_2.0-4                BiocParallel_1.42.1        
# [133] lazyeval_0.2.2              spatstat.geom_3.4-1         Matrix_1.7-3               
# [136] RcppHNSW_0.6.0              hms_1.1.3                   patchwork_1.3.0            
# [139] statmod_1.5.0               shiny_1.10.0                SummarizedExperiment_1.38.1
# [142] rbibutils_2.3               ROCR_1.0-11                 igraph_2.1.4    

