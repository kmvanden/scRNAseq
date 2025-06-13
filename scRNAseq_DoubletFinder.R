# DoubletFinder: detection of doublets in scRNAseq data to filter them out
# https://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_doublets.R

# homotypic doublets: doublets derived from transcriptionally similar cells
# heterotypic doublets: doublets derived from transcriptionally distinct cells

# DoubletFinder needs 3 parameters
   # pN = the number of artificial doublets   
   # pK = the neighborhood size (pK) used to compute the number of artificial nearest neighbors
   # Exp = the number of expected real doublets (can be determined based on the number of cells loaded and the number of cells recovered)

# 1. stimulates artificial doublets from existing scRNAseq data by averaging the gene expression profile from random pairs of cells
# 2. merges artificial and real data and performs Seurat preprocessing steps
# 3. performs dimensionality reduction (PCA) to determine the similarity between real and artificial cells
# 4. detects the k nearest neighbors for every real cell in the PC space and computes each cell's proportion of artificial nearest neighbors (pANN)
   # highly dependent on the pK parameter (different methods to select based on whether there is a ground truth)
# 5. pANN that is identified for each cell is thresholded using the total number of expected doublets and is used to predict the real doublets

# DoubletFinder is a QC step and should not be applied to aggregated or integrated data
# should be applied on each sample individually (not merged data)
# should be applied on quality controlled cells


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)


# data (feature/cell matrix (raw)): https://www.10xgenomics.com/datasets/10k-human-pbmcs-3-v3-1-chromium-controller-3-1-high 
# create counts matrix
cts <- ReadMtx(mtx = "raw_feature_bc_matrix/matrix.mtx.gz",
               features = "raw_feature_bc_matrix/features.tsv.gz",
               cells = "raw_feature_bc_matrix/barcodes.tsv.gz")
cts[1:10,1:10] # sparse Matrix of class "dgCMatrix"


# create a Seurat object
pbmc.seurat <- CreateSeuratObject(counts = cts)
str(pbmc.seurat) # Formal class 'Seurat' [package "SeuratObject"] with 13 slots


### quality control and filtering
pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')
pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
                                 nFeature_RNA > 500 &
                                 mitoPercent < 10)

pbmc.seurat # 36601 features across 2099284 samples within 1 assay 
pbmc.seurat.filtered # 36601 features across 10017 samples within 1 assay 


### standard Seurat pre-processing workflow
pbmc.seurat.filtered <- pbmc.seurat.filtered %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(pbmc.seurat.filtered, ndims = 50)

pbmc.seurat.filtered <- pbmc.seurat.filtered %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:20)


### pK identification (no ground-truth)
sweep.res.list_pbmc <- paramSweep(pbmc.seurat.filtered, PCs = 1:20, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
head(sweep.stats_pbmc)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
head(bcmvn_pbmc)

ggplot(bcmvn_pbmc, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line() 

pK <- bcmvn_pbmc %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]])) # stores the optimal pK value in a variable | 0.26


### homotypic doublet proportion estimate
annotations <- pbmc.seurat.filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) # proportion of homotypic doublets
nExp_poi <- round(0.076*nrow(pbmc.seurat.filtered@meta.data)) # expected number of doublets: based on number of cells loaded and number of cells recovered
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # 691


# run doubletFinder 
# packageVersion("Seurat") # 5.3.0
# packageVersion("DoubletFinder") # 2.0.6
# DoubletFinder expects the v4 structure of Seurat

# pbmc.seurat.filtered <- doubletFinder(pbmc.seurat.filtered, PCs = 1:20, 
#                                       pN = 0.25, pK = pK, nExp = nExp_poi.adj, 
#                                       reuse.pANN = FALSE, sct = FALSE)
# 
# # visualize doublets
# DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = "DF.classifications_0.25_0.21_691")
# # number of singlets and doublets
# table(pbmc.seurat.filtered@meta.data$DF.classifications_0.25_0.21_691)


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
#   [1] DoubletFinder_2.0.6 future_1.58.0       lubridate_1.9.4     forcats_1.0.0       stringr_1.5.1      
# [6] dplyr_1.1.4         purrr_1.0.4         readr_2.1.5         tidyr_1.3.1         tibble_3.3.0       
# [11] ggplot2_3.5.2       tidyverse_2.0.0     Seurat_5.3.0        SeuratObject_5.1.0  sp_2.2-0           
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     rstudioapi_0.17.1      jsonlite_2.0.0         magrittr_2.0.3        
# [5] spatstat.utils_3.1-4   ggbeeswarm_0.7.2       farver_2.1.2           fields_16.3.1         
# [9] vctrs_0.6.5            ROCR_1.0-11            spatstat.explore_3.4-3 htmltools_0.5.8.1     
# [13] curl_6.3.0             sctransform_0.4.2      parallelly_1.45.0      KernSmooth_2.23-26    
# [17] desc_1.4.3             htmlwidgets_1.6.4      ica_1.0-3              plyr_1.8.9            
# [21] plotly_4.10.4          zoo_1.8-14             igraph_2.1.4           mime_0.13             
# [25] lifecycle_1.0.4        iterators_1.0.14       pkgconfig_2.0.3        Matrix_1.7-3          
# [29] R6_2.6.1               fastmap_1.2.0          fitdistrplus_1.2-2     shiny_1.10.0          
# [33] digest_0.6.37          colorspace_2.1-1       ps_1.9.1               patchwork_1.3.0       
# [37] tensor_1.5             RSpectra_0.16-2        irlba_2.3.5.1          labeling_0.4.3        
# [41] progressr_0.15.1       spatstat.sparse_3.1-0  timechange_0.3.0       mgcv_1.9-3            
# [45] httr_1.4.7             polyclip_1.10-7        abind_1.4-8            compiler_4.5.0        
# [49] remotes_2.5.0          bit64_4.6.0-1          withr_3.0.2            doParallel_1.0.17     
# [53] fastDummies_1.7.5      pkgbuild_1.4.8         maps_3.4.3             MASS_7.3-65           
# [57] tools_4.5.0            vipor_0.4.7            lmtest_0.9-40          beeswarm_0.4.0        
# [61] httpuv_1.6.16          future.apply_1.20.0    goftest_1.2-3          glue_1.8.0            
# [65] callr_3.7.6            nlme_3.1-168           promises_1.3.3         grid_4.5.0            
# [69] Rtsne_0.17             cluster_2.1.8.1        reshape2_1.4.4         generics_0.1.4        
# [73] hdf5r_1.3.12           gtable_0.3.6           spatstat.data_3.1-6    tzdb_0.5.0            
# [77] data.table_1.17.4      hms_1.1.3              spatstat.geom_3.4-1    RcppAnnoy_0.0.22      
# [81] ggrepel_0.9.6          RANN_2.6.2             foreach_1.5.2          pillar_1.10.2         
# [85] spam_2.11-1            RcppHNSW_0.6.0         later_1.4.2            splines_4.5.0         
# [89] lattice_0.22-7         survival_3.8-3         bit_4.6.0              deldir_2.0-4          
# [93] tidyselect_1.2.1       miniUI_0.1.2           pbapply_1.7-2          gridExtra_2.3         
# [97] scattermore_1.2        spacexr_2.2.1          matrixStats_1.5.0      stringi_1.8.7         
# [101] lazyeval_0.2.2         codetools_0.2-20       cli_3.6.5              uwot_0.2.3            
# [105] xtable_1.8-4           reticulate_1.42.0      processx_3.8.6         dichromat_2.0-0.1     
# [109] Rcpp_1.0.14            globals_0.18.0         spatstat.random_3.4-1  png_0.1-8             
# [113] ggrastr_1.0.2          spatstat.univar_3.1-3  parallel_4.5.0         dotCall64_1.2         
# [117] listenv_0.9.1          viridisLite_0.4.2      scales_1.4.0           ggridges_0.5.6        
# [121] crayon_1.5.3           rlang_1.1.6            cowplot_1.1.3    

