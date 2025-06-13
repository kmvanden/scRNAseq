# Using SCTransform in Seurat
# https://satijalab.org/seurat/articles/sctransform_vignette

# biological heterogeneity in scRNAseq data is often confounded by technical factors (e.g., sequrncing depth)
# the number of molecules detected in each cell can vary significantly between cells, thus interpretation of scRNAseq data required effective pre-processing and normalization to remove this technical variability

# log normalization relies on the assumption that each cell originally contained the same number of molecules
# SCTransform doesn't make this assumption
# SCTransform omits the need for heuristic steps including pseudocount addition


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# load libraries
library(Seurat)
library(ggplot2)
library(sctransform)

# load data
pbmc_data <- Read10X("hg19/")
pbmc <- CreateSeuratObject(counts = pbmc_data)


### apply SCTransform normalization
# SCTransform replaces: NoramlizeData(), FindVariableFeatures() and ScaleData()
# transformed data will be available in the SCT assay (set as the default after running SCTransform)
# confounding sources of variation like high mitochondrial percentages, can be removed during normalization

# store mitochondiral percentage in the object's metadata
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

# run SCTransform
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)


### perform dimensionality reduction by PCA and UMAP embedding
pbmc <- pbmc %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE) %>%
  FindClusters(verbose = FALSE)
DimPlot(pbmc, label = TRUE)

# using more PCs when using SCTransform than you would with the standard Seurat workflow is often beneficial
# after standard log-normalization, variation in sequencing depth is still a confounding factor, which can subtly influence higher PCs
# this effect is substantially mitigated with SCTransform, thus higher PCs are more likely to represent subtle, but biologically relevant, sources of heterogeneity
# thus including higher PCs with SCTransform may improve downstream analysis

# SCTrasnform returns 3000 variable features, unlike the 2000 in the standard Seurat workflow
# additional variable features are less likely to be driven by technical differences across cells, and instead may represent more subtle biological findings

# following code replicates the full end-to-end workflow in a single command
# pbmc <- CreateSeuratObject(pbmc_data) %>%
#   PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
#   SCTransform(vars.to.regress = "percent.mt") %>%
#   RunPCA() %>%
#   FindNeighbors(dims = 1:30) %>%
#   RunUMAP(dims = 1:30) %>%
#   FindClusters()


pbmc[["SCT"]]$scale.data # normalized values
# only values for variable genes are stored
# return.only.var.genes = TRUE is the defaut in SCTransform

# Pearson residual are converted back to corrected UMI counts (the UMI counts you would expect the observe if all cells had been sequenced to the same depth)
pbmc[["SCT"]]$counts # corrected UMI counts
pbmc[["SCT"]]$data # log normalized versions of the corrected counts


### visualize canonical marker genes as violin plots
VlnPlot(pbmc, features = c("CD8A", "GZMK", "CCL5", "S100A4", 
                           "ANXA1", "CCR7", "ISG15", "CD3D"),
        pt.size = 0.2, ncol = 4)

# visualize canonical marker genes on the sctransform embedding.
FeaturePlot(pbmc, features = c("CD8A", "GZMK", "CCL5", "S100A4", 
                               "ANXA1", "CCR7"), 
            pt.size = 0.2, ncol = 3)
FeaturePlot(pbmc, features = c("CD3D", "ISG15", "TCL1A", 
                               "FCER2", "XCL1", "FCGR3A"), 
            pt.size = 0.2, ncol = 3)


### use of SCTransform with Harmonry integration
# https://github.com/satijalab/seurat/issues/4896
# you can do SCTransform(object, return.only.var.genes = FALSE) to return residuals for all genes.

data("pbmc")

pbmc.list <- SplitObject(pbmc, split.by="Method")
pbmc.list <- lapply(X = pbmc.list, 
                       FUN = SCTransform, 
                       method = "glmGamPoi", 
                       return.only.var.genes = FALSE)
var.features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 3000)

pbmc.sct <- merge(x = pbmc.list[[1]], y = pbmc.list[2:length(pbmc.list)], merge.data=TRUE)
VariableFeatures(pbmc.sct) <- var.features
pbmc.sct <- RunPCA(pbmc.sct, verbose = FALSE)
pbmc.sct <- RunHarmony(pbmc.sct, assay.use="SCT", group.by.vars = "Method")
pbmc.sct <- RunUMAP(pbmc.sct, reduction = "harmony", dims = 1:30)
pbmc.sct <- FindNeighbors(pbmc.sct, reduction = "harmony", dims = 1:30) %>% FindClusters()
DimPlot(pbmc.sct, group.by = c("Method", "ident", "CellType"), ncol = 3)


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
#   [1] sctransform_0.4.2  future_1.58.0      Seurat_5.3.0       SeuratObject_5.1.0
# [5] sp_2.2-0           lubridate_1.9.4    forcats_1.0.0      stringr_1.5.1     
# [9] dplyr_1.1.4        purrr_1.0.4        readr_2.1.5        tidyr_1.3.1       
# [13] tibble_3.3.0       ggplot2_3.5.2      tidyverse_2.0.0   
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3          rstudioapi_0.17.1           jsonlite_2.0.0             
# [4] magrittr_2.0.3              spatstat.utils_3.1-4        ggbeeswarm_0.7.2           
# [7] farver_2.1.2                vctrs_0.6.5                 ROCR_1.0-11                
# [10] DelayedMatrixStats_1.30.0   spatstat.explore_3.4-3      S4Arrays_1.8.1             
# [13] htmltools_0.5.8.1           SparseArray_1.8.0           parallelly_1.45.0          
# [16] KernSmooth_2.23-26          htmlwidgets_1.6.4           ica_1.0-3                  
# [19] plyr_1.8.9                  plotly_4.10.4               zoo_1.8-14                 
# [22] igraph_2.1.4                mime_0.13                   lifecycle_1.0.4            
# [25] iterators_1.0.14            pkgconfig_2.0.3             Matrix_1.7-3               
# [28] R6_2.6.1                    fastmap_1.2.0               GenomeInfoDbData_1.2.14    
# [31] MatrixGenerics_1.20.0       fitdistrplus_1.2-2          shiny_1.10.0               
# [34] digest_0.6.37               colorspace_2.1-1            S4Vectors_0.46.0           
# [37] patchwork_1.3.0             tensor_1.5                  RSpectra_0.16-2            
# [40] irlba_2.3.5.1               GenomicRanges_1.60.0        beachmat_2.24.0            
# [43] labeling_0.4.3              progressr_0.15.1            spatstat.sparse_3.1-0      
# [46] timechange_0.3.0            mgcv_1.9-3                  httr_1.4.7                 
# [49] polyclip_1.10-7             abind_1.4-8                 compiler_4.5.0             
# [52] bit64_4.6.0-1               withr_3.0.2                 doParallel_1.0.17          
# [55] fastDummies_1.7.5           R.utils_2.13.0              MASS_7.3-65                
# [58] DelayedArray_0.34.1         tools_4.5.0                 vipor_0.4.7                
# [61] lmtest_0.9-40               beeswarm_0.4.0              httpuv_1.6.16              
# [64] future.apply_1.20.0         goftest_1.2-3               glmGamPoi_1.20.0           
# [67] R.oo_1.27.1                 glue_1.8.0                  nlme_3.1-168               
# [70] promises_1.3.3              grid_4.5.0                  Rtsne_0.17                 
# [73] cluster_2.1.8.1             reshape2_1.4.4              generics_0.1.4             
# [76] hdf5r_1.3.12                gtable_0.3.6                spatstat.data_3.1-6        
# [79] tzdb_0.5.0                  R.methodsS3_1.8.2           data.table_1.17.4          
# [82] hms_1.1.3                   XVector_0.48.0              BiocGenerics_0.54.0        
# [85] spatstat.geom_3.4-1         RcppAnnoy_0.0.22            ggrepel_0.9.6              
# [88] RANN_2.6.2                  foreach_1.5.2               pillar_1.10.2              
# [91] spam_2.11-1                 RcppHNSW_0.6.0              later_1.4.2                
# [94] splines_4.5.0               lattice_0.22-7              survival_3.8-3             
# [97] bit_4.6.0                   deldir_2.0-4                tidyselect_1.2.1           
# [100] miniUI_0.1.2                pbapply_1.7-2               gridExtra_2.3              
# [103] IRanges_2.42.0              SummarizedExperiment_1.38.1 scattermore_1.2            
# [106] stats4_4.5.0                Biobase_2.68.0              spacexr_2.2.1              
# [109] matrixStats_1.5.0           UCSC.utils_1.4.0            stringi_1.8.7              
# [112] lazyeval_0.2.2              codetools_0.2-20            cli_3.6.5                  
# [115] uwot_0.2.3                  xtable_1.8-4                reticulate_1.42.0          
# [118] GenomeInfoDb_1.44.0         dichromat_2.0-0.1           Rcpp_1.0.14                
# [121] globals_0.18.0              spatstat.random_3.4-1       png_0.1-8                  
# [124] ggrastr_1.0.2               spatstat.univar_3.1-3       parallel_4.5.0             
# [127] dotCall64_1.2               sparseMatrixStats_1.20.0    listenv_0.9.1              
# [130] viridisLite_0.4.2           scales_1.4.0                ggridges_0.5.6             
# [133] crayon_1.5.3                rlang_1.1.6                 cowplot_1.1.3 

