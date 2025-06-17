# Mapping and annotating query datasets
# https://satijalab.org/seurat/articles/integration_mapping


# build an integrated reference to annotate query datasets
# generating an integrated reference follows the same workflow used to integrate data to deal with batch effects
# once generated, the integrated reference can be used to analyze additional query datasets
    # using cell type label transfer and projecting query cells onto reference UMAPs


# load libraries
library(Seurat)
library(SeuratData)
library(ggplot2)
library(tidyverse)

# load data
# InstallData("panc8") # human pancreatic islet cell datasets (four technologies)
panc8 <- LoadData("panc8")

# select data from celseq2 and smartseq2 technologies
table(panc8$tech)
pancreas.ref <- subset(panc8, tech %in% c("celseq2", "smartseq2"))
pancreas.ref[["RNA"]] <- split(pancreas.ref[["RNA"]], f = pancreas.ref$tech)


### pre-process the data without integration
pancreas.ref <- pancreas.ref %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30)

# visualize data (batch effects)
DimPlot(pancreas.ref, group.by = c("celltype", "tech"))


### integrate the datasets into a shared reference 
pancreas.ref <- IntegrateLayers(object = pancreas.ref, method = CCAIntegration, 
                                orig.reduction = "pca", new.reduction = "integrated.cca", 
                                verbose = FALSE)

# re-run clustering and run UMAP
pancreas.ref <- pancreas.ref %>%
  FindNeighbors(reduction = "integrated.cca", dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(reduction = "integrated.cca", dims = 1:30)

# visualize data
DimPlot(pancreas.ref, group.by = c("tech", "celltype"))


### cell type classification using an integrated reference 
# projection of reference data onto a query object
# Seurat doesn't correct or modify the query expression data
# Seurat projects the PCA structure of a reference onto the query (default), instead of learning a joint structure with CCA

# select two technologies for the query datasets
pancreas.query <- subset(panc8, tech %in% c("fluidigmc1", "celseq"))
pancreas.query <- NormalizeData(pancreas.query)


# find anchors
pancreas.anchors <- FindTransferAnchors(reference = pancreas.ref, query = pancreas.query, 
                                        dims = 1:30, reference.reduction = "pca")
# classify the query cells based on a reference dataset
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.ref$celltype, 
                            dims = 1:30)

# add the predictions to the metadata
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)


# compare the original label annotations (celltype) to the predicted celltype annotations (predicted.id)
pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)
# FALSE  TRUE 
#    63  1579  # most cells were labelled correctly

# verify predictions further by examining canonical cell type markers
table(pancreas.query$predicted.id)
# acinar activated_stellate              alpha               beta              delta             ductal 
#    262                 39                436                419                 73                330 
# endothelial              gamma         macrophage               mast            schwann 
#          19                 41                 15                  2                  6 

VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")


#### Unimodal UMAP Projection
# project the query onto the reference UMAP structure
# compute the reference UMAP model then call MapQuery() rather than TransferData()
pancreas.ref <- RunUMAP(pancreas.ref, dims = 1:30, reduction = "integrated.cca", 
                        return.model = TRUE)
pancreas.query <- MapQuery(anchorset = pancreas.anchors, reference = pancreas.ref, 
                           query = pancreas.query, refdata = list(celltype = "celltype"), 
                           reference.reduction = "pca", reduction.model = "umap")

### MapQuery() is a wrapper around:
# # TransferData(): transfers cell type label and imputes ADT values
# pancreas.query <- TransferData(anchorset = pancreas.anchors, reference = pancreas.ref, 
#                                query = pancreas.query, refdata = list(celltype = "celltype"))
# # IntegrateEmbeddings(): integrates the reference with the query by correcting hte query's projected low dimensional embeddings
# pancreas.query <- IntegrateEmbeddings(anchorset = pancreas.anchors, reference = pancreas.ref, 
#                                       query = pancreas.query, new.reduction.name = "ref.pca")
# # ProjectUMAP(): projects the query data onto the UMAP structure of the reference
# pancreas.query <- ProjectUMAP(query = pancreas.query, query.reduction = "ref.pca", 
#                               reference = pancreas.ref, reference.reduction = "pca", 
#                               reduction.model = "umap")


### visualize the query cells alongside the reference cells
p1 <- DimPlot(pancreas.ref, reduction = "umap", group.by = "celltype", label = TRUE, 
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", 
              label = TRUE, label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2


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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] panc8.SeuratData_3.0.2      stxBrain.SeuratData_0.1.2   ssHippo.SeuratData_3.1.4   
# [4] pbmcsca.SeuratData_3.0.0    pbmcref.SeuratData_1.0.0    pbmc3k.SeuratData_3.1.4    
# [7] ifnb.SeuratData_3.1.0       SeuratData_0.2.2.9002       future_1.58.0              
# [10] pheatmap_1.0.13             lubridate_1.9.4             forcats_1.0.0              
# [13] stringr_1.5.1               dplyr_1.1.4                 purrr_1.0.4                
# [16] readr_2.1.5                 tidyr_1.3.1                 tibble_3.3.0               
# [19] ggplot2_3.5.2               tidyverse_2.0.0             Seurat_5.3.0               
# [22] SeuratObject_5.1.0          sp_2.2-0                    celldex_1.18.0             
# [25] SingleR_2.10.0              SummarizedExperiment_1.38.1 Biobase_2.68.0             
# [28] GenomicRanges_1.60.0        GenomeInfoDb_1.44.0         IRanges_2.42.0             
# [31] S4Vectors_0.46.0            BiocGenerics_0.54.0         generics_0.1.4             
# [34] MatrixGenerics_1.20.0       matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] spatstat.sparse_3.1-0     httr_1.4.7                RColorBrewer_1.1-3       
# [4] doParallel_1.0.17         tools_4.5.0               sctransform_0.4.2        
# [7] alabaster.base_1.8.0      utf8_1.2.6                R6_2.6.1                 
# [10] HDF5Array_1.36.0          lazyeval_0.2.2            uwot_0.2.3               
# [13] rhdf5filters_1.20.0       withr_3.0.2               gridExtra_2.3            
# [16] progressr_0.15.1          cli_3.6.5                 spatstat.explore_3.4-3   
# [19] fastDummies_1.7.5         alabaster.se_1.8.0        labeling_0.4.3           
# [22] spatstat.data_3.1-6       ggridges_0.5.6            pbapply_1.7-2            
# [25] dichromat_2.0-0.1         parallelly_1.45.0         rstudioapi_0.17.1        
# [28] RSQLite_2.4.1             ica_1.0-3                 spatstat.random_3.4-1    
# [31] Matrix_1.7-3              ggbeeswarm_0.7.2          abind_1.4-8              
# [34] lifecycle_1.0.4           yaml_2.3.10               rhdf5_2.52.1             
# [37] SparseArray_1.8.0         BiocFileCache_2.16.0      Rtsne_0.17               
# [40] grid_4.5.0                blob_1.2.4                promises_1.3.3           
# [43] ExperimentHub_2.16.0      crayon_1.5.3              miniUI_0.1.2             
# [46] lattice_0.22-7            beachmat_2.24.0           cowplot_1.1.3            
# [49] KEGGREST_1.48.0           pillar_1.10.2             future.apply_1.20.0      
# [52] codetools_0.2-20          spacexr_2.2.1             glue_1.8.0               
# [55] spatstat.univar_3.1-3     data.table_1.17.4         vctrs_0.6.5              
# [58] png_0.1-8                 gypsum_1.4.0              spam_2.11-1              
# [61] gtable_0.3.6              cachem_1.1.0              S4Arrays_1.8.1           
# [64] mime_0.13                 survival_3.8-3            iterators_1.0.14         
# [67] fitdistrplus_1.2-2        ROCR_1.0-11               nlme_3.1-168             
# [70] bit64_4.6.0-1             alabaster.ranges_1.8.0    filelock_1.0.3           
# [73] RcppAnnoy_0.0.22          irlba_2.3.5.1             vipor_0.4.7              
# [76] KernSmooth_2.23-26        colorspace_2.1-1          DBI_1.2.3                
# [79] ggrastr_1.0.2             tidyselect_1.2.1          bit_4.6.0                
# [82] compiler_4.5.0            curl_6.3.0                httr2_1.1.2              
# [85] BiocNeighbors_2.2.0       h5mread_1.0.1             hdf5r_1.3.12             
# [88] DelayedArray_0.34.1       plotly_4.10.4             scales_1.4.0             
# [91] lmtest_0.9-40             rappdirs_0.3.3            digest_0.6.37            
# [94] goftest_1.2-3             spatstat.utils_3.1-4      alabaster.matrix_1.8.0   
# [97] XVector_0.48.0            htmltools_0.5.8.1         pkgconfig_2.0.3          
# [100] sparseMatrixStats_1.20.0  dbplyr_2.5.0              fastmap_1.2.0            
# [103] rlang_1.1.6               htmlwidgets_1.6.4         UCSC.utils_1.4.0         
# [106] shiny_1.10.0              DelayedMatrixStats_1.30.0 farver_2.1.2             
# [109] zoo_1.8-14                jsonlite_2.0.0            BiocParallel_1.42.1      
# [112] magrittr_2.0.3            GenomeInfoDbData_1.2.14   dotCall64_1.2            
# [115] patchwork_1.3.0           Rhdf5lib_1.30.0           Rcpp_1.0.14              
# [118] viridis_0.6.5             reticulate_1.42.0         stringi_1.8.7            
# [121] alabaster.schemas_1.8.0   MASS_7.3-65               AnnotationHub_3.16.0     
# [124] plyr_1.8.9                parallel_4.5.0            listenv_0.9.1            
# [127] ggrepel_0.9.6             deldir_2.0-4              Biostrings_2.76.0        
# [130] splines_4.5.0             tensor_1.5                hms_1.1.3                
# [133] igraph_2.1.4              spatstat.geom_3.4-1       RcppHNSW_0.6.0           
# [136] reshape2_1.4.4            BiocVersion_3.21.1        BiocManager_1.30.26      
# [139] tzdb_0.5.0                foreach_1.5.2             httpuv_1.6.16            
# [142] RANN_2.6.2                polyclip_1.10-7           scattermore_1.2          
# [145] xtable_1.8-4              RSpectra_0.16-2           later_1.4.2              
# [148] viridisLite_0.4.2         memoise_2.0.1             beeswarm_0.4.0           
# [151] AnnotationDbi_1.70.0      cluster_2.1.8.1           timechange_0.3.0         
# [154] globals_0.18.0           

