# Integrative analysis in Seurat v5
# https://satijalab.org/seurat/articles/seurat5_integration

# integration of single-cell sequencing datasets across experimental batches, donors or conditions
# can help match shared cell types and states across datasets 
   # boosts statistical power and facilitates accurate comparative analysis across datasets


# load libraries
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
# options(future.globals.maxSize = 1e9) # not large enough for RPCAIntegration
# The total size of the 10 globals exported for future expression ('FUN()') is 1.69 GiB.
options(future.globals.maxSize = 4 * 1024^3)  # helps prevent accidental overuse of memory in parallel jobs

# load data
# InstallData("pbmcsca")
obj <- LoadData("pbmcsca") # dataset of human PBMCs profiles with 7 different technologies
str(obj)
obj[["RNA"]]@layers$counts # raw, un-normalized counts
obj[["RNA"]]@layers$data # normalized counts


### remove low-quality cells
obj <- subset(obj, nFeature_RNA > 1000)

### obtain predicted cell annotations using Azimuth (to be used later for assessing integration)
# Warning: Overwriting miscellaneous data for model
# Warning: Adding a dimensional reduction (refUMAP) without the associated assay being present
# Warning: Adding a dimensional reduction (refUMAP) without the associated assay being present
# detected inputs from HUMAN with id type Gene.name
# reference rownames detected HUMAN with id type Gene.name
# Error in ValidateParams_FindTransferAnchors(reference = reference, query = query,  : 
#                                               Reference assay is SCT, but query assay is RNA. Mixing SCT and non-SCT in FindTransferAnchors is not supported.
# https://github.com/satijalab/seurat/issues/9871
# Hello! This bug has been fixed in the latest development version of Seurat; please try to reinstall Seurat from our GitHub repository (the main branch) and do let us know if the error arises again. SeuratObject does not need to be reinstalled. Thank you!

obj <- RunAzimuth(obj, reference = "pbmcref")
obj
# An object of class Seurat 
# 33789 features across 10434 samples within 4 assays 
# Active assay: RNA (33694 features, 0 variable features)
# 2 layers present: counts, data
# 3 other assays present: prediction.score.celltype.l1, prediction.score.celltype.l2, prediction.score.celltype.l3
# 2 dimensional reductions calculated: integrated_dr, ref.umap


# the object contains data from nine different batches representing seven different technologies
table(obj$Method)

# in previous versions of Seurat, the data would need to be represented as nine different Seurat objects
# in Seurat v5, the data can be kept in one object, simply split into layers (counts and data for all 9 batches)

# a standard scRNAseq analysis (i.e., without integration) can be run
# since the data is split into layers, normalization and variable feature identification is performed for each batch independently (a consensus set of variable features is automatically identified)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$Method)
obj
# An object of class Seurat 
# 33789 features across 10434 samples within 4 assays 
# Active assay: RNA (33694 features, 0 variable features)
# 18 layers present: counts.Smart-seq2, counts.CEL-Seq2, counts.10x_Chromium_v2_A, counts.10x_Chromium_v2_B, counts.10x_Chromium_v3, counts.Drop-seq, counts.Seq-Well, counts.inDrops, counts.10x_Chromium_v2, data.Smart-seq2, data.CEL-Seq2, data.10x_Chromium_v2_A, data.10x_Chromium_v2_B, data.10x_Chromium_v3, data.Drop-seq, data.Seq-Well, data.inDrops, data.10x_Chromium_v2
# 3 other assays present: prediction.score.celltype.l1, prediction.score.celltype.l2, prediction.score.celltype.l3
# 2 dimensional reductions calculated: integrated_dr, ref.umap

# standard Seurat workflow
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# visualize the results of the standard analysis without integration by batch and cell type annotation
   # cell type annotations were previously added by Azimuth
   # cells group by cell type and by method
   # clustering of the dataset just returns predominately batch-specific clusters (downstream analysis of this would be extremely difficult)
DimPlot(obj, reduction = "umap.unintegrated", group.by = c("Method", "predicted.celltype.l2"))


################################################
### PERFORM STREAMLINED INTEGRATIVE ANALYSIS ###
################################################

# Seurat allows streamlined integrative analysis using the IntegrateLayers() function
# Seurat currently supports five integration methods
   # these methods perform integration in low-dimensional space and return a dimensional reduction that aims to co-embed shared cell types across batches
# Anchor-based CCA integration (method=CCAIntegration)
# Anchor-based RPCA integration (method=RPCAIntegration)
# Harmony (method=HarmonyIntegration)
# FastMNN (method= FastMNNIntegration)
# scVI (method=scVIIntegration)

# scVI required reticualte and scvi-tools and its dependencies installed in a conda environment
# bash
# conda create -n scvi-env python=3.12  # any python 3.10 to 3.13
# conda activate scvi-env
# conda install scvi-tools -c conda-forge
# Error in py_module_import(module, convert = convert) : 
#   ModuleNotFoundError: No module named 'scanpy'
# conda install -c conda-forge scanpy
# conda env list # path to environment
library(reticulate)


### perform integration using the five available methods
obj <- IntegrateLayers(object = obj, method = CCAIntegration,
                       orig.reduction = "pca", new.reduction = "integrated.cca", 
                       verbose = FALSE)

obj <- IntegrateLayers(object = obj, method = RPCAIntegration, 
                       orig.reduction = "pca", new.reduction = "integrated.rpca", 
                       verbose = FALSE)

obj <- IntegrateLayers(object = obj, method = HarmonyIntegration, 
                       orig.reduction = "pca", new.reduction = "harmony", 
                       verbose = FALSE)

obj <- IntegrateLayers(object = obj, method = FastMNNIntegration, 
                       new.reduction = "integrated.mnn", verbose = FALSE)


# obj <- IntegrateLayers(object = obj, method = scVIIntegration, 
#                        new.reduction = "integrated.scvi",
#                        conda_env = "/Users/kristinvandenham/miniforge3/envs/scvi-env", 
#                        verbose = FALSE)
# R Session Aborted. R encountered a fatal error. The session was terminated.
# tried running the function using a smaller data set to see if it was a memory issue, but session still aborted
   # obj_small <- subset(obj, cells = sample(colnames(obj), 1000))
# Would need to patch scVIIntegration function in SeuratWrapper package to force PyTorch to run on CPU (accelerator="cpu")


### CCA integration
# cluster and visualize the datasets
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:30)
obj <- FindClusters(obj, resolution = 2, cluster.name = "cca_clusters")

obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
p1 <- DimPlot(obj, reduction = "umap.cca", 
              group.by = c("Method", "predicted.celltype.l2", "cca_clusters"), 
              combine = FALSE, label.size = 2)

### RPCA integration
# cluster and visualize the datasets
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:30)
obj <- FindClusters(obj, resolution = 2, cluster.name = "rpca_clusters")

obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
p2 <- DimPlot(obj, reduction = "umap.rpca", 
              group.by = c("Method", "predicted.celltype.l2", "rpca_clusters"), 
              combine = FALSE, label.size = 2)

### Harmony integration
# cluster and visualize the datasets
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj, resolution = 2, cluster.name = "harmony_clusters")

obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
p3 <- DimPlot(obj, reduction = "umap.harmony", 
              group.by = c("Method", "predicted.celltype.l2", "harmony_clusters"), 
              combine = FALSE, label.size = 2)

### MNN integration
# cluster and visualize the datasets
obj <- FindNeighbors(obj, reduction = "integrated.mnn", dims = 1:30)
obj <- FindClusters(obj, resolution = 2, cluster.name = "mnn_clusters")

obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")
p4 <- DimPlot(obj, reduction = "umap.mnn", 
              group.by = c("Method", "predicted.celltype.l2", "mnn_clusters"), 
              combine = FALSE, label.size = 2) 

wrap_plots(p4, ncol = 3, byrow = F, show.legend = FALSE) 


### compare the expression of biological markers based on different clustering solutions
# violin plots
p1 <- VlnPlot(obj, features = "rna_CD8A", group.by = "unintegrated_clusters") + 
  NoLegend() + ggtitle("CD8A - Unintegrated Clusters")
p2 <- VlnPlot(obj, features = "rna_CD8A", group.by = "cca_clusters") + 
  NoLegend() + ggtitle("CD8A - CCA Clusters")
p3 <- VlnPlot(obj, features = "rna_CD8A", group.by = "rpca_clusters") + 
  NoLegend() + ggtitle("CD8A - RPCA Clusters")

p1 | p2 | p3

# UMAP plot
p1 <- DimPlot(obj, reduction = "umap.unintegrated", group.by = "unintegrated_clusters")
p2 <- DimPlot(obj, reduction = "umap.cca", group.by = "cca_clusters")
p3 <- DimPlot(obj, reduction = "umap.rpca", group.by = "rpca_clusters")

p1 | p2 | p3


### rejoin the layers
# once the integrative analysis is complete, the layers can be rejoined
# this collapses the individual datasets together and recreates the original counts and data layers
# the layers need to be rejoined before you can perform any differential analysis
# layers can be resplit if you want to re-perform integrative analysis
obj <- JoinLayers(obj)
obj
# An object of class Seurat 
# 35789 features across 10434 samples within 5 assays 
# Active assay: RNA (33694 features, 2000 variable features)
# 3 layers present: data, counts, scale.data
# 4 other assays present: prediction.score.celltype.l1, prediction.score.celltype.l2, prediction.score.celltype.l3, mnn.reconstructed
# 12 dimensional reductions calculated: integrated_dr, ref.umap, pca, umap.unintegrated, integrated.cca, integrated.rpca, harmony, integrated.mnn, umap.cca, umap.rpca, umap.harmony, umap.mnn


#########################################
### INTEGRATION FOLLOWING SCTRANSFORM ###
#########################################

# run SCTransform normalization 
library(harmony)

obj <- SCTransform(obj)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- RunHarmony(obj, assay.use="SCT", group.by.vars = "Method")
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30)
DimPlot(obj, group.by = c("Method", "ident", "CellType"), ncol = 3)


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
#   [1] reticulate_1.42.0         harmony_1.2.3             Rcpp_1.0.14               future_1.58.0            
# [5] patchwork_1.3.0           ggplot2_3.5.2             Azimuth_0.5.0             shinyBS_0.61.1           
# [9] SeuratWrappers_0.4.0      stxBrain.SeuratData_0.1.2 ssHippo.SeuratData_3.1.4  pbmcsca.SeuratData_3.0.0 
# [13] pbmcref.SeuratData_1.0.0  pbmc3k.SeuratData_3.1.4   ifnb.SeuratData_3.1.0     SeuratData_0.2.2.9002    
# [17] Seurat_5.0.0              SeuratObject_5.1.0        sp_2.2-0                 
# 
# loaded via a namespace (and not attached):
#   [1] fs_1.6.6                          ProtGenerics_1.40.0               matrixStats_1.5.0                
# [4] spatstat.sparse_3.1-0             bitops_1.0-9                      DirichletMultinomial_1.50.0      
# [7] TFBSTools_1.46.0                  httr_1.4.7                        RColorBrewer_1.1-3               
# [10] doParallel_1.0.17                 tools_4.5.0                       sctransform_0.4.2                
# [13] R6_2.6.1                          DT_0.33                           lazyeval_0.2.2                   
# [16] uwot_0.2.3                        rhdf5filters_1.20.0               withr_3.0.2                      
# [19] gridExtra_2.3                     progressr_0.15.1                  cli_3.6.5                        
# [22] Biobase_2.68.0                    spatstat.explore_3.4-3            fastDummies_1.7.5                
# [25] EnsDb.Hsapiens.v86_2.99.0         shinyjs_2.1.0                     labeling_0.4.3                   
# [28] spatstat.data_3.1-6               ggridges_0.5.6                    pbapply_1.7-2                    
# [31] Rsamtools_2.24.0                  R.utils_2.13.0                    dichromat_2.0-0.1                
# [34] parallelly_1.45.0                 BSgenome_1.76.0                   rstudioapi_0.17.1                
# [37] RSQLite_2.4.1                     generics_0.1.4                    BiocIO_1.18.0                    
# [40] gtools_3.9.5                      ica_1.0-3                         spatstat.random_3.4-1            
# [43] googlesheets4_1.1.1               dplyr_1.1.4                       Matrix_1.7-3                     
# [46] S4Vectors_0.46.0                  abind_1.4-8                       R.methodsS3_1.8.2                
# [49] lifecycle_1.0.4                   yaml_2.3.10                       SummarizedExperiment_1.38.1      
# [52] rhdf5_2.52.1                      SparseArray_1.8.0                 Rtsne_0.17                       
# [55] grid_4.5.0                        blob_1.2.4                        promises_1.3.3                   
# [58] shinydashboard_0.7.3              crayon_1.5.3                      pwalign_1.4.0                    
# [61] miniUI_0.1.2                      lattice_0.22-7                    cowplot_1.1.3                    
# [64] GenomicFeatures_1.60.0            KEGGREST_1.48.0                   pillar_1.10.2                    
# [67] GenomicRanges_1.60.0              rjson_0.2.23                      future.apply_1.20.0              
# [70] codetools_0.2-20                  fastmatch_1.1-6                   spacexr_2.2.1                    
# [73] leiden_0.4.3.1                    glue_1.8.0                        spatstat.univar_3.1-3            
# [76] data.table_1.17.4                 remotes_2.5.0                     vctrs_0.6.5                      
# [79] png_0.1-8                         spam_2.11-1                       cellranger_1.1.0                 
# [82] gtable_0.3.6                      cachem_1.1.0                      Signac_1.14.0                    
# [85] S4Arrays_1.8.1                    mime_0.13                         survival_3.8-3                   
# [88] gargle_1.5.2                      RcppRoll_0.3.1                    iterators_1.0.14                 
# [91] fitdistrplus_1.2-2                ROCR_1.0-11                       nlme_3.1-168                     
# [94] bit64_4.6.0-1                     RcppAnnoy_0.0.22                  GenomeInfoDb_1.44.0              
# [97] irlba_2.3.5.1                     KernSmooth_2.23-26                SeuratDisk_0.0.0.9021            
# [100] colorspace_2.1-1                  seqLogo_1.74.0                    BiocGenerics_0.54.0              
# [103] DBI_1.2.3                         tidyselect_1.2.1                  bit_4.6.0                        
# [106] compiler_4.5.0                    curl_6.3.0                        hdf5r_1.3.12                     
# [109] DelayedArray_0.34.1               plotly_4.10.4                     rtracklayer_1.68.0               
# [112] scales_1.4.0                      caTools_1.18.3                    lmtest_0.9-40                    
# [115] rappdirs_0.3.3                    stringr_1.5.1                     digest_0.6.37                    
# [118] goftest_1.2-3                     presto_1.0.0                      spatstat.utils_3.1-4             
# [121] XVector_0.48.0                    htmltools_0.5.8.1                 pkgconfig_2.0.3                  
# [124] MatrixGenerics_1.20.0             fastmap_1.2.0                     ensembldb_2.32.0                 
# [127] rlang_1.1.6                       htmlwidgets_1.6.4                 UCSC.utils_1.4.0                 
# [130] shiny_1.10.0                      farver_2.1.2                      zoo_1.8-14                       
# [133] jsonlite_2.0.0                    BiocParallel_1.42.1               R.oo_1.27.1                      
# [136] RCurl_1.98-1.17                   magrittr_2.0.3                    GenomeInfoDbData_1.2.14          
# [139] dotCall64_1.2                     Rhdf5lib_1.30.0                   stringi_1.8.7                    
# [142] MASS_7.3-65                       plyr_1.8.9                        parallel_4.5.0                   
# [145] listenv_0.9.1                     ggrepel_0.9.6                     deldir_2.0-4                     
# [148] Biostrings_2.76.0                 splines_4.5.0                     tensor_1.5                       
# [151] BSgenome.Hsapiens.UCSC.hg38_1.4.5 igraph_2.1.4                      spatstat.geom_3.4-1              
# [154] RcppHNSW_0.6.0                    reshape2_1.4.4                    stats4_4.5.0                     
# [157] TFMPvalue_0.0.9                   XML_3.99-0.18                     BiocManager_1.30.26              
# [160] JASPAR2020_0.99.10                foreach_1.5.2                     httpuv_1.6.16                    
# [163] RANN_2.6.2                        tidyr_1.3.1                       purrr_1.0.4                      
# [166] polyclip_1.10-7                   scattermore_1.2                   rsvd_1.0.5                       
# [169] xtable_1.8-4                      restfulr_0.0.15                   AnnotationFilter_1.32.0          
# [172] RSpectra_0.16-2                   later_1.4.2                       googledrive_2.1.1                
# [175] viridisLite_0.4.2                 tibble_3.3.0                      memoise_2.0.1                    
# [178] AnnotationDbi_1.70.0              GenomicAlignments_1.44.0          IRanges_2.42.0                   
# [181] cluster_2.1.8.1                   globals_0.18.0   



################################################################################


# scRNAseq integration: overview of comparative analyses
# integrate control and IFNb-stimulated human PBMCs to jointly identify cell subpopulations across datasets and explore how each group differs across conditions
# in Seurat v5, all data is kept in one object and simply split into multiple layers
# https://satijalab.org/seurat/articles/integration_introduction.html

# load libraries
library(Seurat)
library(SeuratData)
library(patchwork)
library(tidyverse)

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# load data
ifnb <- LoadData("ifnb") # IFNb-stimulated and control human PBMCs 
table(ifnb@meta.data$stim)
# CTRL STIM 
# 6548 7451 

# split the RNA measurements into two layers: one for control cells and one for IFNb-stimulated cells
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb
# An object of class Seurat 
# 14053 features across 13999 samples within 1 assay 
# Active assay: RNA (14053 features, 0 variable features)
# 4 layers present: counts.CTRL, counts.STIM, data.CTRL, data.STIM


############################################
### PERFORM ANALYSIS WITHOUT INTEGRATION ###
############################################

# run standard Seurat analysis workflow
ifnb <- ifnb %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30, reduction = "pca") %>%
  FindClusters(resolution = 2, cluster.name = "unintegrated_clusters")

# plot UMAP of unitegrated data
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
plot_unint <- DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))


#########################################
### PERFORM ANALYSIS WITH INTEGRATION ###
#########################################

# integrate the data from the two conditions so that the same cell types will cluster together
# returns a single dimensional reduction that captures the shared sources of variance across multiple layers (i.e., cells in similar biological states will cluster)
# the method returns a dimensional reduction that can be used fro visualization and unsupervised clustering analysis
# performance can be evaluated using cell type labels that are pre-loaded in the seurat_annotations metadata column
table(ifnb@meta.data$seurat_annotations)

# perform CCA integration
ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", 
                        new.reduction = "integrated.cca", verbose = FALSE)

# re-join layers after integration
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

# perform clustering
ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)

# visualize the integration
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))

# visualize the UMAP of control and stimulated cells side-by-side
DimPlot(ifnb, reduction = "umap", split.by = "stim")


############################################
### IDENTIFY CONSERVED CELL TYPE MARKERS ###
############################################

# FindConservedMarkers() function can be used to identify canonical cell type marker genes that are conserved across conditions
   # differential gene expression testing is performed for each dataset/group and the p-values combined using meta-analysis methods implemented in the MetaDE package

Idents(ifnb) <- "seurat_annotations"
nk.markers <- FindConservedMarkers(ifnb, ident.1 = "NK", grouping.var = "stim", verbose = FALSE)
head(nk.markers)


# the same analysis can be performed on the unsupervised clustering results (stored in seurat_clusters)
# and these conserved markers can be used to annotate the cell types in the dataset
Idents(ifnb) <- factor(Idents(ifnb), levels = c("pDC", "Eryth", "Mk", "DC", "CD14 Mono", 
                                                "CD16 Mono","B Activated", "B", "CD8 T", "NK", 
                                                "T activated", "CD4 Naive T", "CD4 Memory T"))

markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", 
                     "NKG7", "CCL5", "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", 
                     "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", "GPR183", "PPBP", 
                     "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")

# view the conserved cell type markers across conditions
   # showing both expression level and the percentage of cells in a cluster expressing any given gene
# plot 2-3 strong markers fro each of the 13 clusters
DotPlot(ifnb, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, 
        split.by = "stim") + RotatedAxis()


#####################################################################
### IDENTIFY DIFFERENTIALLY EXPRESSED GENES ACROSS THE CONDITIONS ###
#####################################################################

# load libraries
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# aggregate cells of a similar type and condition to create a pseudobulk profile
aggregate_ifnb <- AggregateExpression(ifnb, group.by = c("seurat_annotations", "stim"), 
                                      return.seurat = TRUE)

# compare the pseudobulk profiles of naive CD4 T cells and CD14 monocytes, and their gene expression profiles before and after stimulation
genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", 
                   "IFIT2", "IFIT1", "CXCL10", "CCL8")

p1 <- CellScatter(aggregate_ifnb, "CD14 Mono_CTRL", "CD14 Mono_STIM", 
                  highlight = genes.to.label)
p2 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

p3 <- CellScatter(aggregate_ifnb, "CD4 Naive T_CTRL", "CD4 Naive T_STIM", 
                  highlight = genes.to.label)
p4 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE)

p2 + p4


### determine which genes change in different conditions for cells of the same cell type
# create a column in the metadata slot to hold both cell type and stimulation information
ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")

# switch Idents to that column
Idents(ifnb) <- "celltype.stim"

# use FindMarkers() to find the genes that are different between the stimulated and control B cells
b.interferon.response <- FindMarkers(ifnb, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
head(b.interferon.response, n = 15)
# p-values obtained from this analysis should be interpreted with caution -> these tests treat each cells as an independent replicate and ignore the inherent correlations between cells originating from the same sample
# DE tests across multiple conditions should expressly utilize multiple samples/replicates and be performed after pseudobulking cells from the same sample and subpopulation together
# in this dataset there is only a single replicate 


### changes in gene expression can also be visualized using FeaturePlot() and VlnPlot()

# plot FeaturePlots and VlnPlots using a list of genes/features, split by a grouping variable (stimulation)
# some canonical cell type markers (CD3D) are largely unaffected by IFNb stimulation
# core interferon genes (IFI6) are upregulated after IFNb stimulation
# certain markers are only upregulated on certain cell types (e.g., CXCL10 on monocytes and B cells)
FeaturePlot(ifnb, features = c("CD3D", "IFI6", "CXCL10"), split.by = "stim", 
            max.cutoff = 3, cols = c("grey", "red"), reduction = "umap")

plots <- VlnPlot(ifnb, features = c("CD3D", "IFI6", "CXCL10"), split.by = "stim", 
                 group.by = "seurat_annotations", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)


################################################################
### PERFORM INTEGRATION WITH SCTRANSFORM-NORMALIZED DATASETS ###
################################################################

# SCTransform can be used rather than log-normalization

# load libraries
library(harmony)

# load data
ifnb <- LoadData("ifnb")

# split datasets and process without integration
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)

# perform SCTransform workflow steps
ifnb <- SCTransform(ifnb)
ifnb <- RunPCA(ifnb)
ifnb <- RunUMAP(ifnb, dims = 1:30)
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))

# integrate datasets and re-perform workflow steps
ifnb <- RunHarmony(ifnb, assay.use="SCT", group.by.vars = "stim")
ifnb <- FindNeighbors(ifnb, reduction = "harmony", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 0.6)
ifnb <- RunUMAP(ifnb, reduction = "harmony", dims = 1:30)
DimPlot(ifnb, group.by = c("stim", "seurat_annotations"), ncol = 3)
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))

# perform differential expression
ifnb <- PrepSCTFindMarkers(ifnb)
ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
Idents(ifnb) <- "celltype.stim"
b.interferon.response <- FindMarkers(ifnb, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
head(b.interferon.response)


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
#   [1] harmony_1.2.3             Rcpp_1.0.14               cowplot_1.1.3            
# [4] future_1.58.0             lubridate_1.9.4           forcats_1.0.0            
# [7] stringr_1.5.1             dplyr_1.1.4               purrr_1.0.4              
# [10] readr_2.1.5               tidyr_1.3.1               tibble_3.3.0             
# [13] ggplot2_3.5.2             tidyverse_2.0.0           patchwork_1.3.0          
# [16] stxBrain.SeuratData_0.1.2 ssHippo.SeuratData_3.1.4  pbmcsca.SeuratData_3.0.0 
# [19] pbmcref.SeuratData_1.0.0  pbmc3k.SeuratData_3.1.4   ifnb.SeuratData_3.1.0    
# [22] SeuratData_0.2.2.9002     Seurat_5.0.0              SeuratObject_5.1.0       
# [25] sp_2.2-0                 
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22            splines_4.5.0               later_1.4.2                
# [4] polyclip_1.10-7             fastDummies_1.7.5           lifecycle_1.0.4            
# [7] Rdpack_2.6.4                doParallel_1.0.17           globals_0.18.0             
# [10] lattice_0.22-7              MASS_7.3-65                 magrittr_2.0.3             
# [13] limma_3.64.1                plotly_4.10.4               plotrix_3.8-4              
# [16] qqconf_1.3.2                httpuv_1.6.16               sn_2.1.1                   
# [19] glmGamPoi_1.20.0            sctransform_0.4.2           spam_2.11-1                
# [22] spatstat.sparse_3.1-0       reticulate_1.42.0           pbapply_1.7-2              
# [25] RColorBrewer_1.1-3          multcomp_1.4-28             abind_1.4-8                
# [28] GenomicRanges_1.60.0        Rtsne_0.17                  presto_1.0.0               
# [31] BiocGenerics_0.54.0         TH.data_1.1-3               rappdirs_0.3.3             
# [34] sandwich_3.1-1              GenomeInfoDbData_1.2.14     IRanges_2.42.0             
# [37] S4Vectors_0.46.0            ggrepel_0.9.6               irlba_2.3.5.1              
# [40] listenv_0.9.1               spatstat.utils_3.1-4        TFisher_0.2.0              
# [43] goftest_1.2-3               RSpectra_0.16-2             spatstat.random_3.4-1      
# [46] fitdistrplus_1.2-2          parallelly_1.45.0           DelayedMatrixStats_1.30.0  
# [49] DelayedArray_0.34.1         leiden_0.4.3.1              codetools_0.2-20           
# [52] tidyselect_1.2.1            UCSC.utils_1.4.0            farver_2.1.2               
# [55] matrixStats_1.5.0           stats4_4.5.0                spatstat.explore_3.4-3     
# [58] mathjaxr_1.8-0              jsonlite_2.0.0              multtest_2.64.0            
# [61] progressr_0.15.1            ggridges_0.5.6              survival_3.8-3             
# [64] iterators_1.0.14            foreach_1.5.2               tools_4.5.0                
# [67] ica_1.0-3                   glue_1.8.0                  SparseArray_1.8.0          
# [70] mnormt_2.1.1                gridExtra_2.3               metap_1.12                 
# [73] MatrixGenerics_1.20.0       GenomeInfoDb_1.44.0         withr_3.0.2                
# [76] numDeriv_2016.8-1.1         BiocManager_1.30.26         fastmap_1.2.0              
# [79] digest_0.6.37               timechange_0.3.0            R6_2.6.1                   
# [82] mime_0.13                   scattermore_1.2             tensor_1.5                 
# [85] dichromat_2.0-0.1           spatstat.data_3.1-6         RhpcBLASctl_0.23-42        
# [88] generics_0.1.4              data.table_1.17.4           S4Arrays_1.8.1             
# [91] httr_1.4.7                  htmlwidgets_1.6.4           uwot_0.2.3                 
# [94] pkgconfig_2.0.3             gtable_0.3.6                lmtest_0.9-40              
# [97] XVector_0.48.0              htmltools_0.5.8.1           dotCall64_1.2              
# [100] scales_1.4.0                Biobase_2.68.0              png_0.1-8                  
# [103] spatstat.univar_3.1-3       rstudioapi_0.17.1           tzdb_0.5.0                 
# [106] reshape2_1.4.4              nlme_3.1-168                zoo_1.8-14                 
# [109] spacexr_2.2.1               KernSmooth_2.23-26          parallel_4.5.0             
# [112] miniUI_0.1.2                vipor_0.4.7                 ggrastr_1.0.2              
# [115] pillar_1.10.2               grid_4.5.0                  vctrs_0.6.5                
# [118] RANN_2.6.2                  promises_1.3.3              beachmat_2.24.0            
# [121] xtable_1.8-4                cluster_2.1.8.1             beeswarm_0.4.0             
# [124] mvtnorm_1.3-3               cli_3.6.5                   compiler_4.5.0             
# [127] rlang_1.1.6                 crayon_1.5.3                future.apply_1.20.0        
# [130] mutoss_0.1-13               labeling_0.4.3              plyr_1.8.9                 
# [133] ggbeeswarm_0.7.2            stringi_1.8.7               viridisLite_0.4.2          
# [136] deldir_2.0-4                lazyeval_0.2.2              spatstat.geom_3.4-1        
# [139] Matrix_1.7-3                RcppHNSW_0.6.0              hms_1.1.3                  
# [142] sparseMatrixStats_1.20.0    statmod_1.5.0               shiny_1.10.0               
# [145] SummarizedExperiment_1.38.1 rbibutils_2.3               ROCR_1.0-11                
# [148] igraph_2.1.4  

