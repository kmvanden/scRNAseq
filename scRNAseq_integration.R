# Single-cell RNAseq Integration to Correct for Batch Effects

# Canonical Correlation Analysis 
# https://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_integration.R
# https://satijalab.org/seurat/articles/integration_introduction.html

# Harmony
# https://portals.broadinstitute.org/harmony/articles/quickstart.html
# https://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_integrate_harmony.R
# https://satijalab.org/seurat/articles/seurat5_integration

### integration can be used for
# integrating multiple scRNAseq datasets and correct for batch effects
# cell label transfer: transfer type classification from a reference to a query dataset
# integration of multi-modal cell data (e.g., scRNAseq and scATACseq): integrate signals collected from seperate assays into a single-cell multiomics dataset

# horizontal integration: assays anchored by common gene set (e.g., different patients)
# vertical integration: assays are anchored by the cells (e.g, different modalities collected)


######################################
### Canonical Correlation Analysis ###
######################################


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

### get data - filtered matrices from GSE180665
# for i in *.gz; do tar -xvzf $i; done

# get data location
dirs <- list.dirs(path = "GSE180665/", recursive = F, full.names = F)

# create a Seurat object for each sample
for(x in dirs){
  name <- gsub('_filtered_feature_bc_matrix','', x)
  
  # get the matrix, features and barcode files from each folder
  cts <- ReadMtx(mtx = paste0('GSE180665/',x,'/matrix.mtx.gz'),
                 features = paste0('GSE180665/',x,'/features.tsv.gz'),
                 cells = paste0('GSE180665/',x,'/barcodes.tsv.gz'))
  
  # create Seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}

# merge the Seurat objects
ls() # 3:9 are the cell ids
# [1] "cts"             "dirs"            "HB17_background" "HB17_PDX"       
# [5] "HB17_tumor"      "HB30_PDX"        "HB30_tumor"      "HB53_background"
# [9] "HB53_tumor"      "merged_seurat"   "name"            "x"    
merged_seurat <- merge(HB17_background, 
                       y = c(HB17_PDX, HB17_tumor, HB30_PDX, HB30_tumor, 
                             HB53_background, HB53_tumor),
                       add.cell.ids = ls()[3:9],
                       project = 'HB')
merged_seurat
# An object of class Seurat 
# 33538 features across 77936 samples within 1 assay 
# Active assay: RNA (33538 features, 0 variable features)
# 7 layers present: counts.1, counts.2, counts.3, counts.4, counts.5, counts.6, counts.7


### quality control and filtering 
head(merged_seurat@meta.data)

# create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)
# split sample column into patient id, type and barcode by underscore
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Patient', 'Type', 'Barcode'), 
                                    sep = '_')
table(merged_seurat$Patient)
table(merged_seurat$Type)
length(merged_seurat$Barcode) # 77936


# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')

# filtering (using the thresholds used in the paper)
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
                                   nFeature_RNA > 500 &
                                   mitoPercent < 10)

length(merged_seurat_filtered$Barcode) # 67851


# perform standard workflow steps
merged_seurat_filtered <- merged_seurat_filtered %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>% 
  RunPCA()

ElbowPlot(merged_seurat_filtered, ndims = 50)

merged_seurat_filtered <- merged_seurat_filtered %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:20)

# plot data by Patient and by tissue type
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Patient')
p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Type',
              cols = c('red','green','blue'))

p1/p2 # the cells coming from different patients cluster differently (technical differences), masking the biological variations


### perform integration to correct for batch effects
# split Seurat object by patient (batch effects are coming from patients)
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Patient') # object list for each patient

# normalize and find variable features for each object in the list
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}

# select integration features from object list
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA) to integrate the data across the patients
# CCA method is computationally intensive
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)


# scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- seurat.integrated %>%
  ScaleData() %>%
  RunPCA() %>% 
  RunUMAP(dims = 1:50)

p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Patient')
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Type',
              cols = c('red','green','blue'))

# compare plots from before and after integration
grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)


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
#   [1] gridExtra_2.3      lubridate_1.9.4    forcats_1.0.0      stringr_1.5.1     
# [5] dplyr_1.1.4        purrr_1.0.4        readr_2.1.5        tidyr_1.3.1       
# [9] tibble_3.3.0       tidyverse_2.0.0    ggplot2_3.5.2      Seurat_5.3.0      
# [13] SeuratObject_5.1.0 sp_2.2-0          
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-4           pbapply_1.7-2          rlang_1.1.6            magrittr_2.0.3        
# [5] RcppAnnoy_0.0.22       matrixStats_1.5.0      ggridges_0.5.6         compiler_4.5.0        
# [9] spatstat.geom_3.4-1    png_0.1-8              vctrs_0.6.5            reshape2_1.4.4        
# [13] pkgconfig_2.0.3        fastmap_1.2.0          promises_1.3.3         tzdb_0.5.0            
# [17] jsonlite_2.0.0         goftest_1.2-3          later_1.4.2            spatstat.utils_3.1-4  
# [21] irlba_2.3.5.1          parallel_4.5.0         cluster_2.1.8.1        R6_2.6.1              
# [25] ica_1.0-3              stringi_1.8.7          RColorBrewer_1.1-3     spatstat.data_3.1-6   
# [29] reticulate_1.42.0      parallelly_1.45.0      spatstat.univar_3.1-3  lmtest_0.9-40         
# [33] scattermore_1.2        iterators_1.0.14       Rcpp_1.0.14            tensor_1.5            
# [37] future.apply_1.20.0    zoo_1.8-14             sctransform_0.4.2      timechange_0.3.0      
# [41] httpuv_1.6.16          Matrix_1.7-3           splines_4.5.0          igraph_2.1.4          
# [45] tidyselect_1.2.1       rstudioapi_0.17.1      dichromat_2.0-0.1      abind_1.4-8           
# [49] doParallel_1.0.17      codetools_0.2-20       spatstat.random_3.4-1  miniUI_0.1.2          
# [53] spatstat.explore_3.4-3 listenv_0.9.1          lattice_0.22-7         plyr_1.8.9            
# [57] withr_3.0.2            shiny_1.10.0           ROCR_1.0-11            Rtsne_0.17            
# [61] future_1.58.0          fastDummies_1.7.5      survival_3.8-3         polyclip_1.10-7       
# [65] fitdistrplus_1.2-2     pillar_1.10.2          spacexr_2.2.1          KernSmooth_2.23-26    
# [69] foreach_1.5.2          plotly_4.10.4          generics_0.1.4         RcppHNSW_0.6.0        
# [73] hms_1.1.3              scales_1.4.0           globals_0.18.0         xtable_1.8-4          
# [77] glue_1.8.0             lazyeval_0.2.2         tools_4.5.0            data.table_1.17.4     
# [81] RSpectra_0.16-2        RANN_2.6.2             dotCall64_1.2          cowplot_1.1.3         
# [85] grid_4.5.0             colorspace_2.1-1       nlme_3.1-168           patchwork_1.3.0       
# [89] cli_3.6.5              spatstat.sparse_3.1-0  spam_2.11-1            viridisLite_0.4.2     
# [93] uwot_0.2.3             gtable_0.3.6           digest_0.6.37          progressr_0.15.1      
# [97] ggrepel_0.9.6          htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.8.1     
# [101] lifecycle_1.0.4        httr_1.4.7             mime_0.13              MASS_7.3-65


###########################
### HARMONY INTEGRATION ###
###########################

# integration across conditions using Harmony
# Harmony only calculates the corrected dimensionally reduced values (embeddings)
# raw data, data and scaled data slots are not changed

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# set seed
set.seed(1234)

# load libraries
library(harmony)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(ggplot2)


### load dataset
# InstallData("ifnb") # IFNb stimulated PBMC dataset
ifnb <- LoadData("ifnb")
str(ifnb)


### quality control and filtering
# calculate the percent of mitochondrial reads
ifnb$mito.percent <- PercentageFeatureSet(ifnb, pattern = '^MT-')
table(ifnb@meta.data$mito.percent)
# 0 
# 13999 | data already filtered

# filtering using thresholds from the tutorial
length(ifnb@meta.data$orig.ident) # 13999
ifnb.filtered <- subset(ifnb, subset = nCount_RNA > 800 &
                          nFeature_RNA > 200 & 
                          mito.percent < 5)
length(ifnb@meta.data$orig.ident) # 13999 | data already filtered


# standard Seurat workflow steps
ifnb.filtered <- ifnb.filtered %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(ifnb.filtered, ndims = 50)

ifnb.filtered <- RunUMAP(ifnb.filtered, dims = 1:20, reduction = 'pca')

before <- DimPlot(ifnb.filtered, reduction = 'umap', group.by = 'stim') # before Harmony integration
before # batch effects by condition observed

### run Harmony integration
# returns corrected dimensionality reductions -> embeddings
ifnb.harmony <- ifnb.filtered %>%
  RunHarmony(group.by.vars = 'stim', plot_convergence = FALSE)

ifnb.harmony@reductions # slots for reductions: pca, umap and harmony

# get Harmony embeddings
ifnb.harmony.embed <- Embeddings(ifnb.harmony, "harmony") 
ifnb.harmony.embed[1:10,1:10]


# create UMAP and clustering using Harmony embeddings
ifnb.harmony <- ifnb.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# visualize --> Harmony integration corrected batch effects
after <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'stim')
before/after


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
# [13] ssHippo.SeuratData_3.1.4  pbmc3k.SeuratData_3.1.4   ifnb.SeuratData_3.1.0    
# [16] SeuratData_0.2.2.9002     Seurat_5.3.0              SeuratObject_5.1.0       
# [19] sp_2.2-0                  harmony_1.2.3             Rcpp_1.0.14              
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     rstudioapi_0.17.1      jsonlite_2.0.0        
# [4] magrittr_2.0.3         spatstat.utils_3.1-4   farver_2.1.2          
# [7] vctrs_0.6.5            ROCR_1.0-11            spatstat.explore_3.4-3
# [10] htmltools_0.5.8.1      sctransform_0.4.2      parallelly_1.45.0     
# [13] KernSmooth_2.23-26     htmlwidgets_1.6.4      ica_1.0-3             
# [16] plyr_1.8.9             plotly_4.10.4          zoo_1.8-14            
# [19] igraph_2.1.4           mime_0.13              lifecycle_1.0.4       
# [22] iterators_1.0.14       pkgconfig_2.0.3        Matrix_1.7-3          
# [25] R6_2.6.1               fastmap_1.2.0          fitdistrplus_1.2-2    
# [28] shiny_1.10.0           digest_0.6.37          colorspace_2.1-1      
# [31] patchwork_1.3.0        tensor_1.5             RSpectra_0.16-2       
# [34] irlba_2.3.5.1          labeling_0.4.3         progressr_0.15.1      
# [37] spatstat.sparse_3.1-0  timechange_0.3.0       httr_1.4.7            
# [40] polyclip_1.10-7        abind_1.4-8            compiler_4.5.0        
# [43] withr_3.0.2            doParallel_1.0.17      fastDummies_1.7.5     
# [46] MASS_7.3-65            rappdirs_0.3.3         tools_4.5.0           
# [49] lmtest_0.9-40          httpuv_1.6.16          future.apply_1.20.0   
# [52] goftest_1.2-3          glue_1.8.0             nlme_3.1-168          
# [55] promises_1.3.3         grid_4.5.0             Rtsne_0.17            
# [58] cluster_2.1.8.1        reshape2_1.4.4         generics_0.1.4        
# [61] gtable_0.3.6           spatstat.data_3.1-6    tzdb_0.5.0            
# [64] data.table_1.17.4      hms_1.1.3              spatstat.geom_3.4-1   
# [67] RcppAnnoy_0.0.22       ggrepel_0.9.6          RANN_2.6.2            
# [70] foreach_1.5.2          pillar_1.10.2          spam_2.11-1           
# [73] RcppHNSW_0.6.0         later_1.4.2            splines_4.5.0         
# [76] lattice_0.22-7         survival_3.8-3         deldir_2.0-4          
# [79] tidyselect_1.2.1       miniUI_0.1.2           pbapply_1.7-2         
# [82] gridExtra_2.3          scattermore_1.2        RhpcBLASctl_0.23-42   
# [85] spacexr_2.2.1          matrixStats_1.5.0      stringi_1.8.7         
# [88] lazyeval_0.2.2         codetools_0.2-20       cli_3.6.5             
# [91] uwot_0.2.3             xtable_1.8-4           reticulate_1.42.0     
# [94] dichromat_2.0-0.1      globals_0.18.0         spatstat.random_3.4-1 
# [97] png_0.1-8              spatstat.univar_3.1-3  parallel_4.5.0        
# [100] dotCall64_1.2          listenv_0.9.1          viridisLite_0.4.2     
# [103] scales_1.4.0           ggridges_0.5.6         crayon_1.5.3          
# [106] rlang_1.1.6            cowplot_1.1.3   

