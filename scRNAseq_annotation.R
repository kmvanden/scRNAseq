# Automatic cell type annotation for ssRNAseq data
# https://github.com/kpatel427/YouTubeTutorials/blob/main/annotateSingleR.R
# https://github.com/kpatel427/YouTubeTutorials/blob/main/annotateSingleR_multipleRefs.R


# marker-based annotation: labels cells or cell clusters on the basis of the characteristic expression of known marker genes
   # known relationships between marker genes and cell types can be obtained from databases
   # MSigDB, PanglaoDB, CellMarker

# reference-based annotation: transfers labels from a reference cell or cluster (from an expertly annotated scRNAseq dataset) to a sufficiently similar one in the query (the data to be annotated)
   # GEO, Single Cell Expression Atlas, cell atlas projects


# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# load libraries
library(SingleR) # reference-based annotation approach (correlates to reference gene expression)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)

# load data from 10X CellRanger (20k Human PBMCs from a healthy female donor)
hdf5_obj <- Read10X_h5(filename = "20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",
                       use.names = TRUE, unique.features = TRUE)
pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj) # create Seurat object


### quality control and filtering
pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = "^MT-")
VlnPlot(pbmc.seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 & 
                                 nFeature_RNA > 500 & mitoPercent < 10)

### standard pre-processing workflow
pbmc.seurat.filtered <- pbmc.seurat.filtered %>%
  NormalizeData() %>% 
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(pbmc.seurat.filtered) # dimensionality of the data

pbmc.seurat.filtered <- pbmc.seurat.filtered %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:20)

head(pbmc.seurat.filtered@meta.data) # seurat clusters
DimPlot(pbmc.seurat.filtered, reduction = "umap")


##############################################
### USE SINGLER WITH ONE REFERENCE DATASET ###
##############################################


### get reference dataset
# choose a reference that has similar cell types and that was generated in a similar way as the query dataset
# expression values in the reference datasets are log normalized counts
ref <- celldex::HumanPrimaryCellAtlasData()
unique(ref$label.fine) # cell types included in the reference dataset


# get counts data from the query data set
pbmc_counts <- GetAssayData(pbmc.seurat.filtered, layer = "counts")

# run Single R on default mode (performs annotation on single cells )
# use label.main column of the reference dataset
pred <- SingleR(test = pbmc_counts, ref = ref, labels = ref$label.main)
head(pred)

# save the labels in the Seurat object
pbmc.seurat.filtered$singleR.labels <- pred$labels[match(rownames(pbmc.seurat.filtered@meta.data), rownames(pred))]
# visualize the labels
DimPlot(pbmc.seurat.filtered, reduction = "umap", group.by = "singleR.labels")
table(pbmc.seurat.filtered$singleR.labels)


### annotation diagnostics (how well SingleR performed annotation)
pred$scores # prediction scores (how ambiguous the assignment was)
plotScoreHeatmap(pred)

# identify poor quality or ambiguous cell identity assignments
plotDeltaDistribution(pred)

# compare cell type assignments with the results of unsupervised clustering 
tab <- table(Assigned = pred$labels, Clusters=pbmc.seurat.filtered$seurat_clusters)
table(tab)
pheatmap(log10(tab + 10), color = colorRampPalette(c("white","blue"))(10))
 

##########################################
### USE SINGLER WITH MULTIPLE DATASETS ###
##########################################


### load two reference datasets
hpca <- celldex::HumanPrimaryCellAtlasData()
dice <- celldex::DatabaseImmuneCellExpressionData()


### STRATEGY 1: combining the references, but keeping the cell types  separate between the references
# cell type assigned by the highest scoring label from either reference
# cell type labels
hpca$label.main
dice$label.main

# add which reference database the label is from to the cell type label
hpca$label.main <- paste0("HPCA.", hpca$label.main)
dice$label.main <- paste0("DICE.", dice$label.main)

# create a combined reference database based on shared genes
shared <- intersect(rownames(hpca), rownames(dice))
combined <- cbind(hpca[shared,], dice[shared,])
table(combined$label.main)

# save the counts into a separate object
pbmc_counts <- GetAssayData(pbmc.seurat.filtered, slot = "counts")

# run singleR using the combined reference database on the counts data
com.res1 <- SingleR(test = pbmc_counts, ref = combined, labels = combined$label.main)
table(com.res1$labels) # check the number of cells assigned to each label

# save the labels into a metadata column
pbmc.seurat.filtered$com.res1.labels <- com.res1[match(rownames(pbmc.seurat.filtered@meta.data), rownames(com.res1)), "labels"]
head(pbmc.seurat.filtered@meta.data)

DimPlot(pbmc.seurat.filtered, reduction = "umap", group.by = "com.res1.labels", label = TRUE)


### STRATEGY 2: Combining the scores across the reference databases

# remove the reference dataset label from the cell type label
hpca$label.main <- gsub("HPCA\\.","", hpca$label.main)
dice$label.main <- gsub("DICE\\.","", dice$label.main)

# run singleR using the reference dataset separately on the counts data
com.res2 <- SingleR(test = pbmc_counts, ref = list(HPCA = hpca, DICE = dice),
                    labels = list(hpca$label.main, dice$label.main))

table(com.res2$labels) # examine the label assignment (doesn't deal with inconsistency in the cell type labelling)

# determine which reference dataset scored scored the best for each assignment
grouping <- paste0(com.res2$labels,".", com.res2$reference)
best_ref <- as.data.frame(split(com.res2, grouping))
View(best_ref)

# get the marker genes for each cell type from each reference dataset
metadata(com.res2$orig.results$HPCA)$de.genes
metadata(com.res2$orig.results$DICE)$de.genes

# create heatmap to visualize the scores for each cell type assignment 
# from the individual reference datasets and from the combined reference (re-computed scores)
plotScoreHeatmap(com.res2)

# save the labels into a metadata column
pbmc.seurat.filtered$com.res2.labels <- com.res2[match(rownames(pbmc.seurat.filtered@meta.data), rownames(com.res2)), "labels"]
head(pbmc.seurat.filtered@meta.data)

DimPlot(pbmc.seurat.filtered, reduction = "umap", group.by = "com.res2.labels", label = TRUE)


### STRATEGY 3: using Harmonized labels

# cell.ont = "nonna": all samples without a valid term are discarded
hpca.ont <- celldex::HumanPrimaryCellAtlasData(cell.ont = "nonna")
dice.ont <- celldex::DatabaseImmuneCellExpressionData(cell.ont = "nonna")

# retain only the shared set of genes (in both references separately)
shared <- intersect(rownames(hpca.ont), rownames(dice.ont))
hpca.ont <- hpca.ont[shared,]
dice.ont <- dice.ont[shared,]

shared <- intersect(rownames(hpca), rownames(dice))
combined <- cbind(hpca[shared,], dice[shared,])
table(combined$label.main)

# top 10 most frequent cell ontology terms
tail(sort(table(hpca.ont$label.ont)),10)
tail(sort(table(dice.ont$label.ont)), 10)

# use SingleR with the cell ontology terms (label.ont) instead of cell types (label.main)
com.res3 <- SingleR(test = pbmc_counts, ref = list(HPCA = hpca.ont, DICE = dice.ont),
                    labels = list(hpca.ont$label.ont, dice.ont$label.ont))
table(com.res3$labels) # cell ontology terms


#### convert cell ontology terms into cell types and add to Seurat metadata

### 1. using mapping files
# mapping file for HumanPrimaryCellAtlasData
hpca.fle <- system.file("mapping","hpca.tsv", package = "celldex")
hpca.mapping <- read.delim(hpca.fle, header = F)
# mapping file for DatabaseImmuneCellExpressionData
dice.fle <- system.file("mapping","dice.tsv", package = "celldex")
dice.mapping <- read.delim(dice.fle, header = F)
# combine mapping files, keeping only unique values
combined_map <- unique(rbind(hpca.mapping, dice.mapping))

# create vector for mapping
map_vector <- setNames(combined_map$V1, combined_map$V2)
# use vector to map ontology to cell types and convert to data.frame
converted <- map_vector[com.res3.labels]
converted <- as.data.frame(converted)
sum(is.na(converted)) # 0

# add com.res3.labels as a column to the Seurat metadata
head(pbmc.seurat.filtered)
pbmc.seurat.filtered$com.res3.labels <- converted$converted

# even using the ontologies, there are still "duplicates"
table(converted$converted)
# Monocytes, CD14+ = CL:0002057 (7491)
# Monocyte:CD14+ =  CL:0001054 (27)
# also a group that is just CD4+ T cells, plus more specific subtypes
# also a group that is just CD8+ T cells, plus more specific subtypes


### 2. using col data
# colData for hpca.ont, keeping only the label.main and label.ont columns
col.hpca.ont <- as.data.frame(colData(hpca.ont))[, c(2,3)]
# colData for dice.ont, keeping only the label.main and label.ont columns
col.dice.ont <- as.data.frame(colData(dice.ont))[, c(2,3)]
# combine mapping files, keeping only unique values
combined_map_col <- unique(rbind(col.hpca.ont, col.dice.ont))

# create vector for mapping
map_vector_col <- setNames(combined_map_col$label.fine, combined_map_col$label.ont)
converted_col <- map_vector_col[com.res3.labels]
converted_col <- as.data.frame(converted_col)
sum(is.na(converted_col)) # 0
length(converted_col$converted_col) # 20844
table(converted_col$converted_col)

# add com.res3.ont.labels as a column to the Seurat metadata
head(pbmc.seurat.filtered)
pbmc.seurat.filtered$com.res3.ont.labels <- converted_col$converted_col


# UMAP with mapping data
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'com.res3.labels', label = TRUE)
# UMAP with colData
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'com.res3.ont.labels', label = TRUE)


### with colData using label.main (more general labels)
# colData for hpca.ont, keeping only the label.main and label.ont columns
main.hpca.ont <- as.data.frame(colData(hpca.ont))[, c(1,3)]
# colData for dice.ont, keeping only the label.main and label.ont columns
main.dice.ont <- as.data.frame(colData(dice.ont))[, c(1,3)]
# combine mapping files, keeping only unique values
combined_map_main <- unique(rbind(main.hpca.ont, main.dice.ont))

# create vector for mapping
map_vector_main <- setNames(combined_map_main$label.main, combined_map_main$label.ont)
converted_main <- map_vector_main[com.res3.labels]
converted_main <- as.data.frame(converted_main)
sum(is.na(converted_main)) # 0
length(converted_main$converted_main) # 20844
table(converted_main$converted_main)

# add com.res3.ont.labels as a column to the Seurat metadata
head(pbmc.seurat.filtered)
pbmc.seurat.filtered$com.res3.main.labels <- converted_main$converted_main


# UMAP with mapping data
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'com.res3.main.labels', label = TRUE)


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
#   [1] future_1.58.0               pheatmap_1.0.13             lubridate_1.9.4             forcats_1.0.0              
# [5] stringr_1.5.1               dplyr_1.1.4                 purrr_1.0.4                 readr_2.1.5                
# [9] tidyr_1.3.1                 tibble_3.3.0                ggplot2_3.5.2               tidyverse_2.0.0            
# [13] Seurat_5.3.0                SeuratObject_5.1.0          sp_2.2-0                    celldex_1.18.0             
# [17] SingleR_2.10.0              SummarizedExperiment_1.38.1 Biobase_2.68.0              GenomicRanges_1.60.0       
# [21] GenomeInfoDb_1.44.0         IRanges_2.42.0              S4Vectors_0.46.0            BiocGenerics_0.54.0        
# [25] generics_0.1.4              MatrixGenerics_1.20.0       matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22          splines_4.5.0             later_1.4.2               filelock_1.0.3           
# [5] polyclip_1.10-7           fastDummies_1.7.5         httr2_1.1.2               lifecycle_1.0.4          
# [9] doParallel_1.0.17         hdf5r_1.3.12              globals_0.18.0            lattice_0.22-7           
# [13] MASS_7.3-65               alabaster.base_1.8.0      magrittr_2.0.3            plotly_4.10.4            
# [17] yaml_2.3.10               httpuv_1.6.16             sctransform_0.4.2         spam_2.11-1              
# [21] spatstat.sparse_3.1-0     reticulate_1.42.0         cowplot_1.1.3             pbapply_1.7-2            
# [25] DBI_1.2.3                 RColorBrewer_1.1-3        abind_1.4-8               Rtsne_0.17               
# [29] rappdirs_0.3.3            GenomeInfoDbData_1.2.14   ggrepel_0.9.6             irlba_2.3.5.1            
# [33] listenv_0.9.1             spatstat.utils_3.1-4      goftest_1.2-3             RSpectra_0.16-2          
# [37] spatstat.random_3.4-1     fitdistrplus_1.2-2        parallelly_1.45.0         DelayedMatrixStats_1.30.0
# [41] codetools_0.2-20          DelayedArray_0.34.1       tidyselect_1.2.1          UCSC.utils_1.4.0         
# [45] farver_2.1.2              viridis_0.6.5             BiocFileCache_2.16.0      spatstat.explore_3.4-3   
# [49] jsonlite_2.0.0            BiocNeighbors_2.2.0       progressr_0.15.1          ggridges_0.5.6           
# [53] survival_3.8-3            iterators_1.0.14          foreach_1.5.2             tools_4.5.0              
# [57] ica_1.0-3                 Rcpp_1.0.14               glue_1.8.0                gridExtra_2.3            
# [61] SparseArray_1.8.0         HDF5Array_1.36.0          gypsum_1.4.0              withr_3.0.2              
# [65] BiocManager_1.30.26       fastmap_1.2.0             rhdf5filters_1.20.0       digest_0.6.37            
# [69] timechange_0.3.0          R6_2.6.1                  mime_0.13                 colorspace_2.1-1         
# [73] scattermore_1.2           tensor_1.5                dichromat_2.0-0.1         spatstat.data_3.1-6      
# [77] RSQLite_2.4.1             h5mread_1.0.1             utf8_1.2.6                data.table_1.17.4        
# [81] httr_1.4.7                htmlwidgets_1.6.4         S4Arrays_1.8.1            uwot_0.2.3               
# [85] pkgconfig_2.0.3           gtable_0.3.6              blob_1.2.4                lmtest_0.9-40            
# [89] XVector_0.48.0            htmltools_0.5.8.1         dotCall64_1.2             alabaster.matrix_1.8.0   
# [93] scales_1.4.0              png_0.1-8                 spatstat.univar_3.1-3     rstudioapi_0.17.1        
# [97] tzdb_0.5.0                reshape2_1.4.4            nlme_3.1-168              curl_6.3.0               
# [101] rhdf5_2.52.1              zoo_1.8-14                cachem_1.1.0              spacexr_2.2.1            
# [105] BiocVersion_3.21.1        KernSmooth_2.23-26        vipor_0.4.7               parallel_4.5.0           
# [109] miniUI_0.1.2              AnnotationDbi_1.70.0      ggrastr_1.0.2             pillar_1.10.2            
# [113] grid_4.5.0                alabaster.schemas_1.8.0   vctrs_0.6.5               RANN_2.6.2               
# [117] promises_1.3.3            dbplyr_2.5.0              beachmat_2.24.0           xtable_1.8-4             
# [121] cluster_2.1.8.1           beeswarm_0.4.0            cli_3.6.5                 compiler_4.5.0           
# [125] rlang_1.1.6               crayon_1.5.3              future.apply_1.20.0       labeling_0.4.3           
# [129] ggbeeswarm_0.7.2          plyr_1.8.9                stringi_1.8.7             alabaster.se_1.8.0       
# [133] viridisLite_0.4.2         deldir_2.0-4              BiocParallel_1.42.1       Biostrings_2.76.0        
# [137] lazyeval_0.2.2            spatstat.geom_3.4-1       Matrix_1.7-3              ExperimentHub_2.16.0     
# [141] RcppHNSW_0.6.0            hms_1.1.3                 patchwork_1.3.0           sparseMatrixStats_1.20.0 
# [145] bit64_4.6.0-1             Rhdf5lib_1.30.0           KEGGREST_1.48.0           shiny_1.10.0             
# [149] alabaster.ranges_1.8.0    AnnotationHub_3.16.0      ROCR_1.0-11               igraph_2.1.4             
# [153] memoise_2.0.1             bit_4.6.0 

