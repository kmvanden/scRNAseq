# How to create pseudobulk from scRNAseq data
# to perform differential expression analysis on scRNAseq data
# artificially add up the counts for cells from the same cell type of the same sample
# https://divingintogeneticsandgenomics.com/post/how-to-create-pseudobulk-from-single-cell-rnaseq-data/


### BASH CODE
# access data
# GEO data of a pan cancer single cell transcriptional atlas of tumor infiltrating myeloid cells
# wget -nH --cut-dirs=5  --no-clobber --convert-links --random-wait -r -p -E -e robots=off https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154763/suppl/
# ls *gz | xargs gunzip # unzip files

# look at the file structure
# normalized expression data
# less -S GSE154763_ESCA_metadata.csv
# less -S GSE154763_ESCA_normalized_expression.csv

# convert log normalized counts back to raw counts
# https://github.com/immunitastx/recover-counts
# regular expression to remove the normalized_expression.csv and rename the output to _raw_counts.txt
# ls *expression.csv | parallel --rpl '{%(.+?)} s/$$1$//;' ./recover_counts_from_log_normalized_data.py -m 10000 -d CSV {} -o {%_normalized_expression.csv}_raw_count.txt


### R

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/GSE154763/")

# load libraries
library(stringr)
library(dplyr)
library(ggplot2)
library(Seurat)
library(purrr)
# applies a function to each element of a list or atomic vector and returning an object of the same length as the input
library(readr)
library(harmony)
library(scCustomize)
library(SeuratDisk)
library(ggpointdensity)
library(viridis)
library(grid)


### read in counts files
read_counts <- function(file){
  x <- read_csv(file)
  x <- as.data.frame(x)
  cells <- x$index
  
  mat <- as.matrix(x[,-1])
  
  rownames(mat) <- cells
  mat <- t(mat)
  return(mat)
}

counts_files <- list.files(pattern = "*raw_count.txt")
counts_files

# only the last 3 files are used in the tutorial
counts_files <- counts_files %>%
  tail(n=3)

samples <- map_chr(counts_files, basename) 
samples <- str_replace(samples, "(GSE[0-9]+)_(.+)_raw_count.txt", "\\2")
names(counts_files) <- samples
counts <- purrr::map(counts_files, read_counts)
counts$PAAD[1:5, 1:5]


### read in the metadata files
read_meta <- function(file){
  y <- read_csv(file, guess_max = 100000, 
               col_types = cols(tissue = col_character()))
  y <- as.data.frame(y)
  cells <- y$index
  y <- y[,-1]
  rownames(y) <- cells
  return(y)
}

meta_files <- list.files(pattern = "metadata.csv")

# only the last 3 files are used in the tutorial
meta_files <- meta_files %>% 
  tail(n=3)

meta_names <- map_chr(meta_files, basename)
meta_names <- str_replace(meta_names, "(GSE[0-9]+)_(.+)_metadata.csv", "\\2")
names(meta_files) <- meta_names
meta <- purrr::map(meta_files, read_meta)
head(meta$PAAD)

# the different cancer types have different numbers of genes
map(counts, nrow)
# $PAAD
# [1] 14140
# 
# $THCA
# [1] 15211
# 
# $UCEC
# [1] 15849

# determine how many genes are in common among the three datasets
map(counts, ~rownames(.x)) %>%
  purrr::reduce(function(x,y) {intersect(x,y)}) %>%
  length() # 13762

# make Seurat objects
objs <- purrr::map2(counts, meta,
                    ~ CreateSeuratObject(counts = as(.x, "sparseMatrix"),
                                         meta.data = .y))

# datasets contain both normal (N) and tumor (T) tissue
table(objs$PAAD$tissue)
# N    T 
# 225 2628 

# subset tissue just to tumor cells
objs <- purrr::map(objs, ~subset(.x, subset = tissue =="T"))

objs # a list of 3 Seurat objects
# $PAAD
# An object of class Seurat 
# 14140 features across 2628 samples within 1 assay 
# Active assay: RNA (14140 features, 0 variable features)
# 1 layer present: counts
# 
# $THCA
# An object of class Seurat 
# 15211 features across 4171 samples within 1 assay 
# Active assay: RNA (15211 features, 0 variable features)
# 1 layer present: counts
# 
# $UCEC
# An object of class Seurat 
# 15849 features across 4724 samples within 1 assay 
# Active assay: RNA (15849 features, 0 variable features)
# 1 layer present: counts

preprocessSeurat <- function(obj){
  obj <-
    NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    RunHarmony(group.by.vars = "patient", dims.use = 1:30) %>%
    RunUMAP(reduction = "harmony", dims = 1:30) %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = 0.6)

  return(obj)
}

objs <- purrr::map(objs, preprocessSeurat)

# plot UMAP of UCEC
p1 <- scCustomize::DimPlot_scCustom(objs$UCEC) # colored by cluster id
p2 <- scCustomize::DimPlot_scCustom(objs$UCEC, group.by = "MajorCluster") # colored by annotation from the original author
p1/p2


##################################################
### CREATE PSEUDOBULK OF THE SINGLE CELL DATA ###
##################################################

# ensure that annotation is consistent between the different cancer types
purrr::map(objs, ~.x$MajorCluster %>% unique() %>% sort())

# remove the M## at the start of the names
clean_cluster_name <- function(obj){
  annotation <- obj@meta.data %>%
    mutate(annotation = str_replace(MajorCluster, "^M[0-9]{2}_", "")) %>%
    pull(annotation)
  obj$annotation <- annotation
  return(obj)
}

objs <- purrr::map(objs, clean_cluster_name)
purrr::map(objs, ~.x$annotation %>% unique() %>% sort())


### CREATE PSEUDOBULK
# 1. using presto
library(presto)

# devtools::install_github("immunogenomics/presto")
# /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rcpp/include/Rcpp/platform/compiler.h:100:10: fatal error: 'cmath' file not found
# 100 | #include <cmath>
#   |          ^~~~~~~
# it is caused by Apple: they prompt you to update Xcode, and then you have a mismatch between the Xcode version and the OS version. 
# sudo rm -rf /Library/Developer/CommandLineTools # completely uninstall the xcode command line tools
# xcode-select --install # reinstall (installation dialog should pop up)

# create a pseudobulk per cancer type (cancer), per sample (patient) and per cell type (annotation)

head(objs$PAAD@meta.data)
colnames(objs$PAAD@meta.data)
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent_mito"    "n_counts"       
# [6] "percent_hsp"     "barcode"         "batch"           "library_id"      "cancer"         
# [11] "patient"         "tissue"          "n_genes"         "MajorCluster"    "source"         
# [16] "tech"            "UMAP1"           "UMAP2"           "RNA_snn_res.0.6" "seurat_clusters"
# [21] "annotation"

objs$PAAD@assays$RNA
# Assay (v5) data with 14140 features for 2628 cells
# Top 10 variable features:
#   CCL17, SERPINB2, CXCL10, HSPA6, FSCN1, CXCL1, CCL20, CCL2, CCL22, CCL18 
# Layers:
#   counts, data, scale.data 

create_pseudo_bulk <- function(obj){
  data_collapsed <- presto::collapse_counts(obj@assays$RNA$counts, 
                                            obj@meta.data, 
                                            c("annotation", "patient", "cancer"))
  meta_data <- data_collapsed$meta_data 
  
  mat <- data_collapsed$counts_mat
  
  colnames(mat) <- paste(meta_data$cancer, colnames(mat), sep = "_")
  return(list(mat = mat, meta_data = meta_data))
}

pseudo_bulk_objs <- purrr::map(objs, create_pseudo_bulk)
str(pseudo_bulk_objs) # look at created object


purrr::map(pseudo_bulk_objs, ~nrow(.x$mat)) ### different number of genes per dataset 
# intersect only the common genes so that datasets can be combined
genes_per_data <- purrr::map(pseudo_bulk_objs, ~rownames(.x$mat)) 
common_genes <- purrr::reduce(genes_per_data, intersect) 
length(common_genes) # 13762 | common number of genes

  
# subset the matrix to just the common genes
subset_common_mat <- function(mat){
  mat <- mat[common_genes, ]
  return(mat)
}
mats <- purrr::map(pseudo_bulk_objs, ~subset_common_mat(.x$mat))
purrr::map(mats, dim)
mats$PAAD[1:5, 1:6]
str(mats) # list of 3

# combine the counts data
mats <- purrr::reduce(mats, cbind)
dim(mats) # 13762   255
str(mats)

# combine the metadata
meta_data <- purrr::map(pseudo_bulk_objs, ~.x$meta_data)
str(meta_data) # list of 3 data frames
meta_data <- purrr::reduce(meta_data, rbind) 
final_meta <- meta_data
str(final_meta) # data frame of 255 obs. of 3 variables
dim(final_meta) # 255   3


####################################
##### PCA ANALYSIS AND HEATMAP #####
####################################

### log normalize the counts matrix
total_reads <- colSums(mats)
final_mat <- t(t(mats)/total_reads)
final_mat <- log2(final_mat + 1)


### choose the top 1000 most variable genes
library(genefilter)
top_genes <- genefilter::rowVars(final_mat) %>% 
  sort(decreasing = TRUE) %>%
  names() %>%
  head(1000)

# subset to only the top 1000 most variable genes
expression_mat_sub <- final_mat[top_genes, ]


### calculate the PCA
pca <- prcomp(t(expression_mat_sub), center = TRUE, scale. = TRUE) 

# check the order of the samples are the same
all.equal(rownames(pca$x), colnames(expression_mat_sub))

PC1_and_PC2 <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2]) # extract the first two PCs
PC1_and_PC2 <- cbind(PC1_and_PC2, final_meta) # bind the metadata together

head(PC1_and_PC2)


### plot the PCA plot
library(ggplot2)
library(Polychrome)
set.seed(123)

# color the plot by cancer type
p1 <- ggplot(PC1_and_PC2, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = cancer)) +
  theme_bw(base_size = 14) 

# color the plot by cell type
length(unique(PC1_and_PC2$annotation))
mypal <- kelly.colors(14)[-1] # to get rid of the first color (white)

p2 <- ggplot(PC1_and_PC2, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = annotation)) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = mypal %>% unname()) 

p1/p2 # the same cell types from different cancers cluster together


# visualize in a heatmap
library(RColorBrewer)
library(ComplexHeatmap)

set.seed(123)
cols <- RColorBrewer::brewer.pal(n = 3, name = "Set1") # 3 colors
unique(final_meta$cancer) # "PAAD" "THCA" "UCEC"


length(rownames(final_meta)) # 255
length(colnames(expression_mat_sub)) # 255
rownames(final_meta) <- colnames(expression_mat_sub)

ha <- HeatmapAnnotation(df = final_meta[, -2],
                       col = list(cancer = setNames(cols, unique(final_meta$cancer) %>% sort()),
                                annotation = setNames(unname(mypal), unique(final_meta$annotation) %>% sort())),
                       show_annotation_name = TRUE)

col_fun <- circlize::colorRamp2(c(-2, 0, 2), colors = c("blue", "white", "red"))

Heatmap(t(scale(t(expression_mat_sub))), top_annotation = ha,
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = FALSE,
        col = col_fun,
        name = "expression",
        column_dend_reorder = TRUE,
        raster_quality = 3,
        use_raster = FALSE,
)


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
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ComplexHeatmap_2.24.0 RColorBrewer_1.1-3    Polychrome_1.5.4      genefilter_1.90.0    
# [5] presto_1.0.0          data.table_1.17.4     future_1.58.0         viridis_0.6.5        
# [9] viridisLite_0.4.2     ggpointdensity_0.2.0  SeuratDisk_0.0.0.9021 scCustomize_3.0.1    
# [13] harmony_1.2.3         Rcpp_1.0.14           readr_2.1.5           purrr_1.0.4          
# [17] Seurat_5.3.0          SeuratObject_5.1.0    sp_2.2-0              ggplot2_3.5.2        
# [21] dplyr_1.1.4           stringr_1.5.1        
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22            splines_4.5.0               later_1.4.2                
# [4] prismatic_1.1.2             tibble_3.2.1                polyclip_1.10-7            
# [7] janitor_2.2.1               XML_3.99-0.18               fastDummies_1.7.5          
# [10] lifecycle_1.0.4             doParallel_1.0.17           globals_0.18.0             
# [13] lattice_0.22-7              vroom_1.6.5                 hdf5r_1.3.12               
# [16] MASS_7.3-65                 magrittr_2.0.3              plotly_4.10.4              
# [19] httpuv_1.6.16               sctransform_0.4.2           spam_2.11-1                
# [22] spatstat.sparse_3.1-0       reticulate_1.42.0           cowplot_1.1.3              
# [25] pbapply_1.7-2               DBI_1.2.3                   lubridate_1.9.4            
# [28] abind_1.4-8                 GenomicRanges_1.60.0        Rtsne_0.17                 
# [31] BiocGenerics_0.54.0         GenomeInfoDbData_1.2.14     circlize_0.4.16            
# [34] IRanges_2.42.0              S4Vectors_0.46.0            ggrepel_0.9.6              
# [37] irlba_2.3.5.1               listenv_0.9.1               spatstat.utils_3.1-4       
# [40] goftest_1.2-3               RSpectra_0.16-2             spatstat.random_3.4-1      
# [43] annotate_1.86.0             fitdistrplus_1.2-2          parallelly_1.45.0          
# [46] DelayedArray_0.34.1         codetools_0.2-20            tidyselect_1.2.1           
# [49] shape_1.4.6.1               UCSC.utils_1.4.0            farver_2.1.2               
# [52] matrixStats_1.5.0           stats4_4.5.0                spatstat.explore_3.4-3     
# [55] jsonlite_2.0.0              GetoptLong_1.0.5            progressr_0.15.1           
# [58] ggridges_0.5.6              survival_3.8-3              iterators_1.0.14           
# [61] foreach_1.5.2               tools_4.5.0                 ica_1.0-3                  
# [64] glue_1.8.0                  SparseArray_1.8.0           gridExtra_2.3              
# [67] DESeq2_1.48.1               MatrixGenerics_1.20.0       GenomeInfoDb_1.44.0        
# [70] withr_3.0.2                 fastmap_1.2.0               digest_0.6.37              
# [73] timechange_0.3.0            R6_2.6.1                    mime_0.13                  
# [76] ggprism_1.0.6               colorspace_2.1-1            scattermore_1.2            
# [79] tensor_1.5                  dichromat_2.0-0.1           spatstat.data_3.1-6        
# [82] RSQLite_2.4.0               RhpcBLASctl_0.23-42         tidyr_1.3.1                
# [85] generics_0.1.4              S4Arrays_1.8.1              httr_1.4.7                 
# [88] htmlwidgets_1.6.4           scatterplot3d_0.3-44        uwot_0.2.3                 
# [91] pkgconfig_2.0.3             gtable_0.3.6                blob_1.2.4                 
# [94] lmtest_0.9-40               XVector_0.48.0              htmltools_0.5.8.1          
# [97] dotCall64_1.2               clue_0.3-66                 scales_1.4.0               
# [100] Biobase_2.68.0              png_0.1-8                   spatstat.univar_3.1-3      
# [103] snakecase_0.11.1            rstudioapi_0.17.1           rjson_0.2.23               
# [106] tzdb_0.5.0                  reshape2_1.4.4              nlme_3.1-168               
# [109] zoo_1.8-14                  cachem_1.1.0                GlobalOptions_0.1.2        
# [112] spacexr_2.2.1               KernSmooth_2.23-26          parallel_4.5.0             
# [115] miniUI_0.1.2                vipor_0.4.7                 AnnotationDbi_1.70.0       
# [118] ggrastr_1.0.2               pillar_1.10.2               vctrs_0.6.5                
# [121] RANN_2.6.2                  promises_1.3.3              xtable_1.8-4               
# [124] cluster_2.1.8.1             beeswarm_0.4.0              paletteer_1.6.0            
# [127] magick_2.8.6                locfit_1.5-9.12             cli_3.6.5                  
# [130] compiler_4.5.0              rlang_1.1.6                 crayon_1.5.3               
# [133] future.apply_1.20.0         labeling_0.4.3              rematch2_2.1.2             
# [136] plyr_1.8.9                  forcats_1.0.0               ggbeeswarm_0.7.2           
# [139] stringi_1.8.7               BiocParallel_1.42.1         deldir_2.0-4               
# [142] Biostrings_2.76.0           lazyeval_0.2.2              spatstat.geom_3.4-1        
# [145] Matrix_1.7-3                RcppHNSW_0.6.0              hms_1.1.3                  
# [148] patchwork_1.3.0             bit64_4.6.0-1               KEGGREST_1.48.0            
# [151] shiny_1.10.0                SummarizedExperiment_1.38.1 ROCR_1.0-11                
# [154] igraph_2.1.4                memoise_2.0.1               bit_4.6.0   

