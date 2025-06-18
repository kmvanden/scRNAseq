# Pseudobulk DE analysis for scRNAseq data 
# https://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_pseudoBulk.R
# aggregate the counts and metadata to the sample/replicate level -> leverage existing robust bulk RNAseq DE frameworks (e.g, DESeq2)
# single cell methods treat each cell as a sample -> variation across the population is not truly investigated and p-values are inflated and there are issues with unmodelled correlations between samples (not independent of each other)
# DE testing on pseudobulk: combines the resolution offered by single-cell technologies to define the labels and statistical rigor of exisiting methods for DE analysis  

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# load libraries
library(ExperimentHub)
library(Seurat)
library(DESeq2)
library(tidyverse)

### load data
eh <- ExperimentHub()
query(eh, "Kang") # data associated with Kang
sce <- eh[["EH2259"]] # EH2259 | Kang18_8vs8 dataset
sce # SingleCellExperiment
seu.obj <- as.Seurat(sce, data = NULL) # convert into Seurat object
seu.obj
head(seu.obj@meta.data)

### QC and filtering
# get mito percent
seu.obj$mitoPercent <- PercentageFeatureSet(seu.obj, pattern = '^MT-')
head(seu.obj@meta.data)
table(seu.obj@meta.data$mitoPercent) # no mito -> already filtered

# filter low quality cells
# nFeature (number of genes) > 200 & < 2500
# nCount (number of transcripts) > 800
# get rid 
seu.filtered <- subset(seu.obj, subset = nFeature_originalexp > 200 & nFeature_originalexp < 2500 &
                         nCount_originalexp > 800 & 
                         mitoPercent < 5 &
                         multiplets == 'singlet')
seu.obj
seu.filtered


### run Seurat's standard workflow steps
seu.filtered <- seu.filtered %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()
ElbowPlot(seu.filtered) # variation captured by the PCs
seu.filtered <- RunUMAP(seu.filtered, dims = 1:20)


### visualize - cluster by cell type and condition
cell_plot <- DimPlot(seu.filtered, reduction = 'umap', group.by = 'cell', label = TRUE)
cond_plot <- DimPlot(seu.filtered, reduction = 'umap', group.by = 'stim')

cell_plot/cond_plot


############################
### pseudo-bulk workflow ###
############################

# create a new column containing both sample id (ind) and condition information (stim)
head(seu.filtered@meta.data)
seu.filtered$samples <- paste0(seu.filtered$stim, seu.filtered$ind)
head(seu.filtered@meta.data)

# aggregate counts across the cell type and sample (individual)
DefaultAssay(seu.filtered) # "originalexp"
seu.filtered[["originalexp"]]$counts # slot "counts"

cts <- AggregateExpression(seu.filtered, 
                           group.by = c("cell", "samples"),
                           assays = 'originalexp',
                           slot = "counts",
                           return.seurat = FALSE)
str(cts) # list with a dgCMatrix
cts <- as.matrix(cts$originalexp)
class(cts) # "matrix" "array" 


### split the matrix by cell type
cts.t <- t(cts) # transpose the matrix
cts.t <- as.data.frame(cts.t) # convert to data.frame
splitRows <- gsub('_.*', '', rownames(cts.t)) # get cell type from the rownames
cts.split <- split.data.frame(cts.t, f = factor(splitRows)) # split the data.frame
View(cts.split) # list of length 8 of data.frames for each cell type
cts.split$`B cells`[1:10, 1:10]


# remove B cells prefix from rownames and transpose data.frame
# gsub('.*_(.*)', '\\1', 'B cells_ctrl101')

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
})
cts.split.modified$`B cells`[1:10, 1:10]


# get a counts matrix for the B cells
counts_bcell <- cts.split.modified$`B cells`

# generate sample level metadata
colData <- data.frame(samples = colnames(counts_bcell)) # sample name
colData <- colData %>%
  mutate(condition = ifelse(grepl('stim', samples), 'Stimulated', 'Control')) %>%
  column_to_rownames(var = 'samples') # add condition (control or stimulated) column and make samples into rownames

# create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_bcell,
                              colData = colData,
                              design = ~ condition) # compare between conditions (ctrl and stim)

# filter the dds object
keep <- rowSums(counts(dds)) >= 10 # remove all the genes that have lower than 10 reads
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)

# generate results object
resultsNames(dds) # check the coefficients for the comparison
res <- results(dds, name = "condition_Stimulated_vs_Control")
res # genes that are differentially expressed between ctrl and stim in B cells


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
#   [1] muscData_1.22.0             SingleCellExperiment_1.30.1 lubridate_1.9.4            
# [4] forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4                
# [7] purrr_1.0.4                 readr_2.1.5                 tidyr_1.3.1                
# [10] tibble_3.3.0                ggplot2_3.5.2               tidyverse_2.0.0            
# [13] DESeq2_1.48.1               SummarizedExperiment_1.38.1 Biobase_2.68.0             
# [16] MatrixGenerics_1.20.0       matrixStats_1.5.0           GenomicRanges_1.60.0       
# [19] GenomeInfoDb_1.44.0         IRanges_2.42.0              S4Vectors_0.46.0           
# [22] Seurat_5.3.0                SeuratObject_5.1.0          sp_2.2-0                   
# [25] ExperimentHub_2.16.0        AnnotationHub_3.16.0        BiocFileCache_2.16.0       
# [28] dbplyr_2.5.0                BiocGenerics_0.54.0         generics_0.1.4             
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22        splines_4.5.0           later_1.4.2            
# [4] filelock_1.0.3          polyclip_1.10-7         fastDummies_1.7.5      
# [7] lifecycle_1.0.4         doParallel_1.0.17       globals_0.18.0         
# [10] lattice_0.22-7          MASS_7.3-65             magrittr_2.0.3         
# [13] plotly_4.10.4           yaml_2.3.10             httpuv_1.6.16          
# [16] sctransform_0.4.2       spam_2.11-1             spatstat.sparse_3.1-0  
# [19] reticulate_1.42.0       cowplot_1.1.3           pbapply_1.7-2          
# [22] DBI_1.2.3               RColorBrewer_1.1-3      abind_1.4-8            
# [25] Rtsne_0.17              rappdirs_0.3.3          GenomeInfoDbData_1.2.14
# [28] ggrepel_0.9.6           irlba_2.3.5.1           listenv_0.9.1          
# [31] spatstat.utils_3.1-4    goftest_1.2-3           RSpectra_0.16-2        
# [34] spatstat.random_3.4-1   fitdistrplus_1.2-2      parallelly_1.45.0      
# [37] codetools_0.2-20        DelayedArray_0.34.1     tidyselect_1.2.1       
# [40] UCSC.utils_1.4.0        farver_2.1.2            spatstat.explore_3.4-3 
# [43] jsonlite_2.0.0          progressr_0.15.1        ggridges_0.5.6         
# [46] survival_3.8-3          iterators_1.0.14        foreach_1.5.2          
# [49] tools_4.5.0             ica_1.0-3               Rcpp_1.0.14            
# [52] glue_1.8.0              gridExtra_2.3           SparseArray_1.8.0      
# [55] withr_3.0.2             BiocManager_1.30.26     fastmap_1.2.0          
# [58] digest_0.6.37           timechange_0.3.0        R6_2.6.1               
# [61] mime_0.13               colorspace_2.1-1        scattermore_1.2        
# [64] tensor_1.5              dichromat_2.0-0.1       spatstat.data_3.1-6    
# [67] RSQLite_2.4.1           data.table_1.17.4       httr_1.4.7             
# [70] htmlwidgets_1.6.4       S4Arrays_1.8.1          uwot_0.2.3             
# [73] pkgconfig_2.0.3         gtable_0.3.6            blob_1.2.4             
# [76] lmtest_0.9-40           XVector_0.48.0          htmltools_0.5.8.1      
# [79] dotCall64_1.2           scales_1.4.0            png_0.1-8              
# [82] spatstat.univar_3.1-3   rstudioapi_0.17.1       tzdb_0.5.0             
# [85] reshape2_1.4.4          nlme_3.1-168            curl_6.3.0             
# [88] cachem_1.1.0            zoo_1.8-14              spacexr_2.2.1          
# [91] BiocVersion_3.21.1      KernSmooth_2.23-26      parallel_4.5.0         
# [94] miniUI_0.1.2            AnnotationDbi_1.70.0    pillar_1.10.2          
# [97] grid_4.5.0              vctrs_0.6.5             RANN_2.6.2             
# [100] promises_1.3.3          xtable_1.8-4            cluster_2.1.8.1        
# [103] cli_3.6.5               locfit_1.5-9.12         compiler_4.5.0         
# [106] rlang_1.1.6             crayon_1.5.3            future.apply_1.20.0    
# [109] labeling_0.4.3          plyr_1.8.9              stringi_1.8.7          
# [112] viridisLite_0.4.2       deldir_2.0-4            BiocParallel_1.42.1    
# [115] Biostrings_2.76.0       lazyeval_0.2.2          spatstat.geom_3.4-1    
# [118] Matrix_1.7-3            RcppHNSW_0.6.0          hms_1.1.3              
# [121] patchwork_1.3.0         bit64_4.6.0-1           future_1.58.0          
# [124] KEGGREST_1.48.0         shiny_1.10.0            ROCR_1.0-11            
# [127] igraph_2.1.4            memoise_2.0.1           bit_4.6.0   


################################################################################


# DESeq2 workflow tutorial - differential gene expression analysis
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#standard-workflow
# https://github.com/kpatel427/YouTubeTutorials/blob/main/getData.R
# https://github.com/kpatel427/YouTubeTutorials/blob/main/runDESeq2.R

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# load libraries
library(DESeq2)
library(tidyverse)
library(airway)

# get data
data(airway)
airway

# clean up data
sample_info <- as.data.frame(colData(airway)) # convert colData into a data.frame
sample_info <- sample_info[,c(2,3)] # reduce to cell and treatment columns
sample_info$dex <- gsub('trt', 'treated', sample_info$dex) # convert trt to treated
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex) # convert untrt to untreated
names(sample_info) <- c('cellLine', 'dexamethasone') # change column names
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F) # write to a csv file

countsData <- assay(airway) # get counts table
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F) # write to a csv file


# read in counts data
counts_data <- read.csv('counts_data.csv')
head(counts_data)

# read in metadata
colData <- read.csv('sample_info.csv')
head(colData)

# make sure that the row names in colData match the column names in counts_data
all(colnames(counts_data) %in% rownames(colData))
# make sure that the row names in colData are in the same order as the column names in counts_data
all(colnames(counts_data) == rownames(colData))


# construct DESeq object
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ dexamethasone) # compare untreated versus treated
dim(dds) # 63677     8


# pre-filtering: remove rows with low gene counts (< 10 total reads)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dim(dds) # 22369     8

# set the factor level -> use untreated as the reference level
dds$dexamethasone
# [1] untreated treated   untreated treated   untreated treated   untreated treated  
# Levels: treated untreated
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

### collapse technical replicates using collapseReplicates() before performing DE analysis

# run DESeq
dds <- DESeq(dds)
res <- results(dds)
res


### explore results
summary(res) # uses an adj p-value < 0.1
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01) # uses an adj p-value < 0.01

# contrasts
resultsNames(dds) # "dexamethasone_treated_vs_untreated"
# e.g. if you have more levels: treated_4hrs, treated_8hrs, untreated
# results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated")) 


# MA plot -> scatter plot of log2fold change versus mean of normalized counts
plotMA(res)


################################################################################


# Pseudobulking and Differential Expression Analysis
# https://satijalab.org/seurat/articles/de_vignette

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# load libraries and data
library(Seurat)
library(SeuratData)
library(ggplot2)

# InstallData("ifnb")
ifnb <- LoadData("ifnb") # IFNb stimulated PBMC dataset


### perform default differential expression tests
# Seurat performs DE testing based on the non-parametric Wilcoxon rank sum test (default)

# normalize the data
ifnb <- NormalizeData(ifnb)

# find DE features between CD16+ monocytes and CD14+ monocytes
Idents(ifnb) <- "seurat_annotations"
monocyte.de.markers <- FindMarkers(ifnb, ident.1 = "CD16 Mono", ident.2 = "CD14 Mono")

head(monocyte.de.markers) # view results
# adjusted p-value using Bonferroni correction


# find DE features between CD14+ monocytes and all other cells that are more highly expressed on the CD14+ monocytes (ident.1)
monocyte.de.markers <- FindMarkers(ifnb, ident.1 = "CD16 Mono", ident.2 = NULL, only.pos = TRUE)

head(monocyte.de.markers)# view results


### perform DE analysis within the same cell type across across conditions
# control versus stimulated with IFNb
# what genes change in different conditions for cells of the same type

# create a column in the meta.data slot to hold both the cell type and treatment information
table(ifnb$stim)
table(ifnb$seurat_annotations)
ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
head(ifnb$celltype.stim)

# switch the current Idents to that column
Idents(ifnb) <- "celltype.stim"

# find the genes that are different between the control and the stimulated CD14+ monocytes
mono.de <- FindMarkers(ifnb, ident.1 = "CD14 Mono_STIM", ident.2 = "CD14 Mono_CTRL", verbose = FALSE)
head(mono.de, n = 10)

### p-values obtained from this analysis should be interpreted with caution
# these tests treat each cell as an independent replicate adn ignore inherent correlations between cells originating from the same sample
# use pseudobulking to account for this within-sample correlation


###############################################
### PERFORM DE ANALYSIS AFTER PSEUDOBULKING ###
###############################################

### retrieve the sample information for each cell
# load the inferred sample IDs of each cell
ctrl <- read.table(url("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye1.ctrl.8.10.sm.best"), head = T, stringsAsFactors = F)
stim <- read.table(url("https://raw.githubusercontent.com/yelabucsf/demuxlet_paper_code/master/fig3/ye2.stim.8.10.sm.best"), head = T, stringsAsFactors = F)
info <- rbind(ctrl, stim)

# rename the cell IDs by substituting the '-' into '.'
head(info$BARCODE)
info$BARCODE <- gsub(pattern = "\\-", replacement = "\\.", info$BARCODE)
head(info$BARCODE)

# only keep cells with high-confidence sample IDs
table(info$BEST)
info <- info[grep(pattern = "SNG", x = info$BEST), ]

# remove cells with duplicated IDs in both ctrl and stim groups
sum(duplicated(info$BARCODE)) # 250
info <- info[!duplicated(info$BARCODE) & !duplicated(info$BARCODE, fromLast = T), ]

# add the sample IDs to ifnb 
head(info)
rownames(info) <- info$BARCODE
info <- info[, c("BEST"), drop = F]
names(info) <- c("donor_id")
ifnb <- AddMetaData(ifnb, metadata = info)

# remove cells without donor IDs
sum(is.na(ifnb$donor_id)) # 331
ifnb$donor_id[is.na(ifnb$donor_id)] <- "unknown"
ifnb <- subset(ifnb, subset = donor_id != "unknown")
table(ifnb$donor_id)

# sum together the gene counts of all the cells from the same sample and condition type (by donor IDs)
# this results in one gene expression profile per sample and cell type
pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("stim", "donor_id", "seurat_annotations"))
pseudo_ifnb$celltype.stim <- paste(pseudo_ifnb$seurat_annotations, pseudo_ifnb$stim, sep = "_")
table(pseudo_ifnb$celltype.stim)

# perform DE analysis using DESeq2 at the sample level (samples, not the individual cells, are treated as independent observations)
Idents(pseudo_ifnb) <- "celltype.stim"

bulk.mono.de <- FindMarkers(object = pseudo_ifnb, 
                            ident.1 = "CD14 Mono_STIM", 
                            ident.2 = "CD14 Mono_CTRL",
                            test.use = "DESeq2")
head(bulk.mono.de, n = 15)


### compare the DE p-values between the single-cell level and the pseudobulk level results
names(bulk.mono.de) <- paste0(names(bulk.mono.de), ".bulk")
bulk.mono.de$gene <- rownames(bulk.mono.de)

names(mono.de) <- paste0(names(mono.de), ".sc")
mono.de$gene <- rownames(mono.de)

merge_dat <- merge(mono.de, bulk.mono.de, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val.bulk), ]

# number of genes that are significant
common <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 & 
                                 merge_dat$p_val.sc < 0.05)] 
length(common) # in both = 3519

only_sc <- merge_dat$gene[which(merge_dat$p_val.bulk > 0.05 & 
                                  merge_dat$p_val.sc < 0.05)]
length(only_sc) # only sc = 1649

only_bulk <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 & 
                                    merge_dat$p_val.sc > 0.05)]
length(only_bulk) # only bulk = 204

# p-values are correlated between the single-cll and pseudobulk data
# p-values are often much smaller in the single-cell analysis
# (before multiple hypothesis testing)


### investigate the discrepancies between single-cell and pseudobulk analyses using VlnPlot

# examine top genes that are differentialy expressed in both analyses
# create a new column to annotate sample-condition-celltype in the single-cell dataset
ifnb$donor_id.stim <- paste0(ifnb$stim, "-", ifnb$donor_id)

# generate violin plot 
Idents(ifnb) <- "celltype.stim"
print(merge_dat[merge_dat$gene%in%common[1:2],c('gene','p_val.sc','p_val.bulk')])
VlnPlot(ifnb, features = common[1:2], idents = c("CD14 Mono_CTRL", "CD14 Mono_STIM"), group.by = "stim") # ctrl vs stim
VlnPlot(ifnb, features = common[1:2], idents = c("CD14 Mono_CTRL", "CD14 Mono_STIM"), group.by = "donor_id.stim", ncol = 1) # donor id


# if we examine the examples of genes that are only DE under the sc analysis
print(merge_dat[merge_dat$gene%in%c('SRGN','HLA-DRA'),c('gene','p_val.sc','p_val.bulk')])
VlnPlot(ifnb, features <- c('SRGN','HLA-DRA'), idents = c("CD14 Mono_CTRL", "CD14 Mono_STIM"), group.by = "stim") # ctrl vs stim
VlnPlot(ifnb, features <- c('SRGN','HLA-DRA'), idents = c("CD14 Mono_CTRL", "CD14 Mono_STIM"), group.by = "donor_id.stim", ncol = 1) # donor id
### SGRN and HLA-DRA both have very small p-values in the sc analysis, but much larger p-values in the pseudobulk analysis
### there only appears to be a difference when ignoring sample information


# DE analysis using FindMarkers() can be performed with the following tests
# “wilcox” : Wilcoxon rank sum test (default, using ‘presto’ package)
# “wilcox_limma” : Wilcoxon rank sum test (using ‘limma’ package)
# “bimod” : Likelihood-ratio test for single cell feature expression, (McDavid et al., Bioinformatics, 2013)
# “roc” : Standard AUC classifier
# “t” : Student’s t-test
# “poisson” : Likelihood ratio test assuming an underlying negative binomial distribution. Use only for UMI-based datasets
# “negbinom” : Likelihood ratio test assuming an underlying negative binomial distribution. Use only for UMI-based datasets
# “LR” : Uses a logistic regression framework to determine differentially expressed genes. Constructs a logistic regression model predicting group membership based on each feature individually and compares this to a null model with a likelihood ratio test.
# “MAST” : GLM-framework that treates cellular detection rate as a covariate (Finak et al, Genome Biology, 2015) (Installation instructions)
# “DESeq2” : DE based on a model using the negative binomial distribution (Love et al, Genome Biology, 2014) (Installation instructions) For MAST and DESeq2, please ensure that these packages are installed separately in order to use them as part of Seurat. Once installed, use the test.use parameter can be used to specify which DE test to use.

Idents(ifnb) <- "seurat_annotations"
head(FindMarkers(ifnb, ident.1 = "CD14 Mono", ident.2 = "CD16 Mono", test.use = "MAST"))


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
#   [1] ifnb.SeuratData_3.1.0     stxBrain.SeuratData_0.1.2 ssHippo.SeuratData_3.1.4 
# [4] pbmc3k.SeuratData_3.1.4   SeuratData_0.2.2.9002     ComplexHeatmap_2.24.0    
# [7] RColorBrewer_1.1-3        Polychrome_1.5.4          genefilter_1.90.0        
# [10] presto_1.0.0              data.table_1.17.4         future_1.58.0            
# [13] viridis_0.6.5             viridisLite_0.4.2         ggpointdensity_0.2.0     
# [16] SeuratDisk_0.0.0.9021     scCustomize_3.0.1         harmony_1.2.3            
# [19] Rcpp_1.0.14               readr_2.1.5               purrr_1.0.4              
# [22] Seurat_5.3.0              SeuratObject_5.1.0        sp_2.2-0                 
# [25] ggplot2_3.5.2             dplyr_1.1.4               stringr_1.5.1            
# 
# loaded via a namespace (and not attached):
#   [1] matrixStats_1.5.0           spatstat.sparse_3.1-0       lubridate_1.9.4            
# [4] httr_1.4.7                  doParallel_1.0.17           tools_4.5.0                
# [7] sctransform_0.4.2           R6_2.6.1                    lazyeval_0.2.2             
# [10] uwot_0.2.3                  GetoptLong_1.0.5            withr_3.0.2                
# [13] prettyunits_1.2.0           gridExtra_2.3               progressr_0.15.1           
# [16] cli_3.6.5                   Biobase_2.68.0              spatstat.explore_3.4-3     
# [19] fastDummies_1.7.5           labeling_0.4.3              prismatic_1.1.2            
# [22] spatstat.data_3.1-6         ggridges_0.5.6              pbapply_1.7-2              
# [25] dichromat_2.0-0.1           parallelly_1.45.0           limma_3.64.1               
# [28] rstudioapi_0.17.1           RSQLite_2.4.1               generics_0.1.4             
# [31] shape_1.4.6.1               ica_1.0-3                   spatstat.random_3.4-1      
# [34] vroom_1.6.5                 Matrix_1.7-3                ggbeeswarm_0.7.2           
# [37] S4Vectors_0.46.0            abind_1.4-8                 lifecycle_1.0.4            
# [40] scatterplot3d_0.3-44        snakecase_0.11.1            SummarizedExperiment_1.38.1
# [43] SparseArray_1.8.0           Rtsne_0.17                  paletteer_1.6.0            
# [46] blob_1.2.4                  promises_1.3.3              crayon_1.5.3               
# [49] miniUI_0.1.2                lattice_0.22-7              cowplot_1.1.3              
# [52] annotate_1.86.0             KEGGREST_1.48.0             magick_2.8.7               
# [55] pillar_1.10.2               GenomicRanges_1.60.0        rjson_0.2.23               
# [58] future.apply_1.20.0         codetools_0.2-20            spacexr_2.2.1              
# [61] glue_1.8.0                  spatstat.univar_3.1-3       vctrs_0.6.5                
# [64] png_0.1-8                   spam_2.11-1                 gtable_0.3.6               
# [67] rematch2_2.1.2              cachem_1.1.0                S4Arrays_1.8.1             
# [70] mime_0.13                   survival_3.8-3              SingleCellExperiment_1.30.1
# [73] iterators_1.0.14            statmod_1.5.0               fitdistrplus_1.2-2         
# [76] ROCR_1.0-11                 nlme_3.1-168                bit64_4.6.0-1              
# [79] progress_1.2.3              RcppAnnoy_0.0.22            GenomeInfoDb_1.44.0        
# [82] irlba_2.3.5.1               vipor_0.4.7                 KernSmooth_2.23-26         
# [85] colorspace_2.1-1            BiocGenerics_0.54.0         DBI_1.2.3                  
# [88] ggrastr_1.0.2               DESeq2_1.48.1               tidyselect_1.2.1           
# [91] bit_4.6.0                   compiler_4.5.0              hdf5r_1.3.12               
# [94] DelayedArray_0.34.1         plotly_4.10.4               scales_1.4.0               
# [97] lmtest_0.9-40               rappdirs_0.3.3              digest_0.6.37              
# [100] goftest_1.2-3               spatstat.utils_3.1-4        XVector_0.48.0             
# [103] RhpcBLASctl_0.23-42         htmltools_0.5.8.1           pkgconfig_2.0.3            
# [106] MatrixGenerics_1.20.0       fastmap_1.2.0               rlang_1.1.6                
# [109] GlobalOptions_0.1.2         htmlwidgets_1.6.4           UCSC.utils_1.4.0           
# [112] shiny_1.10.0                farver_2.1.2                zoo_1.8-14                 
# [115] jsonlite_2.0.0              BiocParallel_1.42.1         magrittr_2.0.3             
# [118] GenomeInfoDbData_1.2.14     dotCall64_1.2               patchwork_1.3.0            
# [121] reticulate_1.42.0           stringi_1.8.7               MASS_7.3-65                
# [124] MAST_1.33.0                 plyr_1.8.9                  parallel_4.5.0             
# [127] listenv_0.9.1               ggrepel_0.9.6               forcats_1.0.0              
# [130] deldir_2.0-4                Biostrings_2.76.0           splines_4.5.0              
# [133] tensor_1.5                  hms_1.1.3                   circlize_0.4.16            
# [136] locfit_1.5-9.12             igraph_2.1.4                spatstat.geom_3.4-1        
# [139] RcppHNSW_0.6.0              reshape2_1.4.4              stats4_4.5.0               
# [142] XML_3.99-0.18               BiocManager_1.30.26         ggprism_1.0.6              
# [145] tzdb_0.5.0                  foreach_1.5.2               httpuv_1.6.16              
# [148] RANN_2.6.2                  tidyr_1.3.1                 polyclip_1.10-7            
# [151] clue_0.3-66                 scattermore_1.2             janitor_2.2.1              
# [154] xtable_1.8-4                RSpectra_0.16-2             later_1.4.2                
# [157] tibble_3.3.0                memoise_2.0.1               beeswarm_0.4.0             
# [160] AnnotationDbi_1.70.0        IRanges_2.42.0              cluster_2.1.8.1            
# [163] timechange_0.3.0            globals_0.18.0 


