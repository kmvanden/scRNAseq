# Pseudobulking and Differential Expression Analysis
# https://satijalab.org/seurat/articles/de_vignette
# https://www.melbournebioinformatics.org.au/tutorials/tutorials/seurat-de/seurat-de/#step-3-find-differentially-expressed-genes-degs-between-our-two-conditions-using-cd16-mono-cells-as-an-example

# setwd
setwd("/Users/kristinvandenham/kmvanden/RStudio/")

# load libraries
library(Seurat)
library(tidyverse)
library(pheatmap)

# load the data (use ifb_harmony.rds file generated during scRNAseq_markers_clusterid.R script
ifnb_harmony_de <- readRDS(file = "ifnb_harmony_de.rds")

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

# add the sample IDs to ifnb_harmony_de 
head(info)
rownames(info) <- info$BARCODE
info <- info[, c("BEST"), drop = F]
names(info) <- c("donor_id")
head(info)
ifnb_harmony_de <- AddMetaData(ifnb_harmony_de, metadata = info)

# remove cells without donor IDs
sum(is.na(ifnb_harmony_de$donor_id)) # 331
ifnb_harmony_de$donor_id[is.na(ifnb_harmony_de$donor_id)] <- "unknown"
ifnb_harmony_de <- subset(ifnb_harmony_de, subset = donor_id != "unknown")
table(ifnb_harmony_de$donor_id)

# sum together the gene counts of all the cells from the same sample, condition and cell type
pseudo_ifnb <- AggregateExpression(ifnb_harmony_de, assays = "RNA", return.seurat = TRUE, group.by = c("stim", "donor_id", "seurat_annotations"))
head(pseudo_ifnb@meta.data)

# add a column combining both cell type and stimulation
pseudo_ifnb$celltype.stim <- paste(pseudo_ifnb$seurat_annotations, pseudo_ifnb$stim, sep = "_")
table(pseudo_ifnb$celltype.stim)

### perform DE analysis using DESeq2 at the sample level (samples, not the individual cells, are treated as independent observations)
Idents(pseudo_ifnb) <- "celltype.stim"

bulk.ifnb.CD14.de <- FindMarkers(object = pseudo_ifnb, ident.1 = "CD14 Mono_STIM", 
                            ident.2 = "CD14 Mono_CTRL", test.use = "DESeq2")
head(bulk.ifnb.CD14.de, n = 15)

# number of genes with a adj p value less than 0.05
bulk.CD14.de.genes <- rownames(bulk.ifnb.CD14.de)[which(bulk.ifnb.CD14.de$p_val_adj < 0.05)]
length(bulk.CD14.de.genes) # 1387

head(Cells(pseudo_ifnb)) # cells are stim/donor/cell groups


### examine signficiant DEGs defined by the pseudobulk approach
CD14.sig.markers <- bulk.ifnb.CD14.de %>% 
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::mutate(gene = rownames(.))

# create a column combining cell type, donor and stimulation
ifnb_harmony_de$celltype.stim.donor_id <- paste0(ifnb_harmony_de$seurat_annotations, "-",
                                                 ifnb_harmony_de$stim, "-", ifnb_harmony_de$donor_id)
head(ifnb_harmony_de@meta.data)

# set Idents to celltype.stim.donor_id
Idents(ifnb_harmony_de) <- "celltype.stim.donor_id"

# get average pseudobulk expression
all.sig.avg.Expression.mat <- AverageExpression(ifnb_harmony_de, features = CD14.sig.markers$gene, 
                                                layer = 'scale.data')
# View(all.sig.avg.Expression.mat %>% as.data.frame())

# select only data from the CD14 monocyte cell type is being used
table(pseudo_ifnb@meta.data$seurat_annotations)
CD14.sig.avg.Expression.mat <- all.sig.avg.Expression.mat$RNA %>%
  as.data.frame() %>%
  dplyr::select(starts_with("CD14 Mono"))
# View(CD14.sig.avg.Expression.mat)

# plot heatmap of average pseudobulk scaled expression valeus from the DEGs
pheatmap(CD14.sig.avg.Expression.mat, cluster_rows = TRUE, show_rownames = FALSE, 
         border_color = NA, fontsize = 10, scale = "row", fontsize_row = 10, height = 20)

# customize heatmap with metadata
cluster_metadata <- data.frame(row.names = colnames(CD14.sig.avg.Expression.mat)) %>% 
  dplyr::mutate(Cell_Type = "CD14 Mono", 
                Treatment_Group = ifelse(str_detect(row.names(.), "STIM|CTRL"), 
                                         str_extract(row.names(.), "STIM|CTRL")))

sig.DEG.heatmap <- pheatmap(CD14.sig.avg.Expression.mat, cluster_rows = TRUE, show_rownames = FALSE, 
                            annotation = cluster_metadata[, c("Treatment_Group", "Cell_Type")], 
                            border_color = NA, fontsize = 10, scale = "row", fontsize_row = 10, height = 20,
                            annotation_names_col = FALSE)


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
#   [1] pheatmap_1.0.13    lubridate_1.9.4    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4       
# [6] purrr_1.0.4        readr_2.1.5        tidyr_1.3.1        tibble_3.3.0       ggplot2_3.5.2     
# [11] tidyverse_2.0.0    Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0          
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3          rstudioapi_0.17.1           jsonlite_2.0.0             
# [4] magrittr_2.0.3              spatstat.utils_3.1-4        farver_2.1.2               
# [7] vctrs_0.6.5                 ROCR_1.0-11                 spatstat.explore_3.4-3     
# [10] S4Arrays_1.8.1              htmltools_0.5.8.1           SparseArray_1.8.0          
# [13] sctransform_0.4.2           parallelly_1.45.0           KernSmooth_2.23-26         
# [16] htmlwidgets_1.6.4           ica_1.0-3                   plyr_1.8.9                 
# [19] plotly_4.10.4               zoo_1.8-14                  igraph_2.1.4               
# [22] mime_0.13                   lifecycle_1.0.4             iterators_1.0.14           
# [25] pkgconfig_2.0.3             Matrix_1.7-3                R6_2.6.1                   
# [28] fastmap_1.2.0               MatrixGenerics_1.20.0       GenomeInfoDbData_1.2.14    
# [31] fitdistrplus_1.2-2          future_1.58.0               shiny_1.10.0               
# [34] digest_0.6.37               colorspace_2.1-1            patchwork_1.3.0            
# [37] S4Vectors_0.46.0            DESeq2_1.48.1               tensor_1.5                 
# [40] RSpectra_0.16-2             irlba_2.3.5.1               GenomicRanges_1.60.0       
# [43] progressr_0.15.1            spatstat.sparse_3.1-0       timechange_0.3.0           
# [46] httr_1.4.7                  polyclip_1.10-7             abind_1.4-8                
# [49] compiler_4.5.0              withr_3.0.2                 doParallel_1.0.17          
# [52] BiocParallel_1.42.1         fastDummies_1.7.5           MASS_7.3-65                
# [55] DelayedArray_0.34.1         tools_4.5.0                 lmtest_0.9-40              
# [58] httpuv_1.6.16               future.apply_1.20.0         goftest_1.2-3              
# [61] glue_1.8.0                  nlme_3.1-168                promises_1.3.3             
# [64] grid_4.5.0                  Rtsne_0.17                  cluster_2.1.8.1            
# [67] reshape2_1.4.4              generics_0.1.4              gtable_0.3.6               
# [70] spatstat.data_3.1-6         tzdb_0.5.0                  data.table_1.17.4          
# [73] hms_1.1.3                   XVector_0.48.0              BiocGenerics_0.54.0        
# [76] spatstat.geom_3.4-1         RcppAnnoy_0.0.22            ggrepel_0.9.6              
# [79] RANN_2.6.2                  foreach_1.5.2               pillar_1.10.2              
# [82] spam_2.11-1                 RcppHNSW_0.6.0              later_1.4.2                
# [85] splines_4.5.0               lattice_0.22-7              survival_3.8-3             
# [88] deldir_2.0-4                tidyselect_1.2.1            locfit_1.5-9.12            
# [91] miniUI_0.1.2                pbapply_1.7-2               gridExtra_2.3              
# [94] IRanges_2.42.0              SummarizedExperiment_1.38.1 scattermore_1.2            
# [97] stats4_4.5.0                Biobase_2.68.0              spacexr_2.2.1              
# [100] matrixStats_1.5.0           UCSC.utils_1.4.0            stringi_1.8.7              
# [103] lazyeval_0.2.2              codetools_0.2-20            cli_3.6.5                  
# [106] uwot_0.2.3                  xtable_1.8-4                reticulate_1.42.0          
# [109] GenomeInfoDb_1.44.0         dichromat_2.0-0.1           Rcpp_1.0.14                
# [112] globals_0.18.0              spatstat.random_3.4-1       png_0.1-8                  
# [115] spatstat.univar_3.1-3       parallel_4.5.0              dotCall64_1.2              
# [118] listenv_0.9.1               viridisLite_0.4.2           scales_1.4.0               
# [121] ggridges_0.5.6              crayon_1.5.3                rlang_1.1.6                
# [124] cowplot_1.1.3 

