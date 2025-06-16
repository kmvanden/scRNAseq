### How PCA projection and cell label transfer work in Seurat
# https://crazyhottommy.github.io/single-cell-RNAseq-PCA-CCA-cell-annotation/how-seurat-pca-label-transfer.html

# understand the example datasets
# we will use the PBMC3k and PBMC10k data
# we will project the PBMC3k data to the PBMC10k data and get the labels

# load packages
library(Seurat)
library(SeuratData)
library(ggplot2)
library(Matrix)
library(irlba)  # For PCA
library(RcppAnnoy)  # For fast nearest neighbor search
library(dplyr)
library(ComplexHeatmap)

# assuming the PBMC datasets (3k and 10k) are already normalized and represented as sparse matrices
# AvailableData()
# InstallData("pbmc3k")
 
pbmc3k = UpdateSeuratObject(pbmc3k)
pbmc3k@meta.data %>% head()

# routine processing
pbmc3k = pbmc3k %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:10, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE) %>%
  RunUMAP(dims = 1:10, verbose = FALSE)

# examine the pbmc3k data
p1 = DimPlot(pbmc3k, reduction = "umap", label = TRUE, group.by = 
               "RNA_snn_res.0.5")
p2 = DimPlot(pbmc3k, reduction = "umap", label = TRUE, group.by = "seurat_annotations", label.size = 3)
p1 + p2

# read in the pbmc10k data and examine data
pbmc10k = readRDS("/Users/kristinvandenham/kmvanden/RStudio/pbmc_10k_v3.rds")
pbmc10k = UpdateSeuratObject(pbmc10k)
pbmc10k@meta.data %>% head()

DimPlot(pbmc10k, label = TRUE, repel = TRUE) + NoLegend()

# the pbmc3k and pmbc10k have a different number of gene names
# subset to the common genes
# use the 10k pbmc data to transfer the labels

pbmc3k_genes = rownames(pbmc3k)
pbmc10k_genes = rownames(pbmc10k)

head(pbmc3k_genes)
head(pbmc10k_genes)

# find common genes
common_genes = intersect(pbmc3k_genes, pbmc10k_genes)
head(common_genes)
# length(pbmc3k_genes)

pbmc3k = subset(pbmc3k, features = common_genes)
pbmc10k = subset(pbmc10k, features = common_genes)

all.equal(rownames(pbmc3k), rownames(pbmc10k)) # check numbers


# singular value decomposition (SVD) | X = UDV^T
# the decomposition of matrix X (dimensions n x p) | n = # of cells/samples and p = # of features
# U = n x n orthogonal matrix containing the left singular vectors (cells/samples)
# D = n x p diagonal matrix (non-negative real numbers on the diagonal -> singular values of X, which indicate the variance captured by each component)

# principal components (PCs) | Z = UD
# Z contains the projection of you data onto the principal component space

# calculate PCA from scratch using irlba (for big matrices)
# can use built-in svd if the matrices are small


# scale the matrix 
dim(pbmc10k)
pbmc10k_scaled = pbmc10k@assays$RNA@scale.data
dim(pbmc10k_scaled) # 2068 9432

# transpose the matrix to gene x sample
# keep 100 PCs
pca_10k = irlba(t(pbmc10k_scaled), nv = 100) 


#### source code for RunPCA in Seurat ####
# by default, RunPCA computes the PCA on the cell x gene (n x p) matrix
# after running irlba, the v matrix is the gene loadings and the u matrix is the cell embeddings

# pcs.compute = min(pcs.compute, nrow(x = data.use)-1)
# pca.results = irlba(A = t(x = data.use), nv = pcs.compute, ...)
# 
# gene.loadings = pca.results$v
# 
# sdev = pca.results$d/sqrt(max(1, ncol(data.use) - 1))
# 
# if(weight.by.var){
#   cell.embeddings = pca.results$u %*% diag(pca.results$d)
# } else {
#   cell.embeddings = pca.results$u
# }
# 
# rownames(x = gene.loadings) = rownames(x = data.use)
# colnames(x = gene.loadings) = paste0(reduction.key, 1:pcs.compute)
# rownames(x = cell.embeddings) = colnames(x = data.use)
# colnames(x = cell.embeddings) = colnames(x = gene.loadings)

#### source code for RunPCA in Seurat ####


# get the gene loadings (V matrix)
# gene loadings (feature/genes in rows and PCS in columns)
gene_loadings_10k = pca_10k$v
dim(gene_loadings_10k) # 2068 100 | 2068 features/most variable genes and 100 PCs

rownames(gene_loadings_10k) = rownames(pbmc10k_scaled)
colnames(gene_loadings_10k) = paste0("PC", 1:100)

# colnames(gene_loadings_10k) %>% head()

length(rownames(gene_loadings_10k)) # 2068
length(VariableFeatures(pbmc10k)) # 2068

# get PCA/cell embeddings (U matrix x D matrix)
cell_embeddings_10k = pca_10k$u %*% diag(pca_10k$d) # 10k cells in rows
dim(cell_embeddings_10k)

rownames(cell_embeddings_10k) = colnames(pbmc10k_scaled)
colnames(cell_embeddings_10k) = colnames(gene_loadings_10k)

length(rownames(cell_embeddings_10k)) # 9432
length(colnames(cell_embeddings_10k)) # 100
cell_embeddings_10k[1:5, 1:10]


# center the 3k pbmc data with the 10k gene means and scale
# aligns data distributions: ensures both datasets are aligned, reducing biases from differences in gene expression profiles
# ensures consistency: if the datasets come from difference conditions, centering standardizes them, making them more comparable
# variance representation: ensures the variance is accurately captured by the principle components (PCA is sensitive to the data's mean)
# improves projection accuracy: proper centering improves projection accuracy and enhances the label transfer process by focusing on biological variation instead of technical noise

pbmc3k_normalized = pbmc3k@assays$RNA$data
pbmc3k_scaled = scale(t(pbmc3k_normalized), 
                       center = rowMeans(pbmc10k@assays$RNA$data), 
                       scale = TRUE)
dim(pbmc3k_scaled) # 2700 11774

# subset the same genes for the scaled data
pbmc3k_scaled = pbmc3k_scaled[, rownames(pbmc10k_scaled)]
dim(pbmc3k_scaled) # 2700 2068

# project the 3k cells onto the PCA space of 10k dataset
# the expression matrix is the pbmc3k_scaled matrix and the V matrix (gene loadings) is from the 10k pbmc data
cell_embeddings_3k = as.matrix(pbmc3k_scaled) %*% gene_loadings_10k
cell_embeddings_3k[1:5, 1:5]
all.equal(rownames(cell_embeddings_3k), rownames(pbmc3k@meta.data)) # check numbers

cbind(cell_embeddings_3k, pbmc3k@meta.data) %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = seurat_annotations)) +
  theme_classic(base_size = 14)

# PCA space based on pbmc3k on its own 
DimPlot(pbmc3k, reduction = "pca", group.by = "seurat_annotations", 
        label = TRUE) +
  NoLegend()

# the cells are split roughly into three major islands
# B cells, myeloid cells (CD14+CD16+ monocytes), and T cells and NK cells

# the pbmc3k cell embeddings have been projected onto the pmbc10k PCA space
# find the k nearest neighbors in the 10k dataset for every cell in the 3k dataset

# AnnoyAngular calculates the cosine distance in the PCA space to find the NN between the 3k and 10k datasets
# create Annoy index for the pmbc10k dataset | build index with 10 trees
annoy_index = new(AnnoyAngular, ncol(cell_embeddings_10k)) 
for (i in 1:nrow(cell_embeddings_10k)) {
  annoy_index$addItem(i - 1, cell_embeddings_10k[i, ])
}
annoy_index$build(10)

# find the nearest neighbors for each cell in the 3k dataset
n_neighbors = 30 # number of nearest neighbors to find
nn_indices = t(sapply(1:nrow(cell_embeddings_3k), function(i) {
  annoy_index$getNNsByVector(cell_embeddings_3k[i, ], n_neighbors)
})) # gives the indices of NN in the pmbc10k dataset

# the rows are cells from the pmbc3k dataset and columns are the 30 nearest cells in the 10k dataset
dim(nn_indices)
head(nn_indices)

# transfer labels based on majority vote from nearest neighbors
labels_10k<- as.character(pbmc10k$celltype)

transfer_labels = apply(nn_indices, 1, function(neighbors) {# get labels for the nearest neighbors
  neighbor_labels = labels_10k[neighbors + 1]  # add 1 for R's 1-based index
    most_common_label = names(sort(table(neighbor_labels), decreasing = TRUE))[1]
  return(most_common_label)
})

head(transfer_labels) # transfer_labels contains the predicted labels for the pmbc3k dataset
pbmc3k$predicted = transfer_labels
DimPlot(pbmc3k, reduction = "umap", group.by = "predicted", label = TRUE, repel = TRUE) +
  NoLegend()



#### compare with Seurat's wrapper ####
# step 1: find transfer anchors
anchors = FindTransferAnchors(
  reference = pbmc10k,     # reference dataset
  query = pbmc3k,          # query dataset
  dims = 1:100,            # dimensions to use for anchor finding
  reduction = "pcaproject" # default
) # 2539 anchors

# step 2: transfer labels
predictions = TransferData(
  anchors = anchors,           # anchors identified in previous step
  refdata = pbmc10k$celltype,  # true labels in pmbc10k
  dims = 1:30                  # dimensions to use for transferring
)

# step 3: add predictions to the query dataset
pbmc3k = AddMetaData(pbmc3k, metadata = predictions)

#### compare with Seurat's wrapper ####



# predicted.id = Seurat's wrapper function | predicted = naive implementation
table(pbmc3k$predicted, pbmc3k$predicted.id)

# visualize in a heatmap
table(pbmc3k$predicted, pbmc3k$predicted.id) %>%
  as.matrix() %>%
  scale() %>%
  Heatmap(cluster_rows = FALSE, cluster_columns = FALSE, name = "scaled\ncell number")

# mutual nearest neighbors (MNN) is a key part of anchor identification during label transfer in Seurat
# used to match cells from two datasets (query and reference)

# based on their proximity in the shared feature space (e.g., PCA space)
# find nearest neighbors in dataset b for a cell in dataset a, and the nearest neighbors in dataset a for a cell in dataset b
# if the two cells are nearest neighbors of each other they are MNNs
# MNNs serve as anchors between the two datasets that help align the datasets for further downstream tasks (i.e., label transfer)

# MNN: the nearest neighbor relationship is mutual (designed to be more robust | important if the datasets have batch or other technical differences)

# kNN with PCA projection: the query dataset is projected into the reference dataset's PCA space, and then find the nearest neighbors in that space


# find the nearest neighbors for each dataset seperately
n_neighbors = 30 # number of nearest neighbors to find
### RcppAnnoy in C is 0 based, R is 1 based --> add 1 to the index

# build an annoy index for the 10k dataset using cosine distance (not Euclidean = AnnoyEuclidean)
annoy_index_10k = new(AnnoyAngular, ncol(cell_embeddings_10k))
for (i in 1:nrow(cell_embeddings_10k)) {
  annoy_index_10k$addItem(i - 1, cell_embeddings_10k[i, ])  # 0-based index for Annoy
} # add the PCA embeddings for each cell to the index
annoy_index_10k$build(10) # build index for fast nearest neighbor search


# build an annoy index for the 3k dataset using cosine distance (not Euclidean = AnnoyEuclidean)
annoy_index_3k = new(AnnoyAngular, ncol(cell_embeddings_3k)) 
for (i in 1:nrow(cell_embeddings_3k)) {
  annoy_index_3k$addItem(i - 1, cell_embeddings_3k[i, ])  # 0-based index for Annoy
} # add the PCA embeddings for each cell to the index
annoy_index_3k$build(10) # build index for fast nearest neighbor search


# find the nearest neighbors in the pmbc10k for each cell in the pmbc3k
nn_10k_for_3k = t(sapply(1:nrow(cell_embeddings_3k), function(i) {
  annoy_index_10k$getNNsByVector(cell_embeddings_3k[i, ], n_neighbors)
}))
nn_10k_for_3k = nn_10k_for_3k + 1  # convert to 1-based indexing for R
head(nn_10k_for_3k)


# find the nearest neighbors in the pmbc3k for each cell in the pmbc10k
nn_3k_for_10k = t(sapply(1:nrow(cell_embeddings_10k), function(i) {
  annoy_index_3k$getNNsByVector(cell_embeddings_10k[i, ], n_neighbors)
}))
nn_3k_for_10k = nn_3k_for_10k + 1  # convert to 1-based indexing for R
head(nn_3k_for_10k)


# identify mutual nearest neighbors (MNN)
labels_10k = as.character(labels_10k)

# create empty vectors to store the scores and labels
pbmc3k_transferred_labels = rep(NA, nrow(cell_embeddings_3k))
pbmc3k_transfer_scores = rep(0, nrow(cell_embeddings_3k))

# loop through each cell in the pbmc3k dataset to find the mutual nearest neighbors
for (i in 1:nrow(cell_embeddings_3k)) {
  # Get nearest neighbors of the i-th 3k cell in 10k
  nn_in_10k <- nn_10k_for_3k[i, ]
  
  # initialize count for mutual nearest neighbors
  mutual_count = 0
  
  # check mutual nearest neighbors
  for (nn in nn_in_10k) {
    # check if i-th pbmc3k cell is a nearest neighbor for the nn-th pbmc10k cell
    if (i %in% nn_3k_for_10k[nn, ]) {  # correct for 1-based indexing in R
      mutual_count <- mutual_count + 1
      
      # transfer the label from the pmbc10k cell to the pmbc3k cell
      pbmc3k_transferred_labels[i] <- labels_10k[nn]
    }
  }
  
  # calculate the transfer score (mutual neighbor count / total neighbors)
  pbmc3k_transfer_scores[i] <- mutual_count / n_neighbors
}

# weighting transfer system: how Seurat handles cells that are not MNNs (don't have a direct anchor)
# labels are predicted based on the similarity (distance) to the identified anchors
# the influence of each anchors is weighted according to its distance from the query cells
# the labels for cells that don't have MNN is inferred based on their relative position to the MNN cells

# Seurat's TransferData function takes into account the proximity of these non-MNN cells to the anchor cells 
# and extrapolates the labels based on the information from the anchor set
# Seurat also provides prediction scores that indicate the confidence of the label transfer for each cell

# fill in missing labels for cells without MNN based on nearest neighbor in pmbc10k
for (i in 1:length(pbmc3k_transferred_labels)) {
  if (is.na(pbmc3k_transferred_labels[i])) {
    # assign the label of the nearest pmbc10k cell
    nearest_10k_cell = nn_10k_for_3k[i, 1]  # First nearest neighbor
    pbmc3k_transferred_labels[i] = labels_10k[nearest_10k_cell]
    
    # assign a lower score for non-mutual neighbors
    pbmc3k_transfer_scores[i] = 0.01  # assign a small, non-zero, score for non-mutual neighbors
  }
}

head(pbmc3k_transferred_labels)
head(pbmc3k_transfer_scores)

pbmc3k$pbmc3k_transferred_labels = pbmc3k_transferred_labels # add predictions to the query dataset
# predicted.id = Seurat's wrapper function | pbmc3k_transferred_labels = naive MNN implementation
table(pbmc3k$pbmc3k_transferred_labels, pbmc3k$predicted.id)

# visualize in heatmap
table(pbmc3k$pbmc3k_transferred_labels, pbmc3k$predicted.id) %>%
  as.matrix() %>%
  scale() %>%
  Heatmap(cluster_rows = FALSE, cluster_columns= FALSE, name= "scaled\ncell number")

# Seurat's MNN implementation includes additional optimization (e.g., PC scaling and anchor filtering)


sessionInfo()
# R version 4.5.0 (2025-04-11)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.4.1
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
#   [1] ComplexHeatmap_2.24.0   ggplot2_3.5.2           future_1.40.0           dplyr_1.1.4            
# [5] RcppAnnoy_0.0.22        irlba_2.3.5.1           Matrix_1.7-3            pbmc3k.SeuratData_3.1.4
# [9] SeuratData_0.2.2.9002   Seurat_5.3.0            SeuratObject_5.1.0      sp_2.2-0               
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     shape_1.4.6.1          rstudioapi_0.17.1      jsonlite_2.0.0        
# [5] magrittr_2.0.3         spatstat.utils_3.1-3   farver_2.1.2           GlobalOptions_0.1.2   
# [9] vctrs_0.6.5            ROCR_1.0-11            spatstat.explore_3.4-2 htmltools_0.5.8.1     
# [13] sctransform_0.4.1      parallelly_1.43.0      KernSmooth_2.23-26     htmlwidgets_1.6.4     
# [17] ica_1.0-3              plyr_1.8.9             plotly_4.10.4          zoo_1.8-14            
# [21] igraph_2.1.4           mime_0.13              lifecycle_1.0.4        iterators_1.0.14      
# [25] pkgconfig_2.0.3        R6_2.6.1               fastmap_1.2.0          clue_0.3-66           
# [29] fitdistrplus_1.2-2     shiny_1.10.0           digest_0.6.37          colorspace_2.1-1      
# [33] patchwork_1.3.0        S4Vectors_0.46.0       tensor_1.5             RSpectra_0.16-2       
# [37] labeling_0.4.3         progressr_0.15.1       spatstat.sparse_3.1-0  httr_1.4.7            
# [41] polyclip_1.10-7        abind_1.4-8            compiler_4.5.0         withr_3.0.2           
# [45] doParallel_1.0.17      fastDummies_1.7.5      MASS_7.3-65            rappdirs_0.3.3        
# [49] rjson_0.2.23           tools_4.5.0            lmtest_0.9-40          httpuv_1.6.16         
# [53] future.apply_1.11.3    goftest_1.2-3          glue_1.8.0             nlme_3.1-168          
# [57] promises_1.3.2         Rtsne_0.17             cluster_2.1.8.1        reshape2_1.4.4        
# [61] generics_0.1.3         gtable_0.3.6           spatstat.data_3.1-6    tidyr_1.3.1           
# [65] data.table_1.17.0      BiocGenerics_0.54.0    spatstat.geom_3.3-6    ggrepel_0.9.6         
# [69] RANN_2.6.2             foreach_1.5.2          pillar_1.10.2          stringr_1.5.1         
# [73] spam_2.11-1            RcppHNSW_0.6.0         later_1.4.2            circlize_0.4.16       
# [77] splines_4.5.0          lattice_0.22-7         survival_3.8-3         deldir_2.0-4          
# [81] tidyselect_1.2.1       miniUI_0.1.2           pbapply_1.7-2          gridExtra_2.3         
# [85] IRanges_2.42.0         scattermore_1.2        stats4_4.5.0           matrixStats_1.5.0     
# [89] stringi_1.8.7          lazyeval_0.2.2         codetools_0.2-20       tibble_3.2.1          
# [93] BiocManager_1.30.25    cli_3.6.5              uwot_0.2.3             xtable_1.8-4          
# [97] reticulate_1.42.0      munsell_0.5.1          Rcpp_1.0.14            globals_0.17.0        
# [101] spatstat.random_3.3-3  png_0.1-8              spatstat.univar_3.1-2  parallel_4.5.0        
# [105] dotCall64_1.2          listenv_0.9.1          viridisLite_0.4.2      scales_1.4.0          
# [109] ggridges_0.5.6         purrr_1.0.4            crayon_1.5.3           GetoptLong_1.0.5      
# [113] rlang_1.1.6            cowplot_1.1.3  









