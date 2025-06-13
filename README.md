# 🧬 Single-cell RNA Sequencing Workflow

### Create a Seurat object
#### A Seurat object serves as a container for single-cell dataset data (e.g., count matrix) and analyses (PCA or clustering results).
1. From the output of the CellRanger pipeline from 10X Genomics (raw data)
    - Create a sparse count matrix by reading in the data ```Read10X_h5()```
    - Convert the sparse count matrix in to a Seurat object ```CreateSeuratObject(min.cells = 3, min.features = 200)```
5. From datasets stored in SeuratData
    - Load in data using ```LoadData())```
<!-- Merge the multiple Seurat objects into one object, if needed, for easier downstream manipulation. -->

### Quality control and filtering
1. Mitochondrial genes/features
    - A high percentage of mitochondrial genes = loq quality or dying cell
    - Create new column using ```PercentageFeatureSet(pattern = "^MT-")```
2. Number of unique genes/features detected in each cell (nFeature_RNA) and total number counts/molecules detected in each cell (nCount_RNA).
    - A low number genes or counts = low quality cells or empty droplets
    - A high number of genes or counts = doublets or multiplets
- Create plots to explore data and determine thresholds 
    - ```VlnPlot(nsclc_seu, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"), ncol = 3)``` 
    - ```FeatureScatter(nsclc_seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')```
- Filter data according to the determined thresholds
    - ```subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)```
> [!NOTE]
> Quality control metrics are dependent on cell type.

### Normalize data
#### Eliminates batch effects/technical variation (e.g., sequencing depth).
- Normalizes the feature expression measurements for each cell by the total expression, then multiples this by a scale factor (10,000) and log-transforms the result.
    - ```NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000)```

### Find variable features
#### Cells with low cell-to-cell variation (i.e., housekeeping genes) are not very informative, whereas focusing on variable genes helps to highlight the biological signal in single-cell datasets.
1. Find variable features ```FindVariableFeatures(nfeatures = 2000)``` 
2. Examine the variable features ```VariableFeaturePlot(): look at most variable genes```
* Variable features are used with in downstream analyses (e.g., PCA).

### Scale data
#### This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate.
1. Shifts the expression of each gene, so that the mean expression across cells is 0.
2. Scales the expression of each gene, so that the variance across cells is 1.
    - ```ScaleData(features = default is variable features)```
* Standard pre-processing step prior to dimensional reduction techniques like PCA

### Linear dimensional reduction (PCA)
#### Linear transformation of the original dataset into principal components ranked in decreasing order of variance 
  #### The variance of the data is maximized in the lower dimensional space
- Dimensionality reduction has two components:
    - Feature selection (selection of a smaller subset from the original set of variables)
        - Based on the assumption that genes showing high variability correspond to biological variation
    - Feature extraction (high dimension data is projected to a lower dimension)
1. Run the principal component analysis ```RunPCA()```
  - By default, only previously determined variable features are used as input.
3. Visualize the results
    - ```VizDimLoadings(dims = 1:2, reduction = "pca")```
    - ```DimPlot(reduction = “pca”)```
    - ```DimHeatMap(dims = 1:4, cells = 500, balanced = TRUE)```
        - Setting cells to a number plots the “extreme” cells on both ends of the spectrum and dramatically decreases plotting time for large datasets.
3. Determine the dimensionality of the dataset (number of PCs to use)
    - ```ElbowPlot()```
> [!NOTE]
> scRNAseq has a highly non-linear structure, so PCA alone is not best suited for data visualization.

### Cluster the cells
#### 1. Cells are embedded into a graph structure, edges are drawn between cells with similar feature expression patterns, and the graph is partitioned into highly interconnected communties/clusters.

1. A k-nearest neighbors (kNN) graph is constructed based on the Euclidean distance in the PCA space and the edge weights between any two cells are refined based on the shared overlap in their local neighborhoods (Jaccard similarity).
    - ```FindNeighbors(dims = # of PCs decided upon in the previous step)```
2. Modularity optimization techniques (default = Louvain algorithm) are applied to iteratively group cells together with the goal of optimizing the standard modularity function.
    - ```FindClusters(resolution = c(0.1,0.3, 0.5, 0.7, 1))```
    - Resolution parameter sets the granularity of the downstream clustering = higher values result in greater numbers of clusters
    * 0.4 -1.2 typically has good results for single-cell datasets of around 3000 cells and optimal resolution often increases for larger datasets
    - Visulaize the effect of the different resolution values
      - ```DimPlot(group.by = 'RNA_snn_res.0.1', label = TRUE)```
      - ```Idents() <- 'RNA_snn_res.0.1'```

### Non-linear dimensional reduction (UMAP)
#### Cells that are grouped together within graph-based clusters determined during the clustering step should co-localize on these dimension reduction plots
1. Run non-linear dimensional reduction and examine the results
    - ```RunUMAP(dims = # of PCs determined previously)```
    - ```DimPlot(reduction = “umap”)```
* UMAP aims to preserve local distances in the dataset, but often does not preserve more glocal relationships
> [!NOTE]
> UMAP can be used for visualization, but biological conclusions should not be drawn solely from this visualization technique

### Correct for batch effects 
#### Were batch effects observed in the UMAP = cells clustering for reasons other the biological variability
1. Perform integration
    - IntegrateLayers(method = HarmonyIntegration, orig.reduction = “pca”, new.reduction = “harmony”)
    * Performs integration in low-dimensional space and returns a dimensional object
2. Repeat clustering with integrated data
    - ```FindNeighbors(reduction = “harmony”)``` 
    - ```FindClusters(cluster.name = “harmony_clusters”)```
3. Re-run non-linear dimensional reduction
    - ```RunUMAP(reduction = “harmony”, reduction.name = “umap_harmony”)```
4. Visualize the results to see if batch effects persist
    - ```DimPlot(reduction = “umap_harmony”)```
    - ```VlnPlot(group.by = “harmony_clusters”)```
5. Rejoin the layers if you want to perform differential expression analysis
    - Collapses the individual datasets together and recreates the original counts and data layers
    - ```JoinLayers()```  (this needs to be done before performing any differential expression analysis)
> [!NOTE]
> Seurat supports five integration methods: Anchor-based CCA integration, Anchor-based RPCA integration, Harmony, FastMNN, and scVI.
