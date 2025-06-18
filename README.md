# 🧬 Single-cell RNA Sequencing Workflow

### Create a Seurat object
##### A Seurat object serves as a container for single-cell dataset data (e.g., count matrix) and analyses (e.g., PCA or clustering results).
1. From CellRanger output (folder containing barcodes.tsv, genes.tsv, and matrix.mtx)
    - Create a sparse count matrix by reading in the data: ```Read10X()```
    - Convert the sparse count matrix in to a Seurat object: ```CreateSeuratObject()```
2. From CellRanger output (file ending in: count_raw_feature_bc_matrix.h5)
    - Create a sparse count matrix by reading in the data: ```Read10X_h5()```
    - Convert the sparse count matrix in to a Seurat object: ```CreateSeuratObject()```
3. From GEO datasets with multiple samples
    - Read in count and metadata files using: ```read_tsv()```
    - Convert the files into Seurat objects: ```purr::map2(counts, meta, ~CreateSeuratObject(counts = as(.x, "sparseMatrix"), meta.data = .y))```
    - Merge the objects into one Seurat object: ```purrr::reduce(obj, function(x,y) {merge(x,y)})```
4. From datasets stored in SeuratData
    - Load in data using: ```LoadData()```

### Quality control and filtering
##### Low quality cells can be the result of cell damage or failure in library preparation, and typically manifest as cells with low total counts, a low level of expressed genes ot high mitochondrial proportions. If low quality cells are not removed, the first few principal components and the genes with the largest variances will be driven quality rather than biology.
1. Proportion of mitochondrial genes
    - A high percentage of mitochondrial genes is indicative of low quality or dying cells, due to the loss of cytoplasmic RNA from perforated cells.
    - Create new column using: ```PercentageFeatureSet(pattern = "^MT-")```
2. Number of unique genes detected in each cell (nFeature_RNA)
    - Low numbers of unique genes is likely a low quality cell, because the diverse transcript population has not been successfully captured.
3. Total number counts detected in each cell (nCount_RNA)
    - Low number of counts is indicative of low quality cells, due to RNA having been lost at some point during the library preparation (e.g., cell lysis or inefficient cDNA capture and amplification).
- A high number of genes or counts can be indicative of doublets or multiplets.
- Create plots to explore data and determine thresholds. 
    - ```VlnPlot(nsclc_seu, features = c("nCount_RNA", "nFeature_RNA", "percent_mt"), ncol = 3)``` 
    - ```FeatureScatter(nsclc_seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')```
- Filter data according to the determined thresholds.
    - ```subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)```
> [!NOTE]
> Quality control metrics are dependent on cell type (e.g, hepatocytes are highly metabolically active and therefore have higher mitochondrial percentages than other cells) and data type (e.g., counts-based data vs UMI-based data).

### Normalize data
##### Eliminates batch effects/technical variation (e.g., sequencing depth).
- Normalizes the feature expression measurements for each cell by the total expression, then multiples this by a scale factor (10,000) and log-transforms the result.
    - ```NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000)```
> ```SCTransform()``` can be used in the workflow in place of ```NormalizeData()```, ```FindVariableFeatures()``` and ```ScaleData()```.
> Log normalization relies on the assumption that each cell originally contained the same number of molecules: SCTransform does not make this assumption.
> Confounding sources of variation, like high mitochondrial percentages, can also be removed with SCTransform().
   > ```PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt")``` followed by
   > ```SCTransform(vars.to.regress = "percent.mt", verbose = FALSE)```

### Find variable features
##### Cells with low cell-to-cell variation (i.e., housekeeping genes) are not very informative, whereas focusing on variable genes helps to highlight the biological signal in single-cell datasets.
1. Find variable features ```FindVariableFeatures(nfeatures = 2000)``` 
2. Examine the variable features ```VariableFeaturePlot(): look at most variable genes```
* Variable features are used with in downstream analyses (e.g., PCA).

### Scale data
##### This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate.
1. Shifts the expression of each gene, so that the mean expression across cells is 0.
2. Scales the expression of each gene, so that the variance across cells is 1.
    - ```ScaleData(features = default is variable features)```
* Standard pre-processing step prior to dimensional reduction techniques like PCA

### Linear dimensional reduction (PCA)
##### Linear transformation of the original dataset into principal components ranked in decreasing order of variance. The variance of the data is maximized in the lower dimensional space.
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
##### Cells are embedded into a graph structure, edges are drawn between cells with similar feature expression patterns, and the graph is partitioned into highly interconnected communties/clusters.
> k-means clustering randomly initializes _k_ cluster centers, assigns points to nearest center and then updates the cluster centers and repeats the assignment process until the centers stop changing. Intitial cluster centers are randomly assigned: set a seed ```set.seed(1234)``` for consistent clustering. 

1. A k-nearest neighbors (kNN) graph is constructed based on the Euclidean distance in the PCA space and the edge weights between any two cells are refined based on the shared overlap in their local neighborhoods (Jaccard similarity).
    - ```FindNeighbors(dims = # of PCs decided upon in the previous step)```
2. Modularity optimization techniques (default = Louvain algorithm) are applied to iteratively group cells together with the goal of optimizing the standard modularity function.
    - ```FindClusters(resolution = c(0.1, 0.3, 0.5, 0.7, 1))```
    - Resolution parameter sets the granularity of the downstream clustering = higher values result in greater numbers of clusters
    * 0.4 -1.2 typically has good results for single-cell datasets of around 3000 cells and optimal resolution often increases for larger datasets
    - Visulaize the effect of the different resolution values
      - ```DimPlot(group.by = 'RNA_snn_res.0.1', label = TRUE)```
      - ```Idents() <- 'RNA_snn_res.0.1'```

### Non-linear dimensional reduction (UMAP)
##### Cells that are grouped together within graph-based clusters determined during the clustering step should co-localize on these dimension reduction plots
1. Run non-linear dimensional reduction and examine the results
    - ```RunUMAP(dims = # of PCs determined previously)```
    - ```DimPlot(reduction = “umap”)```
* UMAP aims to preserve local distances in the dataset, but often does not preserve more glocal relationships
> [!NOTE]
> UMAP can be used for visualization, but biological conclusions should not be drawn solely from this visualization technique. More information on UMAP interpretation can be found here: [Is UMAP Accurate?](https://nikolay-oskolkov.medium.com/is-umap-accurate-fad1b3f14fb5)

### Correct for batch effects 
##### Perform integration if batch effects were observed in the UMAP (e.g., cells clustering based on batches, donors or conditions)
1. Perform integration
    - ```IntegrateLayers(method = HarmonyIntegration, orig.reduction = “pca”, new.reduction = “harmony”)```
    - Performs integration in low-dimensional space and returns a dimensional object (in reductions slot: pca, umap and harmony)
> [!NOTE]
> Seurat supports five integration methods: Anchor-based CCA integration, Anchor-based RPCA integration, Harmony, FastMNN, and scVI. More information about batch correction can be found here: [Batch Correction](https://www.biorxiv.org/content/10.1101/2024.03.19.585562v1).
2. Repeat clustering with integrated data
    - ```FindNeighbors(reduction = “harmony”)``` 
    - ```FindClusters(cluster.name = “harmony_clusters”)```
3. Re-run non-linear dimensional reduction
    - ```RunUMAP(reduction = “harmony”, reduction.name = “umap_harmony”)```
4. Visualize the results to see if batch effects persist
    - ```DimPlot(reduction = “umap_harmony”)```
    - ```VlnPlot(group.by = “harmony_clusters”)```
5. Integration can be performed after using ```SCTransform()````
    - ```SCTransform()```
    - ```RunPCA()```
    - ```RunHarmony(assay.use="SCT", group.by.vars = "Method")```
    - ```FindNeighbors(reduction = "harmony")```
    - ```FindClusters()```
    - ```RunUMAP(reduction = "harmony")```
    - ```DimPlot(group.by = "Method")```
6. Rejoin the layers
    - Collapses the individual datasets together and recreates the original counts and data layers.
    - ```JoinLayers()```
> [!IMPORTANT]
> Layers need to be rejoined before performing differential expression analysis.

### Cluster identification
1. Identify marker genes for the clusters and assign the annotations manually.
    - If your clusters are composed of cells from more than one condition (e.g., control and stimulated).
    - ```FindConservedMarkers(ident.1 = 3, grouping.var = "stim")```
        - Finds the markers in cluster three that are concerved between the control and stimulated groups.
    - Not setting ```ident.2``` to a specific cluster compares the cluster in ```ident.1``` to all other clusters.
    - Visulaize the identified markers
        - ```VlnPlot(features = c("FCGR3A"))```
        - ```FeaturePlot(features = c("FCGR3A"), min.cutoff = "q10")```
        - ```min.cutoff = "q10"```: genes in the tenth quartile will appear grey, increasing contrast.
    - Rename the cluster based on the identified markers.
        - ```Idents()```
        - ```RenameIdents("3" = "CD16 Mono")```
2. Annotation using a reference Seurat object
- Find a set of anchors between a reference and a query object.
    - ```anchors <- FindTransferAnchors(reference = reference, query = query, dims = 1:30, reference.reduction = "pca")```
- Use anchors to transfer data from the reference to the query object.
    - ```predictions <- TransferData(anchorset = anchors, refdata = reference$celltype, dims = 1:30)```
- Add the cell type predictions to the query metadata
    - ```query <- AddMetaData(query, metadata = predictions)```
- Verify cell type predictions by examining canonical cell type marker expression in the clusters.
    - ```VlnPlot(query, c("REG1A", "PPY"), group.by = "predicted.id")```
3. Automated annotation using an annotated object from a database
- The reference needs to contain a cell population sufficiently similar to the cell population in your query object. Using multiple references increases the likelihood that the cell types in your query are accounted for.
    - celldex: provides a collection of reference expression datasets with curated cell type labels.
    - SingleR: uses reference transcriptomic datasets of pure cell types to infer the cell of origin of each of the single cells independently.
    - Default method: performs classification seperately with each reference and then combines the results and chooses the highest overall prediction score.
    - Load the reference libraries: ```hpca <- celldex::HumanPrimaryCellAtlasData()``` and ```dice <- celldex::DatabaseImmuneCellExpressionData()```
    - Get the counts from you query dataset: ```counts <- GetAssayData(obj, slot = "counts")```
    - Run SingleR: ```pred <- SingleR(test = counts, ref = list(HPCA = hpca, DICE = dice), labels = list(hpca$label.main, dice$label.main))```
    - Save the cell type labels to the Seurat metadata: ```obj$pred.labels <- pred[match(rownames(obj@meta.data), rownames(pred)), "labels"]```
        - May have to fix inconsistencies with cell labels.
    - Visualize the labels: ```DimPlot(obj, reduction = "umap", group.by = "pred.labels", label = TRUE)```
    - Get the marker genes from each reference dataset for each cell type: ```metadata(pred$orig.results$HPCA)$de.genes``` and ```metadata(com.res2$orig.results$DICE)$de.genes```
> [!NOTE]
> Lack of consistency in cell labels across references can complicate interpretation.

### Differential gene expression
1. Identify genes that are differentially expressed between conditions for a given cluster.
    - Create a column that has information about both cluster and condition.
        - ```obj$cluster.cnd <- paste0(obj$cluster,"_", obj$condition)```
        - ```Idents(obj) <- obj$cluster.cnd```
    - Find differential markers between the groups
        - ```data.frame <- FindMarkers(obj, ident.1 = "cluster.cnd_1", ident.2 = "cluster.cnd_2")```
    - Visualize the differential gene expression
        - ```FeaturePlot(features = c("FCGR3A", "VMO1"), split.by = "stim", min.cutoff = "q10")```
        - ```VlnPlot(features = c("FCGR3A", "VMO1"))```
        - ```DotPlot(features = c("FCGR3A", "VMO1"), split.by = "stim")```
2. Pseudobulking + differential expression analysis
    - Single cell methods treat each cell as a sample: p-values are inflated, variation across the population is not truly investigated, and there are issues with unmodelled correlations between samples (the samples/single cells are not independent of each other).
    - Aggregate counts across sample, condition and cell type.
        - ```pseudo_obj <- AggregateExpression(obj, assays = "RNA", return.seurat = TRUE, group.by = c("stim", "donor_id", "seurat_annotations"))```
    - Create a column that has information about both cell type and condition.
        - ```obj$celltype.stim <- paste(obj$seurat_annotations, obj$stim, sep = "_")```
        - ```Idents(obj) <- "celltype.stim"```
    - Perform differential expression analysis using DESeq2 at the sample level
        - ```FindMarkers(object = obj, ident.1 = "CD14 Mono_STIM", ident.2 = "CD14 Mono_CTRL", test.use = "DESeq2")```
> [!NOTE]
> Differential expression analysis using ```FindMarkers()``` can be performed using alternative tests, including: "wilcox", "wilcox_limma", "bimod", "roc", "t", "poisson", "negbiom", "LR" and "MAST".

### Useful references
- [Basics of Single-Cell Analysis with Bioconductor
](https://bioconductor.org/books/3.21/OSCA.basic/)
