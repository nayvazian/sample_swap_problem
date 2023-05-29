library(Seurat)
seuratObj <- CreateSeuratObject(counts = counts_df, 
                                      project = "UMAP of Swapped Samples",
                                      assay = "RNAseq",
                                      features = rownames(counts_df))

seuratObj <- AddMetaData(seuratObj, metadata = metadata$donor, col.name = "donor")
seuratObj <- AddMetaData(seuratObj, metadata = metadata$tissue, col.name = "tissue")

# Preprocess and normalize the data
seuratObj <- NormalizeData(seuratObj)
seuratObj <- FindVariableFeatures(seuratObj)
seuratObj <- ScaleData(seuratObj)

#Run PCA and UMAP
seuratObj <- RunPCA(seuratObj, features = VariableFeatures(object = seuratObj), npcs = 30)

seuratObj <- RunUMAP(seuratObj, reduction = "pca", dims = 1:10,n.neighbors = 50)