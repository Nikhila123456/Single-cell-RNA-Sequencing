import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

# load sparse matrix:
X = io.mmread("/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/2_export_SeuratInfo/Lsc_cell_counts.mtx")

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)


# load cell metadata:
cell_meta = pd.read_csv("/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/2_export_SeuratInfo/Lsc_cells_metadata_renamed.csv")

# load gene names:
with open("/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/2_export_SeuratInfo/Lsc_gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes (name them cells since the index cannot have the same name as a column name), var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.obs.index.name = "cells"
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/2_export_SeuratInfo/Lsc_pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:

adata.obs['seurat_clusters'] = list(adata.obs['seurat_clusters'])
adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('category')

#set color palette here (use hex codes for the colors used in R
gyr = ['#F8766D', '#E9842C', '#D69100', '#BC9D00', '#9CA700', '#6FB000', '#00B813', '#00BD61', '#00C08E', '#00C0B4']
##00BDD4', '#00B5EE', '#00A7FF', '#7F96FF', '#BC81FF', '#E26EF7', '#F863DF', '#FF62BF', '#FF6A9A'
sc.pl.umap(
        adata,
        color=['seurat_clusters'],
        frameon=False,
        palette= gyr, 
        save="_Lsc_third_clustering_orig")

# save dataset as anndata format
adata.write('/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/3_scvelo/Lsc.h5ad')

# reload dataset
#adata = sc.read_h5ad('macs.h5ad')
