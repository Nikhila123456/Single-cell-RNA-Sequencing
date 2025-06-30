import scvelo as scv
import scanpy as sc
#import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
#cr.settings.verbosity = 2

gyr = ['#AE3121','#BCF5BE', '#755F00','#3D85C6','#0AF412','#F7B903','#07F2FD','#FF9D02','#C90076','#8306FA']

adata = scv.read('Lsc_velocity.h5ad')
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='Lsc_embedding.pdf')
scv.pl.velocity_embedding_grid(adata, basis='umap', color='seurat_clusters', palette=gyr, save='Lsc_embedding_grid.pdf', title='', scale=0.25)
scv.pl.velocity_embedding_stream(adata, basis='umap',color='seurat_clusters', palette=gyr,save='Lsc_embedding_stream.pdf')

#plot velocity length and confidence

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95],save='Lsc_confidence.png')


df = scv.get_df(adata, 'fit*', dropna=True).head()
df.to_csv('Lsc_kinetics.csv')

#visualize latent time

scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save='Lsc_latentTimeUMAP.pdf')

top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='seurat_clusters', n_convolve=100, palette=gyr, save='Lsc_latentTime_heatmap.pdf')
