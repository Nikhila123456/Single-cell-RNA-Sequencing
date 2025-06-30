import scvelo as scv
import scanpy as sc
#import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
#cr.settings.verbosity = 2

adata = sc.read_h5ad("/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/3_scvelo/Lsc.h5ad")
#read in cell IDs
sample_obs = pd.read_csv("/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/2_export_SeuratInfo/Lsc_cellID_obs.csv")
#get cell IDs for each sample
cellID_obs_naiv11 = sample_obs[sample_obs["x"].str.contains("_Naiv-1.1")]
cellID_obs_naiv12 = sample_obs[sample_obs["x"].str.contains("_Naiv-1.2")]
cellID_obs_naiv21 = sample_obs[sample_obs["x"].str.contains("_Naiv-2.1")]
cellID_obs_naiv22 = sample_obs[sample_obs["x"].str.contains("_Naiv-2.2")]
cellID_obs_naiv31 = sample_obs[sample_obs["x"].str.contains("_Naiv-3.1")]
cellID_obs_naiv32 = sample_obs[sample_obs["x"].str.contains("_Naiv-3.2")]
cellID_obs_inf11 = sample_obs[sample_obs["x"].str.contains("_Inf-1.1")]
cellID_obs_inf21 = sample_obs[sample_obs["x"].str.contains("_Inf-2.1")]
cellID_obs_inf22 = sample_obs[sample_obs["x"].str.contains("_Inf-2.2")]
cellID_obs_inf31 = sample_obs[sample_obs["x"].str.contains("_Inf-3.1")]
cellID_obs_inf32 = sample_obs[sample_obs["x"].str.contains("_Inf-3.2")]



# load files for spliced/unspliced matrices for each sample:
ldatanai11 = scv.read('/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/1_quantification/velocyto/063395/possorted_genome_bam_RDKTN.loom')
ldatanai12 = scv.read('/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/1_quantification/velocyto/063396/possorted_genome_bam_15T5V.loom')
ldatanai21 = scv.read('/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/1_quantification/velocyto/Nai2_1/possorted_genome_bam_ISXXE.loom')
ldatanai22 = scv.read('/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/1_quantification/velocyto/Nai2_2/possorted_genome_bam_7KGUR.loom')
ldatanai31 = scv.read('/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/1_quantification/velocyto/Nai3_1/possorted_genome_bam_QBT7J.loom')
ldatanai32 = scv.read('/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/1_quantification/velocyto/Nai3_2/possorted_genome_bam_E2HXL.loom')
ldatainf11 = scv.read('/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/1_quantification/velocyto/063394/possorted_genome_bam_PNBU0.loom')
ldatainf21 = scv.read('/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/1_quantification/velocyto/Inf2_1/possorted_genome_bam_HAA1D.loom')
ldatainf22 = scv.read('/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/1_quantification/velocyto/Inf2_2/possorted_genome_bam_2NDU3.loom')
ldatainf31 = scv.read('/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/1_quantification/velocyto/Inf3_1/possorted_genome_bam_CEWGM.loom')
ldatainf32 = scv.read('/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/1_quantification/velocyto/Inf3_2/possorted_genome_bam_R8T80.loom')



#rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldatanai11.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_Naiv-1.1' for bc in barcodes]
ldatanai11.obs.index = barcodes
ldatanai11 = ldatanai11[np.isin(ldatanai11.obs.index, cellID_obs_naiv11["x"])]

barcodes = [bc.split(':')[1] for bc in ldatanai12.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_Naiv-1.2' for bc in barcodes]
ldatanai12.obs.index = barcodes
ldatanai12 = ldatanai12[np.isin(ldatanai12.obs.index,cellID_obs_naiv12["x"])]


barcodes = [bc.split(':')[1] for bc in ldatanai21.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_Naiv-2.1' for bc in barcodes]
ldatanai21.obs.index = barcodes
ldatanai21 = ldatanai21[np.isin(ldatanai21.obs.index, cellID_obs_naiv21["x"])]

barcodes = [bc.split(':')[1] for bc in ldatanai22.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_Naiv-2.2' for bc in barcodes]
ldatanai22.obs.index = barcodes
ldatanai22 = ldatanai22[np.isin(ldatanai22.obs.index, cellID_obs_naiv22["x"])]

barcodes = [bc.split(':')[1] for bc in ldatanai31.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_Naiv-3.1' for bc in barcodes]
ldatanai31.obs.index = barcodes
ldatanai31 = ldatanai31[np.isin(ldatanai31.obs.index, cellID_obs_naiv31["x"])]

barcodes = [bc.split(':')[1] for bc in ldatanai32.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_Naiv-3.2' for bc in barcodes]
ldatanai32.obs.index = barcodes
ldatanai32 = ldatanai32[np.isin(ldatanai32.obs.index, cellID_obs_naiv32["x"])]

barcodes = [bc.split(':')[1] for bc in ldatainf11.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_Inf-1.1' for bc in barcodes]
ldatainf11.obs.index = barcodes
ldatainf11 = ldatainf11[np.isin(ldatainf11.obs.index, cellID_obs_inf11["x"])]

barcodes = [bc.split(':')[1] for bc in ldatainf21.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_Inf-2.1' for bc in barcodes]
ldatainf21.obs.index = barcodes
ldatainf21 = ldatainf21[np.isin(ldatainf21.obs.index, cellID_obs_inf21["x"])]

barcodes = [bc.split(':')[1] for bc in ldatainf22.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_Inf-2.2' for bc in barcodes]
ldatainf22.obs.index = barcodes
ldatainf22 = ldatainf22[np.isin(ldatainf22.obs.index, cellID_obs_inf22["x"])]

barcodes = [bc.split(':')[1] for bc in ldatainf31.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_Inf-3.1' for bc in barcodes]
ldatainf31.obs.index = barcodes
ldatainf31 = ldatainf31[np.isin(ldatainf31.obs.index, cellID_obs_inf31["x"])]

barcodes = [bc.split(':')[1] for bc in ldatainf32.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_Inf-3.2' for bc in barcodes]
ldatainf32.obs.index = barcodes
ldatainf32 = ldatainf32[np.isin(ldatainf32.obs.index, cellID_obs_inf32["x"])]

# make variable names unique
ldatanai11.var_names_make_unique()
ldatanai12.var_names_make_unique()
ldatanai21.var_names_make_unique()
ldatanai22.var_names_make_unique()
ldatanai31.var_names_make_unique()
ldatanai32.var_names_make_unique()
ldatainf11.var_names_make_unique()
ldatainf21.var_names_make_unique()
ldatainf22.var_names_make_unique()
ldatainf31.var_names_make_unique()
ldatainf32.var_names_make_unique()

# concatenate the three loom
ldata = ldatanai11.concatenate([ldatanai12, ldatanai21,ldatanai22,ldatanai31,ldatanai32,ldatainf11,ldatainf21,ldatainf22,ldatainf31,ldatainf32])

# merge matrices into the original adata object

scv.utils.clean_obs_names(adata)
scv.utils.clean_obs_names(ldata)
# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)

#plot a umap
gyr = ['#AE3121','#BCF5BE', '#755F00','#3D85C6','#0AF412','#F7B903','#07F2FD','#FF9D02','#C90076','#8306FA']
sc.pl.umap(adata, color='seurat_clusters', frameon=False, legend_loc='on data', palette= gyr, title='', save='lsc_celltypes.pdf')

#inspect spliced and unspliced reads
scv.pl.proportions(adata, groupby='seurat_clusters', save='/depot/tlratlif/data/LSC_scRNAseq_2023/3_velocity/3_scvelo/lsc_clusters_celltypes.pdf')

# pre-processing
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

# compute velocity
#scv.tl.velocity(adata, mode='stochastic')
scv.tl.recover_dynamics(adata, var_names='all')
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

adata.write('Lsc_velocity.h5ad',compression='gzip')


scv.pl.velocity_embedding(adata, basis='umap',frameon=False, save='Lsc_embedding.pdf')
scv.pl.velocity_embedding_grid(adata, basis='umap',color='seurat_clusters', palette=gyr, save='Lsc_embedding_grid.pdf', title='', scale=0.25)


