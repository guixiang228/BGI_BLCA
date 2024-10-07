## infer the differentiation trajectory by RNA velocity
import os,sys
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy.io import mmread
import anndata as ad
import matplotlib.pyplot as plt
import scvelo as scv
import cellrank as cr


msn = sc.read("merge_CAF_5Myo_umap_new.h5ad")

col = ['pANN_0.25_5e.04_1950', 'pANN_0.25_5e.04_1445', 'pANN_0.25_5e.04_2321', 'pANN_0.25_5e.04_1942', 'pANN_0.25_5e.04_3813', 'pANN_0.25_0.001_1079', 'pANN_0.25_0.001_1337', 'pANN_0.25_0.001_1062', 'pANN_0.25_5e.04_3120', 'pANN_0.25_5e.04_3486', 'pANN_0.25_5e.04_2159', 'pANN_0.25_5e.04_1461', 'pANN_0.25_0.001_1379', 'pANN_0.25_0.001_876', 'pANN_0.25_5e.04_2439', 'pANN_0.25_5e.04_2374', 'pANN_0.25_5e.04_2581', 'pANN_0.25_5e.04_2245', 'pANN_0.25_5e.04_1943', 'pANN_0.25_0.04_138', 'pANN_0.25_5e.04_1796', 'pANN_0.25_5e.04_1562', 'pANN_0.25_0.001_1340', 'pANN_0.25_0.001_1397', 'pANN_0.25_0.02_296', 'pANN_0.25_0.001_779']
for i in col:
    msn.obs = msn.obs.drop(i, axis=1)


path = './rawData/'
msn.obs_names.unique()

l = []
if os.path.exists(path):
        file_list1 = ['C230727005','C230727006']
        for file1 in file_list1: 
            a = '/RNAVelocityMatrix' 
            path1 = f'/rawData/{file1}{a}'
            file_list = os.listdir(path1)
            for file in file_list:
                file_dir = os.path.join(path1, file)
                if file == "barcodes.tsv.gz":
                    obs = pd.read_csv(file_dir, header = None, index_col = 0, sep = "\t")
                    obs.index.name = "barcode"
                if file == "features.tsv.gz":
                    var = pd.read_csv(file_dir, header = None, index_col = 0, sep = "\t")
                    var.index.name = "gene"
                if file == "spliced.mtx.gz":
                    splice_mtx = csr_matrix(mmread(file_dir).T)
                if file == "unspliced.mtx.gz":
                    unsplice_mtx = csr_matrix(mmread(file_dir).T)
            mtx = splice_mtx + unsplice_mtx
            adata = sc.AnnData(mtx, obs = obs, var = var)
            adata.obs['sample'] = file1
            adata.var['gene_ids'] = adata.var_names
            adata.layers['spliced'] = splice_mtx
            adata.layers['unspliced'] = unsplice_mtx
            names = adata.obs_names.tolist() 
            samples = adata.obs['sample'].tolist() 
            barcode = [i + '_' + j for i, j in zip(samples, names)]
            adata.obs_names = barcode
            l.append(adata[adata.obs_names.isin(msn.obs_names),adata.var_names.isin(msn.var_names)])


data = ad.concat(l, join = "outer")
print(len(data.obs))

np.random.seed(1223)
num = np.random.choice(len(msn.obs.index), len(data.obs), replace=False)
a = msn.obs.index[num]
msn = msn[msn.obs.index.isin(a)]

msn = msn[:, data.var_names].copy()
msn = msn.raw.to_adata()

msn.layers['spliced'] = data.layers['spliced']
msn.layers['unspliced'] = data.layers['unspliced']
msn.write("add_spliced_unspliced.h5ad")


scv.set_figure_params()
adata = scv.read("add_spliced_unspliced.h5ad", cache = True)

sc.pp.neighbors(adata, n_pcs=10, n_neighbors=10)
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
adata.X = np.nan_to_num(adata.X)
scv.tl.velocity(adata)
adata.layers['velocity'].shape
scv.tl.velocity_graph(adata)
scv.tl.umap(adata)


fig, ax = plt.subplots(frameon=False, figsize=(7.5, 6))
ax.set_facecolor('white')
ax.grid(False)
scv.pl.velocity_embedding_stream(adata, basis='umap', ax=ax, color='celltype', size=60, legend_fontsize=6, alpha=0.9, save='embedding_stream.png')
scv.pl.velocity_embedding_stream(adata, basis='umap', color='celltype', size=60, legend_fontsize=6, alpha=0.9, save='embedding_stream.svg')
scv.pl.velocity_embedding_grid(adata, basis='umap', ax=ax, color='celltype', save = 'embedding_grid.pdf')
adata.write("Scvelo.h5ad")

