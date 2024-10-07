### Extract the tumor cells in the 30um radius of immune cells
import os, sys
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from sklearn.metrics import pairwise_distances
import squidpy as sq

adata = sc.read("immune_subtype.h5ad")
sc.set_figure_params(facecolor="white", figsize=(4, 4))
sc.settings.verbosity = 3
sc.settings.dpi = 300
sc.settings.figdir = "./figures"


sq.gr.spatial_neighbors(adata, coord_type='generic', library_key='batch', radius=60)

existing_categories = ['APOE_Macrophages', 'B_cells', 'CA10_T_cells', 'CD1C_cDC2', 'CD8A_CTL',
       'CLEC9A_cDC1', 'LAMP3_cDC3', 'LYVE1_Macrophages', 'Mast_cells',
       'NCAM1_NK_cells', 'NLRP3_Macrophages', 'Plasma_cells', 'TCF7_T_cells',
       'TNFRSF9_T_cells', 'TOP2A_T_cells']
connectivities = adata.obsp['spatial_connectivities'].tocsr()
cell_indices = np.arange(len(adata))

st_ad_list = []
for immune in existing_categories:
    cells_in_region = cell_indices[(adata.obs['immune_subtype'] == immune)]
    indices = []
    for cell in cells_in_region:
        start = connectivities.indptr[cell]
        end = connectivities.indptr[cell + 1]
        neighbors_of_cell = connectivities.indices[start:end]
        valid_neighbors = neighbors_of_cell[adata.obs['annotated_cluster'].iloc[neighbors_of_cell].isin(['NMIBC', 'MIBC', 'Metastatic_tumor', 'Tumor_Epithelial'])]
        indices.extend(valid_neighbors.tolist())
    indices = list(set(indices))
    adata_select = adata[indices]
    adata_select.obs['tumor_type'] = immune
    st_ad_list.append(adata_select)

adata = sc.AnnData.concatenate(*st_ad_list, join='outer')
adata.write('tumor_immune_nei.h5ad')


### Extract tumor cells in different neighboring radius
adata_30um = sc.read_h5ad('./30um/tumor_immune_nei.h5ad')
adata_50um = sc.read_h5ad('./50um/tumor_immune_nei.h5ad')
adata_100um = sc.read_h5ad('./100um/tumor_immune_nei.h5ad')
adata_200um = sc.read_h5ad('./200um/tumor_immune_nei.h5ad')
adata_300um = sc.read_h5ad('./300um/tumor_immune_nei.h5ad')


for immune in existing_categories:
    st_ad_list = []
    folder_name = f'{immune}'
    folder_path = os.path.join(od, folder_name)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        os.chdir(folder_path)
    adata_30um_immune = adata_30um[adata_30um.obs['tumor_type']==immune]
    adata_50um_immune = adata_50um[adata_50um.obs['tumor_type']==immune]
    adata_100um_immune = adata_100um[adata_100um.obs['tumor_type']==immune]
    adata_200um_immune = adata_200um[adata_200um.obs['tumor_type']==immune]
    adata_300um_immune = adata_300um[adata_300um.obs['tumor_type']==immune]

    cells_30um = set(list(adata_30um_immune.obs.index))
    cells_50um = set(list(adata_50um_immune.obs.index))
    cells_100um = set(list(adata_100um_immune.obs.index))
    cells_200um = set(list(adata_200um_immune.obs.index))
    cells_300um = set(list(adata_300um_immune.obs.index))

    adata_30um_immune.obs['distance'] = '30um'
    st_ad_list.append(adata_30um_immune)

    unique_cells_50um = list(cells_50um-cells_30um)
    adata_50um_unique = adata_50um_immune[adata_50um_immune.obs.index.isin(unique_cells_50um), :]
    adata_50um_unique.obs['distance'] = '50um'
    st_ad_list.append(adata_50um_unique)

    unique_cells_100um =  list(cells_100um - (cells_30um.union(cells_50um)))
    adata_100um_unique = adata_100um_immune[adata_100um_immune.obs.index.isin(unique_cells_100um), :]
    adata_100um_unique.obs['distance'] = '100um'
    st_ad_list.append(adata_100um_unique)

    unique_cells_200um =  list(cells_200um - (cells_30um.union(cells_50um,cells_100um)))
    adata_200um_unique = adata_200um_immune[adata_200um_immune.obs.index.isin(unique_cells_200um), :]
    adata_200um_unique.obs['distance'] = '200um'
    st_ad_list.append(adata_200um_unique)

    unique_cells_300um =  list(cells_300um - (cells_30um.union(cells_50um,cells_100um,cells_200um)))
    adata_300um_unique = adata_300um_immune[adata_300um_immune.obs.index.isin(unique_cells_300um), :]
    adata_300um_unique.obs['distance'] = '300um'
    st_ad_list.append(adata_300um_unique)

    adata = sc.AnnData.concatenate(*st_ad_list, join='outer')
    adata.write('tumor_distance.h5ad')
    