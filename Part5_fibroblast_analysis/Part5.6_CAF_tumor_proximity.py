### calculate the proportion of fibroblast subtype around tumor cells in different domains
import os, sys
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from sklearn.metrics import pairwise_distances
import squidpy as sq

od = './nhood_30um_all/'
os.chdir(od)

adata = sc.read("fibsubtype.h5ad")
sc.set_figure_params(facecolor="white", figsize=(4, 4))
sc.settings.verbosity = 3
sc.settings.dpi = 300
sc.settings.figdir = "./figures"

adata.obs['spatial_domain'] = adata.obs['spatial_domain'].astype('str')
adata.obs['spatial_domain'] = adata.obs['spatial_domain'].astype('category')
adata = adata[adata.obs['spatial_domain'].isin(['11','2','10','0','3','14'])]
cluster_annotation = {'11':'SD11','2':'SD2','10':'SD10','0':'SD0','3':'SD3','14':'SD14'}
adata.obs['region'] = adata.obs['spatial_domain'].map(cluster_annotation).astype('category')
adata = adata[adata.obs['region'].isin(['SD11','SD2','SD10','SD0','SD3','SD14'])]
existing_categories = adata.obs['fibsubtype'].cat.categories
adata.obs['fibsubtype'].cat.add_categories(['other'], inplace=True)
adata.obs['fibsubtype'] = adata.obs['fibsubtype'].fillna('other')

## Calculate neighboring cells
sq.gr.spatial_neighbors(adata, coord_type='generic', library_key='batch', radius=60)

## Calculate the proportion of neighbors of different fibsubtypes around tumor cell
def calculate_fibsubtype_proportion(adata):

    connectivities = adata.obsp['spatial_connectivities'].tocsr()
    result = pd.DataFrame(index=adata.obs['region'].unique(), columns=adata.obs['fibsubtype'].unique())
    result[:] = np.nan
    cell_indices = np.arange(len(adata))
    
    for region in adata.obs['region'].unique():
        for fibsubtype in adata.obs['fibsubtype'].unique():

            cells_in_region = cell_indices[(adata.obs['region'] == region) & (adata.obs['annotated_cluster'].isin(['NMIBC','MIBC','Metastatic_tumor']))]
            different_fibsubtype_count = 0
            total_neighbors = 0
            
            for cell in cells_in_region:
                start = connectivities.indptr[cell]
                end = connectivities.indptr[cell + 1]
                neighbors_of_cell = connectivities.indices[start:end]
                neighbors_fibsubtype = adata.obs['fibsubtype'].iloc[neighbors_of_cell]
                different_fibsubtype_count += (neighbors_fibsubtype == fibsubtype).sum()
                total_neighbors += len(neighbors_of_cell)
                
            if total_neighbors > 0:
                result.at[region, fibsubtype] = different_fibsubtype_count / total_neighbors
    return result

result_df = calculate_fibsubtype_proportion(adata)
print(result_df)

result_df_normalized = result_df.drop('other', axis=1)
result_df_normalized = result_df_normalized.div(result_df_normalized.abs().sum(axis=0), axis=1)

for col in result_df_normalized.columns:
    result_df_normalized[col] = pd.to_numeric(result_df_normalized[col], errors='coerce')
result_df_normalized = result_df_normalized.T
order = ['SD11','SD2','SD10','SD0','SD3','SD14']
result_df_normalized = result_df_normalized[order]

import seaborn as sns
sns.set_style('whitegrid')
fig,ax = plt.subplots(figsize = (3.5,6))
sns.heatmap(result_df_normalized, cmap='coolwarm',center = 0,ax=ax,linewidths=0.3,vmin=0,vmax=result_df_normalized.max().max())
ax.set_facecolor('#dcdcdc')
plt.savefig('region_nhood_heatmap.pdf', dpi=300, bbox_inches='tight')
