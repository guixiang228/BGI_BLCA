### calculate the proportion of fibroblast subtype around margin cells in different domains
import os, sys
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from sklearn.metrics import pairwise_distances
import squidpy as sq

od = './fibsubtype/'
os.chdir(od)

adata = sc.read("all_fibsubtype.h5ad")

sc.set_figure_params(facecolor="white", figsize=(4, 4))
sc.settings.verbosity = 3
sc.settings.dpi = 300
sc.settings.figdir = "./figures"

adata.obs['spatial_domain'] = adata.obs['spatial_domain'].astype('str')
adata.obs['spatial_domain'] = adata.obs['spatial_domain'].astype('category')
cluster_annotation = {'0':'SD0','1':'SD1','2':'SD2','3':'SD3','4':'SD4','5':'SD5','6':'SD6','7':'SD7','8':'SD8','9':'SD9','10':'SD10','11':'SD11','12':'SD12','13':'SD13','14':'SD14'}
adata.obs['region'] = adata.obs['spatial_domain'].map(cluster_annotation).astype('category')

adata.obs['fibsubtype'] = adata.obs['fibsubtype'].astype('category')
adata.obs['fibsubtype'] = adata.obs['fibsubtype'].cat.add_categories(['other'])
adata.obs['fibsubtype'] = adata.obs['fibsubtype'].fillna('other')

sq.gr.spatial_neighbors(adata, coord_type='generic', library_key='batch', radius=60)

def calculate_fibsubtype_proportion(adata):
    connectivities = adata.obsp['spatial_connectivities'].tocsr() 
    result = pd.DataFrame(index=adata.obs['region'].unique(), columns=adata.obs['fibsubtype'].unique())
    result[:] = np.nan
    
    cell_indices = np.arange(len(adata))
    
    for region in adata.obs['region'].unique():
        for fibsubtype in adata.obs['fibsubtype'].unique():
            
            cells_in_region = cell_indices[adata.obs['region'] == region]
            different_fibsubtype_count = 0
            total_neighbors = 0
            
            for cell in cells_in_region:
                
                start = connectivities.indptr[cell]
                end = connectivities.indptr[cell + 1]
                neighbors_of_cell = connectivities.indices[start:end]
                ## select margin cells
                neighbors_region = adata.obs['region'].iloc[neighbors_of_cell]
                proportion = (neighbors_region == region).sum() / len(neighbors_of_cell)
                if proportion < 0.5:
                    neighbors_fibsubtype = adata.obs['fibsubtype'].iloc[neighbors_of_cell]
                    different_fibsubtype_count += (neighbors_fibsubtype == fibsubtype).sum()
                    total_neighbors += len(neighbors_of_cell)
                
            if total_neighbors > 0:
                result.at[region, fibsubtype] = different_fibsubtype_count / total_neighbors
    return result

result_df = calculate_fibsubtype_proportion(adata)
result_df_normalized = result_df.drop('other', axis=1)
result_df_normalized = result_df_normalized.div(result_df_normalized.abs().sum(axis=0), axis=1)

for col in result_df_normalized.columns:
    result_df_normalized[col] = pd.to_numeric(result_df_normalized[col], errors='coerce')
result_df_normalized = result_df_normalized.T
order = ['SD11','SD2','SD10','SD0','SD3','SD14','SD4','SD8','SD9','SD5','SD1','SD6','SD7','SD12','SD13']
result_df_normalized = result_df_normalized[order]

import seaborn as sns
sns.set_style('whitegrid')
fig,ax = plt.subplots(figsize = (7,7))
sns.heatmap(result_df_normalized, cmap='coolwarm',center = 0,ax=ax,linewidths=0.3,vmin=0,vmax=result_df_normalized.max().max())
ax.set_facecolor('#dcdcdc')
plt.savefig('region_nhood_heatmap_all.pdf', dpi=300, bbox_inches='tight')
