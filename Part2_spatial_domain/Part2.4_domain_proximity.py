## Calculate the proportion of domains in the neighbors of each domain
import os, sys
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from sklearn.metrics import pairwise_distances
import squidpy as sq

od = './SD/'
os.chdir(od)

adata = sc.read("adata.h5ad")
sc.set_figure_params(facecolor="white", figsize=(4, 4))
sc.settings.verbosity = 3
sc.settings.dpi = 300
sc.settings.figdir = "./figures"

adata.obs['spatial_domain'] = adata.obs['spatial_domain'].astype('str')
adata.obs['spatial_domain'] = adata.obs['spatial_domain'].astype('category')
cluster_annotation = {'0':'SD0','1':'SD1','2':'SD2','3':'SD3','4':'SD4','5':'SD5','6':'SD6','7':'SD7','8':'SD8','9':'SD9','10':'SD10','11':'SD11','12':'SD12','13':'SD13','14':'SD14'}
adata.obs['region'] = adata.obs['spatial_domain'].map(cluster_annotation).astype('category')
print(adata.obs['region'].value_counts())


### define the neighboring cells in specific radius
sq.gr.spatial_neighbors(adata, coord_type='generic', library_key='batch', radius=60)
### calculate proportion
def calculate_proportion(adata):
    connectivities = adata.obsp['spatial_connectivities'].tocsr()
    result = pd.DataFrame(index=adata.obs['region'].unique(), columns=adata.obs['region'].unique())
    result[:] = np.nan
    
    cell_indices = np.arange(len(adata))
    
    for region in adata.obs['region'].unique():
        for subtype in adata.obs['region'].unique():
            cells_in_region = cell_indices[adata.obs['region'] == region]
            different_subtype_count = 0
            total_neighbors = 0
            
            for cell in cells_in_region:
                start = connectivities.indptr[cell]
                end = connectivities.indptr[cell + 1]
                neighbors_of_cell = connectivities.indices[start:end]
                neighbors_subtype = adata.obs['region'].iloc[neighbors_of_cell]
                proportion = (neighbors_subtype == region).sum() / len(neighbors_of_cell)
                if proportion < 0.5:
                    different_subtype_count += (neighbors_subtype == subtype).sum()
                    total_neighbors += len(neighbors_of_cell)
                
            if total_neighbors > 0:
                result.at[region, subtype] = different_subtype_count / total_neighbors
    return result

result_df = calculate_proportion(adata)
print(result_df)


result_df_normalized = result_df
for col in result_df_normalized.columns:
    result_df_normalized[col] = pd.to_numeric(result_df_normalized[col], errors='coerce')
order = ['SD0','SD1','SD2','SD3','SD4','SD5','SD6','SD7','SD8','SD9','SD10','SD11','SD12','SD13','SD14']
result_df_normalized = result_df_normalized[order]
result_df_normalized = result_df_normalized.reindex(order)

import seaborn as sns
sns.set_style('whitegrid')
fig,ax = plt.subplots(figsize = (7,7))
sns.heatmap(result_df_normalized, cmap='coolwarm',center = 0,ax=ax,linewidths=0.3,vmin=0,vmax=result_df_normalized.max().max())
ax.set_facecolor('#dcdcdc')
plt.savefig('region_nhood_heatmap.pdf', dpi=300, bbox_inches='tight')
