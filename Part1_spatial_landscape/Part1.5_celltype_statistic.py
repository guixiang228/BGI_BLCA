## cell type proportion in samples of different stage
import os,sys
import scanpy as sc
import numpy as np
import pandas as pd
from anndata import AnnData

import matplotlib.pyplot as plt
from collections import defaultdict
import seaborn as sns
from PIL import Image, ImageDraw, ImageFont
import random

os.chdir("./figures")

adata = sc.read_h5ad('allsample.h5ad')
adata.obs.loc[adata.obs['batch'].isin(['HBCP1','HBCP7','HBCP9','HBCP13','HBCP14','HBCP15','HBCP17','HBCP5A','HBCP5B','HBCP3C','HBCP11']),'stage'] = 'T1'
adata.obs.loc[adata.obs['batch'].isin(['HBCP8','HBCP3B','HBCP6','HBCP18','HBCP10','HBCP3A','HBCP12']),'stage'] = 'T2+'
adata = adata[adata.obs['stage'].notna()]
adata.obs['batch'] = pd.Categorical(adata.obs['batch'], categories=['HBCP1', 'HBCP3C', 'HBCP5A', 'HBCP5B','HBCP7', 'HBCP9', 'HBCP11', 'HBCP13', 'HBCP14', 'HBCP15', 'HBCP17',  'HBCP12','HBCP18','HBCP3B','HBCP6','HBCP3A','HBCP10','HBCP8'], ordered=True)

### Box plot for all samples
cell_colors = {
    "Bcell": "#FFBBFF",
    "DCS": "#81d8ae",
    "Endothelial": "#607d3a",
    "Fibroblast": "#FFD700",
    "Macrophage": "#87CEFF",
    "Mastcell": "#AB82FF",
    "Myocyte": "#3399FF",
    "Pericyte": "#FFD39B",
    "Plasmocyte": "#6bd155",
    "Tcell": "#B23AEE",
    "Tumor_cell": "#CD3278",
    "Normal_Epithelial": "#663300"
}
cell_type_counts = adata.obs.groupby(['batch', 'annotated_cluster']).size()
group_totals = adata.obs.groupby('batch').size()
cell_type_proportions = cell_type_counts.divide(group_totals, level='batch')
cell_type_proportions_df = cell_type_proportions.reset_index()
cell_type_proportions_df.columns = ['batch', 'annotated_cluster', 'Proportion']
pivot_df = cell_type_proportions_df.pivot_table(index='batch', columns='annotated_cluster', values='Proportion', aggfunc='sum')

fig, ax = plt.subplots(figsize=(8,5))
bottom = np.zeros(len(pivot_df))
for column in pivot_df.columns:
    color = cell_colors.get(column, 'gray')
    ax.bar(pivot_df.index, pivot_df[column], bottom=bottom, color=color, label=column)
    bottom += pivot_df[column]
ax.set_title('Cell Type Proportions by Batch', fontsize=14)
ax.set_xlabel('Batch', fontsize=12)
ax.set_ylabel('Proportion', fontsize=12)
ax.set_xticks(range(len(pivot_df.index)))
ax.set_xticklabels(pivot_df.index, rotation=90, fontsize=10)
ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
plt.savefig('batch_type.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.close()



#### boxplot of celltype proportion with stage seperation
df = adata.obs.copy()
proportions = df.groupby(['annotated_cluster', 'stage', 'batch']).size().unstack()
proportions = proportions.replace(0, np.nan)
proportions = proportions.reset_index()
proportions_long = proportions.melt(id_vars=['stage', 'annotated_cluster'], var_name='batch', value_name='count')


proportions_long['proportion'] = proportions_long.groupby('batch').apply(
    lambda x: x['count'] / x['count'].sum()
).reset_index(level=0, drop=True)

g = sns.FacetGrid(proportions_long, col='annotated_cluster', col_wrap=12, height=4, aspect=0.4, sharey=False)
g.map_dataframe(sns.boxplot, x='stage', y='proportion', palette=['#4169E1','#FF7F24'], zorder=1)

proportions_long['color'] = proportions_long.apply(lambda row: '#4169E1' if row['stage'] == 'T1' else '#FF7F24', axis=1)

def scatter_white(data, x, y, c, **kwargs):
    plt.scatter(data[x], data[y], facecolors='white', edgecolors=data[c], s=80, linewidths=1.5, zorder=3, **kwargs)

g.map_dataframe(scatter_white, x='stage', y='proportion', c='color')
g.set_titles(col_template="{col_name}", fontsize=14)
g.set_axis_labels("Stage", "Proportion", fontsize=14)
g.set_xticklabels(rotation=30, fontsize=12)
g.set_yticklabels(fontsize=12)
g.tight_layout()
plt.savefig('boxplot_cellType_2stage_by_cluster.pdf', dpi=300, bbox_inches='tight')
plt.close()
