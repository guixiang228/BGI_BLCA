# annotation for BCa_immune_cell
import os
import gc
import glob
import pandas as pd
import numpy as np
import scanpy as sc
import scanpy.external as sce
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from anndata import read_csv
from loompy import make_row_attrs_from_gene_annotations
from matplotlib.pyplot import rc_context
import warnings

warnings.filterwarnings ('ignore')
np.random.seed(42)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Arial']
plt.rcParams['font.size'] = 6  
sc.set_figure_params(dpi=300, color_map="gist_heat_r")
os.chdir('/Users/zhangqian/Desktop/Immune/')
keys = sorted([f"leiden_res_{res:4.2f}" for res in np.arange(0.1, 2.1, 0.1)])

# STEP01: annotation for seprarate myeloid cell and lymphocytes
sam = 'Immune'
stage = 'BCa'
adata = sc.read_h5ad('BCa_human_immune.h5ad') # load data
adata.obs_names_make_unique()
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata)
sce.pp.harmony_integrate(adata, 'orig.ident')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
sc.pp.neighbors(adata)
for res in np.arange(0.1, 2.1, 0.1):
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res
    )
sc.tl.umap(adata)

res = 0.30
adata.obs["CellTypeL1"] = adata.obs[f"leiden_res_{res:4.2f}"].map(
    {"0" : "Myeloid_cells",
     "1" : "Lymphocytes",
     "2":  "Lymphocytes",
     "3": "Myeloid_cells",
     "4": "Myeloid_cells"}
)
adata.obs['CellTypeL1'] = adata.obs['CellTypeL1'].cat.add_categories('unknown')
adata.obs['CellTypeL1'] = adata.obs['CellTypeL1'].fillna('unknown')
adata.write("BCa_immune_L1_annotation.h5ad") # cluster 5 have mix signal
adata.obs.to_csv('BCa_immune_L1_annotation.csv')


# S2-01, Lymphocytes annotataion
adata = adata[adata.obs['CellTypeL1'] == "Lymphocytes",:]
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata)
sce.pp.harmony_integrate(adata, 'orig.ident')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
sc.pp.neighbors(adata)
for res in np.arange(0.1, 2.1, 0.1):
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res
    )
sc.tl.umap(adata)
res = 0.90
adata.obs["CellTypeL2"] = adata.obs[f"leiden_res_{res:4.2f}"].map(
    {
        "0": "CA10_Th",
        "1": "CCL5_NK",  # "CCL5", "ITGAE", "KLRD1"
        "2": "TCF7_Tex",
        "3": "TNFRSF9_Treg",
        "4": "Plasma_B",
        "5": "B_cell",
        "6": "TOP2A_Tpro",
        "7": "unknown",        
        "8": "Plasma_B",
        "9": "Plasma_B",
        "10": "unknown",
        "11": "Plasma_B",
    }
)
adata.obs.to_csv('lymphocytes.csv')
adata.write_h5ad('myeloid_cell.h5ad')


# S2-01, myeloid cell annotataion
adata = adata[adata.obs['CellTypeL1'] == "Myeloid_cells",:]
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata)
sce.pp.harmony_integrate(adata, 'orig.ident')
adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
sc.pp.neighbors(adata)
for res in np.arange(0.1, 2.1, 0.1):
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res
    )
sc.tl.umap(adata)
res = 1.20
adata.obs["CellTypeL2"] = adata.obs[f"leiden_res_{res:4.2f}"].map(
    {"0" : "NLRP3_Macrophages",
     "1" : "LYVE1_Macrophages",
     "2":  "NLRP3_Macrophages",
     "3": "unknown",
     "4": "LYVE1_Macrophages",
     "5": "APOE_Macrophages",
     "6": "CD1C_cDC2",
     "7": "APOE_Macrophages",
     "8": "LYVE1_Macrophages",
     "9": "APOE_Macrophages",
     "10": "CLEC9A_cDC1",
     "11": "Mast_cells",
     "12": "LAMP3_cDC3",
     "13": "unknown",
     }
)
adata.obs.to_csv('Myeloid_cells.csv')
adata.write_h5ad('lymphocytes.h5ad')

# S03 join together
adata = sc.read_h5ad('BCa_immune_L1_annotation.h5ad')
csv1 = pd.read_csv('lymphocytes.csv', index_col=0)
csv2 = pd.read_csv('myeloid_cell.csv', index_col=0)
csv = pd.concat([csv1, csv2])

adata = adata[adata.obs["CellTypeL1"] != "unknown" ]
adata.obs['CellTypeL2'] = csv['CellTypeL2']
adata = adata[adata.obs["CellTypeL2"] != "unknown" ]
marker_genes = {
    "Macrophages_APOE": ["C1QC", "APOE", "FCGR3A"],
    "B_Cells"   : ["CD79B", "MS4A1","FCRL1"],
    "CA10_Th": ["SLC4A10", "CA10", "GNG4"],
    "CCL5_NK": ["CCL5", "KLRD1", "NKG7"],
    "cDC2": ["CD1C", "FCER1A", "GPAT3"],
    "cDC1": ["CLEC9A", "IRF8", "XCR1"],
    "cDC3": ["LAMP3", "CCR7", "FSCN1"],
    "Macrophages_LYVE1": ["LYVE1", "PLTP","CD163"],
    "Mast_cells": ["CPA3", "KIT", "IL1RL1"],
    "Macrophages_NLRP3": ["NLRP3","MMP19","EMP1"],
    "Plasma_cells":["IGHA1", "IGHG1", "MZB1",],
    "TCF7_Tex": ["IL7R", "TC2N", "TCF7"],
    "TNFRSF9_Treg": ["TNFRSF9", "CTLA4", "IL2RA"],
    "TOP2A_Tpro": ["TOP2A", "MKI67", "SOX4"],
}
sc.pl.dotplot(adata,
              marker_genes,
              standard_scale="var",
              groupby="CellTypeL2",
              save='CellTypeL2_Marker.pdf')
sc.pl.umap(adata,
           color=["CellTypeL1", "CellTypeL2", "nFeature_RNA"],
           ncols=5,
           legend_loc="on data",
           legend_fontsize=4,
           save='UMAP_Immune.pdf')
adata.write_h5ad('BCa_Immune_L2_annotation.h5ad')


