### identify marker genes of fibroblast subtype
import sys
import os
import scanpy as sc
import pandas as pd
import numpy as np
import re
from matplotlib import pyplot as plt
import matplotlib
import anndata as ad
matplotlib.use("Agg")
from anndata import AnnData
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.serif'] = ['Arial']
sc.settings.set_figure_params(dpi=100,dpi_save=300,facecolor='white',frameon=False,fontsize=10,vector_friendly=True,figsize=(3,3))

sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
def pic(pdf):
    searchObj = re.search( r'(.*).pdf', pdf)
    png = f"{searchObj.group(1)}.png"
    plt.savefig(pdf, bbox_inches="tight")
    plt.savefig(png, bbox_inches="tight", dpi=100)
    plt.close()
def removeBiasGenes(adata):
    malat1 = adata.var_names.str.startswith('MALAT1')
    MTgenes = adata.var_names.str.startswith('MT')
    hb_genes = adata.var_names.str.contains('^HB[^(P)]')
    RPgenes = adata.var_names.str.startswith('RP') & adata.var_names.str.contains('-')
    RPgenes2 = adata.var_names.str.contains('^RP[SL]')
    CTCgenes = adata.var_names.str.startswith('CTC') & adata.var_names.str.contains('-')
    MIRgenes = adata.var_names.str.startswith('MIR')
    ACgenes = adata.var_names.str.contains('^AC[0-9]') & adata.var_names.str.contains('.')
    CTgenes = adata.var_names.str.startswith('CT') & adata.var_names.str.contains('-')
    LINCgenes = adata.var_names.str.contains('^LINC[0-9]')
    ALgenes = adata.var_names.str.contains('^AL') & adata.var_names.str.contains('.')
    remove_genes = malat1 | MTgenes | hb_genes | RPgenes | RPgenes2 | CTCgenes | MIRgenes | ACgenes | CTgenes | LINCgenes | ALgenes
    keep = np.invert(remove_genes)
    res = adata[:,keep]
    return res

#### calculate DEGs of fibroblast subtypes
adata = sc.read_h5ad('fib_type.h5ad')
group_counts = adata.obs['fib_type'].value_counts()
valid_groups = group_counts[group_counts > 1].index
adata = adata[adata.obs['fib_type'].isin(valid_groups)]

adata = removeBiasGenes(adata)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
del adata.raw
adata.obs['fib_type'] = adata.obs['fib_type'].astype('category')
sc.tl.rank_genes_groups(adata, 'fib_type', method='wilcoxon')

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
de_genes_df = pd.DataFrame()
for group in groups:
    gene_names = result['names'][group]
    group_cells = adata.obs['fib_type'] == group
   if group_cells.sum() > 0:
        expression_values = adata[group_cells, gene_names].X.mean(axis=0)
        if isinstance(expression_values, np.matrix) or len(expression_values.shape) > 1:
            expression_values = np.array(expression_values).flatten()
        print(f"Group: {group}, Expression Values Shape: {expression_values.shape}")
        group_df = pd.DataFrame({
            group + '_gene': gene_names,
            group + '_logfoldchange': result['logfoldchanges'][group],
            group + '_pvals': result['pvals'][group],
            group + '_pvals_adj': result['pvals_adj'][group],
            group + '_expression': expression_values
        })
            
        de_genes_df = pd.concat([de_genes_df, group_df],axis=1)

de_genes_df.to_csv('all_degs_fib.csv', index=False)


### plot gene dot plot
genes_of_interest = ["ITGA8", "PLAT", "NR4A1",
                     "NAV3", "KCNN3", "MATN2",
                     "FBLN1", "CD34", "C3",
                     "MYL9", "CSRP1", "SYNM",
                     "ZBTB16", "ERBB2", "FKBP5",
                     "IL6", "LIF", "CCL2",               
                     "ABL2","KAZN", "RORA",
                     "DEPTOR", "PRICKLE1", "BICC1",
                     "SPINT2", "SIPA1L3", "ACTN4",
                     "MYOCD", "CARMN", "LMOD1",
                     "FBXO32", "KIF26B", "KIAA1217",
                     "ABCC4", "FAM155A", "MCTP2",
                     "SLC14A1", "B2M", "CXCL14", 
                     "POSTN", "COL1A1", "SPARC"]
groupby_order=['ITGA8+myCAF','NAV3+CAF','FBLN1+CAF','MYL9+myCAF','ZBTB16+CAF','IL6+iCAF','ABL2+CAF','DEPTOR+CAF','SPINT2+CAF','MYOCD+myCAF','FBXO32+myCAF','ABCC4+myCAF','SLC14A1+myCAF','POSTN+myCAF'] 
adata.obs['fib_type'] = pd.Categorical(adata.obs['fib_type'], categories=groupby_order, ordered=True)
sc.pl.dotplot(
    adata,
    var_names=genes_of_interest,
    groupby='fib_type',
    dendrogram=False,
    swap_axes=False,
    standard_scale="var",
    save='DEG_fibtype.pdf',
    figsize=(12,4)
)