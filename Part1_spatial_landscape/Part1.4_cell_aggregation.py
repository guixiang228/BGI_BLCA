#### Algorithm for cell aggregation
# %%
import os, sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import BisectingKMeans
from scipy.sparse import csr_matrix

def merge_big_cell(adata, resolution, n=30):
    adata_list = []
    for category in adata.obs[resolution].unique():
        sub_adata = adata[adata.obs[resolution] == category].copy()
        adata_list.append(sub_adata)
    merged_adatas = []
    merged_cluster_labels = pd.DataFrame(columns=["cell_id", "merged_cluster"])
    for i, st_adata in enumerate(adata_list):
        t = st_adata.shape[0] // n
        if t == 0:
            t = 1
        X = st_adata.obs[["x", "y"]].values
        spectral_clustering = BisectingKMeans(n_clusters=t, bisecting_strategy="largest_cluster", random_state=0)
        cluster_labels = spectral_clustering.fit_predict(X)
        # cluster_sizes = np.bincount(cluster_labels)
        # sorted_clusters = np.argsort(cluster_sizes)[::-1]
        # max_cluster_index = sorted_clusters[0]
        # while max_cluster_index
        cluster = st_adata.obs[resolution].iloc[0]
        adata_idx = adata.obs[resolution].unique().to_list().index(cluster)
        merged_cluster_labels_i = pd.DataFrame(
            list(zip(st_adata.obs_names, ["%s_%s" % (i - 1, adata_idx) for i in cluster_labels])),
            columns=["cell_id", "merged_cluster"],
        )
        merged_cluster_labels = pd.concat([merged_cluster_labels, merged_cluster_labels_i], axis=0)
        
        cluster_centers = []
        cluster_gene_counts = []
        cluster_cell_counts = []
        for cluster_id in np.unique(cluster_labels):
            cluster_spots = st_adata.obs.index[cluster_labels == cluster_id]
            cluster_center = np.mean(st_adata.obs.loc[cluster_spots, ["x", "y"]], axis=0)
            cluster_gene_count = np.sum(st_adata[cluster_spots].X, axis=0)
            cluster_centers.append(cluster_center)
            cluster_gene_counts.append(cluster_gene_count)
            cluster_cell_counts.append(cluster_spots.shape[0])
        
        adata_cluster = sc.AnnData(
            X=np.array(cluster_gene_counts)[:, 0, :],
            obsm={"spatial": np.array(cluster_centers)},
            obs=pd.DataFrame(
                cluster_cell_counts, index=np.arange(len(cluster_centers)), columns=["cluster_cell_counts"]
            ),
            var=adata.var,
        )
        adata_cluster.obs["merged_cluster"] = st_adata.obs[resolution].iloc[0]
        merged_adatas.append(adata_cluster)
    merged_cluster_labels.to_csv("merged_cluster_labels.txt", sep="\t", index=None)
    
    merged_adata = sc.AnnData.concatenate(*merged_adatas, join="outer")
    return merged_adata

samplelist = [
    "HBCP1", "HBCP3A", "HBCP3B", "HBCP3C","HBCP5A","HBCP5B","HBCP6","HBCP7","HBCP8","HBCP9","HBCP10","HBCP11","HBCP12","HBCP13","HBCP14","HBCP15","HBCP17","HBCP18"]

for sample in samplelist:    
    adata = sc.read_h5ad(f'{sample}.h5ad')
    adata.obs['annotated_cluster'] = adata.obs['annotated_cluster'].astype('str')
    merged_adata = merge_big_cell(adata, 'annotated_cluster', 30)
    sparse_X = csr_matrix(merged_adata.X)
    merged_adata.X = sparse_X
    merged_adata.obs['x'] = merged_adata.obsm['spatial'][:, 0]
    merged_adata.obs['y'] = merged_adata.obsm['spatial'][:, 1]
    merged_adata.write_h5ad(f'{sample}_bigcell.h5ad')



