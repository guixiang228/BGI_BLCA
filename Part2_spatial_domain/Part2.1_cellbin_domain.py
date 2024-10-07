#### Algorithm for domain partitioning
import sys
import os
from SPACEL.setting import set_environ_seed
set_environ_seed(42)
from SPACEL import Splane
import scanpy as sc
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.serif'] = ['Arial']
sc.settings.set_figure_params(dpi=100,
                              dpi_save=300,
                              facecolor='white',
                              fontsize=10,
                              vector_friendly=True,
                              figsize=(3, 3))

sc.settings.set_figure_params(
    dpi=50, dpi_save=300, facecolor="white", fontsize=10, vector_friendly=True, figsize=(4, 4)
)
sc.settings.verbosity = 3
use_gpu = "cuda:6"


args = sys.argv
n_cluster = int(args[1])
n_neighbors = int(args[2])
k = int(args[3])
gnn_dropout = float(args[4])

od = './n_cluster_{}_n_neighbors_{}_k_{}_gnn_dropout_{}'.format(n_cluster, n_neighbors, k, gnn_dropout)
os.system(f'mkdir -p {od}')
os.chdir(od)
save_path = f'{od}/Splane_models/'

# %%
import yaml

samplelist = ["HBCP1","HBCP3A","HBCP3B","HBCP3C","HBCP5A","HBCP5B","HBCP6","HBCP7","HBCP8","HBCP9","HBCP10","HBCP11","HBCP12","HBCP13","HBCP14","HBCP15","HBCP17","HBCP18"]

celltypelist = ['Bcell', 'DCS', 'Endothelial', 'Fibroblast', 'MIBC', 'Macrophage','Mastcell', 'Metastatic_tumor', 'Myocyte', 'NMIBC', 'Pericyte','Plasmocyte', 'Tcell', 'Tumor_Epithelial','Normal_Epithelial']


st_ad_list = []
for sample in samplelist:
    adata = sc.read(f'./{sample}/image_anno.h5ad')
    adata.uns['celltypes'] = celltypelist
    for celltype in celltypelist:
        if celltype not in adata.obs.columns:
            adata.obs[celltype] = 0
    adata.obs["celltype"] =  adata.obs['annotated_cluster'].str.extract(r'^([^_]+)')[0]
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.filter_cells(adata, min_genes=5)
    st_ad_list.append(adata)


# %%
splane = Splane.init_model(st_ad_list,
                           n_clusters = n_cluster,
                           use_gpu = use_gpu,
                           n_neighbors=n_neighbors,
                           k=k,
                           gnn_dropout=gnn_dropout)

# %% [markdown]
# Here, we train the model to obtain latent feature of each spots/cells. The parameter `d_l` affects the level of batch effect correction between slices. By default, `d_l` is `0.2`.

# %%
splane.train(d_l=0.1, save_path=save_path)

# %% [markdown]
# Then, we can identify the spatial domain to which each spot/cell belongs. By default, the results will be saved in `spatial_domain` column in `.obs`. If the `key` parameter is provided, the results will be saved in `.obs[key]`.

# %%
splane.identify_spatial_domain()
merged_adata = sc.AnnData.concatenate(*st_ad_list,
                                      join='outer',
                                      batch_categories = samplelist,
                                      batch_key='batch')
merged_adata.write("cellbin_merged_adata.h5ad")

## Plot spatial domain results
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

def getDefaultColors(n, type = 1):
    if type == 1:
        colors = ["#ff1a1a", "#1aff1a", "#1a1aff", "#ffff1a", "#ff1aff",
                "#ff8d1a", "#7cd5c8", "#c49a3f", "#5d8d9c", "#90353b",
                "#507d41", "#502e71", "#1B9E77", "#c5383c", "#0081d1",
                "#674c2a", "#c8b693", "#aed688", "#f6a97a", "#c6a5cc",
                "#798234", "#6b42c8", "#cf4c8b", "#666666", "#ffd900",
                "#feb308", "#cb7c77", "#68d359", "#6a7dc9", "#c9d73d"]
    elif type == 2:
        if n <= 15:
            colors = ["#87b3d4", "#d5492f", "#6bd155", "#683ec2", "#c9d754",
                  "#d04dc7", "#81d8ae", "#d34a76", "#607d3a", "#6d76cb",
                  "#ce9d3f", "#81357a", "#d3c3a4", "#3c2f5a", "#b96f49"]
        elif n <= 20:
            colors = ["#87b3d4", "#d5492f", "#6bd155", "#683ec2", "#c9d754",
                  "#d04dc7", "#81d8ae", "#d34a76", "#607d3a", "#6d76cb",
                  "#ce9d3f", "#81357a", "#d3c3a4", "#3c2f5a", "#b96f49",
                  "#4e857e", "#6e282c", "#d293c8", "#393a2a", "#997579"]
        elif n <= 30:
            colors = ["#628bac", "#ceda3f", "#7e39c9", "#72d852", "#d849cc",
                  "#5e8f37", "#5956c8", "#cfa53f", "#392766", "#c7da8b",
                  "#8d378c", "#68d9a3", "#dd3e34", "#8ed4d5", "#d84787",
                  "#498770", "#c581d3", "#d27333", "#6680cb", "#83662e",
                  "#cab7da", "#364627", "#d16263", "#2d384d", "#e0b495",
                  "#4b272a", "#919071", "#7b3860", "#843028", "#bb7d91"]
        else:
            colors = ["#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#4836be",
                  "#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
                  "#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
                  "#e08930", "#6a8bd3", "#da4f1e", "#83e6d6", "#df4341",
                  "#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
                  "#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
                  "#d0cdb8", "#421b28", "#5eae99", "#a03259", "#406024",
                  "#e598d7", "#343b20", "#bbb5d9", "#975223", "#576e8b",
                  "#d97f5e", "#253e44", "#de959b", "#417265", "#712b5b",
                  "#8c6d30", "#a56c95", "#5f3121", "#8f846e", "#8f5b5c"]
    elif type == 3:
        colors = ["#588dd5", "#c05050", "#07a2a4", "#f5994e",
                "#9a7fd1", "#59678c", "#c9ab00", "#7eb00a"]
    elif type == 4:
        colors = ["#FC8D62", "#66C2A5", "#8DA0CB", "#E78AC3",
                "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"]    
    elif type == 5:
        colors = ["#c14089", "#6f5553", "#E5C494", "#738f4c",
                "#bb6240", "#66C2A5", "#2dfd29", "#0c0fdc"]
    if n:
        if n <= len(colors):
            colors = colors[0:n]
        else:
            step = 16777200 // (n - len(colors)) - 2
            add_colors = []
            tmp = random.sample(range(step),1)[0]
            for i in range(n-len(colors)):
                hextmp = str(hex(tmp)).replace("0x",'')
                if len(hextmp) == 5:
                    hextmp = "0"+hextmp
                add_colors.append("#" + hextmp)
                tmp = tmp + step
            colors = colors + add_colors
    return colors


adata = merged_adata[merged_adata.obs['batch']=='HBCP1']
adata.obs['spatial_domain'] = adata.obs['spatial_domain'].astype('str')
adata.obs['spatial_domain'] = adata.obs['spatial_domain'].astype('category')

res = pd.DataFrame(adata.obs, columns = ["x", "y", 'spatial_domain','annotated_cluster'], index = adata.obs.index)
res.to_csv("bin1clu_celltype.txt",sep = '\t',index =False)
clusters = adata.obs['spatial_domain'].cat.categories
cluster_number = clusters.shape[0]
colors = getDefaultColors(cluster_number, type = 2)
flout = open("color_celltype.list",'w')
flout.write('0\t#FF8C00\n')
flout.write('1\t#d5492f\n')
flout.write('2\t#683ec2\n')
flout.write('3\t#ede737\n')
flout.write('4\t#6d76cb\n')
flout.write('5\t#ce9d3f\n')
flout.write('6\t#81357a\n')
flout.write('7\t#d3c3a4\n')
flout.write('8\t#3c2f5a\n')
flout.write('9\t#b96f49\n')
flout.write('10\t#6bd155\n')
flout.write('11\t#d34a76\n')
flout.write('12\t#c9d754\n')
flout.write('13\t#d04dc7\n')
flout.write('14\t#00FFFF\n')
flout.close()
os.system('cell_bin_plot bin1clu_celltype.txt %s color_celltype.list celltype.tif'%(mask))

Image.MAX_IMAGE_PIXELS = None
categories = []
colors = []
with open('color_celltype.list', 'r') as f:
    for line in f:
        category, color = line.strip().split('\t')
        categories.append(category)
        colors.append(color)
im = Image.open('celltype.tif')
draw = ImageDraw.Draw(im)
font = ImageFont.truetype("arial.ttf", 500)
text_height = font.getsize('A')[1]
text_widths = [font.getsize(category)[0] for category in categories]
max_text_width = max(text_widths)
rect_width = max_text_width + 60
rect_height = (text_height + 10) * len(categories) + 10
rect_x = im.width - rect_width - 10
rect_y = im.height - rect_height - 10
for i in range(len(categories)):
    y = rect_y + i*(text_height+10) + 10
    draw.ellipse((rect_x+10-text_height, y, rect_x, y+text_height), fill=colors[i])
    draw.text((rect_x+40, y), categories[i], font=font, fill='white')
im.save('celltype_add_legend.tif')
