## spacel algorithm for spatial data annotation
import pandas as pd
import scanpy as sc
import anndata
import os
from tqdm import tqdm 
import scanpy as sc
import numpy as np
import sys

import SPACEL
from SPACEL.setting import set_environ_seed
set_environ_seed()
from SPACEL import Spoint

import matplotlib
import seaborn as sns
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.serif'] = ['Arial']
import torch
sc.settings.set_figure_params(
    dpi=50, dpi_save=300, facecolor="white", fontsize=10, vector_friendly=True, figsize=(4, 4)
)
sc.settings.verbosity = 3
use_gpu = "cuda:1"
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:128"

data_path = './spacel/'
os.chdir('data_path')
sc_ad = sc.read_h5ad('mouse_sc_celltypel1_no_gland.h5ad')

i = sys.argv[1]
print("Arguments: sample_name_is_", sys.argv[1])
st_ad = sc.read_h5ad(f'{data_path}mouse_st_{i}_bin50_1.h5ad')
spoint_model = Spoint.init_model(sc_ad, st_ad, celltype_key='CellTypeL1', sm_size=500000, use_gpu=True, n_threads = 8)
spoint_model.train(max_steps=5000, batch_size=1024)
pre = spoint_model.deconv_spatial()
st_ad = spoint_model.st_ad
st_ad.write(f'spacel_mouse_{i}_no_gland_bulk_qsub_L1.h5ad')

