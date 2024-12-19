import os
import scanpy as sc
import numpy as np
import pandas as pd
import celltypist
import warnings

import matplotlib.pyplot as plt
import mplscience

# 1. Train gbmap data
adata = sc.read_h5ad('../data/gbmap/extended_gbmap.h5ad')

# 2. Prepare data : normalize & feature selection
sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)

## Feature selection (1) : highly variable genes
model_fs = celltypist.train(adata, 'cell_type', n_jobs = 70, max_iter = 5, use_SGD = True)

gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -100, axis = 1)[:, -100:]
gene_index = np.unique(gene_index)
print(f"Number of genes selected: {len(gene_index)}")

model = celltypist.train(adata[:, gene_index], 'cell_type', check_expression = False, n_jobs = 70, max_iter = 200)
model.write('../output/gbmap_ct_fs_model.pkl')

