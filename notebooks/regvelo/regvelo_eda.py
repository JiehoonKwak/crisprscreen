# RevVelo env
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import torch
import regvelo
from regvelo import REGVELOVI,sanity_check,prior_GRN,abundance_test,TFscreening
# from velovi import preprocess_data
import cellrank as cr

import mplscience
import anndata
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.spatial.distance import cdist


# Inspect required data
## 1. adata
adata = regvelo.datasets.zebrafish_nc()
adata.var
adata.var.is_tf.sum()

## 2. GRN - gene x gene matrix ?
prior_net = regvelo.datasets.zebrafish_grn()
prior_net

## 3. skeleton
W = adata.uns["skeleton"].copy(); W
W = torch.tensor(np.array(W))