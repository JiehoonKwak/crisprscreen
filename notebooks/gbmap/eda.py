import os
os.chdir('../')
import warnings
warnings.filterwarnings('ignore')
import scanpy as sc
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)

# Subset data to look for
adata = sc.read('data/gbmap/extended_gbmap.h5ad', cache = True)
adata

adata.obs['author'].value_counts()
adata.obs['assay'].value_counts()

bdata = adata[adata.obs['assay'] == "10x 3' v3"].copy()
bdata.write('output/extended_gbmap_filtered.h5ad')


# EDA
bdata.obs.method.value_counts() # cell vs nuclei
bdata.obs.stage.value_counts() # primary
bdata.obs