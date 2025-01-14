import os
os.chdir('../../')
import warnings
import pickle
warnings.filterwarnings('ignore')
import math

import scanpy as sc
import numpy as np
import pandas as pd
import SEACells

def run(file, folder = 'tmp'):
    sample_name = file.split('.h5ad')[0]
    # 1. Load sample
    print('1. Load sample')
    adata = sc.read_h5ad(os.path.join(folder, file))

    # 2. Configure model
    print('2. Configure model')
    n_SEACells = math.ceil(adata.shape[0]/75)

    print(f'n_SEACells: {n_SEACells}')
    build_kernel_on = "X_pca"
    n_waypoint_eigs = min(10, n_SEACells - 1)
    convergence_epsilon = 1e-5
    model = SEACells.core.SEACells(adata, 
                      build_kernel_on=build_kernel_on, 
                      n_SEACells=n_SEACells, 
                      n_waypoint_eigs=n_waypoint_eigs,
                      convergence_epsilon = convergence_epsilon)

    # 3. Run SEACells
    print('3. Run SEACells')
    model.construct_kernel_matrix()
    model.initialize_archetypes()
    model.fit(min_iter=10, max_iter=100)

    # 4. Summarize by SEACell
    print('4. Summarize by SEACell')
    seacell = SEACells.core.summarize_by_SEACell(adata, SEACells_label='SEACell', celltype_label = "tumor", summarize_layer='counts')
    seacell.obs['ncells_per_SEACell'] = adata.obs.groupby('SEACell').size().loc[seacell.obs_names]

    # 5. Transfer metadata
    print('5. Transfer metadata')
    seacell.obs['author'] = adata.obs['author'].iloc[0]
    seacell.obs['donor_id'] = adata.obs['donor_id'].iloc[0]
    seacell.obs['method'] = adata.obs['method'].iloc[0]
    seacell.obs['assay'] = adata.obs['assay'].iloc[0]
    seacell.obs.index = sample_name + "-" + seacell.obs.index


    # 6. Save
    adata.write(f'tmp/seacell/{sample_name}_full.h5ad')
    seacell.write(f'tmp/seacell/{sample_name}_seacell.h5ad')
    with open(f'tmp/seacell/{sample_name}_model.pkl', 'wb') as f:
        pickle.dump(model, f)
    print('6. Done')
    print('-------------------------------------')
    return seacell
    
if __name__ == '__main__':
    i = 0
    seacells = []
    files = os.listdir('tmp')
    files = [f for f in files if f.endswith('.h5ad')]
    for file in files:
        i += 1
        print(f'running {i}/{len(files)}')
        seacell = run(file, folder = 'tmp')    
        seacells.append(seacell)

    seacells = pd.concat(seacells)
    seacells.write('tmp/seacell.h5ad')

