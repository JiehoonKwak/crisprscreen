import os
os.chdir('../../')
import warnings
import pickle
import traceback
warnings.filterwarnings('ignore')
import math

import scanpy as sc
import numpy as np
import pandas as pd
import SEACells

def run(file, folder='tmp'):
    sample_name = file.split('.h5ad')[0]
    try:
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
                          convergence_epsilon=convergence_epsilon)

        # 3. Run SEACells
        print('3. Run SEACells')
        model.construct_kernel_matrix()
        model.initialize_archetypes()
        model.fit(min_iter=10, max_iter=100)

        # 4. Summarize by SEACell
        print('4. Summarize by SEACell')
        seacell = SEACells.core.summarize_by_SEACell(adata, SEACells_label='SEACell', celltype_label="tumor", summarize_layer='counts')
        seacell.obs['ncells_per_SEACell'] = adata.obs.groupby('SEACell').size().loc[seacell.obs_names]

        # 5. Transfer metadata
        print('5. Transfer metadata')
        seacell.obs['author'] = adata.obs['author'].iloc[0]
        seacell.obs['donor_id'] = adata.obs['donor_id'].iloc[0]
        seacell.obs['method'] = adata.obs['method'].iloc[0]
        seacell.obs['assay'] = adata.obs['assay'].iloc[0]
        seacell.obs.index = sample_name + "-" + seacell.obs.index

        # 6. Save
        print('6. Save')
        adata.write(f'tmp/seacell/{sample_name}_full.h5ad')
        seacell.write(f'tmp/seacell/{sample_name}_seacell.h5ad')
        with open(f'tmp/seacell/{sample_name}_model.pkl', 'wb') as f:
            pickle.dump(model, f)
        print('Done')
        print('-------------------------------------')
        return seacell

    except Exception as e:
        error_msg = traceback.format_exc()
        print(f'Error processing {sample_name}: {str(e)}')
        with open('tmp/seacell/error_samples.txt', 'a') as f:
            f.write(f'{sample_name}\t{str(e)}\t{error_msg}\n')
        return None

if __name__ == '__main__':
    # Create seacell directory if it doesn't exist
    os.makedirs('tmp/seacell', exist_ok=True)
    
    # Create or clear error log file
    with open('tmp/seacell/error_samples.txt', 'w') as f:
        f.write('sample_name\terror\tfull_traceback\n')
    
    # Get list of completed samples
    completed_samples = set(f.replace('_seacell.h5ad', '.h5ad') 
                          for f in os.listdir('tmp/seacell') 
                          if f.endswith('_seacell.h5ad'))
    
    # Load existing combined results if available
    combined_file = 'tmp/seacell.h5ad'
    seacells = []
    if os.path.exists(combined_file):
        print("Loading existing combined results...")
        existing_data = sc.read_h5ad(combined_file)
        seacells.append(existing_data)
        # Get sample names from existing data to skip
        completed_samples.update(
            set(name.split('-')[0] + '.h5ad' 
                for name in existing_data.obs.index)
        )
    
    # Filter files to process
    files = os.listdir('tmp')
    files = [f for f in files if f.endswith('.h5ad') and f not in completed_samples]
    
    # Process remaining files
    for i, file in enumerate(files, start=1):
        print(f'Running {i}/{len(files)} ({file})')
        seacell = run(file, folder='tmp')
        if seacell is not None:
            seacells.append(seacell)
    
    # Combine and save results only if we have new data
    if seacells:
        print("Combining all results...")
        combined_seacells = sc.concat(seacells, join='outer')
        combined_seacells.write(combined_file)
        print(f"Results saved to {combined_file}")
    else:
        print("No new data to process")