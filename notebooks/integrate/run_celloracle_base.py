import os
os.chdir('../../')
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import scanpy as sc
import celloracle as co
import matplotlib
matplotlib.use('Agg')

def run(file, folder='tmp', base_GRN=None, n_jobs=64):
    sample_name = file.split('.h5ad')[0]
    
    # 1. Load sample
    print('1. Load sample')
    adata = sc.read_h5ad(os.path.join(folder, file))
    
    # Check if tumor cells exist
    if 'tumor' not in adata.obs['tumor'].unique():
        print(f'Skipping {file}: No tumor cells found')
        with open('tmp/celloracle_base/skipped_samples.txt', 'a') as f:
            f.write(f'{sample_name}\n')
        return None
        
    adata.X = adata.layers['counts'].copy()
    adata = adata[:, adata.var.highly_variable].copy()

    # Rest of your original run function remains the same
    print('2. Configure model')
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="tumor", embedding_name="X_umap")
    oracle.import_TF_data(TF_info_matrix=base_GRN)
    oracle.perform_PCA()

    n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
    n_comps = min(n_comps, 50)

    n_cell = oracle.adata.shape[0]
    k = int(0.025*n_cell)
    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs=n_jobs)

    print('3. Calculate network')
    links = oracle.get_links(cluster_name_for_GRN_unit="tumor", alpha=10, verbose_level=0, n_jobs=n_jobs)

    print('4. Filter link')
    links.filter_links(p=0.05, weight="coef_abs", threshold_number=len(links.links_dict['tumor']))
    links.get_network_score()
    
    df = links.filtered_links['tumor'].reset_index(drop=True)
    df['doner_id'] = adata.obs['donor_id'].iloc[0]
    df['n_cells'] = adata.shape[0]
    df['n_tumor_cells'] = adata.obs['tumor'].value_counts().loc['tumor']
    
    print('5. Save')
    df.to_csv(f'tmp/celloracle_base/{sample_name}_grn.csv', index=False)
    links.to_hdf5(file_path=f"tmp/celloracle_base/{sample_name}_link.celloracle.links")
    
    print('6. Done')
    print('-------------------------------------')
    return df

if __name__ == '__main__':
    # Create or clear the skipped samples file
    with open('tmp/celloracle_base/skipped_samples.txt', 'w') as f:
        f.write('sample_name\treason\n')
    
    folder_path = 'tmp'
    n_jobs = 66
    
    base_GRN = co.data.load_human_promoter_base_GRN()
    
    # Get list of completed samples
    completed_samples = set(f.replace('_grn.csv', '.h5ad') 
                          for f in os.listdir('tmp/celloracle_base') 
                          if f.endswith('_grn.csv'))
    
    # Filter files to process
    files = os.listdir(folder_path)
    files = [f for f in files if f.endswith('.h5ad') and f not in completed_samples]
    
    # Start from where we left off
    dfs = []
    
    # Load existing results if available
    if os.path.exists('tmp/celloracle_base/all_grn.csv'):
        existing_results = pd.read_csv('tmp/celloracle_base/all_grn.csv')
        dfs.append(existing_results)
    
    for i, file in enumerate(files, start=1):
        print(f'Running {i}/{len(files)} ({file})')
        df = run(file, folder=folder_path, base_GRN=base_GRN, n_jobs=n_jobs)
        if df is not None:  # Only append if tumor cells were found
            dfs.append(df)
    
    res = pd.concat(dfs, axis=0, ignore_index=True)
    res.to_csv('tmp/celloracle_base/all_grn.csv', index=False)