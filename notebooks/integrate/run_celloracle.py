import os
os.chdir('../../')
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import scanpy as sc
import celloracle as co # biothings_client==0.2.6
import matplotlib
matplotlib.use('Agg')


def run(file, folder = 'tmp', base_GRN = None, n_jobs = 64):
    sample_name = file.split('.h5ad')[0]
    
    # 1. Load sample
    print('1. Load sample')
    adata = sc.read_h5ad(os.path.join(folder, file))
    adata.X = adata.layers['counts'].copy()
    adata = adata[:, adata.var.highly_variable].copy()

    # 2. Configure model
    print('2. Configure model')
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name="tumor", embedding_name="X_umap")
    oracle.import_TF_data(TF_info_matrix=base_GRN)
    oracle.perform_PCA()

    n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
    n_comps = min(n_comps, 50)

    n_cell = oracle.adata.shape[0]
    k = int(0.025*n_cell)
    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs= n_jobs)

    # 3. Calculate network
    print('3. Calculate network')
    links = oracle.get_links(cluster_name_for_GRN_unit="tumor", alpha=10, verbose_level=0, n_jobs = n_jobs)

    # 4. Filter link
    ### don't filter edges by strength or p-value, later consensus aggreagation will be used
    print('4. Filter link')
    links.filter_links(p=0.05, weight="coef_abs", threshold_number=len(links.links_dict['tumor']))
    links.get_network_score()
    
    df = links.filtered_links['tumor'].reset_index(drop=True)
    df['doner_id'] = adata.obs['donor_id'].iloc[0]
    df['n_cells'] = adata.shape[0]
    df['n_tumor_cells'] = adata.obs['tumor'].value_counts().loc['tumor']
    
    # 5. Save
    print('5. Save')
    df.to_csv(f'tmp/celloracle/{sample_name}_grn.csv', index=False)
    links.to_hdf5(file_path= f"tmp/celloracle/{sample_name}_link.celloracle.links")
    
    print('6. Done')
    print('-------------------------------------')
    return df


if __name__ == '__main__':
    i = 0
    dfs = []
    folder_path = 'tmp'
    n_jobs = 56
    base_GRN = pd.read_parquet('output/base_GRN_from_bulkATAC.parquet')
    
    files = os.listdir(folder_path)
    files = [f for f in files if f.endswith('.h5ad')]
    for file in files:
        i += 1
        print(f'running {i}/{len(files)}')
        df = run(file, folder=folder_path,  base_GRN = base_GRN, n_jobs = n_jobs)
        dfs.append(df)
        
    res = pd.concat(dfs, axis = 0, ignore_index=True)
    res.to_csv('tmp/celloracle/all_grn.csv', index=False)