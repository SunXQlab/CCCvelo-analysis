import scvelo as scv
import multiprocessing
import os
import torch
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scipy
import time
import psutil
import matplotlib
# matplotlib.use('AGG')

if __name__=="__main__":

    process = psutil.Process(os.getpid())
    before_memory = process.memory_info().rss / 1024 ** 2  
    start_time = time.time()

    data_path = "F:/STdataset_yll/2023_Biorxiv_stereoseq_mousebrain"
    adata = sc.read_h5ad(data_path + '/MouseBrain_bin60_clustered.h5ad')
    print(adata)

    adata_cpoy = adata.copy()
    del adata_cpoy.uns['neighbors']
    del adata_cpoy.uns['spatial_neighbors']
    del adata_cpoy.obsp['connectivities']
    del adata_cpoy.obsp['distances']
    del adata_cpoy.obsp['spatial_distances']
    del adata_cpoy.obsp['spatial_connectivities']
    
    # filter cells
    celltype_ls = adata_cpoy.obs['scc_anno'].to_list()
    celltype_name = ['Isocortex L6', 'Isocortex L5', 'Isocortex L4', 'Isocortex L2/3']
    ct_index_ls = []
    for i in celltype_name:
        celltype = get_index1(celltype_ls, i)
        ct_index_ls.extend(celltype)
    adata_ct = adata_cpoy[ct_index_ls,]
    
    # remove outliers
    loc = pd.DataFrame(adata_ct.obsm['spatial'])
    loc_sele = loc[loc[0] < 5500]
    adata = adata_ct[loc_sele.index,]
    
    # sc.pl.embedding(adata, basis="spatial", color="scc_anno", size=8)
    # scv.pl.proportions(adata)
    # proprecess
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    print(adata.uns['neighbors']['indices'])
    # run dynmaical model
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode='dynamical')
    print(adata)

    scv.tl.velocity_graph(adata, basis='spatial',vkey='velocity',xkey='Imputate')
    scv.pl.velocity_embedding_stream(adata,basis='spatial',vkey='velocity',smooth=0.5,cutoff_perc=0,color='Cluster',legend_loc='right margin')
    scv.tl.velocity_pseudotime(adata)
    scv.pl.scatter(adata,basis='spatial',color='velocity_pseudotime',cmap='gnuplot')

    # creat new file to save results
    file_path = "D:/AA-luluyan-phd/code/03_velocity_analysis/apply_in_cortex/OtherMethod/"
    results_path = file_path + "scVelo/"
    if not os.path.exists(results_path):
        os.makedirs(results_path)
    adata.write(results_path+ "/scvelo_results_new.h5ad")
    after_memory = process.memory_info().rss / 1024 ** 2  
    print(f"The memory usage is: {after_memory - before_memory} MB")
    
    end_time = time.time()
    run_time = (end_time - start_time) / 60
    print(f"Running time is: {run_time} mins")

  
    




