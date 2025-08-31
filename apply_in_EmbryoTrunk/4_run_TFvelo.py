import pandas as pd    
import TFvelo as TFv
import anndata as ad
import numpy as np
import scanpy as sc
import scvelo as scv
import matplotlib
import scipy
import os, sys
import time
import psutil

def check_data_type(adata):
    for key in list(adata.var):
        if adata.var[key][0] in ['True', 'False']:
            adata.var[key] = adata.var[key].map({'True': True, 'False': False})
    return          

def data_type_tostr(adata, key):
    if key in adata.var.keys():
        if adata.var[key][0] in [True, False]:
            adata.var[key] = adata.var[key].map({True: 'True', False:'False'})
    return  

# 'cortical glutamatergic neuroblast',
def get_index1(lst=None, item=''):
    return [index for (index, value) in enumerate(lst) if value == item]

def get_pseudotime(adata):
    TFv.tl.velocity_graph(adata, basis='spatial', vkey='velocity', xkey='M_total')
    TFv.tl.velocity_pseudotime(adata, vkey='velocity', modality='M_total') 
    TFv.pl.scatter(adata, basis=basis, color='velocity_pseudotime', cmap='gnuplot', fontsize=20, save='pseudotime')
    return adata

def get_sort_positions(arr):
    positions = np.argsort(np.argsort(arr))
    positions_normed = positions/(len(arr)-1)
    return positions_normed

def get_sort_t(adata):
    t = adata.layers['fit_t_raw'].copy()
    normed_t = adata.layers['fit_t_raw'].copy()
    n_bins = 20
    n_cells, n_genes = adata.shape
    sort_t = np.zeros([n_cells, n_genes])
    non_blank_gene = np.zeros(n_genes, dtype=int)
    hist_all, bins_all = np.zeros([n_genes, n_bins]), np.zeros([n_genes, n_bins+1])
    for i in range(n_genes):
        gene_name = adata.var_names[i]
        tmp = t[:,i].copy()
        if np.isnan(tmp).sum():
            non_blank_gene[i] = 1 
            continue
        hist, bins = np.histogram(tmp, bins=n_bins)
        hist_all[i], bins_all[i] = hist, bins
        if not (0 in list(hist)):
            if (tmp.min() < 0.1) and (tmp.max() > 0.8):
                blank_start_bin_id = np.argmin(hist)
                blank_end_bin_id = blank_start_bin_id
                non_blank_gene[i] = 1
                blank_start_bin = bins[blank_start_bin_id]
                blank_end_bin = bins[blank_end_bin_id]
                tmp = (tmp < blank_start_bin)*1 + tmp 
            else:
                blank_end_bin = tmp.min()
        else:
            blank_start_bin_id = list(hist).index(0)
            for j in range(blank_start_bin_id+1, len(hist)):
                if hist[j] > 0:
                    blank_end_bin_id = j
                    break
            blank_start_bin = bins[blank_start_bin_id]
            blank_end_bin = bins[blank_end_bin_id]
            tmp = (tmp < blank_start_bin)*1 + tmp 
            
        t[:,i] = tmp
        tmp = tmp - blank_end_bin
        tmp = tmp/tmp.max()
        normed_t[:,i] = tmp
        sort_t[:,i] = get_sort_positions(tmp)

    adata.layers['latent_t'] = sort_t.copy() 
    adata.var['non_blank_gene'] = non_blank_gene.copy()
    return adata


data_path = "/home/yll/velocity_methods/01_analysis/apply_in_embryo/input_trunk/"
df_count = pd.read_csv(data_path + "raw_expression_mtx.csv") 
df_meta = pd.read_csv(data_path + "cell_meta.csv")  # meta data info
df_loca = pd.read_csv(data_path + "cell_location.csv")  # cell location

adata = sc.AnnData(X=df_count.values.astype(np.float64))  
adata.obs_names = df_count.index  
adata.var_names = df_count.columns  
adata.obs['Cluster'] = df_meta['Cluster'].values
adata.obsm['spatial'] = df_loca.values.astype(np.float64)
adata.uns['Cluster_colors'] = ["#CE3D32E5",  "#5050FFE5","#749B58E5", "#F0E685E5"]

filtered_cells = (adata.obsm["spatial"][:, 0] <= 2790) & (adata.obsm["spatial"][:, 1] >= 1800)
adata_filtered = adata[filtered_cells].copy()
sc.pl.embedding(adata_filtered, basis="spatial", color="Cluster", s=30,
                palette=["#CE3D32E5", "#5050FFE5", "#749B58E5", "#F0E685E5"], show=True)
print(adata_filtered)
adata = adata_filtered.copy()

adata.var_names_make_unique()
adata.obs_names_make_unique()
adata.uns['genes_all'] = np.array(adata.var_names)

if "spliced" in adata.layers:
    adata.layers["total"] = adata.layers["spliced"].todense() + adata.layers["unspliced"].todense()
elif "new" in adata.layers:
    adata.layers["total"] = np.array(adata.layers["total"].todense())
else:
    adata.layers["total"] = adata.X
    
adata.layers["total_raw"] = adata.layers["total"].copy()
n_cells, n_genes = adata.X.shape
sc.pp.filter_genes(adata, min_cells=int(n_cells/50))
sc.pp.filter_cells(adata, min_genes=int(n_genes/500))
TFv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000, log=True) #include the following steps
adata.X = adata.layers["total"].copy()

result_path = "/home/yll/velocity_methods/01_analysis/apply_in_trunk/TFvelo_result/"

gene_names = []
for tmp in adata.var_names:
    gene_names.append(tmp.upper())
adata.var_names = gene_names
adata.var_names_make_unique()
adata.obs_names_make_unique()

TFv.pp.moments(adata, n_pcs=30, n_neighbors=30)

TFv.pp.get_TFs(adata, databases="ChEA")
print(adata)
adata.uns['genes_pp'] = np.array(adata.var_names)
adata.write(result_path + 'pp_new.h5ad')

n_jobs = 60
n_jobs_max = np.max([int(os.cpu_count()/2), 1])
if n_jobs >= 1:
    n_jobs = np.min([n_jobs, n_jobs_max])
else:
    n_jobs = n_jobs_max
print('n_jobs:', n_jobs)
flag = TFv.tl.recover_dynamics(adata, n_jobs=n_jobs, max_iter=20, var_names="all",
    WX_method = "lsq_linear", WX_thres=20, max_n_TF=99, n_top_genes=2000,
    fit_scaling=True, use_raw=0, init_weight_method="correlation", 
    n_time_points=1000) 

adata.write(result_path + 'rc.h5ad')

if 'highly_variable_genes' in adata.var.keys():
    data_type_tostr(adata, key='highly_variable_genes')
adata.write(result_path + 'rc.h5ad')

losses = adata.varm['loss'].copy()
losses[np.isnan(losses)] = 1e6
adata.var['min_loss'] = losses.min(1)

n_cells = adata.shape[0]
expanded_scaling_y = np.expand_dims(np.array(adata.var['fit_scaling_y']),0).repeat(n_cells,axis=0)
adata.layers['velocity'] = adata.layers['velo_hat'] / expanded_scaling_y  

basis = 'spatial'
dataset_name = 'stereoseq'
if 'X_pca' not in adata.obsm.keys():
    print('PCA ing')
    sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
if (basis=='umap') and ('X_umap' not in adata.obsm.keys()):
    print('Umap ing')
    if dataset_name == 'hesc1':
        sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
        sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=30, n_pcs=5)
        sc.tl.umap(adata)
    else:
        sc.tl.umap(adata)  
        sc.pl.umap(adata, color='clusters', save=True)  

adata = get_pseudotime(adata)
print(adata)

adata_copy = adata.copy()
adata_copy = get_sort_t(adata_copy) 

adata_copy_1 = adata_copy.copy()
data_type_tostr(adata_copy_1,key = None)
print(adata_copy_1)
adata_copy_1.write(result_path + 'rc_new.h5ad')

# parameters setting
loss_percent_thres = 50
spearmanr_thres = 0.1
adata_rc = adata_copy_1.copy()

adata_copy = adata_rc.copy()

thres_loss = np.percentile(adata_copy.var['min_loss'], loss_percent_thres)
adata_copy = adata_copy[:, adata_copy.var['min_loss'] < thres_loss]

thres_n_cells = adata_copy.X.shape[0] * 0.1
adata_copy = adata_copy[:, adata_copy.var['n_cells'] > thres_n_cells]

adata_copy = adata_copy[:, adata_copy.var['non_blank_gene'] == 0]

adata_copy = get_metric_pseudotime(adata_copy)
adata_copy = adata_copy[:, adata_copy.var['spearmanr_pseudotime'] > spearmanr_thres]
print(adata_copy)

plt_path = data_path + "figure_new/"
if not os.path.exists(plt_path):
    os.makedirs(plt_path)

adata_sp = adata_copy.copy()
scv.pl.scatter(adata_sp, basis="spatial", color='velocity_pseudotime', cmap='gnuplot',
               save=plt_path + "spatial" + '_velocity_pseudotime')

scv.tl.latent_time(adata_sp)
scv.pl.scatter(adata_sp, basis='spatial', color='latent_time', color_map='gnuplot', size=80,
                save=plt_path + "spatial" + '_latent_time')

scv.tl.velocity_graph(adata_sp, basis="spatial", xkey='M_total',n_jobs=1)
scv.pl.velocity_embedding_stream(
    adata_sp,
    basis="spatial",
    color="Cluster",
    cutoff_perc=0,
    legend_fontsize=4,
    dpi=300,           # increase dpi for higher resolution
    save=plt_path + 'spatial_embedding_stream'
)


