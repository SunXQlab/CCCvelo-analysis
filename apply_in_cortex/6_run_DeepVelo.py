import anndata as ann
import deepvelo as dv
import scvelo as scv
import pandas as pd
import time 
import psutil
import os

from deepvelo import Constants
from deepvelo.utils import update_dict

data_path = "/home/yll/velocity_methods/dataset/scSTData/2023_Biorxiv_stereoseq_mousebrain/"
adata = ann.read_h5ad(data_path + 'mousebrain_bin60_clustered.h5ad')
adata_copy=adata.copy()
print(adata_copy)

del adata_copy.uns['neighbors']
del adata_copy.uns['spatial_neighbors']
del adata_copy.obsp['connectivities']
del adata_copy.obsp['distances']
del adata_copy.obsp['spatial_distances']
del adata_copy.obsp['spatial_connectivities']

# 'cortical glutamatergic neuroblast',
def get_index1(lst=None, item=''):
    return [index for (index, value) in enumerate(lst) if value == item]

# filter cells
celltype_ls = adata_copy.obs['scc_anno'].to_list()
celltype_name = ['Isocortex L6', 'Isocortex L5', 'Isocortex L4', 'Isocortex L2/3']
ct_index_ls = []
for i in celltype_name:
    celltype = get_index1(celltype_ls, i)
    ct_index_ls.extend(celltype)
adata_ct = adata_copy[ct_index_ls,]

# remove outliers
loc = pd.DataFrame(adata_ct.obsm['spatial'])
loc_sele = loc[loc[0] < 5500]
adata = adata_ct[loc_sele.index,]
print(adata)

# preprocess the data
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_neighbors=30, n_pcs=30)

# specific configs to overide the default configs
configs = {
    "name": "DeepVelo", # name of the experiment
    "n_gpu": 0
}
configs = update_dict(Constants.default_configs, configs)

# run DeepVelo using the default configs
trainer = dv.train(adata, configs)

result_path = "/home/yll/velocity_methods/01_analysis/apply_in_stereo_cortex/DeepVelo_result/"
adata.write(result_path + 'deepvelo_velo.h5ad')

adata.obs['Cluster'] = adata.obs['scc_anno']
adata.uns['Cluster_colors'] = ['#d62728', '#2ca02c', '#1f77b4', '#ff7f0e']  
plt_path = data_path + "figure_new/"
if not os.path.exists(plt_path):
    os.makedirs(plt_path)

adata_sp = adata.copy()
scv.tl.velocity_graph(adata_sp, basis='spatial', vkey='velocity')
scv.pl.velocity_embedding_stream(adata_sp, basis='spatial', color='Cluster',
                                 cutoff_perc=0,
                                 dpi=150, title='',
                                 save=plt_path + 'spatial_velocity_stream')

latent_time(adata_sp)
scv.pl.scatter(adata_sp, basis='spatial', color='latent_time', color_map='gnuplot', size=80,
               save=plt_path + "spatial" + '_latent_time')

scv.tl.velocity_pseudotime(adata_sp)  
scv.pl.scatter(adata_sp, basis="spatial", color='velocity_pseudotime', cmap='gnuplot',
               save=plt_path + "spatial" + '_velocity_pseudotime')

correlation, _ = scipy.stats.spearmanr(adata_sp.obs['latent_time'], adata_sp.obs['velocity_pseudotime'])
print('the correlation between fit_t and velocity_pesudotime is:', correlation)
adata_sp.write(data_path + 'deepvelo_spatial_results.h5ad')

fig, ax = plt.subplots(figsize=(6, 5.5))  
order = ["Isocortex L6", "Isocortex L5", "Isocortex L4","Isocortex L2/3"]
sc.pl.violin(adata_sp, keys=["velocity_pseudotime"], groupby="Cluster",order=order, ax=ax,alpha=0.9)
plt.savefig(os.path.join(plt_path, "deepvelo_violin_plot_psd.pdf"), format='pdf')
plt.close()




