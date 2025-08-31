import pandas as pd   
import anndata as ad
import scanpy as sc
import scvelo as scv
scv.settings.verbosity = 0
import unitvelo as utv
import os 
import time
import psutil

data_path = "/home/yll/velocity_methods/dataset/scSTData/2023_Biorxiv_stereoseq_mousebrain/"
adata = ad.read_h5ad(data_path + 'mousebrain_bin60_clustered.h5ad')
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
# loc_sele = loc_sele[loc[1] > -7500]
print(loc_sele)
adata = adata_ct[loc_sele.index,]
print(adata)

adata.uns['louvain_colors'] = adata.uns['louvain_colors'].astype(str)
adata.uns['scc_anno_colors'] = adata.uns['scc_anno_colors'].astype(str)
adata.uns['scc_colors'] = adata.uns['scc_colors'].astype(str)

# dataset = '../data/BoneHuman/adata.h5ad'
dataset = './data/MouseCortex/cortex_neurous.h5ad'
label = 'scc_anno'
exp_metrics = {}

velo_config.R2_ADJUST = True 
velo_config.IROOT = None
velo_config.FIT_OPTION = '1'
velo_config.GPU = 0

adata = utv.run_model(dataset, label, config_file=velo_config)

result_path = "/home/yll/velocity_methods/01_analysis/apply_in_stereo_cortex/UniTVelo_result/"
adata.write(result_path+"unitvelo_velo.h5ad")

plt_path = result_path + "figure_new/"
if not os.path.exists(plt_path):
    os.makedirs(plt_path)

scv.tl.velocity_graph(adata, basis='spatial', vkey='velocity')
scv.pl.velocity_embedding_stream(adata, basis = 'spatial', color=adata.uns['label'], dpi=150, title='',save=plt_path + 'spaVelo_embedding_stream')
scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, basis='spatial', color='velocity_pseudotime', cmap='gnuplot',
               save=plt_path + 'spatial' + '_velocity_pseudotime')

adata.write(result_path + 'unitvelo_results.h5ad')

fig, ax = plt.subplots(figsize=(6, 5.5))  
order = ["Isocortex L6", "Isocortex L5", "Isocortex L4", "Isocortex L2/3"]
sc.pl.violin(adata_sp, keys=["velocity_pseudotime"], groupby="Cluster", order=order, ax=ax,alpha=0.9)
plt.savefig(os.path.join(plt_path, "UniTvelo_violin_plot_psd.pdf"), format='pdf')
plt.close()

