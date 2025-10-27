import scvelo as scv
import os
import scanpy as sc
import numpy as np
import pandas as pd 
import scipy.sparse as sp
import anndata as ad
scv.logging.print_version()

print("The current work path is：", os.getcwd())
# Load annotated binning data
data_path = "/home/yll/velocity_methods/dataset/scSTData/2023_Biorxiv_stereoseq_mousebrain/"

fname_bin60 = data_path + "mousebrain_bin60_clustered.h5ad"
adata_bin60 = sc.read_h5ad(fname_bin60)

sc.pl.embedding(adata_bin60, basis="spatial", color="scc_anno",size=5)

save_path = "/home/yll/velocity_methods/01_analysis/apply_in_stereo_cortex/data/"

obs = pd.DataFrame(adata_bin60.obs)
obs.to_csv(save_path + "bin60_clustered_with_meta.csv")

loc = pd.DataFrame(adata_bin60.obsm['X_spatial'])
loc.to_csv(save_path + "bin60_clustered_with_loc.csv")

count = pd.DataFrame(adata_bin60.X.todense())
count.to_csv(save_path + "bin60_clustered_with_count.csv")

gene_id = pd.DataFrame(adata_bin60.var_names)
gene_id.to_csv(save_path + "bin60_clustered_with_gene.csv")
