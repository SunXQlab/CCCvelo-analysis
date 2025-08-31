import warnings
import scanpy as sc
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="tkinter")

from models.plot_CCCvelo import *
from models.utils import *

VISUALIZE_DIR = './Result/visualize/'

if __name__=="__main__":

    current_path = os.getcwd()
    print("Current Path:", current_path)

    results_path = "./Result/trained_model/"

    adata_velo = torch.load(results_path + "CCCvelo.pt")
    print(adata_velo)

    adata_velo.obs["Cluster"] = adata_velo.obs["Cluster"].astype("category")
    print("Cluster categories:", adata_velo.obs["Cluster"].cat.categories)
    adata_velo.uns['Cluster_colors'] = ["#CE3D32E5",  "#5050FFE5","#749B58E5", "#F0E685E5"]

    model = torch.load(results_path + "model_spa_velo.pth")
    
    # ====================================plot spatial velocity streamline and pseudotime=================================
    adata_spa = adata_velo.copy()
    scv.tl.velocity_graph(adata_spa, basis='spatial',vkey='velocity',xkey='Imputate')
    scv.pl.velocity_embedding_stream(adata_spa,basis='spatial',vkey='velocity',smooth=0.5,cutoff_perc=0,color='Cluster',legend_loc='right margin',
                                     save=VISUALIZE_DIR + 'spatial_velocity_streamline.pdf')
    scv.pl.velocity_embedding(adata_spa, basis='spatial', color='Cluster', dpi=300, arrow_size=2, arrow_length=2,
                              save=VISUALIZE_DIR + 'spatial_velocity_embedding')
    
    scv.tl.velocity_pseudotime(adata_spa)
    scv.pl.scatter(adata_spa,basis='spatial',color='velocity_pseudotime',cmap='gnuplot',
                   save=VISUALIZE_DIR + 'spatial_velocity_psd')
    
    scv.pl.scatter(adata_spa,basis='spatial',color='fit_t',cmap='gnuplot',
                   save=VISUALIZE_DIR + 'spatial_fit_t')
    
    scv.pl.heatmap(adata_spa, var_names=adata_spa.var_names, sortby='fit_t', col_color='Cluster', n_convolve=150,
                   save=VISUALIZE_DIR + 'heatmap_TGs_with_fit_t_all')
    
    adata_spa.write_h5ad(os.path.join('E:/CCCvelo/Model3_result/EmbryoTrunkResult2/trained_model/adata_spa_result.h5ad'))

    correlation, _ = scipy.stats.spearmanr(adata_spa.obs['fit_t'], adata_spa.obs['velocity_pseudotime'])
    print('the correlation between fit_t and velocity_pesudotime is:', correlation)

    fig, ax = plt.subplots(figsize=(6, 5.5))  
    order = ['NMPs', 'PSM', 'Somites',  'Neural tube']
    sc.pl.violin(adata_spa, keys=["velocity_pseudotime"], groupby="Cluster",order=order, ax=ax,alpha=0.9)
    plt.savefig(os.path.join(VISUALIZE_DIR, "CCCvelo_violin_plot_psd.pdf"), format='pdf')
    plt.close()

    #====================================plot umap velocity streamline and pseudotime================================= 
    adata_umap = adata_velo.copy()
    del adata_umap.obsm['X_umap']
    adata_umap.obsm['X_umap'] = adata_cccvelo.obsm['X_umap']
    scv.tl.velocity_graph(adata_umap, basis='umap',vkey='velocity',xkey='Imputate')
    scv.pl.velocity_embedding_stream(adata_umap,basis='umap',vkey='velocity',smooth=0.5,cutoff_perc=0,color='Cluster',
                                     save=VISUALIZE_DIR + 'umap_velocity_streamline.pdf')
    
    scv.tl.velocity_pseudotime(adata_umap)
    scv.pl.scatter(adata_umap,basis='umap',color='velocity_pseudotime',cmap='gnuplot',
                   save=VISUALIZE_DIR + 'umap_velocity_psd')
    
    cccvelo_psd = adata_spa.obs[['Cluster', 'velocity_pseudotime']]
    cluster_value_map = {
        'NMPs': 1,
        'PSM': 2,
        'Somites': 3,
        'Neural tube':3
    }
    value_series = adata_spa.obs['Cluster'].map(cluster_value_map)
    cccvelo_gdt_ord = pd.DataFrame({'Cluster': adata_spa.obs['Cluster'], 'value': value_series})
    correlation, p_value = spearmanr(cccvelo_gdt_ord['value'], cccvelo_psd['velocity_pseudotime'])
    print(f"cccvelo: The Spearman correlation between gdt_ord and psd {correlation}:")
    
    psd = adata_spa.obs['velocity_pseudotime']
    cluster = adata_spa.obs['Cluster']
    data_psd = pd.DataFrame({
        'cluster': cluster,
        'psd': psd
    })

    palette = {
        'NMPs': "#CE3D32E5",
        'PSM': "#749B58E5",
        'Somites': "#F0E685E5",
        'Neural tube': "#5050FFE5"
    }
    fig1, ax = plt.subplots(figsize=(6, 5.5))
    sns.boxplot(data=data_psd, x='cluster', y='psd', hue='cluster', orient='v', order=order,  ax=ax, width=0.5,fliersize=1.5,
                palette=palette)

    plt.yticks(np.arange(0, 1, 0.2))
    plt.title('Velocity pseudotime')
    plt.xlabel('Cell type')
    plt.savefig(os.path.join(VISUALIZE_DIR, "CCCvelo_boxplot_psd.pdf"), format='pdf')
    plt.close()

    traj1 = ['NMPs','PSM', 'Somites']
    mask1 = adata_velo.obs['Cluster'].isin(traj1)
    traj1_ind = np.where(mask1)[0]
    adata_traj1 = adata_velo[traj1_ind, ]
    adata_traj1_spa = adata_velo[traj1_ind, ]
    adata_traj1_spa.uns['Cluster_colors'] = ["#CE3D32E5","#749B58E5", "#F0E685E5"]
    scv.tl.velocity_graph(adata_traj1_spa, basis='spatial', vkey='velocity', xkey='Imputate')
    scv.pl.heatmap(adata_traj1_spa, var_names=adata_traj1_spa.var_names, sortby='psd', col_color='Cluster', n_convolve=150,
                   save=VISUALIZE_DIR + 'heatmap_TGs_with_psd_traj1')

    traj2 = ['NMPs', 'Neural tube']
    mask2 = adata_velo.obs['Cluster'].isin(traj2)
    traj2_ind = np.where(mask2)[0]
    adata_traj2 = adata_velo[traj2_ind, ]
    adata_traj2_spa = adata_velo[traj2_ind, ]
    adata_traj2_spa.uns['Cluster_colors'] = ["#CE3D32E5", "#5050FFE5"]
    scv.tl.velocity_graph(adata_traj2_spa, basis='spatial', vkey='velocity', xkey='Imputate')
    scv.pl.heatmap(adata_traj2_spa, var_names=adata_traj2_spa.var_names, sortby='psd', col_color='Cluster', n_convolve=150,
                   save=VISUALIZE_DIR + 'heatmap_TGs_with_psd_traj2')






    






