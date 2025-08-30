import warnings
import scanpy as sc
import pandas as pd
from scipy.stats import spearmanr
import seaborn as sns

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="tkinter")


from models2.plot_CCCvelo import *
from models2.utils import *

VISUALIZE_DIR = './results/visualize/'

if __name__=="__main__":

    current_path = os.getcwd()
    print("Current Path:", current_path)
    result_path = "./results/trained_model/"
    
    adata_velo = torch.load(results_path + "CCCvelo.pt")
    print(adata_velo)

    adata_velo.obs["Cluster"] = adata_velo.obs["Cluster"].astype("category")
    print("Cluster categories:", adata_velo.obs["Cluster"].cat.categories)
    adata_velo.uns['Cluster_colors'] = ['#2ca02c', '#1f77b4',  '#d62728','#ff7f0e']
    
    model = torch.load(results_path + "model_spa_velo.pth")
    
    # ====================================plot spatial velocity streamline and pseudotime=================================
    adata_spa = adata_velo.copy()
    scv.tl.velocity_graph(adata_spa, basis='spatial',vkey='velocity',xkey='Imputate')
    scv.pl.velocity_embedding_stream(adata_spa,basis='spatial',vkey='velocity',smooth=0.5,cutoff_perc=0,color='Cluster',legend_loc='right margin',
                                     save=VISUALIZE_DIR + 'spatial_velocity_streamline.pdf')
    
    scv.tl.velocity_pseudotime(adata_spa)
    scv.pl.scatter(adata_spa,basis='spatial',color='velocity_pseudotime',cmap='gnuplot',
                   save=VISUALIZE_DIR + 'spatial_velocity_psd')
    
    scv.pl.scatter(adata_spa,basis='spatial',color='fit_t',cmap='gnuplot',
                   save=VISUALIZE_DIR + 'spatial_fit_t')
    
    correlation, _ = scipy.stats.spearmanr(adata_spa.obs['fit_t'], adata_spa.obs['velocity_pseudotime'])
    print('the correlation between fit_t and velocity_pesudotime is:', correlation)
    
    fig, ax = plt.subplots(figsize=(6, 5.5))  
    order = ["Isocortex L6", "Isocortex L5", "Isocortex L4","Isocortex L23"]
    sc.pl.violin(adata_spa, keys=["velocity_pseudotime"], groupby="Cluster",order=order, ax=ax,alpha=0.9)
    plt.savefig(os.path.join(VISUALIZE_DIR, "CCCvelo_violin_plot_psd.pdf"), format='pdf')
    plt.close()

    psd = adata_spa.obs['velocity_pseudotime']
    cluster = adata_spa.obs['Cluster']
    data_psd = pd.DataFrame({
        'cluster': cluster,
        'psd': psd,
        'velocity_pseudotime':adata_spa.obs['velocity_pseudotime']
    })
    palette = {
        'Isocortex L6': "#2ca02c",
        'Isocortex L5': "#1f77b4",
        'Isocortex L4': "#d62728",
        'Isocortex L23': "#ff7f0e"
    }
    fig1, ax = plt.subplots(figsize=(6.5, 6))
    order = ["Isocortex L6", "Isocortex L5", "Isocortex L4","Isocortex L23"]
    sns.boxplot(data=data_psd, x='cluster', y='velocity_pseudotime', hue='cluster', orient='v', order=order,  ax=ax, width=0.5,fliersize=1.5,
                palette=palette)
    plt.yticks(np.arange(0, 1, 0.2))
    plt.title('Velocity pseudotime')
    plt.xlabel('Cell type')
    plt.savefig(os.path.join(VISUALIZE_DIR, "CCCvelo_boxplot_psd_on_cortex.pdf"), format='pdf')
    plt.close()
  
    cccvelo_psd = adata_spa.obs[['Cluster', 'velocity_pseudotime']]
    cluster_value_map = {
        'Isocortex L23': 4,
        'Isocortex L4': 3,
        'Isocortex L5': 2,
        'Isocortex L6': 1
    }
    value_series = adata_spa.obs['Cluster'].map(cluster_value_map)
    cccvelo_gdt_ord = pd.DataFrame({'Cluster': adata_spa.obs['Cluster'], 'value': value_series})
    correlation, p_value = spearmanr(cccvelo_gdt_ord['value'], cccvelo_psd['velocity_pseudotime'])
    print(f"cccvelo: The Spearman correlation between gdt_ord and psd {correlation}:")
    
    






