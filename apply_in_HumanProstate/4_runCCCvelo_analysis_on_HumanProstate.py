import warnings
import scanpy as sc
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="tkinter")


from models2.plot_CCCvelo import *
from models2.utils import *

VISUALIZE_DIR = './Result/visualize/'

if __name__=="__main__":

    current_path = os.getcwd()
    print("Current Path:", current_path)

    results_path = "./Result/trained_model/"

    adata_velo = torch.load(results_path + "CCCvelo.pt")
    print(adata_velo)
    
    adata_velo.obs["Cluster"] = adata_velo.obs["Cluster"].astype("category")
    print("Cluster categories:", adata_velo.obs["Cluster"].cat.categories)
    adata_velo.uns['Cluster_colors'] = ["#FF7F0EE5", "#2CA02CE5","#1F77B4E5"]
    
    model = torch.load(results_path + "model_spa_velo.pth")
    
     # ====================================plot spatial velocity streamline and pseudotime=================================
    adata_spa = adata_velo.copy()
    scv.tl.velocity_graph(adata_spa, basis='spatial',vkey='velocity',xkey='Imputate')
    scv.pl.velocity_embedding_stream(adata_spa,basis='spatial',vkey='velocity',smooth=0.5,cutoff_perc=0,color='Cluster',legend_loc='right margin',
                                     save=VISUALIZE_DIR + 'spatial_velocity_streamline_v0.pdf')
    scv.pl.velocity_embedding(adata_spa, basis='spatial', color='Cluster', dpi=300, arrow_size=2, arrow_length=2,
                              save=VISUALIZE_DIR + 'spatial_velocity_embedding')
    
    scv.tl.velocity_pseudotime(adata_spa)
    scv.pl.scatter(adata_spa,basis='spatial',color='velocity_pseudotime',cmap='gnuplot',
                   save=VISUALIZE_DIR + 'spatial_velocity_psd')
    
    scv.pl.scatter(adata_spa,basis='spatial',color='fit_t',cmap='gnuplot',
                   save=VISUALIZE_DIR + 'spatial_fit_t')
    
    adata_spa.write_h5ad(os.path.join('E:/CCCvelo/Model3_result/HumanProstateResult/trained_model/adata_spa_result.h5ad'))
    # adata_spa = sc.read_h5ad(os.path.join('E:/CCCvelo/Model3_result/HumanProstateResult/trained_model/adata_spa_result.h5ad'))
    print(adata_spa)
    correlation, _ = scipy.stats.spearmanr(adata_spa.obs['fit_t'], adata_spa.obs['velocity_pseudotime'])
    print('the correlation between fit_t and velocity_pesudotime is:', correlation)

    psd = adata_spa.obs['velocity_pseudotime']
    cluster = adata_spa.obs['Cluster']
    data_psd = pd.DataFrame({
        'cluster': cluster,
        'psd': psd
    })
    fit_t = adata_spa.obs['fit_t']
    corr, _ = scipy.stats.spearmanr(adata_spa.obs['fit_t'], adata_spa.obs['velocity_pseudotime'])
    print('the correlation between fit_t and velocity_pesudotime is:', corr)

    cccvelo_psd = adata_spa.obs[['Cluster','fit_t', 'velocity_pseudotime']]

    cluster_value_map = {
        'E.state tumor': 1,
        'ICS.state tumor': 2,
        'M.state tumor': 3
    }
    value_series = adata_spa.obs['Cluster'].map(cluster_value_map)
    gdt_ord = pd.DataFrame({'Cluster': adata_spa.obs['Cluster'], 'value': value_series})
    cccvelo_gdt_ord = pd.DataFrame({'Cluster': adata_spa.obs['Cluster'], 'value': value_series})
    correlation, p_value = spearmanr(cccvelo_gdt_ord['value'], cccvelo_psd['velocity_pseudotime'])
    correlation, p_value = spearmanr(cccvelo_gdt_ord['value'], cccvelo_psd['velocity_pseudotime'])
    print(f"cccvelo: The Spearman correlation between gdt_ord and psd {correlation}")


    fig1, ax = plt.subplots(figsize=(6.5, 6))
    sns.boxplot(data=data_psd, x='cluster', y='psd', hue='cluster', orient='v', ax=ax, width=0.5,fliersize=1.5,
                palette=["#FF7F0EE5", "#2CA02CE5", "#1F77B4E5"])  #  ["#DAA0B0", "#908899","#9D5A38"]
    
    plt.yticks(np.arange(0, 1, 0.2))
    plt.title('Velocity pseudotime')
    plt.xlabel('Cell type')
    plt.savefig(os.path.join(VISUALIZE_DIR, "boxplot_psd_v0.pdf"), format='pdf')
    plt.close()

  
