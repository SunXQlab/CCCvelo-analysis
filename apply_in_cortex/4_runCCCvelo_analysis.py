import warnings
import scanpy as sc
import pandas as pd
from scipy.stats import spearmanr
import seaborn as sns

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="tkinter")


from models2.plot_CCCvelo import *
from models2.utils import *

VISUALIZE_DIR = './Model3_result/MouseCortexResult/visualize/'

if __name__=="__main__":

    current_path = os.getcwd()
    print("Current Path:", current_path)
    result_path = 'E:/CCCvelo/apply_in_cortex/OtherMethod/CCCvelo/seed3_train2/'
    # results_path = "E:/CCCvelo/Model3_result/MouseCortexResult/trained_model/"

    # adata_velo = torch.load(results_path + "CCCvelo.pt")
    # print(adata_velo)

    # adata_velo.obs["Cluster"] = adata_velo.obs["Cluster"].astype("category")
    # print("Cluster categories:", adata_velo.obs["Cluster"].cat.categories)
    # adata_velo.uns['Cluster_colors'] = ['#2ca02c', '#1f77b4',  '#d62728','#ff7f0e']
    #
    # model = torch.load(results_path + "model_spa_velo.pth")
    #

    # file_path = "E:/CCCvelo/apply_in_cortex/OtherMethod/scVelo/"
    # adata_scvelo = sc.read_h5ad(file_path + "/scvelo_umap_results.h5ad")
    # spatial_spa = adata_velo.obsm['spatial']
    # spatial_scvelo = adata_scvelo.obsm['spatial']
    #
    # scvelo_coord_dict = {tuple(coord): idx for idx, coord in enumerate(spatial_spa)}
    # reordered_indices = []
    #
    # for i, coord_spa in enumerate(spatial_scvelo):
    #     coord_tuple = tuple(coord_spa)  # 转化为元组，方便查找
    #     # print('the i truple is:',coord_tuple)
    #     if coord_tuple in scvelo_coord_dict:  # 检查是否存在匹配的坐标
    #         j = scvelo_coord_dict[coord_tuple]  # 获取匹配的索引
    #         # print('the j truple is:', coord_tuple)
    #         # print(f"================Row {i} in coord_spa matches with row {j} in coord_scv.")
    #         reordered_indices.append(j)  # 保存匹配的索引
    #     else:
    #         print(f"Row {i} in coord_spa has no exact match in coord_scv.")
    #

    # print('the length of reordered_indices is:',len(reordered_indices))
    # adata_velo2 = adata_velo[reordered_indices,]
    #
    # # sc.pl.embedding(adata_velo, basis="spatial", color="Cluster", s=30, show=True)
    # # plt.savefig("mousecortex_spatial_plot.png", dpi=300, bbox_inches='tight')  # 可选格式：.pdf, .svg, .jpg 等
    # # plt.close()
    #
    # # ====================================plot spatial velocity streamline and pseudotime=================================
    # adata_spa = adata_velo2.copy()
    # scv.tl.velocity_graph(adata_spa, basis='spatial',vkey='velocity',xkey='Imputate')
    # scv.pl.velocity_embedding_stream(adata_spa,basis='spatial',vkey='velocity',smooth=0.5,cutoff_perc=0,color='Cluster',legend_loc='right margin',
    #                                  save=VISUALIZE_DIR + 'spatial_velocity_streamline_v0.pdf')
    #
    # scv.tl.velocity_pseudotime(adata_spa)
    # scv.pl.scatter(adata_spa,basis='spatial',color='velocity_pseudotime',cmap='gnuplot',
    #                save=VISUALIZE_DIR + 'spatial_velocity_psd_v0')
    #
    # scv.pl.scatter(adata_spa,basis='spatial',color='fit_t',cmap='gnuplot',
    #                save=VISUALIZE_DIR + 'spatial_fit_t_v0')
    #
    # adata_spa.write_h5ad(os.path.join('E:/CCCvelo/Model3_result/MouseCortexResult/trained_model/adata_spa_result.h5ad'))
    #
    # correlation, _ = scipy.stats.spearmanr(adata_spa.obs['fit_t'], adata_spa.obs['velocity_pseudotime'])
    # print('the correlation between fit_t and velocity_pesudotime is:', correlation)
    #
    # fig, ax = plt.subplots(figsize=(6, 5.5))  # 设定图像大小
    # order = ["Isocortex L6", "Isocortex L5", "Isocortex L4","Isocortex L23"]
    # sc.pl.violin(adata_spa, keys=["velocity_pseudotime"], groupby="Cluster",order=order, ax=ax,alpha=0.9)
    # plt.savefig(os.path.join(VISUALIZE_DIR, "CCCvelo_violin_plot_psd.pdf"), format='pdf')
    # plt.close()

    # adata_spa = sc.read_h5ad('E:/CCCvelo/Model3_result/MouseCortexResult/trained_model/adata_spa_result.h5ad')
    adata_spa = torch.load(result_path+'CCCvelo_results.pt')
    # 1.画箱型图
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
    plt.savefig(os.path.join(VISUALIZE_DIR, "CCCvelo_boxplot_psd_on_cortex_v1.pdf"), format='pdf')
    plt.close()
    exit()



    #====================================plot umap velocity streamline and pseudotime================================= 

    adata_umap = adata_velo2.copy()
    del adata_umap.obsm['X_umap']
    adata_umap.obsm['X_umap'] = adata_scvelo.obsm['X_umap']
    scv.tl.velocity_graph(adata_umap, basis='umap',vkey='velocity',xkey='Imputate')
    scv.pl.velocity_embedding_stream(adata_umap,basis='umap',vkey='velocity',smooth=0.5,cutoff_perc=0,color='Cluster',
                                     save=VISUALIZE_DIR + 'umap_velocity_streamline_v0.pdf')
    
    scv.tl.velocity_pseudotime(adata_umap)
    scv.pl.scatter(adata_umap,basis='umap',color='velocity_pseudotime',cmap='gnuplot',
                   save=VISUALIZE_DIR + 'umap_velocity_psd_v0')
    
    scv.pl.scatter(adata_umap,basis='umap',color='fit_t',cmap='gnuplot',
                   save=VISUALIZE_DIR + 'umap_fit_t_v0')
    
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
    
    






