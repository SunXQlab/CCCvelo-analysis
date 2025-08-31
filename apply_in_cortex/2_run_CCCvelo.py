import warnings
import psutil
import random
import time

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="tkinter")

from models.runMLnet import *
from models.Input_prepare import *
from models.calculateLRscore import *
from models.plot_CCCvelo import *
from models.utils import *

# setting global name and path

DATA_DIR = "./data/processed/"
MLNET_DIR = "./results/mlnet/"
MODEL_DIR = "./results/trained_model/"
TG_PRED_DIR = "./results/tg_prediction/"
LOSS_DIR = "./results/loss_curves/"
VISUALIZE_DIR = './results/visualize/'

# create folders
for dir_path in [DATA_DIR, MODEL_DIR, MLNET_DIR, TG_PRED_DIR, LOSS_DIR, VISUALIZE_DIR]:
    os.makedirs(dir_path, exist_ok=True)

# ========== main ==========

def main(
    seed,
    rec_clusters=None,
    hidden_dims=[200, 200, 200],
    batch_size=1500,
    learning_rate=0.001,
    lambda_reg=0.01,
    n_epochs=20
):

    input_dir = DATA_DIR

    print("Loading data...")
    data_files = {
        'count_file': 'raw_expression_mtx.csv',
        'imput_file': 'imputation_expression_mtx.csv',
        'meta_file': 'cell_meta.csv',
        'loca_file': 'cell_location.csv'
    }
    paths = {key: os.path.join(input_dir, fname) for key, fname in data_files.items()}
    adata = ReadData(**paths)

    print("Loading database...")
    TGs_list = load_json(os.path.join(input_dir, "TGs_list.json"))
    Ligs_list = load_json(os.path.join(input_dir, "Ligs_list.json"))
    Recs_list = load_json(os.path.join(input_dir, "Recs_list.json"))
  
    print("Building multilayer network...")
    resMLnet = runMLnet(
        adata=adata,
        LigClus=None,
        RecClus=rec_clusters,
        OutputDir=MLNET_DIR,
        Databases=None,
        RecTF_method="Search",
        TFTG_method="Search",
        TGList=TGs_list,
        LigList=Ligs_list,
        RecList=Recs_list
    )
 
    ex_mulnetlist = {
        name: mlnet
        for receiver, sender_dict in resMLnet["mlnets"].items()
        for name, mlnet in sender_dict.items()
        if not mlnet["LigRec"].empty
    }
    # print(ex_mulnetlist.items())
  
    print("Multilayer network nodes summary:")
    print(summarize_multilayer_network(ex_mulnetlist))
   
    loop_calculate_LRTF_allscore(
        adata=adata,
        ex_mulnetlist=ex_mulnetlist,
        receiver_celltype=rec_clusters,
        diff_LigRecDB_path='E:/CCCvelo/data/Database/diff_LigRecDB.csv', 
        cont_LigRecDB_path='E:/CCCvelo/data/Database/cont_LigRecDB.csv', 
        OutputDir=MLNET_DIR
    )

    TFLR_all_score = get_TFLR_allactivity(
        mulNetList=ex_mulnetlist,
        OutputDir=MLNET_DIR
    )
    
    save_LRscore_and_MLnet(
        adata,
        mulNetList=ex_mulnetlist,
        TFLR_all_score=TFLR_all_score,
        save_path=MLNET_DIR
    )
    
    print("Selecting receiver cells...")
    celltype_ls = adata.obs['Cluster'].to_list()
    ct_index_ls = []
    for name in rec_clusters:
        ct_index_ls.extend(get_index1(celltype_ls, name))

    adata = adata[ct_index_ls, :].copy()
    x_coord = adata.obsm['spatial'][:, 0]
    selected_cells = adata.obs.index[x_coord < 5500]
    adata = adata[selected_cells].copy()
    print(adata)

    link_files = {
        'LR_link_file': 'LR_links.csv',
        'TFTG_link_file': 'TFTG_links.csv',
        'LRTF_score_file': 'TFLR_score/'
    }
    paths = {key: os.path.join(MLNET_DIR, fname) for key, fname in link_files.items()}
    print('Loading link files from:', paths)

    adata = PrepareInputData(adata, **paths)
    adata.uns['Cluster_colors'] = ['#ff7f0e', '#2ca02c', '#1f77b4',  '#d62728']

    torch.save(adata, os.path.join(MLNET_DIR, "pp_adata.pt"))

    adata = root_cell(adata, select_root='UMAP')
    print('Root cell cluster is:', adata.obs['Cluster'][adata.uns['iroot']])

    print("Training spatial velocity model...")

    n_cells = adata.n_obs
    print(f"Number of receiver cells: {n_cells}")

    if n_cells <= 10000:
        print("Training with standard SpatialVelocity (full batch)...")
        
        from models.train_CCCvelo import SpatialVelocity
        from models.plot_CCCvelo import plot_gene_dynamic

        data = PrepareData(adata, hidden_dims=hidden_dims)
        model = SpatialVelocity(*data, lr=learning_rate, Lambda=lambda_reg)
        iteration_adam, loss_adam = model.train(200)
    
    else:
        print("Training with batch SpatialVelocity (mini-batch mode)...")

        from models.train_CCCvelo_batchs import SpatialVelocity
        from models.plot_CCCvelo_batch import plot_gene_dynamic

        data = PrepareData(adata, hidden_dims=hidden_dims)
        model = SpatialVelocity(*data, lr=learning_rate, Lambda=lambda_reg, batch_size=batch_size)
        iteration_adam, loss_adam = model.train(n_epochs)
   

    adata.write_h5ad(os.path.join(MODEL_DIR, 'adata_pyinput.h5ad'))

    adata_copy = adata[:, adata.var['TGs'].astype(bool)]
    adata_velo = get_raw_velo(adata_copy, model)

    plot_gene_dynamic(adata_velo, model, VISUALIZE_DIR)
    save_model_and_data(model, adata_velo, MODEL_DIR)

    adata_umap = adata_velo.copy()
    del adata_umap.obsm['X_umap']
    adata_umap.obsm['X_umap'] = adata_velo.obsm['X_umap']

    # print('===============plot umap stream line and calculate velocity_psd without setting root cell=================')
    # plot_velocity_streamline(adata_umap, basis='umap', vkey='velocity', xkey='Imputate', root_key=False, density=2,
    #                          smooth=0.5, cutoff_perc=0, plt_path=VISUALIZE_DIR, save_name='v1')
    # scv.pl.scatter(adata_umap, basis='umap', color='velocity_pseudotime', color_map='gnuplot',
    #                save=VISUALIZE_DIR + "umap" + '_velocity_psd_v2')
    
    print('===============plot spatial stream line and calculate velocity_psd without setting root cell==============')
    adata_spa = adata_velo.copy()
    plot_velocity_streamline(adata_spa, basis='spatial', vkey='velocity', xkey='Imputate', root_key=False, density=2,
                             smooth=0.5, cutoff_perc=0, plt_path=VISUALIZE_DIR, save_name='v1')

    # scv.pl.scatter(adata_spa, basis='spatial', color='fit_t', color_map='gnuplot',save=VISUALIZE_DIR + "spatial" + '_latent_time')

    print("Pipeline finished successfully!")

if __name__ == "__main__": 

    seed = 1 # Replace with your seed value 

    torch.manual_seed(seed) 
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

    random.seed(seed)  
    np.random.seed(seed)  

    main(
        seed,
        rec_clusters=['Isocortex L6', 'Isocortex L5', 'Isocortex L4', 'Isocortex L23'],
        hidden_dims=[250, 250, 250],
        batch_size=None,
        learning_rate=0.005,
        lambda_reg=0.01,
        n_epochs=200)

