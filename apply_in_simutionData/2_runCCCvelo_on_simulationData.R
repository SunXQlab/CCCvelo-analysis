import warnings
import time
import random

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="tkinter")

from models.train_CCCvelo import *
# from train_spVelo_test import *
from models.plot_CCCvelo import *
from models.preprocess_CCCvelo import *
from models.evaluation_Metric import *
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import roc_curve, precision_recall_curve

# CUDA support
if torch.cuda.is_available():
    device = torch.device('cuda')
else:
    device = torch.device('cpu')
warnings.filterwarnings("ignore")

def main(seed, slide_num,train_num):
    # Step 1: Load data
    # start_time = time.time()
    base_path = "E:/CCCvelo/apply_in_simu/"
    # results_path = os.path.join(MODEL_DIR, f"slide_{slide_num}/")
    # create_directory(results_path)

    # input_dir = os.path.join(base_path, f"newInput/fan/slide_{slide_num}/")
    input_dir = os.path.join(base_path, f"newInput/20240721_fan/slide_{slide_num}/")
    files = {
        'count_file': 'raw_expression_mtx.csv',
        'meta_file': 'pseudotime_value.csv',
        'loca_file': 'cell_location.csv',
        'LR_link_file': 'LR_links.csv',
        'TFTG_link_file': 'TFTG_links.csv',
        'LRTF_score_file': 'TFLR_score/',
        'TF_activity_file': 'TF_activity_groundtruth_mtx.csv',
        'LRTF_para_file': 'LRTF_paras.csv',
        'TFTG_para_file': 'TFTG_paras.csv'
    }

    paths = {key: os.path.join(input_dir, fname) for key, fname in files.items()}
    adata = PrerocessSimuData(**paths, using_low_emdding=False)  # n_obs × n_vars = 1976 × 14832
    iroot = np.argmin(adata.obs['groundTruth_psd'])
    print('The index of root cell is:', iroot)
    adata = root_cell(adata, select_root=int(iroot))  #
    # adata = root_cell(adata, select_root='UMAP')
    torch.save(adata, os.path.join(MODEL_DIR, "pp_adata.pt"))
    print(adata)

    # # Step 2: Train model
    data = PrepareData(adata, hidden_dims=[6, 6, 6])

    model = SpatialVelocity(*data, lr=0.1, Lambda=0.1)
    iteration_adam, loss_adam = model.train(200)
    plotLossAdam(loss_adam, MODEL_DIR)

    # Extract and save velocity data
    adata_copy = adata[:, adata.var['TGs'].astype(bool)]
    adata_velo = get_raw_velo(adata_copy, model)
    save_model_and_data(model, adata_velo, MODEL_DIR)

    # Step 3: Visualize
    # adata = torch.load(results_path + "pp_adata.pt")
    # adata_velo = torch.load(results_path + "CCCvelo.pt")
    # model = torch.load(results_path + "model_spa_velo.pth")

    # plt_dir = os.path.join(VISUALIZE_DIR, f"slide_{slide_num}/")
    # create_directory(plt_dir)

    plot_gene_dynamic_v2(adata_velo, model, VISUALIZE_DIR)

    adata_spa = adata_velo.copy()
    scv.tl.velocity_graph(adata_spa, basis='spatial', vkey='velocity', xkey='Imputate', n_jobs=1)
    scv.tl.velocity_pseudotime(adata_spa)  # root_key=adata_spa.uns['iroot']
    scv.pl.velocity_embedding_stream(adata_spa, basis='spatial', color='velocity_pseudotime', cmap='gnuplot',
                                     density=2, smooth=0.5, cutoff_perc=0,
                                     save=VISUALIZE_DIR + 'spatial_embedding_stream')

    adata_test1 = calculate_groundTruth_velo_v0(adata)
    adata_test1.layers['velocity'] = adata_test1.layers['groundTruth_velo']
    scv.tl.velocity_graph(adata_test1, basis='spatial', vkey='velocity', xkey='Imputate')
    scv.tl.velocity_pseudotime(adata_test1)  # root_key=adata_test.uns['iroot']
    scv.pl.velocity_embedding_stream(adata_test1, basis='spatial', color='velocity_pseudotime', cmap='gnuplot',
                                     density=2, smooth=0.5, cutoff_perc=0
                                     , save=VISUALIZE_DIR + 'spatial_embedding_stream_with_gd_velo')

    # Step 4: calculate metric
    adata_spa.layers['groundTruth_velo'] = adata_test1.layers['velocity']
    calculate_Spearman(adata_spa)
    calculate_mse(adata_spa)

if __name__ == "__main__":
  
    for i in range(10):

        seed = 1  # Replace with your seed value

        torch.manual_seed(seed)  
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)

        random.seed(seed)  
        np.random.seed(seed)  

        # setting global name and path
        # DATA_DIR = "./data/processed/"
        # MLNET_DIR = "./results3_m3/mlnet/"
        MODEL_DIR = f"./results_simu/slide_{i}/trained_model/"
        TG_PRED_DIR = f"./results_simu/slide_{i}/tg_prediction/"
        LOSS_DIR = f"./results_simu/slide_{i}/loss_curves/"
        VISUALIZE_DIR = f'./results_simu/slide_{i}/visualize/'

        # create folders
        for dir_path in [MODEL_DIR,TG_PRED_DIR, LOSS_DIR, VISUALIZE_DIR]:
            os.makedirs(dir_path, exist_ok=True)

        print(f'========================================Training slide{i+1}=====================================')
        main(seed, slide_num=i+1, train_num=1)

    end_time = time.time()
    run_time = (end_time - start_time) / 60
    print(f"Running time is: {run_time} mins")

