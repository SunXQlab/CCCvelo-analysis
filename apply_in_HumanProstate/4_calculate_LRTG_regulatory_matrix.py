import warnings
import torch
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="tkinter")

from models.plot_CCCvelo_batch import *
from models.train_CCCvelo_batchs import *
from models.plot_CCCvelo import *
from models.preprocess_CCCvelo import *
from models.evaluation_Metric import *
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import roc_curve, precision_recall_curve

if __name__ == "__main__":

    seed = 1
    train_num = 6

    base_path = "E:/CCCvelo/apply_in_prostate/"
    results_path = os.path.join(base_path,
                                f"results/trained_model/")

    save_path = results_path + "regulate_result/"
    create_directory(save_path)

    adata = torch.load(results_path + "pp_adata.pt")
    adata_velo = torch.load(results_path + "CCCvelo.pt")

    model = torch.load(results_path + "model_spa_velo.pth")

    # calculate regulatory matrix between LR and TG
    all_jac_TFTG = torch.zeros(16203,93,40)
    all_jac_LRTF = torch.zeros(16203,40,206)
    all_jac_LRTG = torch.zeros(16203,93,206)
    all_matched_indices = []
    for batch_idx, batch in enumerate(model.dataloader):
        TGs_expr_batch, TFs_expr_batch, TFLR_allscore_batch = [x.to(device) for x in batch]
    
        jac_TFTG = Jacobian_TFTG_batch(model,batch)  
        jac_LRTF = Jacobian_LRTF_batch(model,batch)  
        jac_TFTG = jac_TFTG.float()  
        jac_LRTF = jac_LRTF.float()  
        jac_LRTG = torch.bmm(jac_TFTG, jac_LRTF)  
    
        print('the shape of jac_TFTG:', jac_TFTG.shape)
        print('the shape of jac_LRTF:', jac_LRTF.shape)
        print('the shape of jac_LRTG:', jac_LRTG.shape)
    
        matched_indices = []
    
        for i, tf_batch_row in enumerate(TFs_expr_batch):
            match_found = False
            tf_batch_row = tf_batch_row.float()
            for j, tf_original_row in enumerate(model.TFs_expr):
                if torch.all(torch.eq(tf_batch_row, tf_original_row)): 
                    matched_indices.append((i, j))  
                    match_found = True
                    break
    
            if not match_found:
                print(f"Row {i} in TFs_expr_batch has no exact match in TFs_expr.")
    
            if match_found:
                all_jac_TFTG[j,:,:] = jac_TFTG[i,:,:]
                all_jac_LRTF[j,:,:] = jac_LRTF[i,:,:]
                all_jac_LRTG[j,:,:] = jac_LRTG[i,:,:]
    
        all_matched_indices.append(matched_indices)
    
    print('the shape of all_jac_TFTG:', all_jac_TFTG.shape)
    print('the shape of all_jac_LRTF:', all_jac_LRTF.shape)
    print('the shape of all_jac_LRTG:', all_jac_LRTG.shape)
    
    torch.save(all_jac_TFTG, os.path.join(results_path, "Jacobian_TFTG.pt"))
    torch.save(all_jac_LRTF, os.path.join(results_path, "Jacobian_LRTF.pt"))
    torch.save(all_jac_LRTG, os.path.join(results_path, "Jacobian_LRTG.pt"))
    torch.save(all_matched_indices, os.path.join(results_path, "all_matched_indices.pt"))

    jac_TFTG = torch.load(os.path.join(save_path, "Jacobian_TFTG.pt"))
    jac_LRTF = torch.load(os.path.join(save_path, "Jacobian_LRTF.pt"))
    jac_LRTG = torch.load(os.path.join(save_path, "Jacobian_LRTG.pt"))

    print('the shape of jac_TFTG:', jac_TFTG.shape)
    print('the shape of jac_LRTF:', jac_LRTF.shape)
    print('the shape of jac_LRTG:', jac_LRTG.shape)

    # calculate TAM-each cell type regulate matrix
    base_path = "E:/CCCvelo/apply_in_prostate/"
    LRTF_score_file = os.path.join(base_path, f"Input/TFLR_score/")
    TFTG_link_file = os.path.join(base_path,f"Input/TFTG_links.csv")

    TFTG_link = pd.read_csv(TFTG_link_file)  # TFTG link
    folder_path = LRTF_score_file
    file_names = os.listdir(folder_path)
    obs_names = adata.obs_names
    i = obs_names[0]
    obs_name = i + "_"
    index = [index for index, name in enumerate(file_names) if obs_name in name]

    file_name = file_names[index[0]]
    data_tmp = pd.read_csv(folder_path + file_name)
    LR_pairs = data_tmp.columns.tolist()
    TGs = adata_velo.var_names.tolist()
    TFs = list(np.unique(TFTG_link['TF'].values))

    celltypes = adata_velo.obs['Cluster'].unique()
    for i, ct in enumerate(celltypes):
        print('Current cell type is:', ct)
        ct_indices = np.where(adata_velo.obs['Cluster'] == ct)[0]
        # jac_LRTG_ct = jac_LRTG[ct_indices, :, :]
        # LRTG_regu_mtx_ct = torch.mean(jac_LRTG_ct, dim=0)

        jac_LRTF_ct = jac_LRTF[ct_indices, :, :]
        LRTF_regu_mtx_ct = torch.mean(jac_LRTF_ct, dim=0)

        jac_TFTG_ct = jac_TFTG[ct_indices, :, :]
        TFTG_regu_mtx_ct = torch.mean(jac_TFTG_ct, dim=0)
        print('the shape of TFTG_regu_mtx_ct is:', TFTG_regu_mtx_ct.shape)

        # save results
        # df_LRTG = pd.DataFrame(LRTG_regu_mtx_ct.numpy(), columns=LR_pairs, index=TGs)
        # ct = ct.replace(":", " ")
        # file_path = os.path.join(save_path, f'LRTG_regu_mtx_ct_{ct}.csv')
        # df_LRTG.to_csv(file_path)

        df_LRTF = pd.DataFrame(LRTF_regu_mtx_ct.numpy(), columns=LR_pairs, index=TFs)
        ct = ct.replace(":", " ")
        file_path = os.path.join(save_path, f'LRTF_regu_mtx_ct_{ct}.csv')
        df_LRTF.to_csv(file_path)

        df_TFTG = pd.DataFrame(TFTG_regu_mtx_ct.numpy(), columns=TFs, index=TGs)
        ct = ct.replace(":", " ")
        file_path = os.path.join(save_path, f'TFTG_regu_mtx_ct_{ct}.csv')
        df_TFTG.to_csv(file_path)


    jac_TFTG = torch.load(os.path.join(save_path, "Jacobian_TFTG.pt"))
    jac_LRTF = torch.load(os.path.join(save_path, "Jacobian_LRTF.pt"))
    jac_LRTG = torch.load(os.path.join(save_path, "Jacobian_LRTG.pt"))

    TFTG_regu_mtx = torch.mean(jac_TFTG, dim=0)
    LRTF_regu_mtx = torch.mean(jac_LRTF, dim=0)
    LRTG_regu_mtx = torch.mean(jac_LRTG, dim=0)

    df_TFTG_all = pd.DataFrame(TFTG_regu_mtx.numpy(), columns=TFs, index=TGs)
    file_path = os.path.join(save_path, f'TFTG_regu_mtx_all.csv')
    df_TFTG_all.to_csv(file_path)

    df_LRTF_all = pd.DataFrame(LRTF_regu_mtx.numpy(), columns=LR_pairs, index=TFs)
    file_path = os.path.join(save_path, f'LRTF_regu_mtx_all.csv')
    df_LRTF_all.to_csv(file_path)

    df_LRTG_all = pd.DataFrame(LRTG_regu_mtx.numpy(), columns=LR_pairs, index=TGs)
    file_path = os.path.join(save_path, f'LRTG_regu_mtx_all.csv')
    df_LRTG_all.to_csv(file_path)
