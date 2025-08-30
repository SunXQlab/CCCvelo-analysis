from models.preprocess_CCCvelo import *
from models.evaluation_Metric import *

from sklearn.preprocessing import MinMaxScaler

def pl_GAM_fit(psd,y,n_splines,plt_path,name):

    gam = LinearGAM(s(0), n_splines=n_splines, lam=0.1).fit(psd, y) 
    XX = np.linspace(0, 1, 2000)
    plt.figure(figsize=(8, 6))
    # plt.plot(psd, y, 'o', label='Data', alpha=0.5)
    plt.plot(XX, gam.predict(XX), 'r-', label='GAM Fit')
    confidence_intervals = gam.confidence_intervals(XX, width=0.95)  
    plt.fill_between(XX, confidence_intervals[:, 0], confidence_intervals[:, 1], color='red', alpha=0.3, label='95% CI')
    plt.xlabel('Pseudotime')
    # plt.ylabel('%s regulate score' % key_lr)
    plt.legend()
    # plt.title(key_tg)
    plt.savefig(os.path.join(plt_path, name+"v2_norm_GAM_fit_y_vs_psd_nsplines%s.pdf" % n_splines), format='pdf')

if __name__ == "__main__":

    base_path = "E:/CCCvelo/apply_in_prostate/"
    results_path = "E:/CCCvelo/results/trained_model/"

    adata_velo = torch.load(results_path + "CCCvelo.pt")
    model = torch.load(results_path + "model_spa_velo.pth")

    N_cells, N_TFs, N_LRs = model.TFLR_allscore.shape
    N_TGs = adata_velo.shape[1]

    plt_path = results_path + "figure/GAM_fit"
    if not os.path.exists(plt_path):
        os.makedirs(plt_path)

    scaler = MinMaxScaler()
    adata_pp = torch.load(results_path + 'pp_adata.pt')
    print(adata_pp)
    adata = torch.load(results_path + 'adata_spa_analysis.pt')

    psd = adata.obs['velocity_pseudotime']

    LRTF_score_file = os.path.join(base_path,
                                   f"Input/TFLR_score/")
    folder_path = LRTF_score_file
    file_names = os.listdir(folder_path)
    obs_names = adata.obs_names
    i = obs_names[0]
    obs_name = i + "_"
    index = [index for index, name in enumerate(file_names) if obs_name in name]
    file_name = file_names[index[0]]
    data_tmp = pd.read_csv(folder_path + file_name)
    LR_pairs = data_tmp.columns
    print('the LR_pairs is:\n', LR_pairs)

    # load LR-TG regulate matrx
    save_path = results_path + "regulate_result/"
    jac_LRTG = torch.load(os.path.join(save_path, "Jacobian_LRTG.pt"))
    print('the shape of jac_LRTG is:', jac_LRTG.shape)

    plt_path = results_path + "figure/Ligand_GAM_fit"
    if not os.path.exists(plt_path):
        os.makedirs(plt_path)

    ligand = ["TGFB1","COL1A1","ANGPT1"]

    for lig in ligand:
        lig_name = lig+"_"
        key_lrs = LR_pairs[LR_pairs.str.contains(lig_name)]
        print('the length of key_lrs is:', key_lrs)

        mask = LR_pairs.isin(key_lrs)
        key_lr_ind = np.where(mask)[0]
        key_jac_LRTG = jac_LRTG[:, :, key_lr_ind]  
        print('the shape of key_jac_LRTG is:', key_jac_LRTG.shape)
        tmp1 = torch.norm(key_jac_LRTG, p=1, dim=1)
        tmp1_norm = scaler.fit_transform(tmp1.numpy())
        tmp1_norm = torch.tensor(tmp1_norm, dtype=tmp1.dtype)
        # Lig_score = scaler.fit_transform(tmp1.reshape(-1, 1))
        Lig_score = torch.sum(tmp1_norm, dim=1).numpy().flatten()
        pl_GAM_fit(psd, Lig_score, n_splines=5, plt_path=plt_path, name=lig_name)
