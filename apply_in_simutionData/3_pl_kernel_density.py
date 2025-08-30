import matplotlib.pyplot as plt
import os
import seaborn as sns
import torch
import numpy as np
import pandas as pd
from models.train_CCCvelo import *
from models.plot_CCCvelo import *
from models.preprocess_CCCvelo import *
from models.evaluation_Metric import *
from scipy import stats

def f(x,intercept,trend):
    y=trend*x+intercept
    return y

def calculate_SpearmanMtx(X,Y):
    coef = []
    for i in range(X.shape[0]):
        # X_i = X[:,i]
        # Y_i = Y[:,i]
        X_i = X[i,:]
        Y_i = Y[i,:]
        velo_spear = scipy.stats.spearmanr(X_i, Y_i)[0]
        coef.append(velo_spear)
    correlation_mean = np.mean(coef)
    return correlation_mean

if __name__ == "__main__":

    for i in range(10):

        slide_num = i+1

        base_path = "E:/CCCvelo/apply_in_simu/"
        results_path = os.path.join(base_path, f"results_simu/slide_{slide_num}/")

        adata = torch.load(results_path + "pp_adata.pt")
        adata_gdt = calculate_groundTruth_velo(adata)
        print(adata_gdt)
        adata_velo = torch.load(results_path + "CCCvelo.pt")
        model = torch.load(results_path + "model_spa_velo.pth")

        adata_spa = adata_velo.copy()
        scv.tl.velocity_graph(adata_spa, basis='spatial', vkey='velocity', xkey='Imputate', n_jobs=1)
        scv.tl.velocity_pseudotime(adata_spa)
        # torch.save(adata_spa, os.path.join(results_path, "adata_spa_analysis.pt"))

        psd = adata_spa.obs['velocity_pseudotime']
        gdt_t = adata_spa.obs['groundTruth_psd']
        mse_psd = np.mean((psd - gdt_t) ** 2)
        corr_t, p_val_t = scipy.stats.spearmanr(adata_spa.obs['groundTruth_psd'], adata_spa.obs['velocity_pseudotime'])
        print('the correlation between fit_t and velocity_pesudotime is:', corr_t)
        print('the p_value between fit_t and velocity_pesudotime is:', p_val_t)
        print('the MSE between groundtruth time and velocity_pesudotime is:', mse_psd)

        # calculate the mse between inferred velocity and groundtruth velocity
        gdt_velo = adata_gdt.layers['groundTruth_velo']
        pred_velo = adata_spa.layers['velocity']
        mse_velo = np.mean((gdt_velo - pred_velo) ** 2)
        corr_velo = calculate_SpearmanMtx(adata_gdt.layers['groundTruth_velo'], adata_spa.layers['velocity'])
        print('the correlation between groundtruth velocity and preficted velocit is:', corr_velo)
        # print('the p_value between groundtruth velocity and preficted velocit is:', p_val_velo)
        print('the MSE between groundtruth velocity and preficted velocit is:', mse_velo)

        plt_path = results_path + "figure/"
        if not os.path.exists(plt_path):
            os.makedirs(plt_path)

        gdt_velo_1d = gdt_velo.flatten()
        pred_velo_1d = pred_velo.flatten()

        fig1, ax = plt.subplots(figsize=(6.5, 6))
        sns.kdeplot(x=gdt_t, y=psd, fill=True, ax=ax, cmap="Blues", thresh=0.01, levels=5)  
        plt.plot([-1.1, 1.1], [-1.1, 1.1], 'r--')  
        plt.xlabel("ground truth time")  
        plt.ylabel("velocity pseudotime")
        plt.text(-1.1, 1.2, f'Corr: {corr_t:.2f}', fontsize=12) 
        plt.xlim([-1.1, 1.1])  # plt.xlim([0, 1.2])
        plt.ylim([-1.1, 1.1])  # plt.ylim([0, 1.2])
        plt.savefig(os.path.join(plt_path, "2D_kerneldensity_psd_vs_gdt_t.pdf"), format='pdf')
        plt.close()

        fig2, ax = plt.subplots(figsize=(6.5, 6))
        sns.kdeplot(x=gdt_velo_1d, y=pred_velo_1d, fill=True, ax=ax, cmap="Blues", thresh=0.01, levels=5) 
        plt.plot([-1.1, 1.1], [-1.1, 1.1], 'r--')  
        plt.xlabel("ground truth velocity")
        plt.ylabel("predicted velocity")
        plt.text(-1.1, 1.2, f'Corr: {corr_velo:.2f}', fontsize=12)  
        plt.xlim([-1.1, 1.1])
        plt.ylim([-1.1, 1.1])
        plt.savefig(os.path.join(plt_path, "2D_kerneldensity_velo_vs_gdt_velo.pdf"), format='pdf')
        plt.close()























