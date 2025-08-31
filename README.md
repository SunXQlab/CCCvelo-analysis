# Decoding dynamic cell-cell communication-driven cell state transitions in spatial transcriptomics

## Introduction
In multicellular systems, cell fate determination emerges from the integration of intracellular signaling and intercellular communication. Spatial transcriptomics (ST) provides unprecedented opportunities to elucidate these regulatory processes.

Here, we introduce CCCvelo, a computational framework designed to reconstruct CCC-driven CST dynamics by jointly optimizing a dynamic CCC signaling network and a latent CST clock. To achieve this, CCCvelo formulates a unified multiscale nonlinear kinetic model that integrates intercellular ligand-receptor signaling gradients with intracellular transcription factor activation cascades to capture gene expression dynamics encoding CSTs. Moreover, we devise PINN-CELL, a physics-informed neural network-based coevolution learning algorithm, which simultaneously optimizes model parameters and pseudotemporal ordering. Benchmarking and application of CCCvelo to synthetic data and high-resolution ST datasets, including mouse cortex, embryonic trunk development, and human prostate cancer datasets, demonstrate its ability to successfully recover known morphogenetic trajectories while uncovering dynamic CCC signaling rewiring that orchestrates CST progression. By integrating mechanistic CCC modeling with trajectory learning in spatially resolved omics data, CCCvelo provides a powerful framework for decoding the multiscale and dynamic regulatory principles that govern cell fate decisions.

## Workflow
1. **apply_in_EmbryoTrunk** contains the code to reproduce plots and detailed analysis of the Slide-seq V2 data of mouse embroy trunk<br>
   - 1_select_LRTG.R: the prepration of CCCvelo's inputs, including condiante raw expression matrix, cell meta information, cell spatial coordinates, condiante ligands, condiante receptors and feature genes.
   - 2_run_CCCvelo.py: velocity inference of CCCvelo on mouse embryo trunk dataset.
   - 3_run_CCCvelo_analysis.R: various visualizations of the velocity streamline, velocity-based pseodutime, and gene dynamtic, corresponding to Fig.5a-c and g-l.
   - 5_run_TFvelo.R: velocity inference of TFvelo on mouse embryo trunk dataset and visualizations of the velocity streamline, velocity-based pseodutime, corresponding to Fig.5d-f.
2. **apply_in_HumanProstate** contains the code to reproduce the plot and detailed analysis of the MERFISH data of human prostate.<br>
   - 0_proprecess_data.R: the proprecess of the MERFISH data of human prostate, including regional division.
   - 1_select_LRTG.R: the prepration of CCCvelo's inputs, including condiante raw expression matrix, cell meta information, cell spatial coordinates, condiante ligands, condiante receptors and feature genes.
   - 2_run_CCCvelo.py: velocity inference of CCCvelo on human prostate dataset.
   - 3_run_CCCvelo_analysis.py: various visualizations of the velocity streamline, velocity-based pseodutime, and gene dynamtic, corresponding to Fig.6a-c and g.
   - 4_calculate_LRTG_regulatory_matrix.py: the calculation of LR-TG regulation matrix.
   - 5_pl_LigandPotentialAbility.py: various visualizations of the ligand potential ability, corresponding to Fig7c and Supplementary Fig.S19.
   - 6_pl_MLnet_circro.R: visualizations of the communication between tumor microenvironment and receiver cells (E-state, I-state, and M-state tumor cells), corresponding to Fig.7a.
   - 7_cellLine_analysis.R: visualization of differential expression analysis of TGs downstream of TGFB1 within the multilayer signaling network inferred by CCCvelo, corresponding to Fig.7b.
   - 8_pl_MLnet.R: visualizations of multilayer signaling subnetworks, corresponding to Fig7d-f and Supplementary Fig.S19.
   - 9_run_TFvelo.py: velocity inference of TFvelo on human prostate dataset and visualizations of the velocity streamline, velocity-based pseodutime, corresponding to Fig.6d-f.
3. **apply_in_cortex** contains the code to reproduce the plot and detailed analysis of the Stereo-seq data of mouse cortex. <br>
   - 0_proprecess_data.R: the proprecess of the Stereo-seq data of mouse cortex.
   - 1_select_LRTG.R: the prepration of CCCvelo's inputs, including condiante raw expression matrix, cell meta information, cell spatial coordinates, condiante ligands, condiante receptors and feature genes.
   - 2_run_CCCvelo.py: the velocity inference of CCCvelo.
   - 3_run_CCCvelo_analysis.py: various visualizations of the velocity streamline, velocity-based pseodutime, and gene dynamtic, corresponding to Fig.4 and Supplementary Fig.S12.
   - 4_run_scVelo.py: velocity inference of scvelo on mouse cortex dataset and visualizations of the velocity streamline, velocity-based pseodutime, corresponding to Fig.4.
   - 5_run_TFvelo.py: velocity inference of TFvelo on mouse cortex trunk dataset and visualizations of the velocity streamline, velocity-based pseodutime, corresponding to Fig.4.
   - 6_run_DeepVelo.py: velocity inference of DeepVelo on mouse cortex trunk dataset and visualizations of the velocity streamline, velocity-based pseodutime, corresponding to Fig.4.
   - 7_run_UniTVelo.py: velocity inference of UniTVelo on mouse cortex trunk dataset and visualizations of the velocity streamline, velocity-based pseodutime, corresponding to Fig.4.
4. **apply_in_simu** contains the code to reproduce the simulation study of CCCvelo, corresponding to Fig.3. <br>
   - 0_generate_Fan_data.m: the generation of the Fan-shaped simulated datasets.
   - 0_generate_Spiral_data.m: the generation of the Spiral-shaped simulated datasets.
   - 1_proprecess_Fan_data: the prepration of CCCvelo's inputs.
   - 1_proprecess_Spiral_data: the prepration of CCCvelo's inputs.
   - 2_runCCCvelo_on_Fan_data.py: the velocity inference of CCCvelo on the Fan-shaped simulated datasets, corresponding to Fig.2b-c.
   - 2_runCCCvelo_on_Spiral_data.py: the velocity inference of CCCvelo on the Spiral-shaped simulated datasets to Fig.2e-f.

