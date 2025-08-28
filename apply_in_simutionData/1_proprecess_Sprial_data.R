library(readxl)
library(dplyr)
library(Seurat)

rm(list = ls())
gc()

setwd("/home/yll/velocity_methods/01_analysis/apply_in_simulation/")
source("/home/yll/velocity_methods/01_analysis/apply_in_simulation/R/code.R")

# ground truth 

LR_link <- data.frame(ligand = c('Lig1','Lig3','Lig3','Lig4','Lig5'), 
                      receptor = c('Rec1','Rec1','Rec2','Rec2','Rec1'))
TFTG_link <- data.frame(TF = c('TF1','TF1','TF2','TF2','TF3','TF3'), 
                        TG = c('TG1','TG2','TG1','TG2','TG3','TG4'))
groundTruth_TG <- data.frame(row.names = unlist(lapply(paste0('Lig',1:5), function(x) paste0(x,"_", paste0('Rec',1:2)))),
                             TG1 = c(1,0,0,0,1,1,0,1,1,0),
                             TG2 = c(1,0,0,0,1,1,0,1,1,0),
                             TG3 = c(0,0,0,0,0,1,0,1,0,0),
                             TG4 = c(0,0,0,0,0,1,0,1,0,0))

TGs = colnames(groundTruth_TG)
LRpairs <- lapply(1:ncol(groundTruth_TG), function(x){rownames(groundTruth_TG)})
names(LRpairs) <- TGs
groundTruth_TF <- data.frame(row.names = unlist(lapply(paste0('Lig',1:5), function(x) paste0(x,"_", paste0('Rec',1:2)))),
                             TF1 = c(1,0,0,0,1,0,0,0,1,0),
                             TF2 = c(1,0,0,0,1,1,0,1,1,0),
                             TF3 = c(0,0,0,0,0,1,0,1,0,0))
TFs = colnames(groundTruth_TF)
LRpairs <- lapply(1:ncol(groundTruth_TF), function(x){rownames(groundTruth_TF)[(groundTruth_TF[,x]!=0)]})
names(LRpairs) <- TFs

# load data

# fan_path <- paste0(getwd(),"/newData/fan_modified_initial_y0/")
fan_path <- paste0(getwd(),"/newData/sprial_modified_initial_y0/")
files <- list.files(fan_path)

for (i in files){
  
  print(i)
  
  input_ls <- prepare_input_newdata(fan_path,i)
  
  exprMat <- input_ls$exprMat
  locaMat <- input_ls$locaMat
  TF_act <- input_ls$TF_act
  psdVal <- input_ls$psdVal
  LRTF_para <- input_ls$LRTF_para
  TFTG_para <- input_ls$TFTG_para
  
  distMat <- as.matrix(dist(locaMat))
  diag(distMat) <- 1
  
  # save result for python input
  f_path <- paste0(getwd(),"/newInput/sprial_modified_initial_y0/")
  output_fpath <- paste0(f_path, i)
  if(!dir.exists(output_fpath)){
    dir.create(output_fpath)
  }
  
  wd_score <- paste0(output_fpath,"/TFLR_score/")
  dir.create(wd_score,recursive = T)
  
  for (i in 1:dim(exprMat)[2]){
    
    cellID = colnames(exprMat)[i]
    Receiver = cellID
    Sender = NULL
    LRTF_score_reci <- calculate_LRTF_score_sc(exprMat, distMat,group=NULL,
                                               LRpairs, TFs, Receiver, Sender,
                                               far.ct, close.ct, downsample)
    
    LRs_score <- LRTF_score_reci$LRs_score
    TFs_expr <- LRTF_score_reci$TFs_expr
    
    tflr_score <- matrix(data = 0, nrow = dim(groundTruth_TF)[2], ncol = dim(groundTruth_TF)[1])
    rownames(tflr_score) <- colnames(groundTruth_TF)
    colnames(tflr_score) <- rownames(groundTruth_TF)
    
    for (tf in TFs){
      
      cell_score = LRs_score[[tf]]
      ind = which(colnames(tflr_score) %in% colnames(cell_score))
      tflr_score[tf,ind] = cell_score
      
    }
    
    # write.table(tflr_score, file=paste0(wd_score, cellID, '_TFLR_score.csv'), append= T, sep=','
    #             ,row.names = FALSE,col.names = FALSE)
    write.table(tflr_score, file=paste0(wd_score, cellID, '_TFLR_score.csv'), append= T, sep=','
                ,row.names = FALSE,col.names = TRUE)
    
  }
  
  write.table(exprMat,file=paste0(output_fpath, '/raw_expression_mtx.csv'),
              sep = ",",row.names = TRUE,col.names = TRUE)
  # write.table(normMat,file=paste0(output_fpath, '/norm_expression_mtx.csv'),
  #             sep = ",",row.names = TRUE,col.names = TRUE)
  write.table(TF_act,file=paste0(output_fpath, '/TF_activity_groundtruth_mtx.csv'),
              sep = ",",row.names = TRUE,col.names = TRUE)
  write.table(locaMat,file=paste0(output_fpath, '/cell_location.csv'),
              sep = ",",row.names = TRUE,col.names = TRUE)
  write.table(psdVal,file=paste0(output_fpath, '/pseudotime_value.csv'),
              sep = ",",row.names = TRUE,col.names = TRUE)
  write.table(LR_link,file=paste0(output_fpath, '/LR_links.csv'),
              sep = ",",row.names = FALSE,col.names = TRUE)
  write.table(TFTG_link,file=paste0(output_fpath, '/TFTG_links.csv'),
              sep = ",",row.names = FALSE,col.names = TRUE)
  write.table(LRTF_para,file=paste0(output_fpath, '/LRTF_paras.csv'),
              sep = ",",row.names = FALSE,col.names = TRUE)
  write.table(TFTG_para,file=paste0(output_fpath, '/TFTG_paras.csv'),
              sep = ",",row.names = FALSE,col.names = TRUE)
  
  
}
