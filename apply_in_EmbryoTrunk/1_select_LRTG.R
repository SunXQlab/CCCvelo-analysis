library(Seurat)
library(tidyverse)
library(dbscan)
library(limma)
library(SeuratWrappers)
library(jsonlite)

rm(list = ls())
gc()

setwd('/home/yll/velocity_methods/01_analysis/apply_in_trunk/')

source('/home/yll/velocity_methods/01_analysis/apply_in_stereo_cortex/R/preprocess_code.R')
source('/home/yll/velocity_methods/01_analysis/apply_in_stereo_cortex/R/create_multilayer_network.R')

# load data
data_path <- '/home/yll/velocity_methods/01_analysis/apply_in_trunk/data/'
ser_obj <- readRDS(paste0(data_path,'E8.5_rep1_201104_12_precessed.rds'))
ser_obj$Cluster <- ser_obj@meta.data$cell_state
ser_obj$Cluster <- gsub("Presomitic mesoderm \\(PSM\\)", 'PSM', ser_obj$Cluster)
ser_obj$Cluster <- gsub('Neuromesodermal Progenitors ', 'NMPs',ser_obj$Cluster)
ser_obj$Cluster <- gsub('Secondary heart field (SHF)', 'SHF',ser_obj$Cluster)
Idents(ser_obj) <- ser_obj@meta.data$Cluster

# filter
ct_num <- table(ser_obj$Cluster)
keep_ct <- names(ct_num[ct_num >= 20])
ser_obj <- ser_obj[,which(ser_obj$Cluster %in% keep_ct)]

ser_obj <- SCTransform(ser_obj, assay = "Spatial", verbose = FALSE)

## imputation
seed <- 4321
norm.matrix <- as.matrix(GetAssayData(ser_obj, "data", "SCT"))
exprMat.Impute <- run_Imputation(exprMat = norm.matrix,use.seed = T,seed = seed)

sub_anno <- data.frame(Barcode=colnames(ser_obj),Cluster=ser_obj$Cluster)
sub_loca <- data.frame(x=ser_obj$x,y=ser_obj$y)

## calculate distant
DistMat <- as.matrix(dist(sub_loca))
colnames(DistMat) <- colnames(ser_obj)
rownames(DistMat) <- colnames(ser_obj)
diag(DistMat) <- 1

# load prior databse
Databases <- readRDS('/home/yll/velocity_methods/01_analysis/prior_knowledge/Databases.rds')
quan.cutoff <- 0.98
Databases <- Databases
Databases$RecTF.DB <- Databases$RecTF.DB %>%
  .[.$score > quantile(.$score, quan.cutoff),] %>%
  dplyr::distinct(source, target)
Databases$LigRec.DB <- Databases$LigRec.DB %>%
  dplyr::distinct(source, target) %>%
  dplyr::filter(target %in% Databases$RecTF.DB$source)
Databases$TFTG.DB <- Databases$TFTG.DB %>%
  dplyr::distinct(source, target) %>%
  dplyr::filter(source %in% Databases$RecTF.DB$target)


cts <- unique(ser_obj$Cluster)
save_path <- '/mnt/Newdisk/yll/CCCvelo/apply_in_trunk/input'

# step1 create mulitlayer network
output_fpath <- paste0(save_path, '/MLnet_para/')
if(!dir.exists(output_fpath)){
  dir.create(output_fpath)
}

LRTG_list <- select_LRTG(ser_obj, Databases, log.gc = 0.4, p_val_adj=0.05,
                         pct.ct=0.05, expr.ct = 0.05)

TGs_list <- LRTG_list[["TGs_list"]]
Ligs_expr_list <- LRTG_list[["Ligs_expr_list"]]
Recs_expr_list <- LRTG_list[["Recs_expr_list"]]

## save results

output_fpath <- paste0(getwd(), '/data/processed/')

write_json(TGs_list, path=paste0(output_fpath,"TGs_list.json"), pretty = TRUE, auto_unbox = TRUE)
write_json(Ligs_expr_list, path=paste0(output_fpath,"Ligs_list.json"), pretty = TRUE, auto_unbox = TRUE)
write_json(Recs_expr_list, path=paste0(output_fpath,"Recs_list.json"), pretty = TRUE, auto_unbox = TRUE)
write_json(Databases, path=paste0(output_fpath,"Databases.json"), pretty = TRUE, auto_unbox = TRUE)

# save(TGs_list, Ligs_expr_list, Recs_expr_list, file=paste0(output_fpath,"conditate_signal_g",i,".rda"))

df_count <- as.matrix(GetAssayData(ser_obj, "data", "Spatial"))
rownames(exprMat.Impute) = rownames(df_count)
df_count = df_count[rownames(exprMat.Impute),]
df_count = t(df_count)
exprMat.Impute <- as.matrix(exprMat.Impute)
exprMat.Impute = t(exprMat.Impute)

write.table(df_count,file=paste0(output_fpath, 'raw_expression_mtx.csv'),sep = ",",row.names = TRUE,col.names = TRUE)
write.table(exprMat.Impute,file=paste0(output_fpath, 'imputation_expression_mtx.csv'),sep = ",",row.names = TRUE,col.names = TRUE)
write.table(sub_anno,file=paste0(output_fpath, 'cell_meta.csv'),sep = ",",row.names = FALSE,col.names = TRUE)
write.table(sub_loca,file=paste0(output_fpath, 'cell_location.csv'),sep = ",",row.names = TRUE,col.names = TRUE)



