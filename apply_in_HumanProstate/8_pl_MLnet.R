library(dplyr)
library(ggsci)
library(ggplot2)
library(igraph)
library(plotrix)
library(ggraph)
library(ggalluvial)
library(ggforce)

rm(list=ls())
gc()

setwd('/home/yll/velocity_methods/01_analysis/apply_in_prostate/area_4000x6000_5000x7000_input/')
source('/home/yll/velocity_methods/01_analysis/code/code.R')

###########
## color ##
###########
# color
show_col(pal_d3("category20")(20))
mycolors_nejm <- pal_d3("category20",alpha = 0.9)(16)
mycolor_ct <- mycolors_nejm[1:12]

names(mycolor_ct) <- c("M.state tumor","E.state tumor","ICS.state tumor","Endothelial",
                       "Myeloid","T cells",'B cells','Mast','Epithelial',"Tumor cells",
                       'Smooth muscle cells','Fibroblast')
scales::show_col(mycolor_ct)

mycolor_ct1 <- mycolor_ct[4:12]
mycolor_ct2 <- mycolor_ct[1:3]

# nodekey

scales::show_col(pal_locuszoom(palette = "default", alpha = 0.9)(7))
mycolors_locus <- pal_locuszoom(palette = "default", alpha = 0.9)(7)

nodekey <- c("Ligand","Receptor","TF","Target")
mycolor_key <- mycolors_locus[1:4]
names(mycolor_key) <- nodekey
scales::show_col(mycolor_key)

#############
## workdir ##
#############
plotdir = './visualize/'
dir.create(plotdir,recursive = T)

wd <- "./runscMLnet/"
files = list.files(wd)
files_cps = files[grep(".csv",files,invert = T)]
files_cps = files_cps[grep("_Isocortex L4|Isocortex L23|Isocortex L5|Isocortex L6",files_cps)]

mulNetList = lapply(files_cps, function(files_cp){
  
  readRDS(paste0(wd,files_cp,"/scMLnet.rds"))
  
})
names(mulNetList) = files_cps
mulNetList = mulNetList[!unlist(lapply(mulNetList, function(mulnet){nrow(mulnet$LigRec)==0}))]

mulNet_tab = lapply(mulNetList, function(mlnet){
  
  ligrec = data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target) 
  rectf = data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
  tftg = data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)
  
  res = ligrec %>% merge(., rectf, by = 'Receptor') %>% 
    merge(., tftg, by = 'TF') %>% 
    dplyr::select(Ligand, Receptor, TF, Target) %>% 
    arrange(Ligand, Receptor)
  
}) %>% do.call('rbind',.) %>% as.data.frame()
mulNet_tab = rownames_to_column(mulNet_tab,var='CellPair')
mulNet_tab = mulNet_tab %>% separate(CellPair, into = c("Sender", "Receiver"), sep = "_")
mulNet_tab$Receiver <- gsub("\\.[0-9]+", "",mulNet_tab$Receiver)
mulNet_tab$LRpair <- paste0(mulNet_tab$Ligand,"_",mulNet_tab$Receptor)

plotdir = './visualize/'
dir.create(plotdir,recursive = T)

wd <- "./runscMLnet/"
files = list.files(wd)
files_cps = files[grep(".csv",files,invert = T)]

receivers <- c('E.state tumor','M.state tumor','ICS.state tumor')
for (receiver in receivers){
  
  files_Etumor = files_cps[grep(paste0("_.*",receiver),files_cps)]
  mulNetAllList = lapply(files_Etumor, function(file_tumor){
    
    readRDS(paste0(wd,file_tumor,"/scMLnet.rds"))
    
  })
  names(mulNetAllList) = files_Etumor
  mulNetAllList = mulNetAllList[!unlist(lapply(mulNetAllList, function(mulnet){nrow(mulnet$LigRec)==0}))]
  
  mulNet_tab = lapply(mulNetAllList, function(mlnet){
    
    ligrec = data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target) 
    rectf = data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
    tftg = data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)
    
    res = ligrec %>% merge(., rectf, by = 'Receptor') %>% 
      merge(., tftg, by = 'TF') %>% 
      dplyr::select(Ligand, Receptor, TF, Target) %>% 
      arrange(Ligand, Receptor)
    
  }) %>% do.call('rbind',.) %>% as.data.frame()
  mulNet_tab = rownames_to_column(mulNet_tab,var='CellPair')
  mulNet_tab = mulNet_tab %>% separate(CellPair, into = c("Sender", "Receiver"), sep = "_")
  mulNet_tab$Receiver <- gsub("\\.[0-9]+", "",mulNet_tab$Receiver)
  mulNet_tab$LRpair <- paste0(mulNet_tab$Ligand,"_",mulNet_tab$Receptor)
  # saveRDS(mulNet_tab,file='./TAM_to_Mstatetumor_MLnet.rds')
  
  senders <- lapply(files_Etumor, function(cp){stringr::str_split(cp,"_",simplify = T)[,1]}) %>% do.call('rbind',.) 
  receiver <- lapply(files_Etumor, function(cp){stringr::str_split(cp,"_",simplify = T)[,2]}) %>% do.call('rbind',.) %>% unique(.)
  
  MLnet_merge <- list(
    LigRec = as.data.frame(matrix(ncol = 2,dimnames = list(c(),c('source','target')))),
    RecTF =  as.data.frame(matrix(ncol = 2,dimnames = list(c(),c('source','target')))),
    TFTar =  as.data.frame(matrix(ncol = 2,dimnames = list(c(),c('source','target'))))
  )
  for (sender in senders) {
    
    MLnet <- readRDS(paste0('./runscMLnet/',sender,'_',receiver,'/scMLnet.rds'))
    MLnet_merge$LigRec <- rbind(MLnet_merge$LigRec,MLnet$LigRec) %>% as.data.frame() %>% na.omit()
    MLnet_merge$RecTF <- rbind(MLnet_merge$RecTF,MLnet$RecTF) %>% as.data.frame() %>% na.omit()
    MLnet_merge$TFTar <- rbind(MLnet_merge$TFTar,MLnet$TFTar) %>% as.data.frame() %>% na.omit()
    
  }
  
  # LRTG_im
  LRTG_im <- read_csv(paste0("regulate_result/LRTG_regu_mtx_ct_",receiver,".csv"))
  colnames(LRTG_im)[1] <- "TGs"
  
  LRTG_im_long <- LRTG_im %>%
    pivot_longer(
      cols = 2:207,          
      names_to = "LRpair",       
      values_to = "score"       
    )
  
  LRTG_im_long <- tidyr::separate(LRTG_im_long, LRpair, c('Ligand', 'Receptor'), sep = '_')
  
  #filter cancer target
  LRTG_im_new <- LRTG_im_long[which(LRTG_im_long$TGs %in% unique(mulNet_tab$Target)),]
  # LRTG_im_new$score <- abs(LRTG_im_new$score)
  LRTG_im_new <- LRTG_im_new[order(-LRTG_im_new$score),]
  colnames(LRTG_im_new) <- c('Target','Ligand','Receptor','im_norm')
  LRTG_im_new <- LRTG_im_new[LRTG_im_new$im_norm!=0,]
  
  # LRTG_im_new$score <- (LRTG_im_new$score - min(LRTG_im_new$score))/(max(LRTG_im_new$score)-min(LRTG_im_new$score))
  
  ############################
  ## MLnetPlot: directedCCI ##
  ############################

  Key <- c('TGFB1')  
  Type <- 'Ligand'
  mlnet = MLnet_merge
  lrtg_im = LRTG_im_new
  MLnet_key <- prepareMLnetworkPlotData(mlnet=MLnet_merge,lrtg_im=lrtg_im,Key=Key,
                                        Type=Type,do.check = T)
  length(unique(MLnet_key[["TFTar"]][["target"]]))
  
  MLnet_key$LigRec$target <- paste0(MLnet_key$LigRec$target,'.r')
  MLnet_key$RecTF$source <- paste0(MLnet_key$RecTF$source,'.r')
  MLnet_key$RecTF$target <- paste0(MLnet_key$RecTF$target,'.tf')
  MLnet_key$TFTar$source <- paste0(MLnet_key$TFTar$source,'.tf')
  MLnet_key$TFTar$target <- paste0(MLnet_key$TFTar$target,'.tg')
  
  colodb = pal_locuszoom(palette = "default", alpha = 0.8)(4)
  names(colodb) <- nodekey
  scales::show_col(colodb)
  
  ## plot
  
  downstream <- 'Target'
  gtitle <- paste0('Lignd_tgfb1_Sender-',receiver,'_v1')
  plotdir <- "./visualize/"
  drawMLnetworkPlot(mlnet=MLnet_key,colodb=colodb,downstream = downstream,
                    gtitle=gtitle,wd=plotdir,p_height = 6,p_width = 14)

}



















