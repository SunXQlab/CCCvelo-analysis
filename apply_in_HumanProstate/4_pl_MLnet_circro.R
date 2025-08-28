library(CellChat)
library(ggplot2)
library(ggsci)
library(tidyverse)

rm(list = ls())
gc()

setwd('/home/yll/velocity_methods/01_analysis/apply_in_prostate/area_4000x6000_5000x7000_input_v2')

# color
scales::show_col(pal_igv(palette = "default", alpha = 0.9)(20))
mycolors_nejm <- pal_igv(palette = "default", alpha = 0.9)(20)

mycolor_ct <- mycolors_nejm[1:12]
names(mycolor_ct) <- c("Smooth muscle cells","Tumor cells","Fibroblast","Endothelial",
                       "Myeloid","T cells",'B cells','Mast','Epithelial',"E.state tumor",
                       'M.state tumor','ICS.state tumor')
scales::show_col(mycolor_ct)

mycolor_ct1 <- mycolor_ct[1:9]
mycolor_ct2 <- mycolor_ct[10:12]

#############
## workdir ##
#############

plotdir = './MLnet_para1_46/visualize/'
dir.create(plotdir,recursive = T)

#################
## NetworkPlot ##
#################

wd <- "./MLnet_para1_46/runscMLnet/"
files = list.files(wd)
files_cps = files[grep(".csv",files,invert = T)]

# files_Etumor = files_cps[grep("_.*E.state tumor",files_cps)]
mulNetAllList = lapply(files_cps, function(file_tumor){
  
  readRDS(paste0(wd,file_tumor,"/scMLnet.rds"))
  
})
names(mulNetAllList) = files_cps
mulNetAllList = mulNetAllList[!unlist(lapply(mulNetAllList, function(mulnet){nrow(mulnet$LigRec)==0}))]

mulNet_tab = lapply(mulNetAllList, function(mlnet){
  
  ligrec = data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target) 
  rectf = data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
  tftg = data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)
  
  res = ligrec %>% merge(., rectf, by = 'Receptor') %>% 
    merge(., tftg, by = 'TF') %>% 
    dplyr::select(Ligand, Receptor, TF, Target) %>% 
    arrange(Ligand, Receptor)
  
  res$LRpair = paste0(res$Ligand,"_",res$Receptor)
  c(length(unique(res$LRpair)),length(unique(res$Target)))
  
}) %>% do.call('rbind',.) %>% as.data.frame()

df_cellpair <- names(mulNetAllList) %>% strsplit(.,"_") %>% do.call('rbind',.) %>% as.data.frame()

mulNet_tab <- cbind(df_cellpair,mulNet_tab)
mulNet_tab <- na.omit(mulNet_tab)
colnames(mulNet_tab) <- c('cell_from','cell_to','n_LRs','n_TGs')

cts_less = unique(mulNet_tab$cell_from)[which(!unique(mulNet_tab$cell_from) %in% unique(mulNet_tab$cell_to))]
a = data.frame(cell_from = cts_less,
               cell_to = cts_less,
               n_LRs = rep(0,length(cts_less)),
               n_TGs = rep(0,length(cts_less)))

mulNet_tab <- rbind(mulNet_tab,a)
mulNet_tab$n_LRs <- as.numeric(mulNet_tab$n_LRs)
mulNet_tab$n_TGs <- as.numeric(mulNet_tab$n_TGs)

for (key in colnames(mulNet_tab)[3:4]) {
  
  tmeTab <- mulNet_tab[,c('cell_from','cell_to',key)] %>% spread(cell_to, key) %>% 
    column_to_rownames(.,var = "cell_from")
  tmeTab[is.na(tmeTab)] <- 0 
  tmeTab <- as.matrix(tmeTab)
  
  colordb <- mycolor_ct[rownames(tmeTab)]
  
  pdf(paste0(plotdir,"cellchat_networkPlot_AllCTs_",key,".pdf"),height =6,width = 6)
  netVisual_circle(tmeTab, color.use = colordb,vertex.weight = rowSums(tmeTab),alpha.edge = 0.6, 
                   weight.scale = T, label.edge= F, title.name = "Number of interactions",
                   arrow.width = 1,arrow.size = 0.3,
                   text.x = 15,text.y = 1.5)
  dev.off()
  
}
















