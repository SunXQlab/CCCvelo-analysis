library(Seurat)
library(readr)
library(ggsci)
library(scales)
library(ggplot2)

rm(list = ls())
gc()

setwd('/home/yll/velocity_methods/01_analysis/apply_in_prostate/data/area_4000x6000_5000x7000')

# color 4
show_col(pal_d3("category20")(20))
mycolors_nejm <- pal_d3("category20",alpha = 0.9)(16)
mycolor_ct <- mycolors_nejm[1:12]

names(mycolor_ct) <- c("M.state tumor","E.state tumor","ICS.state tumor","Endothelial",
                       "Myeloid","T cells",'B cells','Mast','Epithelial',"Tumor cells",
                       'Smooth muscle cells','Fibroblast')
scales::show_col(mycolor_ct)

mycolor_ct1 <- mycolor_ct[4:12]
mycolor_ct2 <- mycolor_ct[1:3]

# read data
ser <- readRDS("./ser_obj_with_newannotation.rds")

patch_size <- 2000
x_edges <- seq(-200, max(ser$center_x), patch_size)
inner_x_edges <- x_edges[1:length(x_edges)]
y_edges <- seq(-200, max(ser$center_y), patch_size)
inner_y_edges <- y_edges[1:length(y_edges)]

inner_x_edges <- c(1000,3000,5000,7000,9000)
inner_y_edges <- c(1000,3000,5000,7000,9000)

spa_loc <- data.frame(center_x = ser$center_x,center_y=ser$center_y)
spa_loc$patch_id <- character(length(rownames(spa_loc)))
for (x in inner_x_edges) {
  for (y in inner_y_edges) {
    patch_id <- paste0(as.character(x),"x",as.character(x+patch_size), "_", as.character(y),"x",as.character(y+patch_size))
    patch <- spa_loc[which((spa_loc$center_x > x) & (spa_loc$center_x < x + patch_size) & 
                              (spa_loc$center_y > y) & (spa_loc$center_y < y + patch_size)), ]
    if (length(rownames(patch)) > 0) {
      spa_loc[rownames(patch), ]$patch_id <- patch_id
    }
  }
}

Idents(ser) <- spa_loc$patch_id
patch_id <- unique(spa_loc$patch_id)
for (i in patch_id){
  
  # step1 create newfile
  save_path <- paste0(getwd(), '/area_',i,'/')
  if(!dir.exists(save_path)){
    dir.create(save_path)
  }
  
  sub_ser <- subset(ser, idents = i)
  saveRDS(sub_ser,file=paste0(save_path,'sub_ser_obj.rds'))
  
  cell_freq <- as.data.frame(table(sub_ser$new_celltype)) 
  colnames(cell_freq) <- c('celltype','number')
  write.table(cell_freq,file=paste0(save_path,'cell_freq.csv'),sep = ",",row.names = FALSE,col.names = TRUE)
  
  # 1. plot spatial distribution of all celltype
  df_plot1 <- data.frame(x = sub_ser$center_x,y=sub_ser$center_y,celltype = sub_ser$orig.annotation)
  p1 <- ggplot(df_plot1,aes(x=x,y=y,color=celltype))+
    geom_point(size=0.8)+
    scale_color_manual(values = mycolor_ct1)+
    theme_bw()+
    theme(panel.border = element_blank(),   
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_blank(), 
          axis.text = element_blank(),  
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text = element_text(size = 15,margin = margin(r=10))) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  pdf(paste0(save_path,"spatial_distribution_of_all_celltype_v2.pdf"),
      width = 10,height = 7)
  print(p1)
  dev.off()
  
  # 2. plot spatial distribution of new cell types with three state tumor cells
  df_plot2 <- data.frame(x = sub_ser$center_x,y=sub_ser$center_y,new_celltype = sub_ser$new_celltype)
  p2 <- ggplot(df_plot2,aes(x=x,y=y,color=new_celltype))+
    geom_point(size=0.8)+
    scale_color_manual(values = mycolor_ct)+
    theme_bw()+
    theme(panel.border = element_blank(),   
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_blank(),
          axis.text = element_blank(),  
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text = element_text(size = 15,margin = margin(r=10))) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  pdf(paste0(save_path,"spatial_distribution_of_all_celltype.pdf"),
      width = 10,height = 7)
  print(p2)
  dev.off()
  
  # 2. plot spatial distribution of three state tumor cells
  
  tumor_ct <- c("E.state tumor","ICS.state tumor","M.state tumor") 
  df_plot2 <- data.frame(x = sub_ser$center_x,y=sub_ser$center_y,new_celltype = sub_ser$new_celltype)
  df_plot3 <- df_plot2[which(df_plot2$new_celltype %in% tumor_ct),]
  
  p3 <- ggplot(df_plot3,aes(x=x,y=y,color=new_celltype))+
    geom_point(size=0.8)+
    scale_color_manual(values = mycolor_ct)+
    theme_bw()+
    theme(panel.border = element_blank(),   
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_blank(), 
          axis.text = element_blank(),  
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text = element_text(size = 15,margin = margin(r=10))) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  pdf(paste0(save_path,"spatial_distribution_of_tumor_celltype.pdf"),
      width = 10,height = 7)
  print(p3)
  dev.off()

}





