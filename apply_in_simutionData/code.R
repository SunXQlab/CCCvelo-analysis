## load simulation data
prepare_input_newdata <- function(data_path,slideID)
{
  
  expr_path <- paste0(data_path,slideID,'/ExpressionMatrix.xlsx')
  spa_path <- paste0(data_path,slideID,'/SpatialCoordinates.xlsx')
  LRTF_para_path <- paste0(data_path,slideID,'/LR-TF_parameters.xlsx')
  TFTG_para_path <- paste0(data_path,slideID,'/TF-TG_parameters.xlsx')
  # GWnosie_path <- paste0(data_path,slideID,'/GWnoise.xlsx')
  
  # read the current sheet data
  exprMat <- read_excel(expr_path,col_names = FALSE) %>% as.data.frame(.)
  locaMat <- read_excel(spa_path, col_names = FALSE) %>% as.data.frame(.)
  LRTF_para <- read_excel(LRTF_para_path,col_names = FALSE) %>% as.data.frame(.)
  TFTG_para <- read_excel(TFTG_para_path, col_names = FALSE) %>% as.data.frame(.)
  # GWnosie <- read_excel(GWnosie_path, col_names = FALSE) %>% as.data.frame(.)
  
  # expression matrix and spatial coordinate
  TF_act <- exprMat[15:17,]
  psdVal <- exprMat[18,]
  exprMat <- exprMat[1:14,]
  
  colnames(exprMat) <- paste0('Cell_',1:ncol(exprMat))
  rownames(exprMat) <- c(paste0('Lig',1:5),paste0('Rec',1:2),paste0('TF',1:3),paste0('TG',1:4))
  
  rownames(locaMat) <- colnames(exprMat)
  colnames(locaMat) <- c('dim_x','dim_y')
  
  colnames(TF_act) <- colnames(exprMat)
  rownames(TF_act) <- c(paste0('TF',1:3))
  
  colnames(psdVal) <- colnames(exprMat)
  psdVal <- rbind(colnames(exprMat),psdVal) %>% t() %>% as.data.frame()
  colnames(psdVal) <- c('Barcode','pseudotime')
  rownames(psdVal) <- NULL
  
  # colnames(GWnosie) <- colnames(exprMat)
  # rownames(GWnosie) <- c(paste0('TF',1:3),paste0('TG',1:4))
  
  locaMat <- locaMat[!duplicated(locaMat),]
  exprMat <- exprMat[colnames(exprMat) %in% rownames(locaMat)]
  TF_act <- TF_act[colnames(TF_act) %in% rownames(locaMat)]
  psdVal <-  psdVal[psdVal$Barcode %in% rownames(locaMat),]
  # GWnosie <- GWnosie[colnames(exprMat) %in% rownames(locaMat)]
  
  res = list(exprMat=exprMat, locaMat=locaMat, TF_act=TF_act, psdVal=psdVal,
             LRTF_para=LRTF_para,TFTG_para=TFTG_para)
  
  # res = list(exprMat=exprMat, locaMat=locaMat, TF_act=TF_act, psdVal=psdVal,
  #            LRTF_para=LRTF_para,TFTG_para=TFTG_para,GWnosie=GWnosie)
  return(res)
  
}

## load simulation data
prepare_input_data <- function(expr_path,spa_path,sheetID)
{
  
  # read the current sheet data
  exprMat <- read_excel(expr_path, sheet = sheetID,col_names = FALSE) %>% as.data.frame(.)
  locaMat <- read_excel(spa_path, sheet = sheetID, col_names = FALSE) %>% as.data.frame(.)
  
  TF_act <- exprMat[15:17,]
  psdVal <- exprMat[18,]
  exprMat <- exprMat[1:14,]
  
  colnames(exprMat) <- paste0('Cell_',1:ncol(exprMat))
  rownames(exprMat) <- c(paste0('Lig',1:5),paste0('Rec',1:2),paste0('TF',1:3),paste0('TG',1:4))
  
  rownames(locaMat) <- colnames(exprMat)
  colnames(locaMat) <- c('dim_x','dim_y')
  
  colnames(TF_act) <- colnames(exprMat)
  rownames(TF_act) <- c(paste0('TF',1:3))
  
  colnames(psdVal) <- colnames(exprMat)
  psdVal <- rbind(colnames(exprMat),psdVal) %>% t() %>% as.data.frame()
  colnames(psdVal) <- c('Barcode','pseudotime')
  rownames(psdVal) <- NULL
  
  locaMat <- locaMat[!duplicated(locaMat),]
  exprMat <- exprMat[colnames(exprMat) %in% rownames(locaMat)]
  TF_act <- TF_act[colnames(TF_act) %in% rownames(locaMat)]
  psdVal <-  psdVal[psdVal$Barcode %in% rownames(locaMat),]
  
  res = list(exprMat=exprMat, locaMat=locaMat, TF_act=TF_act, psdVal=psdVal)
  return(res)
  
}

# The distance weights are the reciprocal function 
# Receiver is a single cell ID
# The distance weights are the reciprocal function
calculate_LRTF_score_sc <- function(exprMat, distMat, group=NULL,
                                 LRpairs, TFs, Receiver, Sender = NULL,
                                 far.ct = 0.75, close.ct = 0.25,
                                 downsample = FALSE)
{
  
  receBars = Receiver
  
  # sendBars = colnames(exprMat)[colnames(exprMat) != receBars]
  sendBars = colnames(exprMat)
  
  Receptors = lapply(LRpairs, function(lr){stringr::str_split(lr,"_",simplify = T)[,2]})
  Ligands = lapply(LRpairs, function(lr){stringr::str_split(lr,"_",simplify = T)[,1]})
  
  # get exprMat of Ligand
  LigMats = lapply(TFs, function(tf){   #TG：筛选出的LRpairs对应的Target gene
    # print(tg)
    ligands = Ligands[[tf]]
    if(length(ligands)==1){
      lig_count = exprMat[ligands, sendBars]
      lig_count = matrix(lig_count,nrow = 1)
    }else{
      lig_count = exprMat[ligands, sendBars] %>% as.matrix()
    }
    rownames(lig_count) = LRpairs[[tf]]
    colnames(lig_count) = sendBars
    lig_count
  })
  names(LigMats) = TFs
  
  # get exprMat of Receptor
  RecMats = lapply(TFs, function(tf){
    receptors = Receptors[[tf]]
    if(length(receptors)==1){
      rec_count = exprMat[receptors, receBars]
      rec_count = matrix(rec_count,nrow = 1)
    }else{
      rec_count = exprMat[receptors, receBars] %>% as.matrix()
    }
    rownames(rec_count) = LRpairs[[tf]]
    colnames(rec_count) = receBars
    rec_count
  })
  names(RecMats) = TFs
  
  distMat = distMat[sendBars, receBars]
  
  if(!is.null(group)){
    cpMat <- get_cell_pairs(group, distMat, far.ct, close.ct)
  }else{
    cpMat <- NULL
  }
  
  distMat = 1/distMat
  
  # t1 <- Sys.time(); message(paste0('Start at: ',as.character(t1)))
  LRs_score = lapply(TFs, function(tf){
    
    # print(tg)
    LigMat = LigMats[[tf]]
    RecMat = RecMats[[tf]]
    lr = LRpairs[[tf]]
    
    # 若是 data.frame / 带因子列
    LigMat  <- data.matrix(LigMat)   # 等价于 as.matrix + as.numeric，得到 double 矩阵
    # 若 distMat 是 stats::dist 对象
    if (inherits(distMat, "dist")) {
      distMat <- as.matrix(distMat)
    } else {
      distMat <- data.matrix(distMat)
    }
    storage.mode(LigMat)  <- "double"
    storage.mode(distMat) <- "double"
    
    
    if(is.null(cpMat)){
      
      LR_score = RecMat*(LigMat%*%distMat)
      LR_score = t(LR_score) #Receptor cells * LR pairs
      colnames(LR_score) = lr
      rownames(LR_score) = receBars
      
    }else{
      
      LR_score = lapply(unique(cpMat$Receiver), function(j){
        # j = unique(cpMat$Receiver)[1]
        is <- cpMat$Sender[cpMat$Receiver == j] %>% unique()
        if(length(is)==1){
          RecMat[,j]*(LigMat[,is]*distMat[is,j])
        }else{
          RecMat[,j]*(LigMat[,is]%*%distMat[is,j])
        }
      }) %>% do.call('cbind',.) %>% t()
      colnames(LR_score) = lr
      rownames(LR_score) = unique(cpMat$Receiver)
      
    }
    LR_score
    
  })
  names(LRs_score) = TFs
  # t2 <- Sys.time(); message(paste0('End at: ',as.character(t2)))
  # t2-t1
  
  if(is.null(cpMat)){
    
    TFs_expr = exprMat[,receBars]
    TFs_expr = lapply(TFs, function(tf){ exprMat[tf, receBars] })
    
  }else{
    
    TFs_expr = lapply(TFs, function(tf){ exprMat[tf, unique(cpMat$Receiver)] })
    
  }
  names(TFs_expr) = TFs
  
  # downsample
  if(length(receBars)>500){
    if(isTRUE(downsample)){
      set.seed(2021)
      
      if(is.null(cpMat)){
        keep_cell = sample(receBars, size = 500, replace = F)
      }else{
        keep_cell = sample(unique(cpMat$Receiver), size = 500, replace = F)
      }
      
      LRs_score = lapply(LRs_score, function(LR_score){ LR_score = LR_score[keep_cell,] })
      TFs_expr = lapply(TFs_expr, function(TF_count){ TF_count = TF_count[keep_cell] })
    }
  }
  
  LRTF_score = list(LRs_score = LRs_score, TFs_expr = TFs_expr)
  
  return(LRTF_score)
  
}


