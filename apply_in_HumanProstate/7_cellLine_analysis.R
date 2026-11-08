library(GEOquery)
library(limma)
library(ggpubr)
library(ggthemes)
library(igraph)
library(plotrix)
library(ggraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ROCR)

rm(list = ls())
gc()

setwd('/home/yll/velocity_methods/01_analysis/apply_in_prostate/cellLine/')
source('/home/yll/velocity_methods/01_analysis/code/code.R')
source('/home/yll/cell_cell_interaction/apply_in_stGBM/R/function.R')

# data path
files = list.files(path = '/home/yll/velocity_methods/01_analysis/apply_in_prostate/cellLine/data/GSE230009_RAW/')
expr_mat <- NULL 
for (file in files){
  
  data <- read.delim(
    paste0("./data/GSE230009_RAW/",file),
    sep = "\t",
    header = TRUE,
    comment.char = "#"
  )
  dim(data)
  
  sample_id <- sub("_.*", "", file)

  gene_counts <- data[, c(1, 7)]
  print(dim(gene_counts))
  colnames(gene_counts)[2] <- sample_id
  rownames(gene_counts) <- gene_counts$Geneid
  gene_counts$Geneid <- NULL
  
  if (is.null(expr_mat)) {
    expr_mat <- gene_counts
  } else {
    expr_mat <- cbind(expr_mat, gene_counts)
  }
}

ids <- data.frame(ID = rownames(expr_mat))

ENTREZID <- bitr(ids$ID, 
                 fromType = "ENSEMBL",
                 toType = c("SYMBOL"),
                 OrgDb = org.Hs.eg.db)

ids <- merge(ids, ENTREZID, by.x = "ID", by.y = "ENSEMBL")

expr <- expr_mat[rownames(expr_mat) %in% ids$ID, ]
ids <- ids[match(rownames(expr),ids$ID),]

tmp <- by(expr, ids$SYMBOL, function(x) rownames(x)[which.max(rowMeans(x))])
probes <- as.character(tmp)

expr <- expr[rownames(expr) %in% probes, ]
ids <- ids[ids$ID %in% probes, ]
rownames(expr) <- ids$SYMBOL
saveRDS(expr,"./GSE230009_expr_mat.rds")

count <- readRDS("~/velocity_methods/01_analysis/apply_in_prostate/cellLine/GSE230009_expr_mat.rds")

do.log2 = TRUE
if(do.log2){
  expr = cpm(count, log=TRUE)
}else{
  expr = count$counts
}

# norm

do.norm = FALSE
if(do.norm){
  count <- calcNormFactors(count)
  expr = cpm(count, log=TRUE)
  boxplot(expr, col = group$color, las=2, main="")
}else{
  if(do.log2){
    expr = cpm(count, log=TRUE)
  }else{
    expr = count$counts
  }
}

# filter

cutoff <- 0
drop.gene <- rowSums(expr<=cutoff)==ncol(expr)
table(drop.gene)
expr <- expr[!drop.gene,]

# filter2

genes_0_sc <- which(! rowSums(as.matrix(expr) == 0) == ncol(expr))
expr <- expr[genes_0_sc, ] 
keep_sc <- which(Matrix::rowSums( expr > 0) >= round(0.05 * ncol( expr)))
expr <-  expr[keep_sc,] 

sample_info <- data.frame(geo_accession = colnames(expr), 
                          condition=c('inhibitor',"inhibitor","inhibitor",
                                      "inhibitor","inhibitor","control",
                                      "control","control","control","control"))

compare_group <- c("inhibitor", "control")
sample_info_TGFB1 <- sample_info[sample_info$condition %in% compare_group, ] 

colnames(expr) <- sample_info$geo_accession
expr <- expr[,sample_info_TGFB1$geo_accession]

design <- model.matrix(~0 + sample_info_TGFB1$condition)
colnames(design) <- c("inhibitor","control")
fit <- lmFit(expr, design)

contrast_matrix <- makeContrasts(inhibitor_vs_control = inhibitor - control, levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

diff_expr_results <- topTable(fit2, adjust.method = "fdr", number = Inf)
saveRDS(diff_expr_results,file = "./GSE230009_limma_diff_expr_results.rds")

head(diff_expr_results)

sig_genes <- subset(diff_expr_results, P.Value < 0.05)
head(sig_genes)

# load MLnet result
mlnet_path <- '/home/yll/velocity_methods/01_analysis/apply_in_prostate/area_4000x6000_5000x7000_input_v2/MLnet_para1_46/'
wd <- paste0(mlnet_path, '/runscMLnet/')
files = list.files(wd)
files_tumor = files[grep(".csv",files,invert = T)]
mulNetList = lapply(files_tumor, function(file_tumor){
  
  readRDS(paste0(wd,file_tumor,"/scMLnet.rds"))
  
})
names(mulNetList) = files_tumor
mulNetList_tumor <- mulNetList[names(mulNetList)[grep("_.*tumor",names(mulNetList))]]
mulNet_tab <- lapply(mulNetList_tumor, function(mlnet){
  
  ligrec <- data.frame(Ligand = mlnet$LigRec$source, Receptor = mlnet$LigRec$target)
  rectf <- data.frame(Receptor = mlnet$RecTF$source, TF = mlnet$RecTF$target)
  tftg <- data.frame(TF = mlnet$TFTar$source, Target = mlnet$TFTar$target)
  
  res <- ligrec %>% merge(., rectf, by = 'Receptor') %>% 
    merge(., tftg, by = 'TF') %>% 
    dplyr::select(Ligand, Receptor, TF, Target) %>% 
    arrange(Ligand, Receptor)
})
mulNet_tab <- do.call("rbind", mulNet_tab)

lig_TGFB1_mlnet <- mulNet_tab[mulNet_tab$Ligand=='TGFB1',]
TGs_TGFB1 <- unique(lig_TGFB1_mlnet$Target)
DEG_TGs <- intersect(rownames(sig_genes),TGs_TGFB1) 

df_plot <- diff_expr_results[which(rownames(diff_expr_results)%in% TGs_TGFB1),]
df_plot$logFDR <- -log10(df_plot$P.Value)
df_plot$change <- as.factor(ifelse(df_plot$P.Value < 0.05,'DEG','NOT'))
label_TGs <- df_plot[abs(df_plot$logFC)> 1,]
label_TGs <- df_plot[df_plot$P.Value < 0.05,]

df_plot$label <- NA
df_plot$label[match(rownames(label_TGs), rownames(df_plot))] <- rownames(label_TGs)
write.xlsx(df_plot, file = "ExtendedDataFigure1b.xlsx", sheetName = "Sheet1", rowNames = FALSE)

p2 <- ggplot(data=df_plot, aes(x=logFC, y=logFDR,color=change)) + 
  geom_point(alpha=1, size=3.5) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black")+
  theme_set(theme_set(theme_bw(base_size=20)))+ 
  xlab("log2 fold change") + 
  ylab("-log10 p-value") + 
  theme(plot.title = element_text(size=15,hjust = 0.5)) + 
  scale_colour_manual(values = c('#ff7f0e','#1f77b4'))+   # 
  geom_text(data=label_TGs, aes(label=rownames(label_TGs)), 
            vjust=1, hjust=1, size=4, color="black")+
  theme_base()
p2

pdf("./visualize/volcanoPlot_GSE230009_cellLine_Ligand_TGFB1_v1.pdf",width = 8,height = 7)
print(p2)
dev.off() 

p3 <- ggscatter(df_plot, x = "logFC", y = "logFDR",
                color = "change", 
                palette = c('#ff7f0e','#1f77b4'),
                size = 2.5,
                label = df_plot$label, 
                font.label = c(14,'black'),
                repel = T) + theme_base()+xlab("log2(FC)") + ylab("-log10(FDR)")+
  geom_hline(yintercept = -log10(0.05), linetype="dashed") 
p3
pdf("./visualize/volcanoPlot_GSE230009_cellLine_Ligand_TGFB1_v2.pdf",width = 8,height = 7)
print(p3)
dev.off() 



# calculate AUROC and AUPRC
res_path <- '/home/yll/velocity_methods/01_analysis/apply_in_prostate/area_4000x6000_5000x7000_input_v2/MLnet_para1_46'
LRTG_im <- read_csv(paste0(res_path,"/regulate_result/LRTG_regu_mtx_all.csv"))
colnames(LRTG_im)[1] <- "TGs"

LRTG_im_long <- LRTG_im %>%
  pivot_longer(
    cols = 2:207,          
    names_to = "LRpair",      
    values_to = "score"      
  )

LRTG_im_long <- tidyr::separate(LRTG_im_long, LRpair, c('Ligand', 'Receptor'), sep = '_')
colnames(LRTG_im_long) <- c('Target','Ligand','Receptor','im_norm')
LRTG_im_long <- LRTG_im_long[LRTG_im_long$im_norm!=0,]
LRTG_im_long$im_norm <- (LRTG_im_long$im_norm-min(LRTG_im_long$im_norm))/(max(LRTG_im_long$im_norm)-min(LRTG_im_long$im_norm))

score <- LRTG_im_long %>% filter(Ligand == "TGFB1")

gene <- intersect(score$Target,rownames(diff_expr_results))
pred <- score$im_norm[match(gene,score$Target)]
label <- gene %in% rownames(sig_genes)
res.tmp <- get_evaluate_metrics(pred,label)
res.tmp <- res.tmp$perf_metrics


find_optimal_cutoff <- function(TPR, FPR, threshold){
  
  y = TPR - FPR
  Youden_index = which.max(y)
  # optimal_threshold = threshold[Youden_index]
  return(Youden_index)
  
}
pred <- prediction(pred, label)

perf_ROC <- performance(pred, measure = "tpr", x.measure = "fpr")
ind_ROC <- find_optimal_cutoff(perf_ROC@y.values[[1]],perf_ROC@x.values[[1]],perf_ROC@alpha.values[[1]])
cutoff_ROC <- perf_ROC@alpha.values[[1]][ind_ROC]

perf_PRC <- performance(pred, measure = "prec", x.measure = "rec")
ind_PRC <- find_optimal_cutoff(perf_PRC@y.values[[1]],perf_PRC@x.values[[1]],perf_PRC@alpha.values[[1]])
cutoff_PRC <- perf_PRC@alpha.values[[1]][ind_PRC]

AUC <- performance(pred, measure = "auc")
AUCPR <- performance(pred, measure = "aucpr")

pdf("./visualize/GSE230009_AUROC.pdf", width = 6, height = 6)  
plot(perf_ROC, col = "red", lwd = 2,
     main = "ROC Curve",
     xlab = "FPR",
     ylab = "TPR")
abline(a = 0, b = 1, lty = 2, col = "gray")  
dev.off()


plot(perf_PRC, col = "blue", lwd = 2,
     main = "Precision-Recall Curve",
     xlab = "Recall",
     ylab = "Precision",
     ylim = c(0, 1))

aucpr_value <- AUCPR@y.values[[1]]

fpr <- perf_ROC@x.values[[1]]
tpr <- perf_ROC@y.values[[1]]
df_roc <- data.frame(FPR = fpr, TPR = tpr)

auc_value <- round(AUC@y.values[[1]], 3)
p1 <- ggplot(df_roc, aes(x = FPR, y = TPR)) +
  geom_line(color = "red", size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(title = "ROC Curve", x = "False Positive Rate", y = "True Positive Rate") +
  annotate("text", x = 0.6, y = 0.2, label = paste("AUC =", auc_value), size = 5, color = "black")+
  theme_base() 

pdf("./visualize/GSE230009_AUROC.pdf",
    width = 5.5,height = 5)
p1
dev.off()


recall <- perf_PRC@x.values[[1]]
precision <- perf_PRC@y.values[[1]]
pr_data <- data.frame(Recall = recall, Precision = precision)

aucpr_value <- AUCPR@y.values[[1]]
aucpr_text <- paste0("AUPRC = ", round(aucpr_value, 3))

p2 <- ggplot(pr_data, aes(x = Recall, y = Precision)) +
  geom_line(color = "blue", size = 1) +
  labs(title = "Precision-Recall Curve",
       x = "Recall", y = "Precision") +
  ylim(0, 1) +
  theme_base()  +
  annotate("text", x = 0.6, y = 0.2, label = paste("AUC =", auc_value), size = 5, color = "black")
pdf("./visualize/GSE230009_AUPRC.pdf",
    width = 5.5,height = 5)
p2
dev.off()

