library(GEOquery)
library(limma)
library(ggpubr)
library(ggthemes)
library(igraph)
library(plotrix)
library(ggraph)
library(readxl)
library(tidyverse)
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)

rm(list = ls())
gc()

setwd('/home/yll/velocity_methods/01_analysis/apply_in_prostate/cellLine/')
source('/home/yll/velocity_methods/01_analysis/code/code.R')

# load data
cellline_path <- "/home/yll/velocity_methods/01_analysis/apply_in_prostate/cellLine/"

data <- read_excel("./data/GSE274287/GSE274287_Processed_RNA-seq_cancer_cell_line.xlsx")
data <- column_to_rownames(data,var='geneID')

GPL=getGEO(filename = "./data/GSE274287/GSE274287_family.soft.gz")

ids <- data.frame(ID = rownames(data))
ENTREZID <- bitr(ids$ID, 
                 fromType = "ENSEMBL",
                 toType = c("SYMBOL"),
                 OrgDb = org.Hs.eg.db)
ids <- merge(ids, ENTREZID, by.x = "ID", by.y = "ENSEMBL")

expr <- data[rownames(data) %in% ids$ID, ]
ids <- ids[match(rownames(expr),ids$ID),]

tmp <- by(expr, ids$SYMBOL, function(x) rownames(x)[which.max(rowMeans(x))])
probes <- as.character(tmp)

expr <- expr[rownames(expr) %in% probes, ]
ids <- ids[ids$ID %in% probes, ]
rownames(expr) <- ids$SYMBOL

saveRDS(expr,"./GSE274287_expr_mat.rds")

count <- readRDS("~/velocity_methods/01_analysis/apply_in_prostate/cellLine/GSE274287_expr_mat.rds")

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

sample_info <- data.frame(geo_accession = names(GPL@gsms), 
                          condition=c('control',"treatment ","TGFb",
                                      "treatment","control",
                                      "TGFb"))

colnames(expr) <- sample_info$geo_accession

compare_group <- c("GSM8446946","GSM8446948","GSM8446950", "GSM8446951")
sample_info_TGFB1 <- sample_info[sample_info$geo_accession %in% compare_group, ] 

expr1 <- expr[,compare_group]

# limma analysis
design <- model.matrix(~0 + sample_info_TGFB1$condition)
colnames(design) <- c("control","treatment")
fit <- lmFit(expr1, design)

contrast_matrix <- makeContrasts(control_vs_treatment = control - treatment, levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

diff_expr_results <- topTable(fit2, adjust.method = "fdr", number = Inf)
head(diff_expr_results)


sig_genes <- subset(diff_expr_results, P.Value < 0.05)


head(sig_genes)
saveRDS(sig_genes,'cellLine_ligand_TGFB1_sigGenes.rds')

# load MLnet result
mlnet_path <- '/home/yll/velocity_methods/01_analysis/apply_in_prostate/area_4000x6000_5000x7000_input/'
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

p2 <- ggplot(data=df_plot, aes(x=logFC, y=logFDR,color=change)) + 
  geom_point(alpha=1, size=3.5) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black")+
  theme_set(theme_set(theme_bw(base_size=20)))+ 
  xlab("log2 fold change") + 
  ylab("-log10 p-value") + 
  theme(plot.title = element_text(size=15,hjust = 0.5)) + 
  scale_colour_manual(values = c('#ff7f0e','#1f77b4'))+   
  geom_text(data=label_TGs, aes(label=rownames(label_TGs)), 
            vjust=1, hjust=1, size=4, color="black")+
  theme_base()
p2













