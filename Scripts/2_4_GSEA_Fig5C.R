library(Seurat)
library(SeuratObject)
library(pheatmap)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(Matrix)
library(R6)
library(dplyr)
library(tidyverse)
library(MetBrewer)
library(purrr)
library(harmony)
library(plyr)
library(treedataverse)
library(BiocParallel)
library(DESeq2)
library(scran)
library(glmGamPoi)
library(presto)
library(fgsea)
library(msigdbr)
library(limma)
library(qs)
register(MulticoreParam(50))


base_folder <- "/omics/odcf/analysis/OE0574_projects/brainbreaks/single_cell_BB_E17-5_Alex/result_all_102025/R_files/Results"
annotation_dir <- file.path(base_folder, "new_annotation", "rds")
annotated_cells_path <- file.path(annotation_dir, "single_cell_annotations.rds")
if (!file.exists(annotated_cells_path)) {
  stop("Annotated single-cell object not found at: ", annotated_cells_path)
}
out <- file.path(base_folder, "new_annotation")
gsea_out <- file.path(out, "4_1_Pathway")
dir.create(file.path(gsea_out, "files"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(gsea_out, "plots"), recursive = TRUE, showWarnings = FALSE)

data <- readRDS(annotated_cells_path)
data=SplitObject(data,split.by = "annotations")


#1. single cell gsea-----------------------------

## gsea Hallmark ----------------



run_gsea=function(x){
  out_plot=gsea_out
  n=x$split_groups[1]
  print(n)
  Idents(x)=as.factor(x$Condition)
  wilcox_res=wilcoxauc(x,group_by = "Condition")
  
  print(head(wilcox_res))
  m_db=msigdbr(species="Mus musculus",db_species="MM",category="MH")
  #print(head(m_db))
  fgsea_sets=m_db %>% split(x=.$gene_symbol,f=.$gs_name)
  
  # Control
  control_sets=wilcox_res %>% filter(group=="control")
  wilcox_res_genes=control_sets %>% arrange(desc(logFC))%>% dplyr::select(feature,logFC)
  ranks=deframe(wilcox_res_genes)
  fgsea_res=fgseaMultilevel(fgsea_sets,stats=ranks,eps=0)
  write_csv(fgsea_res,file.path(out_plot,"files",paste0(n,"_control_H.csv")))
  quan_90=quantile(fgsea_res$NES,c(0.1,0.2,0.5,0.9,0.95,1),na.rm = T)["90%"]
  quan_10=quantile(fgsea_res$NES,c(0.1,0.2,0.5,0.9,0.95,1),na.rm = T)["10%"]
  res=fgsea_res %>%arrange(desc(NES)) %>%arrange(padj)%>%head(50)
  ggplot(res %>% dplyr::filter(padj<=0.01) %>% head(15),aes(x=reorder(pathway,NES),y=NES))+
    geom_col(fill="lightgrey")+
    coord_flip()+
    theme_classic()+
    xlab("Pathways")
  ggsave(file.path(out_plot,"plots",paste0(n,"_control_H.pdf")),width = 16,height = 12)
  
  # Complemented
  
  complemented_sets=wilcox_res %>% filter(group=="complemented")
  wilcox_res_genes=complemented_sets %>% arrange(desc(logFC))%>% dplyr::select(feature,logFC)
  ranks=deframe(wilcox_res_genes)
  fgsea_res=fgseaMultilevel(fgsea_sets,stats=ranks,eps=0)
  write_csv(fgsea_res,file.path(out_plot,"files",paste0(n,"_complemented_H.csv")))
  quan_90=quantile(fgsea_res$NES,c(0.1,0.2,0.5,0.9,0.95,1),na.rm = T)["90%"]
  quan_10=quantile(fgsea_res$NES,c(0.1,0.2,0.5,0.9,0.95,1),na.rm = T)["10%"]
  res=fgsea_res %>%arrange(desc(NES)) %>%arrange(padj)%>%head(50)
  ggplot(res %>% dplyr::filter(padj<=0.01) %>% head(15),aes(x=reorder(pathway,NES),y=NES))+
    geom_col(fill="lightgrey")+
    coord_flip()+
    theme_classic()+
    xlab("Pathways")
  ggsave(file.path(out_plot,"plots",paste0(n,"_complemented_H.pdf")),width = 16,height = 12)
}






# 5. Plot -----------------------------
### 1. gsea res -----------------------------------------------
path=file.path(gsea_out,"files")
f=list.files(path,pattern = "^True")
f=f[which(grepl("complemented",f))]

gsea_res=list()

for(i in 1:length(f)){
  print(f[i])
  n=gsub("True_","",f[i])
  n=gsub("_H.csv","",n)
  t=read_csv(file.path(path,f[i]),trim_ws = T) %>% mutate(cells=n) %>% filter(padj <=0.01)
  if(nrow(t) >0){
    gsea_res[[i]]=t
    names(gsea_res)[i]=gsub(".csv","",f[i])
  }
}

gsea_res=Filter(Negate(is.null),gsea_res)

gsea_df=do.call(rbind,gsea_res)


gsea_df$shape=ifelse(gsea_df$NES>0,24,25)
gsea_df$groups=gsub("_complemented","",gsea_df$cells)
library(circlize)
col_fun = colorRamp2(c(0,4,16), c("#EEEEEE","#edc775","#b75347"))

gsea_df$pathway=gsub("HALLMARK_","",gsea_df$pathway)
gsea_df$pathway=factor(gsea_df$pathway,levels=c(sort(unique(gsea_df$pathway))[c(6,13:17,25,26,2:5,9,12,18,19,22,24,7,27,8,11,20,21,1,28,10,23)]))
gsea_df$groups=factor(gsea_df$groups,levels=c(unique(gsea_df$groups)[c(15,8,12,13,10,1,4,14,2,6,7,9,5,11,3)]))



ggplot(gsea_df%>% filter(groups %in% c("Radial_glia","Glioblast","Neuroblast","OPC","Intermediate_progenitor","Microglia")),aes(x=groups,y=pathway,fill=-log(padj,10)))+
  geom_point(aes(shape=shape,size=abs(NES)))+
  scale_shape_identity()+
  scale_fill_gradient(low="#EEEEEE",high="darkred")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=90))+
  labs(size="NES")+ylab("Pathways")

ggsave(file.path(gsea_out,"plots","gsea_res_nc.pdf"),units = "in",width = 8,height = 6)

write.csv(gsea_df,file.path(gsea_out,"files","gsea_df.csv"),row.names = F)

# Summary
# The script performs cell-type–resolved gene-set enrichment analysis (GSEA) on a single-cell Seurat object.
# For every annotated cell type it:

# Ranks genes by Wilcoxon AUC between control and complemented cells.
# Runs fgsea on Hallmark gene sets (mouse) separately for the control-biased and complemented-biased ranks.
# Saves a CSV of enrichment statistics and a bar-plot PDF for each cell type × condition.
# After all cell types are processed, it reloads the “complemented” CSV files, concatenates them, and draws a facetted dot plot that summarises significant Hallmark pathways (padj ≤ 0.01) across selected progenitor / immune populations, with dot size = |NES| and shape = direction of change.
# CSV tables and a final plot (gsea_res_nc.pdf) are written to ~/path_to_output/.

# Step-by-step detail
# Stage	Code / action	Purpose
# Load packages & parallelism	register(MulticoreParam(50))	Enables fgsea multicore execution.
# Read data & split by cell type	```r data <- readRDS(".../new_annotation/rds/single_cell_annotations.rds")	
# data <- SplitObject(data, split.by = "annotations")```	Produces a list of Seurat objects, one per cell-type annotation.	
# Define run_gsea()	Input = one of the split Seurat objects (single cell type).
# Steps inside:
# 1. Set identities to Condition (control vs complemented).
# 2. Wilcoxon AUC (presto::wilcoxauc) → table of logFC per gene.
# 3. Obtain Hallmark gene sets for mouse (msigdbr).
# 4. For “control” ranks: sort genes by logFC, run fgseaMultilevel, write CSV, make bar plot of top NES pathways (padj ≤ 0.01).
# 5. Repeat for “complemented” ranks.	
# (Function not called here)	In this snippet run_gsea() isn’t executed, implying it was run previously and its CSV outputs already exist in ~/path_to_output/.	
# Collect complemented results	```r f <- list.files(out, pattern = "^True*_")	
# f <- f[grepl("complemented", f)]```	Loads every *_complemented_H.csv (Hallmark results for Red = True cells) into a list, filters to padj ≤ 0.01, binds rows into gsea_df.	
# Re-format for plotting	Adds columns:
# • shape = 24 (NES > 0) or 25 (NES < 0).
# • groups = cell-type name.
# • Cleans Hallmark names.
# • Orders factors so facets & axes are readable.	
# Dot-plot heatmap	```r ggplot(gsea_df %>%	
#     filter(groups %in% c("Radial_glia","Glioblast","Neuroblast","OPC",
#                          "Intermediate_progenitor","Microglia")),  
#     aes(x = groups, y = pathway, fill = -log10(padj))) +  
#     geom_point(aes(shape = shape, size = abs(NES))) …``` | Produces a grid where each dot encodes: <br>• *x* = cell type, *y* = Hallmark pathway. <br>• Fill colour = significance (−log10 padj). <br>• Dot size = |NES|, shape indicates up- (▲) or down-regulation (▼). |
# | Save outputs | • Plot saved as gsea_res_nc.pdf. <br>• Full combined table saved as gsea_df.csv. |

# File outputs

# File	Contents
# out/plots/<celltype>_control_H.pdf	Top 15 Hallmark pathways enriched in control-biased genes for that cell type.
# out/plots/<celltype>_complemented_H.pdf	Same for complemented-biased genes.
# out/files/<celltype>_<condition>_H.csv	Full fgsea statistics (NES, padj, leading-edge genes) for each cell type × condition.
# out/plots/gsea_res_nc.pdf	Multi-cell-type dot plot of significant complemented enrichments for selected progenitor/immune groups.
# out/files/gsea_df.csv	Combined table of significant fgsea results (padj ≤ 0.01).
# This pipeline therefore translates single-cell differential-expression signals into pathway-level insights for every annotated population and highlights which biological programs are selectively altered in the complemented condition.
