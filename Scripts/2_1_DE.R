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
library(ggridges)
#BiocManager::install("YuLab-SMU/treedataverse")
library(treedataverse)
library(ComplexHeatmap)
library(circlize)
library(BiocParallel)
library(DESeq2)
library(scran)
library(glmGamPoi)
library(presto)

base_folder <- "/omics/odcf/analysis/OE0574_projects/brainbreaks/single_cell_BB_E17-5_Alex/result_all_102025/R_files/Results"
annotation_dir <- file.path(base_folder, "new_annotation", "rds")
annotated_cells_path <- file.path(annotation_dir, "single_cell_annotations.rds")
if (!file.exists(annotated_cells_path)) {
  stop("Annotated single-cell object not found at: ", annotated_cells_path)
}

register(MulticoreParam(30))


out <- file.path(base_folder, "new_annotation")


# 1. prepare data ---------------------------

data <- readRDS(annotated_cells_path)

Pb=AggregateExpression(data,return.seurat = TRUE,assays="RNA",slot="counts",group.by = c("Mouse","Condition","annotations","Red"))


Pb$meta=rownames(Pb@meta.data)
Pb@meta.data=Pb@meta.data %>% mutate(Red=case_when(grepl(pattern = "_False",Pb$meta)~"False",T~"True"))
Pb@meta.data$meta=sub(pattern="_False",replacement="",Pb$meta)
Pb@meta.data$meta=sub(pattern="_True",replacement="",Pb$meta)

Pb@meta.data[c("Mouse","Condition","Annotation")]=str_split_fixed(Pb$meta,pattern = "_",n=3)
Pb$groups=paste0(Pb$Red,"-",Pb$Annotation)


Cell_list=SplitObject(Pb,split.by = "groups")


# 2. Deseq2 -------------------------

process_deseq2=function(x){
  count_df=data.frame(x@assays[["RNA"]]@counts)
  coldata=data.frame(condition=colnames(count_df))
  rownames(coldata) <- coldata$condition
  
  coldata$condition=sub(pattern="_False",replacement="",coldata$condition)
  coldata$condition=sub(pattern="_True",replacement="",coldata$condition)
  
  coldata[c("mouse","condition","annotations")]=str_split_fixed(coldata$condition,pattern="_",n=3)
  
  dds <- DESeqDataSetFromMatrix(countData =count_df,colData = coldata,design = ~ condition)
  dds = DESeq2::estimateSizeFactors(dds,type="poscount")
  scr=computeSumFactors(dds)
  DESeq2::sizeFactors(dds)=sizeFactors(scr)
  
  dds$condition <- factor(dds$condition, levels = c("control","complemented"))
  dds$condition <- relevel(dds$condition, ref = "control")
  dds$condition
  dds <- DESeq(dds,test="LRT",useT = T,fitType = "glmGamPoi",minmu=1e-06,minReplicatesForReplace = Inf,reduced = ~1,parallel = T)
  res <- results(dds)
  res <- results(dds, name="condition_complemented_vs_control")
  res <- results(dds, contrast=c("condition","complemented","control"))
  res
}

deseq2_res=lapply(Cell_list,process_deseq2)
saveRDS(deseq2_res,"~/Pseudobulk_deseq2_res_list.rds")

deseq2_res_df=list()
for (i in 1:length(deseq2_res)){
  s=names(deseq2_res)[i]
  print(s)
  x=data.frame(deseq2_res[[s]])
  x$gene=rownames(x)
  deseq2_res_df[[i]]=x
  names(deseq2_res_df)[i]=s
  if(grepl("False",s)){
    name=sub(pattern="False-",replacement = "",s)
    write_csv(x,file=paste0("/home/l538g/workingf/brainbreaks/single_cell_BB_E17-5/Giulia_data_analysis/result_all_05082024/R_files/Results/2_2_DE/files/","N_Red_",name,"_DEseq2_nonshrink.csv",sep=""))
  }else{
    name=sub(pattern="True-",replacement = "",s)
    write_csv(x,file=paste0("/home/l538g/workingf/brainbreaks/single_cell_BB_E17-5/Giulia_data_analysis/result_all_05082024/R_files/Results/2_2_DE/files/","Red_",name,"_DEseq2_nonshrink.csv",sep=""))
  }
}

filter_significant=function(x){
  x=x%>% dplyr::filter(padj<=0.05)
  x
}

deseq2_res_df_sig=lapply(deseq2_res_df,filter_significant)


# Summary
# The script performs pseudobulk differential-expression analysis on a single-cell dataset that has already been annotated.
# For every cell-type (“Annotation”) and every Red = True/False subset, it:

# Aggregates counts per mouse × condition × cell-type × Red combination (pseudobulk samples).
# Splits the resulting Seurat object into one list element per Red–Annotation group.
# Runs DESeq2 (LRT with glmGamPoi fit) comparing complemented vs control within that group, using scran’s computeSumFactors for library-size normalisation.
# Writes CSV files with full and FDR-filtered results and saves a list of all DESeq2 result tables.
# Stage	Key code / actions	Purpose / effect
# Load packages & set parallelism	register(MulticoreParam(30))	Loads Seurat, DESeq2, scran, etc.; enables 30-core parallelism for DESeq2.
# 1 Prepare data	data <- readRDS(".../new_annotation/rds/single_cell_annotations.rds")	Loads the annotated Seurat object containing single cells.
# AggregateExpression(…, group.by = c("Mouse","Condition","annotations","Red"), slot = "counts") → Pb	Collapses raw counts into pseudobulk columns for every unique combination of mouse, experimental condition, cell-type annotation, and Red status.
# Post-processing of Pb@meta.data	Parses the column names to extract Mouse, Condition, Annotation; builds a new categorical label groups = Red-Annotation.
# Cell_list <- SplitObject(Pb, split.by = "groups")	Produces a list where each element contains only the pseudobulk samples of one Red/Annotation combination.
# 2 DESeq2 per group	process_deseq2() (applied via lapply)	For each list element:
# • Converts counts to a DESeqDataSet.
# • Uses scran’s computeSumFactors for size-factor normalisation (type = "poscount" fallback).
# • Sets condition factor with “control” as reference.
# • Runs likelihood-ratio test (DESeq(test = "LRT", …, reduced = ~ 1, fitType = "glmGamPoi")).
# • Returns the Wald contrast complemented vs control result table.
# Results saved	• All result objects collected in deseq2_res.
# • Each contrast written to CSV; filename prefix “Red_” or “N_Red_” depending on Red status.
# Significance filter	deseq2_res_df_sig = lapply(deseq2_res_df, filter(padj <= 0.05))
# 3 Save	saveRDS(deseq2_res, "~/Pseudobulk_deseq2_res_list.rds")	Keeps the full list of DESeq2 results for later use.
# In essence, the script converts single-cell data into pseudobulk replicates, performs DESeq2-based differential expression between complemented and control conditions inside every Red/Annotation subgroup, and exports both the full and significant gene lists.


