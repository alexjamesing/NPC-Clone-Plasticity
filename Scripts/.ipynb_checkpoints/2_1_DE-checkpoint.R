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

register(MulticoreParam(30))


out="~/path_to_output/"


# 1. prepare data ---------------------------

data=readRDS("~/data_annotatios.rds")

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



