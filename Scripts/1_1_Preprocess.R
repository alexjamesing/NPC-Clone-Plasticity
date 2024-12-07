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

# 1. read data --------

out="~/path_to_output/"

files=list.files("~/path_to_cellranger_outputs")
files_path="~/path_to_cellranger_outputs"


seurat_sample_list=list()
doub_df=read_csv("~/0_Doublet_prediction/Doublets_corrected_scores.csv")



for ( i in 1:length(files)){
  print(files[i])
  data_path=paste0(files_path[i],"/",files[i],"/outs/filtered_feature_bc_matrix",sep="")
  name=sapply(str_split(files[i],pattern="_result_new",n=2),'[',n=1)
  counts=Read10X(data.dir=data_path)
  seurat <- CreateSeuratObject(counts, project=files[i])
  m_df=doub_df %>% dplyr::filter(samples==name)
  m_df=m_df %>% drop_na()
  m_df$corrected_doublets=ifelse(m_df$scrublet_Scores>m_df$X_pos,"doublet","singlet")
  m_df=m_df %>% select(c("samples","raw_barcodes","corrected_doublets"))
  rownames(m_df)=m_df$raw_barcodes
  seurat=AddMetaData(seurat,metadata = m_df )
  sing=SplitObject(seurat,split.by = "corrected_doublets")
  seurat_sample_list[[i]]=sing[["singlet"]]
  names(seurat_sample_list)[i]=files[i]
  gc()
  
}


raw_process=function(x){
  x[["percent.mt"]]=PercentageFeatureSet(x,pattern="^mt[-\\.]")
  x[["percent.rb"]]=PercentageFeatureSet(x,pattern="^Rp[-\\.]")
  x[["percent.mCherry"]]=PercentageFeatureSet(x,pattern="^mCherry")
  x[["percent.Dta"]]=PercentageFeatureSet(x,pattern="Dta")
  x[["percent.Cre"]]=PercentageFeatureSet(x,pattern="Cre")
  x
}


seurat_sample_list=lapply(seurat_sample_list,raw_process)



qc_process=function(x){
  ncount=quantile(x@meta.data$nCount_RNA,seq(0,1,0.1))
  nfeature=quantile(x@meta.data$nFeature_RNA,seq(0,1,0.1))
  percent_mt=quantile(x@meta.data$percent.mt,seq(0,1,0.1))
  q_df=data.frame(ncount=ncount,
                  nfeature=nfeature,
                  percent_mt=percent_mt)
  q_df$quantile=rownames(q_df)
  q_df$sample_id=x@meta.data$orig.ident[1]
  print(q_df)
  write_csv(q_df,paste0(out,"files/",x@meta.data$orig.ident[1],"_","quantile.csv"))
  x <- subset(x, subset = nFeature_RNA > 500 & percent.mt < 5)
  x
}

seurat_sample_list=lapply(seurat_sample_list,qc_process)



norma_process=function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 3000)
  top_features <- head(VariableFeatures(x), 30)
  plot1 <- VariableFeaturePlot(x,cols = c("black", "grey"))
  plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
  plot2+labs(title=x@meta.data$orig.ident[1])
  gc()
  ggsave(paste0(out,"plots/",x@meta.data$orig.ident[1],"_","VF.pdf"),units = "in",width = 9,height = 6)
  x
  
}

norm_list=lapply(seurat_sample_list,norma_process)

datam=merge(norm_list[[1]],y=c(unlist(norm_list)[2:24]),add.cell.ids = c(names(norm_list)))

saveRDS(datam,paste0(out,"rds/samples_merge.rds"))


