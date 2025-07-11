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

# What the script does -- step by step
# Step	Code block	Purpose
# Load packages	library(Seurat) … library(plyr)	Brings in Seurat (single-cell RNA-seq), plotting and tidy-verse utilities, Harmony for batch correction (though not used later).
# Set I/O paths	out="~/path_to_output/"
# files=list.files("~/path_to_cellranger_outputs")	Defines where to find Cell Ranger “filtered_feature_bc_matrix” folders (files) and where to write results (out).
# Read doublet calls	doub_df = read_csv("~/0_Doublet_prediction/Doublets_corrected_scores.csv")	Loads an external table with Scrublet scores and empirically set X-pos cut-offs.
# Loop over each Cell Ranger run	for (i in 1:length(files)) { … }	• Builds data_path for sample i.
# • Reads the count matrix with Read10X.
# • Creates a Seurat object (CreateSeuratObject).
# • Looks up that sample in doub_df, labels each barcode as “singlet” or “doublet”, and stores the label in the object’s metadata (AddMetaData).
# • Removes doublets by splitting on the new metadata column and keeping only the “singlet” subset.
# • Adds the singlet object to seurat_sample_list.
# Add basic QC metrics	raw_process <- function(x) { … }
# seurat_sample_list = lapply(seurat_sample_list, raw_process)	For every sample, calculates percentage of reads that are: mitochondrial genes (^mt-), ribosomal genes (^Rp-), transgene markers (mCherry, DTA, Cre). Each percentage is stored in @meta.data.
# Per-sample QC filtering	qc_process <- function(x) { … }
# lapply(seurat_sample_list, qc_process)	• Prints and writes a CSV of decile (0–100 %) distributions for nCount, nFeature, and mitochondrial percentage.
# • Filters cells: keep only those with > 500 detected genes and < 5 % mitochondrial content.
# Returns the filtered object.
# Normalisation + variable-feature selection	norma_process <- function(x) { … }
# norm_list = lapply(seurat_sample_list, norma_process)	• NormalizeData (log-normalisation).
# • FindVariableFeatures (top 3 000 genes).
# • Plots the variable-feature scatter, highlights the top 30 features, saves as PDF under out/plots/.
# Merge all samples	datam = merge(norm_list[[1]], y = c(unlist(norm_list)[2:24]), …)	Combines the (up to 24) individual, pre-filtered Seurat objects into one multi-sample object, adding sample-specific prefixes to cell barcodes.
# Save result	saveRDS(datam, paste0(out,"rds/samples_merge.rds"))	Writes the merged, QC-filtered, normalised dataset to disk for downstream analysis.
# In plain language

# Read every Cell Ranger folder, build a Seurat object, and remove doublets using an external Scrublet score sheet.
# Compute standard QC metrics (mitochondrial, ribosomal, transgene reads), write a decile report, and discard low-quality cells.
# Log-normalise each sample and pick 3 000 highly variable genes; save variable-feature plots.
# Merge all cleaned samples into a single Seurat object and store it as samples_merge.rds.
# That RDS can be loaded later for scaling, dimensional reduction, Harmony/CCA integration, clustering, etc.







