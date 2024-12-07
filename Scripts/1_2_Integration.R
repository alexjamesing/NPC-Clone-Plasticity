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
library(purrr)
library(harmony)
library(plyr)

#0. read data-------------------------------------------
out="~/path_to_output/"
data=readRDS(paste0(out,"rds/samples_merge.rds"))

# 1.scaling---------------------------------------

data=data %>%
  NormalizeData()%>%
  FindVariableFeatures(nfeatures=3000) %>%
  ScaleData(vars.to.regress = c("nFeature_RNA", "percent.mt")) %>%
  RunPCA(npcs=50) %>%
  RunUMAP(dims=1:30)


# 2. Integration-------------------------------

data=RunHarmony(data, group.by.vars="orig.ident",dim.use=1:30,max.iter.harmony=50)
data=RunUMAP(data,reduction="harmony",dim=1:30)


saveRDS(data,paste0(out,"rds/sample_merge_integrtion.rds"))



