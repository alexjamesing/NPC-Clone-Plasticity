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

# 0. read samples------------------------------------

data=readRDS(paste0(out,"rds/sample_merge_integrtion.rds"))
out="~/path_to_output/"

# 1. checking resolutions---------------
## manally check
data=FindNeighbors(data,reduction="harmony",dims=1:30) %>% 
  FindClusters(resolution=0.5) %>%
  FindClusters(resolution=1) %>%
  FindClusters(resolution=1.5) %>%
  FindClusters(resolution=2) %>%
  FindClusters(resolution=3) %>%
  FindClusters(resolution=4) 

saveRDS(data,paste0(out,"rds/sample_clusters.rds"))

resolution=c(0.5,1,1.5,2,3,4)
groups=paste(rep("RNA_snn_res.",6),resolution,sep="")

for(i in 1:length(groups)){
  g=groups[i]
  print(g)
  print(table(data@meta.data %>%select(g)))
  h=hist(as.numeric(data@meta.data[[g]]))
  print(h)
  p=DimPlot(data,group.by=groups[i],raster = F)
  print(p)
  ggsave(paste0("/home/l538g/workingf/brainbreaks/single_cell_BB_E17-5/Giulia_data_analysis/result_all_05082024/R_files/Results/1_3_Metacells/plots/",g,"_UMAP.pdf"),units = "in",width = 24,height = 8)
  
}

subset_clusters=Seurat::SplitObject(data,split.by = "RNA_snn_res.4")

saveRDS(subset_clusters,paste0(out,"rds/subset_clusters.rds"))




# 2. checking resolutions---------------

## check cluster size------

check_cell_number=function(x){
  n=nrow(x@meta.data)
  print(n)
  n
} 

cluster_size=unlist(lapply(subset_clusters,check_cell_number))



do_Hclu=function(x,n){
  number=check_cell_number(x)
  if(number>=50){
    x=FindNeighbors(x,dim=1:30)
    x=FindClusters(x,resolution=n)
    cluster_n=paste("RNA_snn_res.",n,sep="")
    gc()
    d=x@meta.data%>% dplyr::select(c("RNA_snn_res.4",paste(cluster_n)))
    d$barcodes=rownames(d)
    d
  }
}


Hclu=list()
filter_group=list()


sub_resoluion=seq(1,10,0.4)

raw_res="RNA_snn_res.4"
do_Hclu=function(x,raw_res){
  sub_resoluion=seq(1,10,0.4)
  features_df=x@meta.data %>% select(c("orig.ident","nCount_RNA","nFeature_RNA"))
  features_df$orig.ident="raw_cluster"
  high_res_df=d=x@meta.data%>% dplyr::select(raw_res)
  colnames(high_res_df)="raw_cluster"
  number=check_cell_number(x)
  for (i in 1:length(sub_resoluion)){
    if(number>2){
      res=sub_resoluion[i]
      x=FindNeighbors(x,dim=1:30)
      x=FindClusters(x,resolution=res)
      cluster_n=paste("RNA_snn_res.",res,sep="")
      print(cluster_n)
      gc()
      d=x@meta.data%>% dplyr::select(paste(cluster_n))
      high_res_df=cbind(high_res_df,d)
      colnames(d)="meta_group"
    
      # check features 
      x=AddMetaData(x,d)
      meta_cell=AggregateExpression(object =x, group.by ="meta_group",return.seurat = TRUE)
      f=meta_cell@meta.data
      f$orig.ident=cluster_n
      features_df=rbind(features_df,f)
    }else{
      break
    }
    number=check_cell_number(meta_cell)
    
  }
  result=list(features_df,high_res_df)
  names(result)=c("nFeatures_df","high_res_df")
  result
}
  


Hclu=lapply(subset_clusters,do_Hclu,raw_res="RNA_snn_res.4")

saveRDS(Hclu,paste0(out,"rds/Hclu.rds"))


plot_nFeatures=function(x,name){
  d=x[["nFeatures_df"]]
  dd=d %>% dplyr::group_by(orig.ident)%>% dplyr::summarise(median_v=median(nFeature_RNA),stand_d=sd(nFeature_RNA))
  t=dd %>% filter(orig.ident!="raw_cluster")
  ggplot()+
    geom_boxplot(data=d,aes(x=as.factor(orig.ident),y=nFeature_RNA))+
    geom_line(data=dd,aes(x=as.factor(orig.ident),y=median_v,group=1),color="#d8443c")+
    geom_point(data= d %>% filter(orig.ident==t[which(t$stand_d==min(t$stand_d)),]$orig.ident),aes(x=as.factor(orig.ident),y=nFeature_RNA),color="#f4c40f")+
    xlab("resolution")+
    scale_x_discrete(guide=guide_axis(angle = 90))+
    theme_classic()+
    labs(title=paste("Cluster",name))
  ggsave(paste0(out,"plots/",name,"_Hclu_nFeatures.pdf"),units = "in",width = 12,height = 8)
  t
}


eval_metacell_size=function(x){
  d=x[["high_res_df"]]
  res_test=colnames(d)
  df=data.frame(raw_cluster=numeric(),
                group=character(),
                min_size=numeric(),
                max_size=numeric(),
                median_size=numeric(),
                mean_size=numeric())
  for (i in 2:length(res_test)){
    resolution=res_test[i]
    print(paste("process ",resolution))
    t=d %>% mutate(ID=rownames(d))%>%select(c("ID","raw_cluster",resolution)) %>% group_by_(resolution) %>%dplyr::mutate(size=n()) 
    dd=data.frame(raw_cluster=t$raw_cluster[1],
                  group=resolution,
                  min_size=min(t$size),
                  max_size=max(t$size),
                  median_size=median(t$size),
                  mean_size=mean(t$size))
    df=rbind(df,dd)
  }
 df$dis=(abs(50-df$min_size)+abs(150-df$max_size))/(150-50)
 df
}

data_size=lapply(Hclu,eval_metacell_size)


res_all=data.frame(
              group=character(),
              median_v=numeric(),
              stand_d=numeric(),
              raw_cluster=numeric(),
              min_size=numeric(),
              max_size=numeric(),
              median_size=numeric(),
              mean_size=numeric(),
              dis=numeric(),
              meandist_125=numeric())


for(i in 1:length(names(Hclu))){
  name=names(Hclu)[i]
  print(name)
  d1=plot_nFeatures(Hclu[[name]],name=name)
  colnames(d1)[1]="group"
  d2=data_size[[name]]
  d_m=merge(d1,d2,by="group")
  d_m$meandist_125=abs(d_m$mean_size-125)
  res=d_m %>% slice_min(meandist_125,n=3) %>% slice_min(stand_d,n=1)
  res_all=rbind(res_all,res)
}

meta_info=res_all[!duplicated(res_all$raw_cluster),]


saveRDS(data_size,paste0(out,"rds/data_size.rds"))
saveRDS(meta_info,paste0(out,"rds/meta_info.rds"))




add_meta=data.frame(
  raw_cluster=character(),
  hclu_res=character(),
  meta_group=character()
  )



for (i in 1:length(unique(meta_info$raw_cluster))){
  c=paste(unique(meta_info$raw_cluster)[i])
  sub_meta=meta_info %>% filter(raw_cluster==c)
  print(c)
  sub_df=Hclu[[c]][["high_res_df"]]%>%select(c("raw_cluster",sub_meta$group[1]))
  colnames(sub_df)[2]="hclu_res"
  sub_df$meta_group=paste0(sub_df$raw_cluster,"_",sub_df$hclu_res)
  add_meta=rbind(add_meta,sub_df)
}

saveRDS(add_meta,paste0(out,"rds/add_meta.rds"))


## Please also check manually for each clusters



# 3. aggregate---------------

data=readRDS(paste0(out,"rds/sample_clusters.rds"))
add_meta=readRDS(paste0(out,"rds/add_meta.rds"))
data<- AddMetaData(data, add_meta)


meta_cell=AggregateExpression(object =data, group.by ="meta_group",return.seurat = TRUE)

## map original umap coordinate--------------


umap_cord <- data[['umap']]@cell.embeddings
data.meta <- data@meta.data %>% dplyr::select(c("orig.ident","meta_group"))
umap_cord <- cbind(umap_cord, data.meta)
umap_cord <- as_tibble(umap_cord, rownames = NA)
colnames(umap_cord)[4] <- "high_res_clusters"
colnames(umap_cord)[3] <- "samples"
umap_cord.aggr <- umap_cord %>% 
  dplyr::group_by(high_res_clusters) %>%
  dplyr::summarise(UMAP_1 = median(as.numeric(UMAP_1)),
                   UMAP_2 = median(as.numeric(UMAP_2))) 


name_col <- c("UMAP_1", "UMAP_2")
name_row <- umap_cord.aggr$high_res_clusters
umap_cord.aggr.mtx <- as.matrix(umap_cord.aggr[,2:3])
rownames(umap_cord.aggr.mtx) <- paste(umap_cord.aggr$high_res_clusters)
colnames(umap_cord.aggr.mtx) <- name_col


meta_cell@reductions[["umap"]] <- CreateDimReducObject(
  embeddings = umap_cord.aggr.mtx, key = "UMAP_", assay = "RNA", global = TRUE)





umap_cord.aggr$ori_cluster=sapply(str_split(umap_cord.aggr$high_res_clusters,pattern="_"),"[",n=1)
umap_cord$ori_cluster=sapply(str_split(umap_cord$high_res_clusters,pattern = "_"),"[",n=1)


saveRDS(meta_cell,paste0("rds/meta_cell.rds"))


# What this (third) script does — step by step
# Stage	Key actions	Purpose / effect
# Load packages	library(Seurat) … library(plyr)	Same Seurat + plotting + Harmony stack as before.
# 0 Read input	r out <- "~/path_to_output/" # (path defined *after* read) data <- readRDS(paste0(out,"rds/sample_merge_integrtion.rds"))	Loads the Harmony-integrated, UMAP-embedded object produced by the previous script.
# 1 Generate a family of clusterings	r data %>% FindNeighbors(reduction="harmony") %>% FindClusters(res=0.5) … FindClusters(res=4)	Computes an SNN graph on the Harmony space and writes six community detections at resolutions 0.5, 1, 1.5, 2, 3, 4, storing each as RNA_snn_res.* metadata columns. Saves the object to sample_clusters.rds.
# • Quick diagnostics	Loop over the six RNA_snn_res.* columns:
# – prints cell counts and a histogram
# – draws a UMAP coloured by that resolution
# – saves each plot to PDF	Lets the user eyeball which global resolution best balances cluster separation vs over-splitting.
# Split by the highest resolution (4.0)	subset_clusters <- SplitObject(data, split.by="RNA_snn_res.4")	Creates a list in which each element is one “raw cluster” from resolution 4. This becomes the starting point for per-cluster refinement.
# 2 Per-cluster sub-clustering (“meta-cells”)	Custom functions check_cell_number, do_Hclu, eval_metacell_size, etc.	Goal: for each raw cluster, find a higher internal resolution (1.0 – 10.0 in 0.4 steps) that yields subclusters (“meta-cells”) of ≈ 50–150 cells and minimal variability in gene counts.

# Workflow per raw cluster:
# 1. Iterate over seq(1,10,0.4). For each res:
# – Re-cluster the cells.
# – Aggregate expression per subcluster (AggregateExpression).
# – Record QC metrics (nFeature_RNA) and cluster sizes.
# 2. Assemble a dataframe of candidate resolutions vs:
# – median nFeature_RNA & its SD
# – min / max / mean cluster size
# – distance from the ideal mean size 125 (custom score).
# 3. Pick the best three by mean-size distance, then the single best by lowest SD.
# 4. Store that chosen resolution for the raw cluster.
# • Outputs of step 2	* Hclu.rds – list of per-cluster QC tables.
# * data_size.rds – cluster-size statistics.
# * meta_info.rds – final chosen resolution for every raw cluster.
# * PDF plots of nFeature_RNA vs resolution.	Gives a reproducible record of how each raw cluster was sub-clustered.
# Map every cell to its chosen meta-cell	Constructs add_meta, a dataframe containing for every barcode:
# – raw_cluster (from res 4)
# – chosen high-resolution label (hclu_res)
# – concatenated meta_group = raw_cluster_hcluRes	Creates a single categorical label (meta_group) that identifies each fine-grained meta-cell across the whole dataset.
# 3 Generate aggregated (“meta-cell”) object	r data <- readRDS(sample_clusters.rds) data <- AddMetaData(data, add_meta) meta_cell <- AggregateExpression(data, group.by="meta_group", return.seurat=TRUE)	Adds the meta_group column to the original object, then aggregates counts across all cells in each meta-cell (sums the raw counts and stores them in a new Seurat object meta_cell).
# • Attach median UMAP coordinates**	Calculates the median (UMAP_1, UMAP_2) of constituent cells per meta-cell, builds a new DimReduc slot in meta_cell, and annotates both meta-cells and single cells with their parent raw_cluster (prefix before the underscore).	Allows immediate plotting of meta-cells in the same UMAP space as single cells.
# Save final object	saveRDS(meta_cell, "rds/meta_cell.rds")	Stores the fully QC’d, resolution-tuned, aggregated meta-cell dataset for downstream analyses such as marker discovery, trajectory inference, or cross-modal integration.
# In plain language

# Creates six global clusterings (0.5 → 4.0 resolution) and lets the user eyeball them on UMAP.
# Takes the finest global clustering (res 4) and, for each cluster individually, searches for an internal resolution that yields meta-cells of ~50–150 cells and stable gene counts.
# Assigns each cell to its chosen meta-cell, aggregates counts, and computes a representative UMAP coordinate for every meta-cell.
# Outputs:
# sample_clusters.rds – object with global clusters
# Hclu.rds, data_size.rds, meta_info.rds – audit trail of the resolution search
# add_meta.rds – per-cell mapping to meta-cells
# meta_cell.rds – final aggregated Seurat object ready for high-throughput analyses.








