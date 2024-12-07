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








