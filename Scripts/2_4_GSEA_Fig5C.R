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


out="~/path_to_output/"
data=readRDS("~/data.rds")
data=SplitObject(data,split.by = "annotations")


#1. single cell gsea-----------------------------

## gsea Hallmark ----------------



run_gsea=function(x){
  out_plot="/home/l538g/workingf/brainbreaks/single_cell_BB_E17-5/Giulia_data_analysis/result_all_05082024/R_files/Results/4_1_Pathway/"
  n=x$split_groups[1]
  print(n)
  Idents(x)=as.factor(x$Condition)
  wilcox_res=wilcoxauc(x,group_by = "Condition")
  
  print(head(wilcox_res))
  m_db=msigdbr("Mus musculus",category="H")
  #print(head(m_db))
  fgsea_sets=m_db %>% split(x=.$gene_symbol,f=.$gs_name)
  
  # Control
  control_sets=wilcox_res %>% filter(group=="control")
  wilcox_res_genes=control_sets %>% arrange(desc(logFC))%>% dplyr::select(feature,logFC)
  ranks=deframe(wilcox_res_genes)
  fgsea_res=fgseaMultilevel(fgsea_sets,stats=ranks,eps=0)
  write_csv(fgsea_res,paste0(out_plot,"files/",n,"_control_H",".csv"))
  quan_90=quantile(fgsea_res$NES,c(0.1,0.2,0.5,0.9,0.95,1),na.rm = T)["90%"]
  quan_10=quantile(fgsea_res$NES,c(0.1,0.2,0.5,0.9,0.95,1),na.rm = T)["10%"]
  res=fgsea_res %>%arrange(desc(NES)) %>%arrange(padj)%>%head(50)
  ggplot(res %>% dplyr::filter(padj<=0.01) %>% head(15),aes(x=reorder(pathway,NES),y=NES))+
    geom_col(fill="lightgrey")+
    coord_flip()+
    theme_classic()+
    xlab("Pathways")
  ggsave(paste0(out_plot,"plots/",n,"_control_H",".pdf"),width = 16,height = 12)
  
  # Complemented
  
  complemented_sets=wilcox_res %>% filter(group=="complemented")
  wilcox_res_genes=complemented_sets %>% arrange(desc(logFC))%>% dplyr::select(feature,logFC)
  ranks=deframe(wilcox_res_genes)
  fgsea_res=fgseaMultilevel(fgsea_sets,stats=ranks,eps=0)
  write_csv(fgsea_res,paste0(out_plot,"files/",n,"_complemented_H",".csv"))
  quan_90=quantile(fgsea_res$NES,c(0.1,0.2,0.5,0.9,0.95,1),na.rm = T)["90%"]
  quan_10=quantile(fgsea_res$NES,c(0.1,0.2,0.5,0.9,0.95,1),na.rm = T)["10%"]
  res=fgsea_res %>%arrange(desc(NES)) %>%arrange(padj)%>%head(50)
  ggplot(res %>% dplyr::filter(padj<=0.01) %>% head(15),aes(x=reorder(pathway,NES),y=NES))+
    geom_col(fill="lightgrey")+
    coord_flip()+
    theme_classic()+
    xlab("Pathways")
  ggsave(paste0(out_plot,"plots/",n,"_complemented_H",".pdf"),width = 16,height = 12)
}






# 5. Plot -----------------------------
### 1. gsea res -----------------------------------------------
path=out
f=list.files(out,pattern = "^True*_")
f=f[which(grepl("complemented",f))]

gsea_res=list()

for(i in 1:length(f)){
  print(f[i])
  n=gsub("True_","",f[i])
  n=gsub("_H.csv","",n)
  t=read_csv(paste0(path,f[i]),trim_ws = T) %>% mutate(cells=n) %>% filter(padj <=0.01)
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

ggsave(paste0(out,"plots/gsea_res_nc.pdf"),units = "in",width = 8,height = 6)

write.csv(gsea_df,paste0(out,"files/gsea_df.csv"),row.names = F)
