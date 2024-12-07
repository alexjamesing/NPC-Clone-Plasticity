library(Seurat)
library(pheatmap)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(Matrix)
library(R6)
library(dplyr)
library(tidyverse)
library(MetBrewer)
library(purrr)
library(forcats)
library(SingleCellExperiment)
library(stats)
library(ggrepel)
library(ggpubr)
library(ggalluvial)
library(ggfittext)
library(patchwork)
library(scales)
library(qs)
library(coin)
library(msigdbr)
library(ggdendro)


out="~/path_to_output/"

#0. Load data --------

## From SCENIC output

n_exp_df=read_csv(paste0(out,"files/n_exp_genes_cells.csv"))
shared_genes_df=read_csv("~/regulons2hallmark.csv")
d=read_csv("~/Mean_values_aucrss.csv")
df_TFs=read_csv("~/df_TFs_significant.csv") %>% filter(p_val <=0.01)
regulon_size=read_csv("~/regulon_size.csv")
regulon_all=read_csv("~/Regulons.csv")


#1. Odds Ratio-------


N=max(n_exp_df$n_exp_genes)

odd_table=data.frame(TFs=character(),Pathways=character(),odds_ratio=numeric(),
                     g=numeric(),p=numeric(),i=numeric(),p_val=numeric())

db=msigdbr("Mus musculus",category="H")

for( i in 1:nrow(regulon_size)){
  tf=regulon_size$TFs[i]
  print(tf)
  sub_g=regulon_all$Targets[which(regulon_all$TFs==tf)]
  g=length(sub_g)
  pathways=split(db$gene_symbol,db$gs_name)
  for(j in 1:length(pathways)){
    n=names(pathways)[j]
    sub_p=pathways[[n]]
    p=length(sub_p)
    
    I=length(intersect(sub_p,sub_g))
    
    odd=(I*(N-g-p+I))/((g-I)*(p-I))
    
    contingency_table <- matrix(c(I, g - I, p - I, N - g - p + I), nrow = 2)
    
    # to calculate the p-value for the association between the GRN and pathway gene sets
    fisher_test <- fisher.test(contingency_table)
    
    d=data.frame(TFs=tf,Pathways=n,odds_ratio=odd,g=g,p=p,i=I,p_val=fisher_test$p.value)
    odd_table=rbind(odd_table,d)
    
  }
  
}

write.csv(odd_table,paste0(out,"files/odd2hallmark_df.csv"),row.names = F)
odd_table=read_csv(paste0(out,"files/odd2hallmark_df.csv"))

#2. Filter non-associated and significant ----------------------------------------------------------------------
## filter out not significant and odds ratio <1
## filter out the intersecting gene sets <=2 

odd_table=odd_table %>% filter(p_val<=0.01 & odds_ratio>=1) %>% filter(i >2)


#3. Check shared gene sets----------------------------------------------------

data=qread("~/rds/Red_neurlin.qs")


## add modulescore -------------------------------------------------
for (i in 1:length(unique(odd_table$Pathways))){
  p=unique(odd_table$Pathways)[i]
  print(p)
  p_sub=odd_table[which(odd_table$Pathways==p),]
  
  for(j in 1:length(unique(p_sub$TFs))){
    tf=unique(p_sub$TFs)[j]
    tf_sub=shared_genes_df%>%dplyr::filter(TFs==tf & Pathways==p)
    f=paste(tf_sub$gene_sets[1])
    data=AddModuleScore(data,features=list(c(gsub(" ","",str_split(f,",")[[1]]))),name=paste0(tf,"_",p))
  }
  
}

module_plot=data@meta.data %>% select(c("annotations","Condition",37:160))
module_plot=pivot_longer(module_plot,cols=colnames(module_plot)[3:126],names_to = "Groups",values_to = "Score")

module_plot[c("TFs","Pathways")]=str_split_fixed(module_plot$Groups,pattern="_",n=2)



## tidy data  ---------------------

temp=list()

for(i in 1:length(unique(df_TFs$annotations))){
  n=unique(df_TFs$annotations)[i]
  print(n)
  targets=df_TFs %>%filter(annotations==n)
  targets=targets %>% select(c("TFs","Condition","auc_mean"))
  temp_d=targets %>% pivot_wider(names_from = Condition,values_from = auc_mean)
  temp_d$diff=temp_d$complemented-temp_d$control
  temp_d=temp_d %>% mutate(activity=ifelse(diff>0,"complemented","control"))
  meta_infor=unique(df_TFs %>%filter(annotations==n)%>%select(c("TFs","p_val","Regulon_size")))
  temp_d=merge(temp_d,meta_infor,all.x=T)
  temp_d$annotations=n
  temp[[i]]=temp_d
}

meta_infor=do.call(rbind,temp)

## ks-test prepare------------------------------------------------------------

## intersect genes >= 3 or >2

pdata3_lists=list()
#pdf(paste0(out,"plots/Significant_regulons2pathways.pdf"), width = 14, height = 8.5)
for(i in 1:length(unique(module_plot$Pathways))){
  pth=unique(module_plot$Pathways)[i]
  print(pth)
  pdata1=module_plot %>% filter(Pathways==pth)
  pdata1$Pathways=gsub("1","",pdata1$Pathways)
  pdata2=odd_table %>% filter(Pathways==gsub("1","",pth)) %>% filter(i>=3)
  pdata3=merge(pdata1,pdata2,by=c("TFs","Pathways"))
  if(nrow(pdata3)!=0){
    pdata3$groups=paste0(pdata3$TFs," (",pdata3$i,"g)")
    plot=ggplot(pdata3,aes(x=Score,fill=Condition))+
      geom_density(alpha=0.6)+
      scale_fill_manual(values=c("#430A5D","grey"))+
      theme_classic()+
      facet_grid(annotations~groups)+
      labs(title=gsub("1","",pth))
    #print(plot)
    pdata3_lists[[i]]=pdata3
  }
}
#dev.off()


pdata3_lists=Filter(negate(is.null),pdata3_lists)

data_all=do.call(rbind,pdata3_lists)


# shared_geneslist-------------------

kegg_db=msigdbr("Mus musculus",category = "H")
#kegg_db=kegg_db %>% filter(gs_subcat=="CP:KEGG" )
##unique(kegg_db$gs_name)
pathways=split(kegg_db$gene_symbol,kegg_db$gs_name)

reguon_list=split(regulon_all$Targets,regulon_all$TFs)
shared_genes_df=data.frame(TFs=character(),Pathways=character(),gene_sets=character(),size=numeric())

for( i in 1:length(reguon_list)){
  reg=reguon_list[[i]]
  TF=names(reguon_list)[i]
  for (j in 1:length(pathways)){
    p=pathways[[j]]
    overlap_genesets=intersect(reg,p)
    size=length(overlap_genesets)
    pathway=names(pathways)[j]
    
    if(is.null(p)){
      inter=data.frame(TFs=TF,Pathways=pathway,gene_sets=paste(as.character(sort(overlap_genesets)), collapse = ", "),size=0)
      shared_genes_df=rbind(shared_genes_df,inter)
    }else{
      inter=data.frame(TFs=TF,Pathways=pathway,gene_sets=paste(as.character(overlap_genesets), collapse = ", "),size=size)
      shared_genes_df=rbind(shared_genes_df,inter)
    }
  }
  
  
}

regulons2pathways=shared_genes_df%>%filter(size!=0)

write.csv(regulons2pathways,paste0(out,"files/regulons2path_df.csv"),row.names = F)



## ks-test------------------------------------------------------------


ks_df=data.frame(annotations=character(),
                 Pathways=character(),
                 TFs=character(),
                 statistics=numeric(),
                 ks_pvals=numeric(),
                 alternative=character())

for (i in 1:length(unique(data_all$Pathways))){
  pathways=unique(data_all$Pathways)[i]
  print(pathways)
  
  data_subp=data_all %>% filter(Pathways==pathways)
  
  for(j in 1:length(unique(data_subp$TFs))){
    tfs=unique(data_subp$TFs)[j]
    data_sub=data_subp%>%filter(TFs==tfs)
    dat_list=split(data_sub,f=~annotations)
    dat_ks=lapply(dat_list,function(x){res=ks.test(x$Score~x$Condition)})
    dat_ks=lapply(dat_ks,function(x){res=as.data.frame.list(x)})
    for(k in 1:length(names(dat_ks))){
      cell=names(dat_ks)[k]
      sub_dat_ks=data.frame(annotations=cell,
                            Pathways=pathways,
                            TFs=tfs,
                            statistics=dat_ks[[cell]]$statistic[1],
                            ks_pvals=dat_ks[[k]]$p.value[1],
                            alternative=dat_ks[[k]]$alternative[1])
      ks_df=rbind(ks_df,sub_dat_ks)
    }
  }
}


## plot ks_test---------------------


library(rstatix)
pdf(paste0(out,"plots/Significant_modules_regulons2pathways.pdf"), width = 16, height = 8.5)
for(i in 1:length(unique(module_plot$Pathways))){
  pth=unique(module_plot$Pathways)[i]
  print(pth)
  pdata1=module_plot %>% filter(Pathways==pth)
  pdata1$Pathways=gsub("1","",pdata1$Pathways)
  pdata2=odd_table %>% filter(Pathways==gsub("1","",pth)) %>% filter(i>=3)
  pdata3=merge(pdata1,pdata2,by=c("TFs","Pathways"))
  pdata3=merge(pdata3,ks_df,by.x=c("TFs","annotations","Pathways"),by.y=c("TFs","annotations","Pathways"))
  # filter based on p_val <=0.05
  pdata3=pdata3 %>% filter(ks_pvals<=0.01)
  
  ## filter out the unsignificant regulon activiy 
  
  pdata3=merge(pdata3,meta_infor,by=c("TFs","annotations"),all.x=T) %>% drop_na()
  
  
  
  
  
  if(nrow(pdata3)!=0){
    pdata3$groups=paste0(pdata3$TFs," (",pdata3$i,"g)")
    plot=ggplot(pdata3,aes(x=Score,fill=Condition,label=paste0("p: ",p_round(ks_pvals))))+
      geom_density(alpha=0.6)+
      scale_fill_manual(values=c("#430A5D","grey"))+
      theme_classic()+
      facet_grid(annotations~groups)+
      labs(title=gsub("HALLMARK_","",pdata3$Pathways[1]))+
      geom_text(x=0,y=0.2,check_overlap = T,size=4)
    print(plot)
    #pdata3_lists[[i]]=pdata3
  }
}
dev.off()


# plot regulons to pathways: Heatmap-------------

d=read_csv("~/Mean_values_aucrss.csv")


d=d %>% group_by(TFs)%>% mutate(z_auc=(auc_mean-mean(auc_mean))/sd(auc_mean))

d=d %>%select(c("annotations","TFs","Condition","z_auc"))%>% pivot_wider(names_from=Condition,values_from = z_auc)
d$auc_dist=d$complemented-d$control



df_merge=read_csv(paste0(out,"files/Genesets_intersect_significant.csv"))

Cells=unique(d$annotations)

for(i in 1:length(Cells)){
  cell=Cells[i]
  print(cell)
  
  active_df=meta_infor%>% filter(annotations==cell) %>%select(c("TFs","activity")) %>% unique()
  cell_df=d %>% filter(annotations==cell)
  
  ks_sub=ks_df %>% mutate(pvals=-log10(ks_pvals)) %>% 
    filter(annotations==cell & ks_pvals<=0.01)%>%
    select(c("TFs","Pathways","pvals")) %>% 
    drop_na() 
  ks_merge=merge(ks_sub,cell_df,by="TFs")
  ks_merge=merge(ks_merge,active_df)
  
  
  ks_mm=ks_df %>% mutate(pvals=-log10(ks_pvals)) %>% 
    filter(annotations==cell & ks_pvals<=0.01)%>%
    select(c("TFs","Pathways","pvals")) %>% 
    drop_na()%>%
    pivot_wider(names_from = TFs,values_from = pvals)%>%
    replace(is.na(.),0) %>% as.matrix()
  row_name=ks_mm[,1]
  row_name=gsub("HALLMARK_","",row_name)
  ks_mm=ks_mm[,2:dim(ks_mm)[2]]
  col_df=as.data.frame(active_df[which(active_df$TFs %in% colnames(ks_mm)),])
  row_colnames=col_df$TFs
  col_df=col_df %>% dplyr::select("activity")
  rownames(col_df)=row_colnames
  hp_color = list(
    activity = c(control="grey",complemented="#430A5D"))
  
  ks_mm=apply(ks_mm, 2, as.numeric)
  rownames(ks_mm)=row_name
  colpalette=colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#FDF0F0")))(49)
  
  ks_mm=ks_mm[,rownames(col_df)]
  ks_mm=ks_mm[apply(ks_mm,1,sum)>0,]
  
  hc_TFs=hclust(dist(t(ks_mm)),"ave")
  dhc1=as.dendrogram(hc_TFs)
  d_ks1=dendro_data(dhc1,typr="rectangle")
  
  label=d_ks1$labels$label
  
  ctrl=rownames(col_df)[which(col_df$activity=="control")]
  comp=base::setdiff(label,ctrl)
  TF_order=c(ctrl,comp)
  
  hc_P=hclust(dist(ks_mm),"ave")
  dhc2=as.dendrogram(hc_P)
  d_ks2=dendro_data(dhc2,typr="rectangle")
  
  P_order=d_ks2$labels$label
  
  
  
  ks_merge$Pathways=gsub("HALLMARK_","",ks_merge$Pathways)
  ks_merge$TFs=factor(ks_merge$TFs,levels=TF_order)
  ks_merge$Pathways=factor(ks_merge$Pathways,levels=rev(P_order))
  p1=ggplot(data=ks_merge,aes(x=TFs,y=Pathways,fill=pvals))+
    geom_point(aes(shape=activity),size=4)+
    scale_shape_manual(values = c(24,25),label=c("Up","Down"))+
    scale_fill_gradientn(colors = colorRampPalette(rev(c("#b2182b","#d6604d","#f4a582","#fddbc7","#FDF0F0")))(10)) +
    theme_minimal()+
    theme(axis.text.x=element_text(angle = 90))+
    labs(fill="-log10(P_vals)",shape="Activity in NesCre::R26-DTA")+xlab("")
  
  p2=ggplot(data=ks_merge,aes(x=TFs,y=annotations,fill=auc_dist))+
    geom_tile()+
    #geom_point(aes(shape=activity))+
    #scale_shape_manual(values = c(22,21))+
    scale_fill_gradient2(low="#091057",mid = "white",high = "#b2182b",midpoint = 0,breaks=seq(-1,1,0.25)) +
    theme_classic()+
    theme(axis.text.x=element_text(angle = 90))+
    labs(fill="Regulon Activity Differrence")+ylab("")+xlab("Regulons")
  
  plots=plot_grid(p1 ,p2,labels = NULL, ncol = 1,rel_heights=c(5,1),align = 'v', axis = 'tblr')
  print(plots)
  ggsave(paste0(out,"plots/",cell,"_pathways_kspval_regulonactiv.pdf"),width = 12,height = 10)
}


