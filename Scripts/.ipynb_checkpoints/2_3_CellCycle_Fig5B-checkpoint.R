library(Seurat)
library(pheatmap)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(Matrix)
library(R6)
library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(forcats)
library(SingleCellExperiment)
library(stats)
library(ggrepel)
library(ggridges)
library(ggpubr)
library(ggalluvial)
library(ggfittext)
library(scales)





# 0.input data --------------------

out="~/path_to_output/"
data=readRDS("~/data.rds")

# 1. cellcycle scores -------------------------

## S phase scores
library(biomaRt)

mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
s.genes <- cc.genes$s.genes
g2m.genes <-cc.genes.updated.2019$g2m.genes

convert_human_to_mouse <- function(gene_list){
  
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="human"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      mouse_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      for(mouse_gene in mouse_genes){
        output = append(output,mouse_gene)
      }
    }
  }
  
  return (output)
}

m.s.genes <- convert_human_to_mouse(s.genes)
m.g2m.genes <- convert_human_to_mouse(g2m.genes)
# Print the first 6 genes found to the screen
data=CellCycleScoring(data, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = FALSE)
data=readRDS(paste0(out,"rds/data_cellycle.rds"))
write_csv(data@meta.data,paste0(out,"files/metadata_S_score.csv"))

# 2. Check distribution --------------------------------


quantile(data$S.Score,prob=seq(0,1,0.05))["95%"]

#0%            5%           10%           15%           20%           25%           30% 
#-0.2322427835 -0.1283919176 -0.1137900537 -0.1037898361 -0.0955874680 -0.0882798131 -0.0815717015 
#35%           40%           45%           50%           55%           60%           65% 
#-0.0752189590 -0.0690021813 -0.0627453180 -0.0563428873 -0.0496116779 -0.0423356887 -0.0343825019 
#70%           75%           80%           85%           90%           95%          100% 
#-0.0251179544 -0.0140694984  0.0006543884  0.0241496606  0.0874513285  0.3374308783  1.2660504646 


quantile(data$G2M.Score,prob=seq(0,1,0.05))



## S phase scores ---------------------


## test the S phase score difference in progenitor cells ------------


target_cells=c("Glioblast","Neuroblast","Radial_glia", "Intermediate_progenitor","OPC", "Cajal_Retzius" )


df=data@meta.data[which(data@meta.data$annotations %in% target_cells),] 



ggplot(data@meta.data[which(data$annotations %in% target_cells),], aes(x =S.Score,fill=Condition,after_stat(scaled))) + 
  #geom_histogram(binwidth = 0.01)+
  #geom_density_ridges(aes(fill=Condition))+
  geom_density(aes(fill=Condition),alpha=0.5)+
  facet_grid(annotations~Red)+
  #scale_fill_manual(values =colorp)+
  theme( legend.position = "none")+
  theme_classic()+
  ylab("S Scores")

ggsave(paste0(out,"/plots/S_Score_histo.pdf"),width = 8,height = 16)



#K-S test ----------------------------------------

res_df=data.frame(annotations=character(),
                  population=character(),
                  statistics_D=numeric(),
                  p_value=numeric(),
                  alternative=character()
                  )

for (i in 1:length(target_cells)){
  n=target_cells[i]
  d=df %>% filter(annotations==n)
  dd=list(red_positive=d[which(d$Red=="True"),],red_negative=d[which(d$Red=="False"),])
  
  for(j in 1:length(dd)){
    class=names(dd)[j]
    temp=dd[[class]]
    res=ks.test(temp$S.Score~temp$Condition)
    s=data.frame(annotations=n,population=class,statistics_D=res$statistic,p_value=res$p.value,alternative=res$alternative)
    res_df=rbind(res_df,s)
  }

}

res_df=res_df %>% mutate(Red=case_when(population=="red_positive"~"True",T~"False"))
res_df[c("X_pos","Y_pos")]=c(rep(0.8,12),rep(0.2,12))
res_df$Condition="control"

p=ggplot(df,aes(x=S.Score,color=Condition))+
  stat_ecdf(geom="step")+
  facet_grid(annotations~Red)+
  theme_classic()
p+ geom_text(data=res_df,mapping=aes(x=X_pos,y=Y_pos,label=paste0("P value: ",round(p_value,2))),color="black")

ggsave(paste0(out,"/plots/S_Score_ecdf_withP.pdf"),width = 8,height = 12)

write_csv(res_df,paste0(out,"files/KStest_res.csv"))






