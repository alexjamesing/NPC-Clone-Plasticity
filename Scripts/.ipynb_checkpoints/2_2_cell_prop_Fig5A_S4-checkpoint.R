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
library(ggplot2)
library(forcats)
library(SingleCellExperiment)
library(stats)
library(ggrepel)
library(ggpubr)
library(ggalluvial)
library(ggfittext)
library(scales)


out="~/path_to_output/"
# 0.input data --------------------

data=readRDS("~/data.rds")

write_csv(data@meta.data,"~/metadata.csv")

mdata=data@meta.data

annot_color=setNames(c(color7[c(3,2,1,4,6)],color1[c(1,2,7,13,12,9,8,3)],color2[c(3,1)],color1[4]),sort(unique(mdata$annotations)))


# 1. GLM ------------------


t= mdata %>% filter(Red=="TRUE") %>%group_by(Condition,Mouse,annotations) %>%dplyr::summarise(counts=n())

prop_d=t %>% group_by(Condition,Mouse)%>% dplyr::add_count(wt=counts)


prop_d$Mouse=factor(prop_d$Mouse,levels=c("mouse3","mouse4","mouse12","mouse16",
                                "mouse1","mouse2","mouse17","mouse18"))



colnames(prop_d)=c("Conditions","Mouse","Cell_Types","cell_counts","Total_cell_counts_per_mouse")


## with data frame 

GLM_binomial=function(df,replicates,conditions,group){
  metad=unique(data.frame(replicate=replicates,condition=conditions))
  metad=setNames(metad$condition,metad$replicate)
  freq <- table(replicates, factor(df$annotations==group, levels=c(TRUE,FALSE)))
  m <- glm(freq ~ metad[rownames(freq)], family = "binomial")
  aov <- anova(m, test = "Chisq")
  res <- setNames(c(coef(m)[2], aov$Pr[2]), c("coef_glm","pval_aov"))
  d=data.frame(t(res)) %>% mutate(cells=group)
  d
}


Cell_groups=as.list(sort(unique(red_d$annotations)))
names(Cell_groups)=sort(unique(red_d$annotations))


red_d=mdata %>% filter(Red=="TRUE")
glm_res=do.call(rbind,lapply(Cell_groups,
                             function(x){GLM_binomial(df=red_d,replicates=red_d$groups,conditions=factor(red_d$Condition,levels=c("control","complemented")),group=x)}))


glm_res$padj_aov <- p.adjust(glm_res$pval_aov)
write_csv(glm_res,"~/glm_res_red_df.csv")


## plot boxplot -----


library(rstatix)

stat_text=glm_res %>% select(c(cells,padj_aov))
colnames(stat_text)[1]="Cell_Types"
stat_text=stat_text %>% mutate(group1="control",group2="complemented",.y.="prop")
stat_text=stat_text%>%mutate(p_sig=case_when(padj_aov <= 0.0001 ~"****",
                                             padj_aov <= 0.001 ~"***",
                                             padj_aov <= 0.01 ~"**",
                                             padj_aov <= 0.05 ~"*",
                                             T~"n.s"))



prop_d$prop=round(prop_d$cell_counts/prop_d$Total_cell_counts_per_mouse,3)
prop_d %>% group_by(Cell_Types)%>% summarise(y.position=max(prop)) -> y_pos
rownames(y_pos)=y_pos$Cell_Types
cbind(stat_text,y_pos[,"y.position"]) -> stat_text

# manually corrected y.pos ----------------

stat_text$y.position[which(stat_text$Cell_Types=="GLUT.")]=0.52
stat_text$y.position[which(stat_text$Cell_Types=="GLYC.")]=0.02
stat_text$y.position[which(stat_text$Cell_Types=="Glioblast")]=0.12

stat_text$y.position=stat_text$y.position+0.02


colorm=setNames(color1[c(1,3,4,5,11:14)],c("mouse1","mouse17","mouse18","mouse2","mouse3","mouse4","mouse12","mouse16"))

prop_d$Conditions=factor(prop_d$Conditions,levels=c("control","complemented" ))
color1 = met.brewer(name="Signac", n=14, type="discrete")
ggplot(data=prop_d,aes(x=Conditions,y=prop))+
  geom_boxplot(outliers = F)+
  geom_point(aes(color=Mouse),size=3,alpha=0.7)+
  #geom_point(aes(color=Mouse),size=3,position = "jitter",alpha=0.6)+
  scale_color_manual(values=colorm)+
  facet_wrap(~Cell_Types,ncol=4,scales = "free_y",shrink = F)+
  scale_y_continuous(labels=scales::number_format(accuracy = 0.01))+
  ylab("Proportion")+
  theme_classic()+
  stat_pvalue_manual(stat_text,label="p_sig",y.position = "y.position")

ggsave(paste0(out,"Cell_prop_GLMBio_res_red",".pdf"),width=12,height = 12)




