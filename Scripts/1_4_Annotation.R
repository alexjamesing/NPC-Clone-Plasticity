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
library(ggridges)
library(treedataverse)
library(ComplexHeatmap)
library(circlize)


# 1. Read data and preprocess  -----------------------------

# ref from La Manno et al., 2021 Nature
ref=readRDS("~/rds/ref_data_rd.rds")
query=readRDS(paste0("rds/meta_cell.rds"))
out="~/path_to_output/"

query[["percent.mt"]]=PercentageFeatureSet(query,pattern="^mt[-\\.]")
query=query %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 3000)%>%
  ScaleData(vars.to.regress = c("nFeature_RNA", "percent.mt")) %>%
  RunPCA(npcs = 50)



# 2. Region  -----------------------------

hv_features=intersect(VariableFeatures(ref), rownames(query))
anchors <- FindTransferAnchors(reference = ref, query = query, dims = 1:30, npcs = 30,reduction="cca",features = hv_features)

predictions_region <- TransferData(anchorset = anchors, refdata = ref$Region, dims = 1:30,weight.reduction = "cca")
write_csv(predictions_region,
          paste0("/home/l538g/workingf/brainbreaks/single_cell_BB_E17-5/Giulia_data_analysis/result_all_05082024/R_files/Results/1_4_Annotation/files/","Linnarsson_prediction_region.csv"))

query$Linnarsson_region_transfer<-predictions_region$predicted.id
query$Linnarsson_region_scores<-predictions_region$prediction.score.max
saveRDS(query,paste0("/home/l538g/workingf/brainbreaks/single_cell_BB_E17-5/Giulia_data_analysis/result_all_05082024/R_files/Results/1_4_Annotation/rds/","query_Linnarsson_region.rds"))





## markers in brains ------------

brain=c("Foxg1","Emx1","Pax6","Dlx2","Dlx5",# Forebrain
  "Otx2","En1","En2","Fgf17", #midbrain,
  "Gbx2", #Thalamus
  "Hoxb2","Hoxb5","Foxa2") # hindbrain

brain=setNames(brain,c(rep("Forebrain",5),rep("Midbrain",5),"Thalamus",rep("Hindbrain",3)))





# 3. SubClass level transfer ------------------------------

predictions_subclass <- TransferData(anchorset = anchors, refdata = ref$Subclass, dims = 1:30,weight.reduction = "cca")
write_csv(predictions_subclass,paste0(out,"Predictions_class.csv"))



# 4. General markers & manually check for whole brains ----------------------------

## Vascular ----------------------------
Endothelia=c("Cldn5","Igfbp7")
Microglia=c("Tmem119","P2ry12","Aif1","C1qb")
Pericyte=c("Pdgfrb","Cspg4")
Vascular=c(Endothelia,Microglia,Pericyte)

## RG Family-------------------------------
Radial_glia=c("Sox2","Pax6","Tnc","Ptprz1","Hopx","Ascl1","Egfr","Olig2","Mki67","Top2a")
Neuroblast=c("Npa3","Adgrb3","Fgf15","Wnt7b","Tnc","Hes6")
Glioblast=c("Aldh1l1","Efna5","Rspo1","Rspo2","Rspo3","Wnt1","Wnt9a","Hes1")
Oligodendrocyte=c("Olig2","Pdgfra")
Choroid_plexus=c("Htr2c","Enpp2")
Intermediate_Progenitor=c("Eomes","Neurog2","Sox2")


RGF=c(Radial_glia,Neuroblast,Glioblast,Oligodendrocyte,Choroid_plexus,Intermediate_Progenitor)
## Neurons -----------------------
Neu=c("Mki67","Nes","Dcx","Neurod2","Neurod6","Map2","Tubb3")
Gabaergic=c("Isl1","Slc32a1","Gad2","Gad1") 
Gluatmatergic=c("Neurod6","Slc17a7","Slc17a6","Slc17a5")
Glycinergic=c("Slc6a5")
Cholinergic=c("Chat","Slc18a3")
Serotonergic=c("Slc6a4","Tph2")
Dopaminergic=c("Slc6a3","Drd2","Th") 
Noradrenergic=c("Slc6a2","Dbh") 

Neurons=c(Neu,Gabaergic,Gluatmatergic,Glycinergic,Cholinergic,Serotonergic,Dopaminergic,Noradrenergic)
### Cortical Layer ------------------
Upper_Layer=c("Satb2","Cux2")
Deep_Layer=c("Bcl11b","Glra2","Grik3")


### Region ---------------
Telen=c("Foxg1")
Mid=c("Otx2","En1","En2","Fgf17")
Midbrain_stem=c("Gata3","Pax5","Sox14")
Hind=c( "Hoxb2","Hoxb5","Foxa2")
MHB=c("Fgf8","Fgf3") #midbrain-hindbrain boundary (MHB)
Cerebellum=c("Esrrb","Car8" )


Ven=c("Dlx1","DLX2","Dlx5","Six3","Sox6","Nr2f2")
Dor=c("Emx1","Pax6","Gli3","Eomes")

Region=c(Upper_Layer,Deep_Layer,Telen,Mid,Midbrain_stem,Hind,MHB,Cerebellum,Ven,Dor)
### Others ---------------------
Oth=str_to_title(c("RSPO3","LHX9","RELN","Col3a1","Rgs5","Tyrobp","Cldn5")) # non-telencephalon related
Microglia=c("P2ry12","Tmem119","Cx3cr1","Cd68")
Hypothalamus=c("Otp", "Vax1", "Nr5a1", "Pmch", "Sim1","Nkx2-1","Lhx6","Lhx8")
#General=str_to_title(c("Map2",Neu,For,Mid,"Gbx2",Hind,Ven,Dor,Oth,Inh,Exc))


Others=c(Oth,Microglia,Hypothalamus)
#FeaturePlot(query,features =Others,reduction = "umap",cols=c("lightgrey","#00ABB3","#497174"),pt.size=1,order = T,raster=F)



Marker_features=tibble::lst(Endothelia, Microglia,Pericyte,
                    Radial_glia,Neuroblast,Glioblast,Oligodendrocyte,Choroid_plexus,Intermediate_Progenitor,
                    Neu,Gabaergic,Gluatmatergic,Glycinergic,Cholinergic,Serotonergic,Dopaminergic,Noradrenergic,
                    Upper_Layer,Deep_Layer,Telen,Mid,Midbrain_stem,Hind,MHB,Cerebellum,Ven,Dor,
                    Oth,Microglia,Hypothalamus)

markers_df=data.frame(gene=character(),class=character())
for(i in 1:length(Marker_features)){
  n=names(Marker_features)[i]
  d=expand.grid(gene=Marker_features[[n]],class=n)
  markers_df=rbind(markers_df,d)
}



## manually annotate after checking marker expression, prediction scores, etc.

glu=as.character(sort(c(42,43,5:7,13,1,9,47,48,29,48,63,50,36,34,77,74,22,78,45,43,15,66,94,66,96,80,101,108,111,16,17,20,25,40,54,8,91)))
gaba=as.character(sort(c(3,1,58,81,12,51,32,27,30,24,38,28,37,19,10,102,73,61,11,75,76,26,97,112)))
mix_neuron=as.character(sort(c(14,18,31,99,106,110,109,60,71,0,4,2,21)))

Pericyte=as.character(sort(c(56,82,107)))
Cajal_Retzius=as.character(sort(c(104,70)))
Microglia=as.character(sort(c(89,98)))
Endothelial=as.character(sort(c(79,49)))

Radial_glia=as.character(sort(c(35,59,93,113)))
Glioblast=as.character(sort(c(68,41,88,39,23,44,100,57)))
Neuroblasts=as.character(sort(c(62,64,65,90,86,105,33,67,72,69,95)))
Choroid_plexus=as.character(sort(c(83)))
Ependymal_cell=as.character(sort(c(84)))
Intermediate_progenitor=as.character(sort(c(46)))
Oligodendrocyte_precursor_cell=as.character(sort(c(103,55,85)))
Erythrocyte=as.character(sort(c(52)))
Glycinergic=as.character(sort(c(53)))
dopaminergic=as.character(sort(c(92)))
Undetermined=as.character(sort(c(87)))
data@meta.data=data@meta.data %>% mutate(annotations=raw_cluster)
data@meta.data=data@meta.data %>% mutate(annotations=case_when(
  annotations %in% glu ~ "Gluta.",
  annotations %in% gaba ~ "GABA.",
  annotations %in% Pericyte ~ "Pericyte",
  annotations %in% mix_neuron ~ "Mixed_region",
  annotations %in% Cajal_Retzius ~ "Cajal_Retzius",
  annotations %in% Microglia ~ "Microglia",
  annotations %in% Endothelial ~ "Endothelial",
  annotations %in% Radial_glia ~ "Radial_glia",
  annotations %in% Glioblast ~ "Glioblast",
  annotations %in% Neuroblasts ~ "Neuroblasts",
  annotations %in% Choroid_plexus ~ "Choroid_plexus",
  annotations %in% Ependymal_cell ~ "Ependymal_cell",
  annotations %in% Intermediate_progenitor ~ "Intermediate_progenitor",
  annotations %in% Oligodendrocyte_precursor_cell ~ "Oligodendrocyte_precursor_cell",
  annotations %in% Erythrocyte ~ "Erythrocyte",
  annotations %in% Glycinergic ~ "Glycinergic",
  annotations %in% dopaminergic ~ "Dopaminergic",
  annotations %in% Undertimined ~ "Undetermined",
  
  
  T~"Others"
))
Idents(data)=as.factor(data$annotations)

saveRDS(data,paste0(out,"rds/data_annotatios.rds"))


# Summary (one-liner)
# This script takes an aggregated “meta-cell” Seurat object, transfers brain-region and subclass labels from a La Manno et al. reference atlas, defines marker-gene panels, hand-assigns each meta-cell to a broad cell-type class, and saves the fully annotated object for downstream analysis.

# Stage	Key code / actions	Detailed purpose
# Load libraries	library(Seurat) … library(circlize)	Brings in Seurat plus plotting and utility packages (some unused here).
# 1 Read and pre-process query	```r ref <- readRDS("~/rds/ref_data_rd.rds") # La Manno reference	
# query <- readRDS("rds/meta_cell.rds") # meta-cells		
# query <- query %>% NormalizeData() %>% FindVariableFeatures(3000) %>%		
#      ScaleData(vars.to.regress = c("nFeature_RNA","percent.mt")) %>%  
#      RunPCA(50)``` | Computes log-normalisation, 3 000 HVGs, scaling (regressing out library size + mito%), and 50 PCs on the meta-cells. |
# | 2 Label transfer – brain region | r hv <- intersect(VariableFeatures(ref),rownames(query)) anchors <- FindTransferAnchors(ref, query, features = hv, reduction = "cca") predictions_region <- TransferData(anchors, ref$Region) | CCA finds anchors between reference and query; transfers “Region” labels (Forebrain, Midbrain, etc.). Results written to CSV and added to query as Linnarsson_region_transfer, with prediction scores. |
# | 3 Label transfer – subclass | predictions_subclass <- TransferData(…, ref$Subclass) | Transfers finer neuronal subclasses; saved to CSV for inspection, but not stored in the object. |
# | 4 Define marker-gene panels | Large blocks that create named vectors for vascular cells, radial-glia family, neuronal sub-classes, cortical layers, regional markers, etc. A tidy tibble markers_df consolidates them. | Provides gene lists for manual FeaturePlots (code commented) and later annotation decisions. |
# | 5 Manual cluster-ID → cell-type mapping | Hard-coded vectors list the cluster IDs (from earlier resolution 4 clustering) that correspond to each major lineage (Glutamatergic, GABAergic, Pericyte, Radial glia, etc.). | Creates a new annotations column by case_when, mapping each cluster ID to its chosen class. |
# | 6 Set identities & save | r Idents(query) <- query$annotations saveRDS(query, paste0(out,"rds/data_annotatios.rds")) | Final object now carries annotations as the active identity—ready for differential expression, plotting, etc. |

# Notes / caveats

# The script writes several intermediate CSV/RDS files (Predictions_class.csv, query_Linnarsson_region.rds) so results can be inspected in a spreadsheet or re-loaded.
# The block that builds markers_df is preparatory; the gene lists are never actively visualised here (FeaturePlot line is commented).
# Undertimined in case_when is misspelled and will cause an error; likewise the final T ~ "Others" should be .default = "Others".
# Some loaded libraries (e.g., ComplexHeatmap, circlize, treedataverse) are not used in this script, suggesting copy-paste from a broader pipeline.
# Overall, the script adds automated atlas labels and curated, marker-based cell-type calls, returning a Seurat object (data_annotatios.rds) that is ready for downstream biological interpretation.


