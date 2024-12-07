library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)



## 0. read files--------------
files=list.files("~/0_Doublet_prediction",pattern="*simulatedDoublet_scores.tsv")

## 1. determine doublets--------------------
columns = c("X_pos","Y_pos","samples") 
res_df = data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(res_df) = columns

for (i in 1:length(files)){
  f=files[i]
  name=sapply(str_split(f,pattern="_result_new_",n=2),'[',n=1)
  print(f)
  path="~/Results/0_Doublet_prediction/"
  d=read_delim(paste0(path,f),delim="\t")
  
  q=quantile(d$sim_scores, probs=seq(0,1,0.2))
  dd=density(d$sim_scores,adjust=0.5)

  thres <- NULL
  for (j in 2:(length(dd$y)-1)) {
    if (dd$y[j-1] >= dd$y[j] & dd$y[j] <= dd$y[j+1]) {
      thres <- cbind(thres, c(dd$x[j], dd$y[j]))
    }
  }
  df=data.frame(t(thres))
  colnames(df)=c("X_pos","Y_pos")
  df$samples=name
  
  res_df=rbind(res_df,df)
  plot(dd)
  points(thres[1,], thres[2,])
}

final_thres=res_df %>% filter(between(X_pos,0.1,0.5))





# manually filter 

final_thres_corrected=final_thres[c(1,7,16,29,42,52,60,68),]
for (i in 1:length(files)){
  f=files[i]
  name=sapply(str_split(f,pattern="_result_new_",n=2),'[',n=1)
  print(f)
  path="~/0_Doublet_prediction/"
  d=read_delim(paste0(path,f),delim="\t")
  
  dd=density(d$sim_scores,adjust=0.5)
  plot(dd)
  points(final_thres_corrected[i,1], final_thres_corrected[i,2])
}


## 2. annotate doublets --------------------------------

files=list.files("~/0_Doublet_prediction",pattern="*scrublet_results.tsv")

read_scrublet=function(x){
  name=sapply(str_split(x,pattern="_result_new_",n=2),'[',n=1)
  path="~/0_Doublet_prediction/"
  d=read_delim(paste0(path,x),delim="\t")
  d$samples=name
  colnames(d)[1]="raw_barcodes"
  d
}

res=lapply(files,read_scrublet)
res=do.call(rbind,res)

df_all=merge(res,final_thres_corrected,all.x=T)

df_all$barcodes=paste0(df_all$samples,"_",df_all$raw_barcodes)
write_csv(df_all,"~/0_Doublet_prediction/Doublets_corrected_scores.csv")



