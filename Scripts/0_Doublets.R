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


# in plain words
# The script tries to decide which cells in each single-cell library are probable doublets (i.e. droplets that contain two cells).
# It does that in two stages:

# Calculate a per-sample threshold on the “simulated‐doublet” scores produced by scrublet (or a similar tool).
# Tag every cell barcode as doublet / singlet by comparing its score with that threshold and write the result to disk.
# Step-by-step walk-through
# Line block	What happens	Why
# library(...)	Load tidyverse packages.	Data wrangling & plotting.
# 0 – read files
# files = list.files(... "*simulatedDoublet_scores.tsv")	Collect all files whose names contain simulatedDoublet_scores.tsv. Each file is one sample’s scrublet simulation output.	
# 1 – determine doublets		
# Loop over files	For each sample:
# • Read the table (sim_scores column).
# • Estimate the kernel density of scores (density(... adjust = 0.5)).
# • Find local minima in that density curve (valleys between the two peaks that represent singlets vs. doublets).
# • Store the x‐coordinate of every minimum together with the sample name in res_df.
# • Plot the density and mark the minima (QC).	The first minimum between the two peaks is commonly used as the score threshold that separates doublets from singlets.
# final_thres = res_df %>% filter(between(X_pos,0.1,0.5))	Keep only minima that fall in a plausible range (0.1–0.5).	
# manual correction
# final_thres_corrected = final_thres[c(1,7,16,29,42,52,60,68),]	Someone eyeballed the plots and hand-picked one “best” minimum per sample (rows 1, 7, 16, …).	
# Second plotting loop	Re-plot each density with its hand-chosen threshold (visual sanity check).	
# 2 – annotate doublets		
# files = list.files(... "*scrublet_results.tsv")	Read the main scrublet output for each sample (one row per barcode).	
# read_scrublet()	Attach sample name, rename the barcode column, return a data frame.	
# res = lapply(...) → do.call(rbind, ...)	Concatenate all samples into one big table.	
# df_all = merge(res, final_thres_corrected, all.x = TRUE)	Add the chosen threshold (X_pos) to every barcode in the corresponding sample.	
# df_all$barcodes = paste0(samples, "_", raw_barcodes)	Make a unique global barcode (sample + raw barcode).	
# write_csv()	Save Doublets_corrected_scores.csv which now contains for every cell:
# • its scrublet score
# • the sample-specific threshold
# • everything needed to flag it as a putative doublet.	
# What you get at the end
# Doublets_corrected_scores.csv
# A master table where you can label a barcode as a doublet if
# sim_score >= X_pos (the stored threshold).

# Caveats / assumptions
# The “true” threshold is chosen manually (final_thres_corrected), so the index list must be updated if file order changes.
# Density‐based valley picking assumes bimodality; if a sample is very clean or very full of doublets the method might fail.
