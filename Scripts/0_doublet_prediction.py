import scrublet as scr
import sys
import os
import scipy.io
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
import numba
import numba.typed
import glob



input_dir="~/path_to_cellranger_outputs"
data_path=glob.glob(input_dir)
for i in range(len(data_path)):
    data_path[i]=data_path[i].split("/",8)[8]
Samples=["Samples"+ str(i) for i in range(1,9)]


data_dict=dict(zip(data_path,Samples))
for i in data_dict:
    data_dict[i]="~/path_to_cellranger_outputs/"+i+"/outs/filtered_feature_bc_matrix"
    print(data_dict[i])
    
    
outdir="~/0_Doublet_prediction/"
for i in data_dict:
    print(i)
    barcode_indir=data_dict[i]+"/barcodes.tsv.gz"
    barcode_outdir=data_dict[i]+"/barcodes.tsv"
    cmd="gunzip -c "+barcode_indir+ " > " + barcode_outdir
    print(cmd)
    os.system(cmd)

    input_dir = data_dict[i]
    counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
    barcodes_df = pd.read_csv(barcode_outdir, delimiter='\t',header=None)
    print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
    print('Number of cells in barcode list: {}'.format(len(barcodes_df)))
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.1)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                              min_cells=3, 
                                                              min_gene_variability_pctl=85, 
                                                              n_prin_comps=30)
    scrub.plot_histogram();
    plt.savefig(outdir+i+'_doublet_score_histogram.png')
    print('Running UMAP...')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    print('Done.')
    scrub.plot_embedding('UMAP', order_points=True);
    plt.savefig(outdir+i+'_UMAP.png')

    results = pd.Series(scrub.predicted_doublets_, name="scrublet_DropletType")
    scores = pd.Series(scrub.doublet_scores_obs_, name="scrublet_Scores")
    dataframe = pd.concat([barcodes_df, results, scores], axis=1)
    dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(True, "doublet")
    dataframe.scrublet_DropletType = dataframe.scrublet_DropletType.replace(False, "singlet")

    dataframe.to_csv(outdir+i+'_scrublet_results.tsv', sep = "\t", index = False)


    ### Make summary of singlets and doublets and write to file ###
    summary = pd.DataFrame(dataframe.scrublet_DropletType.value_counts())
    summary.index.name = 'Classification'
    summary.reset_index(inplace=True)
    summary = summary.rename({'scrublet_DropletType': 'Droplet N'}, axis=1)

    summary.to_csv(outdir+i+'_scrublet_summary.tsv', sep = "\t", index = False) 


# High-level summary
# The script loops over multiple 10x Genomics “filtered_feature_bc_matrix” folders, runs Scrublet to predict doublets, and writes per-sample outputs:

# histogram of Scrublet scores (*_doublet_score_histogram.png)
# UMAP embedding coloured by predicted singlet/doublet (*_UMAP.png)
# a TSV of every barcode with its Scrublet score and call (*_scrublet_results.tsv)
# a short TSV with the singlet/doublet counts (*_scrublet_summary.tsv)
# All results are saved in ~/0_Doublet_prediction/.

# Step-by-step behaviour
# Code section	What it does
# Imports	Loads Scrublet, SciPy, pandas, UMAP, etc.; sets Matplotlib to non-interactive AGG backend so plots render off-screen.
# Collect sample dirs
# input_dir="~/path_to_cellranger_outputs"
# data_path = glob.glob(input_dir)	Suppose ~/path_to_cellranger_outputs contains one sub-directory per 10x sample. glob.glob gathers the directory names. (split("/",8)[8] later keeps only the last folder name.)
# Make name mapping	Creates list Samples = ["Samples1", … "Samples8"] and zips it with folder names to get data_dict, e.g. { "sampleA":"~/…/sampleA/outs/filtered_feature_bc_matrix", … }.
# Decompress barcodes	For each sample path:
# ▪ gunzip barcodes.tsv.gz → plain barcodes.tsv in the same folder (via os.system).
# Load count matrix & barcodes	Reads the sparse matrix.mtx.gz (10x counts; transpose so rows = cells). Reads barcode list into barcodes_df.
# Run Scrublet
# expected_doublet_rate=0.1	Calculates doublet scores and a boolean call (predicted_doublets). Quality/trimming parameters: min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30.
# Plots	• Histogram of doublet scores (*_doublet_score_histogram.png).
# • UMAP embedding of Scrublet manifold (*_UMAP.png).
# Save per-cell results	Concatenates barcodes, predicted class (boolean → “doublet”/“singlet”), and score into one table and writes *_scrublet_results.tsv.
# Save summary counts	Counts singlets vs. doublets and writes *_scrublet_summary.tsv.




