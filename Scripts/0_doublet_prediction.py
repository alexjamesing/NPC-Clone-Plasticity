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
