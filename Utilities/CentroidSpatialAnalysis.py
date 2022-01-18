''' Functions to analyze spatial distribution of centroids/points on the NMJ. 
For single images call analyzeCentroidFile(). For batch processing/aggregation call analyzeCentroidFolder().
Some example plotting code below.
Requires pointpats v2.1 (NOT >2.1)'''


import csv
import numpy as np
import pandas as pd
import os
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import libpysal
import pointpats as pp
import importlib
import glob
import CentroidAnalysisBatch_211213 as ca


def analyzeCentroidFile(file, exp="NA", n=1):
    filename = file.split("\\")[-1]
    table=pd.read_csv(file)
    table.columns = table.columns.str.split("_", expand=True)
    chans = list(set(table.columns.get_level_values(0)))
    # Combine centroids into one list to make alpha shape that describes them cumulatively
    allcents = table.unstack()\
        .unstack(level=1)\
        .reset_index()\
        .drop(["level_0", "level_1"], axis=1)\
        .dropna()
    
    # Get merged alpha shape
    alpha_shape=libpysal.cg.alpha_shape_auto(allcents.to_numpy())
    win = ca.makeWindow(alpha_shape)
    allgs={}
    for c in chans:
        N = table[c].dropna().size//2
        obs = pp.PointPattern(table[c].dropna().to_numpy())
        poisson = pp.PoissonPointProcess(win, N, 1, asPP=True).realizations[0]
        matern = ca.MaternPointProcess(win, N, 1, asPP=True).realizations[0]
        clustered = pp.PoissonClusterPointProcess(win, 
                                              N, 
                                              N/5, 
                                              3, 
                                              1, 
                                              asPP=True
                                              ).realizations[0]
        obsg = pp.G(obs, 40)._stat
        poisg = pp.G(poisson, 40)._stat
        materng = pp.G(matern, 40)._stat
        clustg = pp.G(clustered, 40)._stat
        
        dists = {"obs":obsg, "poisson":poisg, "spaced":materng, "clustered":clustg}
        
        # Make a dictionary where each key is a "EXP_CHAN_TYPE_N" string
        
        gs = {f"{exp}_{c}_{key}_{n}":val for key, val in dists.items()}
        allgs.update(gs)
        
    return pd.DataFrame(allgs)


def analyzeCentroidFolder(root_dir, version=216):
    g_dfs=[]
    # Get all folders containing centroid data
    if(root_dir[-1] != os.path.sep):
        root_dir += os.path.sep
    exp_dirs = glob.glob(root_dir+"*/")
    centroid_dirs = glob.glob(os.path.join(root_dir,"*/", f"*V{version}*/", "*Centroid*/"))
    
    # For each Experiment folder containing centroid data, iterate over csv files and analyze
    print("Running:")
    for c_dir in centroid_dirs:
        files = glob.glob(os.path.join(c_dir, "*.csv"))
        exp = c_dir.split("\\")[-2].split("_")[0]
        for n, file in enumerate(files):
            g_dfs.append(analyzeCentroidFile(file, exp=exp, n=n))
            if (1+n)%25 == 0: print('*')
            else: print('*', end='')
    return pd.concat(g_dfs, axis=1)

def organizeCentroidData(df):
    ''' Take the output dataframe of analyzeCentroidFolder and organize/prep for plotting:
    1. Create multi-index from compound column headings
    2. Average individual images
    3. clean up rogue columns
    4. Melt to d (index), channel (Chan), and data type (Type)'''
       
    df.columns = df.columns.str.split("_", expand=True)
    if ('Unnamed: 0', np.nan, np.nan,  np.nan) in df.columns:
        df = df.drop(('Unnamed: 0', np.nan, np.nan,  np.nan), axis=1).reset_index(col_level=3)
    else: df.reset_index(col_level=3)
    df.columns.names = ["Exp", "Chan", "Type", "N"]
    allCentGroup = df.groupby(axis=1, level=[1,2]).mean()
    allCentGroup = allCentGroup.drop(('', ''), axis=1)
    allCentMelt = allCentGroup.reset_index().melt(id_vars = 'index')
    return allCentMelt


def plotOverallMean(df, savepath=None):
    sns.lineplot(data=df, x='index', y="value", hue="Type", ci="sd")
    if savepath:
        plt.savefig(savepath)


def plotIndividualChannels(df, savepath=None):
    # Plot one graph for each channel
    #chans = list(set(df.Chan))
    chans = ["BRP", "CLC", "DAP160", "NWK", "DYN", "FASII"]
    rows, cols = 3,3
    fig, axes = plt.subplots(rows, cols)
    fig.suptitle("Centroid Spatial Patterning", fontsize=11)
    plt.tight_layout()
    ax_i=0
    for row in range(rows):
        for col in range(cols):
            if ax_i<len(chans):
                ax=axes[row, col]
                show = chans[ax_i]
                sns.lineplot(data=df.loc[df.Chan==show, :], x="index", y="value", hue="Type", 
                                ci="sd", ax=ax)
                ax.set_ylabel("")
                ax.set_xlabel("")
                ax.set_title(chans[ax_i])
                ax.tick_params(axis='x', bottom=False, labelbottom=False)
                ax.tick_params(axis='y', labelsize=10)
                handles, labels = ax.get_legend_handles_labels()
                ax.get_legend().remove()
                ax_i+=1
    plt.figlegend(handles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    while ax_i < rows*cols:
        ax=axes[ax_i//cols, ax_i%cols]
        fig.delaxes(ax)
        ax_i+=1
    
    if savepath:
        plt.savefig(savepath)