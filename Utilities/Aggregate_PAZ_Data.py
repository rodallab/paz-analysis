import pandas as pd
from pathlib import Path
import os
import glob
import tifffile as tif
import PAZ_Processing as paz
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import kstest, kruskal, f_oneway as anova, ttest_ind as ttest, mannwhitneyu as mwu 
import scikit_posthocs as ph

''' Classes to aggregate across experiments the various types of PAZ analysis data, including:
Mesh intensity/localization/colocalization data stored in [X}-ALLTHECELLS.csv files
Radial profile data stored in [X]-ALLTHEPROFILES.csv files
Centroid distribution data stored per [img name]-centroids.csv file in [experiment]_centroids folders'''

# An obvious alternative to these relatively more complex/specific structures would be to simply have a dataframe with experiment as a column value...
# Indeed this is easier. Maybe make a class but just subclass df?

#root_folder = "Z:\\Current members\\DelSignore\\Coding_Projects\\Analyze_PAZ_Data\\TestData"
#version_num = 214
# glob to find all specific version of analysis folders across all experiments
#target_dirs = [path for path in Path(root_folder).rglob(f"Analysis-V{version_num}*")]
#experiment_labels = [dir.name.split("_")[1] for dir in target_dirs]
#for dir in target_dirs:
#    print(dir)

#for lab in experiment_labels:
#    print(lab)

class DataAggregate:
    def __init__(self, datasets=None, data_path=None, version_num=None, file_name=None, *args, **kwargs):
        assert datasets or data_path, "Must provide either a dataset or path etc to construct DataAggregate"
        if(datasets == None):
            self.datasets = self.load_data(data_path, version_num, file_name)
        else: 
            self.datasets = datasets
            
    def load_data(self, data_path, version_num, file_name):
        ''' 
        Load .csv datasets and store as dictionary with experiment number as keys
        Note: Tried using glob and Path.rglob to recursively search dirs but 
        this was *very* slow bc of large number of subdirs.
        '''
        exp_dirs = glob.glob(data_path+"*/")
        analysis_dirs = glob.glob(os.path.join(data_path,"*", f"*V{version_num}*/"))
        exp_labels = [dir.split("_")[0] for dir in exp_dirs]
        datasets = {}
        for exp, exp_dir in zip(exp_labels, exp_dirs):
            datasets[exp] = pd.read_csv(os.path.join(exp_dir, file_name))
        return datasets
    
    def filter(self, filter_cols, filter_vals):
        ''' filter all experiments by indicated columns.
        filter_cols is a list of strings containing column label(s) eg ['AZ_count']
        filter_vals is a list of strings containing filter expressions eg ['>0']'''
        filtered = {}
        filter_strings = [f"['{filter_col}']{filter_val}" for filter_col, filter_val in zip(filter_cols, filter_vals)]
        for experiment in self.datasets:
            filter_expression = [f"(self.datasets['{experiment}']{filter_string})" for filter_string in filter_strings]
            filter_expression = " & ".join(filter_expression) 
            bool_filter = eval(filter_expression)
            filtered[experiment] = self.datasets[experiment][bool_filter]
        return filtered
    
    def average_by_image(self):
        grouped = {}
        for experiment in self.datasets:
            grouped[experiment] = self.datasets[experiment].groupby(['Image']).mean()
        return grouped
    
    def aggregate_experiments(self):
        return pd.concat(self.datasets)
    
    
class MeshData(DataAggregate):
    def __init__(self, datasets=None, data_path=None, version_num=None, file_name=None, *args, **kwargs):
  
        super().__init__(datasets, data_path, version_num, file_name, *args, **kwargs)
    
    def average_by_image(self):
        averaged = super().average_by_image()
        return MeshData(datasets = averaged)
    
    def filter(self, filter_cols, filter_vals):
        filtered = super().filter(filter_cols, filter_vals)
        return MeshData(datasets = filtered)

class ProfileData(DataAggregate):
    def __init__(self, datasets=None, data_path=None, version_num=None, file_name=None, *args, **kwargs):
        super().__init__(datasets, data_path, version_num, file_name, *args, **kwargs)
        
    def average_by_image(self):
        averaged = super().average_by_image()
        return ProfileData(datasets = averaged)
    
    def filter(self, filter_cols, filter_vals):
        filtered = super().filter(filter_cols, filter_vals)
        return ProfileData(datasets = filtered)
    
class CentroidData(DataAggregate):
    def __init__(self, datasets=None, data_path=None, version_num=None, file_name=None, *args, **kwargs):
        super().__init__(datasets, data_path, version_num, file_name, *args, **kwargs)
        
    def average_by_image(self):
        averaged = super().average_by_image()
        return ProfileData(datasets = averaged)
    
    def filter(self, filter_cols, filter_vals):
        filtered = super().filter(filter_cols, filter_vals)
        return ProfileData(datasets = filtered)


def superplot_df(df, measurements):
    ''' Take an aggregated dataframe and organize specific measurements 
    into 'superplot' fashion by channels and measurements
    Row indices are hierarchical - Experiment - image'''
    if type(measurements) is not type([]):
        measurements = [measurements]
    splits = [col.split("_") for col in df.columns.values]
    unique_channels = list(set([colsplit[1] for colsplit in splits if len(colsplit)>1]))
    retrieve_cols = [f'{measure}_{channel}' for measure in measurements for channel in unique_channels]
    # figure out how to handle PCCs. Maybe have those without channels (ie don't underscore)?
    
    superplot = df[[*retrieve_cols]]
    return superplot

def aggregate_experiment_pixels(imgdir):
    ''' 
    Iterate through images in one experiment folder:
      Make mask and signed distance tranform
      Normalize pixel intensities per channel based on masked region
      create aggregated dataframe listing edm and channel intensities per pixel
      
      imdir - path of format [ExperimentLabel]_[c1-label]_[cn-label]...
      returns pd.DataFrame containing Experiment, Image, EDM, and normalized 
        channel intensity columns for all pixels in all images

      '''
    
    if(imgdir[-1] != os.path.sep):
        imgdir += os.path.sep
    exp = (imgdir.split(os.path.sep)[-2]).split("_")[0]
    chans = (imgdir.split(os.path.sep)[-2]).split("_")[1:4]
    alldata=pd.DataFrame(columns = ["Experiment", "Image", "EDM"] + chans)
    imgs = glob.glob(os.path.join(imgdir, "*.tif"))
    
    for img in imgs:
        imgn = tif.imread(img)
        mask = paz.segment(imgn, 1)
        sdm = paz.sdt(mask)
        imgn_norm = paz.norm_by_mask(imgn, mask, 1)
        fg_norm = imgn_norm[mask>0]
        fg_edm = sdm[sdm>-1]
        fg_norm = fg_norm.reshape(-1, imgn_norm.shape[-1])
        fg_edm = fg_edm.flatten()
        imgdf = pd.DataFrame(fg_norm)
        imgdf = imgdf.set_axis(chans, axis=1)
        imgdf["EDM"] = fg_edm
        imgdf["Image"]=os.path.basename(img)
        alldata = pd.concat([alldata, imgdf])
    alldata["Experiment"] = exp
    return alldata

def aggregate_all_pixels(rootpath):
    exp_dirs = glob.glob(rootpath+"*/")
    allexpdata = pd.DataFrame()
    for expdir in exp_dirs:
        print(expdir)
        tempdf = process_folder(expdir)
        allexpdata = pd.concat([allexpdata, tempdf])
    return allexpdata

def aggregate_csvs(root_dir, version, file_name, verbose=False):
    '''
    Open csv file and reorganize into regualrly formatted DataFrame ready to aggregate 
    '''
    
    if(root_dir[-1] != os.path.sep):
        root_dir += os.path.sep
    exp_dirs = glob.glob(root_dir+"*/")
    analysis_dirs = glob.glob(os.path.join(root_dir,"*", f"*V{version}*/"))
    alldata=pd.DataFrame(columns = ["Experiment"])
    for analysis_dir in analysis_dirs:   
        file = os.path.join(analysis_dir, file_name)
        if os.path.exists(file):
            exp = (analysis_dir.split(os.path.sep)[-2]).split("_")[-1]
            df = pd.read_csv(file)
            df["Experiment"] = exp
            alldata = pd.concat([alldata, df])
            if verbose:
                print(f"Loaded {df.shape[0]} meshes from experiment {exp}")
    return alldata

def subset_pairwiseData(df, usecols=None, dropchannels=None):
    ''' 
    Take a PAZ dataframe and return a subset of columns [usecols].
    Optional exclude channels in [dropchannels]
    '''
     
    if usecols:
        if type(usecols) is not type([]):
            usecols=[usecols]
    
        measure_cols = [col for col in df.columns.values if 
                        any([check in col for check in usecols])]
    else:
        measure_cols = df.columns.values
        measure_cols.remove("Experiment")
        measure_cols.remove("Image")
    
    # Start with Experiment and Image columns necessary for sorting later
    dfsub = df.loc[:, ["Experiment"]]
    dfsub[measure_cols] = df.loc[:, measure_cols]
    
    if dropchannels:
        if type(dropchannels) is not type([]):
            dropchannels=[dropchannels]
        
        drop_cols = [col for col in dfsub.columns.values if
                    any([check in col for check in dropchannels])]
        dfsub = dfsub.drop(drop_cols, axis = 1)
    return dfsub

def subset_data(df, usecols=None, dropchannels=None):
    ''' 
    Take a PAZ dataframe and return a subset of columns [usecols].
    Optional exclude channels in [dropchannels]
    '''
     
    if usecols:
        if type(usecols) is not type([]):
            usecols=[usecols]
    
        measure_cols = [col for col in df.columns.values if 
                        any([check in col for check in usecols])]
    else:
        measure_cols = df.columns.values
        measure_cols.remove("Experiment")
        measure_cols.remove("Image")
    
    # Start with Experiment and Image columns necessary for sorting later
    dfsub = df.loc[:, ["Experiment", "Image"]]
    dfsub[measure_cols] = df.loc[:, measure_cols]
    
    if dropchannels:
        if type(dropchannels) is not type([]):
            dropchannels=[dropchannels]
        
        drop_cols = [col for col in dfsub.columns.values if
                    any([check in col for check in dropchannels])]
        dfsub = dfsub.drop(drop_cols, axis = 1)
    return dfsub

def norm_data(df, norm_col=None, melt=True, keep_scale=False):
    '''
    Given dataframe of PAZ data, return normalized dataframe
    Normalization is per row, so if rows are meshes, normalization will be at level of mesh.
    Default behavior is to normalize each column with values against each other column.
    Optionally, can set norm_col to one column against which to normalize all other columns.
    '''
    if norm_col:
        measure_cols = [norm_col]
    else:
        measure_cols = df.columns
        measure_cols = measure_cols.drop(["Experiment", "Image"])
    norm = df.loc[:, ["Experiment", "Image"]]
    for num in measure_cols:
        for den in measure_cols:
            if num is not den:
                norm[f"{num}-{den}"]=df.loc[:, num]/df.loc[:, den]
    if melt:
        melted = norm.melt(id_vars = ["Experiment", "Image"])
        melted[['C1', 'C2']] = melted['variable'].str.split('-', expand=True)
        melted = melted.drop("variable", axis=1)
        return melted
    else:
        return norm
        
def superplot(df, data_order=None, box=True):
    melt = df.melt(id_vars = ["Experiment", "Image"])
    ImgMean = melt.groupby(["Experiment", "Image", "variable"]).mean().reset_index()
    ExpMean = melt.groupby(["Experiment", "variable"]).mean().reset_index()
    sp = sns.swarmplot(data=ImgMean, x="variable", y="value", hue="Experiment", size=5, order=data_order, zorder=0)
    sns.swarmplot(x="variable", y="value", hue="Experiment", size=12, edgecolor="k", linewidth=2, data=ExpMean, order=data_order, zorder=0)
    if box:
        sns.boxplot(data=ImgMean, x="variable", y="value",
            order=data_order, 
            showfliers=False, 
            whis=0, 
            boxprops={'facecolor':'None', 'linewidth':2.5, 'edgecolor':'k'},
            medianprops = {'linewidth':5, 'color':'k'},
            whiskerprops={'linewidth':0}, 
            ax=sp)
    sp.legend_.remove()
    plt.xticks(rotation=45)

def superplot_norm(df, data_order=None):
    ImgMean = df.groupby(["Experiment", "Image", "C1"]).mean().reset_index()
    ExpMean = df.groupby(["Experiment", "C1"]).mean().reset_index()
    sp = sns.swarmplot(data=ImgMean, x="C1", y="value", hue="Experiment", size=5, order=data_order)
    ax = sns.swarmplot(data=ExpMean, x="C1", y="value", hue="Experiment", size=12, edgecolor="k", linewidth=2, order=data_order)
    sp.legend_.remove()
    plt.xticks(rotation=45)

def compareGroups(df, verbose=True):
    ImgPost = df.melt()
    ImgPost = ImgPost.loc[-np.isnan(ImgPost.value), :]

    # Check whether outcomes are normally distributed
    ks = [kstest(df.loc[df[col]>-1, col], 'norm') for col in df]
    kpass = True if np.min(np.array(ks)[:,1])>0.05 else False
    test_name = "ANOVA" if kpass else "Kruskal-Wallis"
    posthoc = pd.DataFrame()
    
    if kpass:
        #Do anova if vars normally distributed
        test = anova(*[df.loc[-np.isnan(df[col]), col] for col in df])
        if test[1] < 0.05:
            # Do posthoc ttest if model is significant
            posthoc = ph.posthoc_ttest(ImgPost, group_col='variable', val_col='value', p_adjust='sidak')  
    else:
        # Do KW test to test if medians different
        test = kruskal(*[df.loc[:, col] for col in df], nan_policy='omit')
        if test[1] < 0.05:
            # Do posthoc dunns if model is significant
            posthoc = ph.posthoc_dunn(ImgPost, group_col='variable', val_col='value', p_adjust='sidak')

    if verbose:
        print(f"Performed {test_name}: test statistic = {test[0]} p-value:{test[1]}")
    
    return {'test':test_name, 'result':test, 'posthoc':posthoc}


def ttestl(df, verbose=True):
    '''Assess distributions of groups to compare. If normal, do ttest. If not, do MWU test.'''
    ImgPost = df.melt()
    ImgPost = ImgPost.loc[-np.isnan(ImgPost.value), :]

    # Check whether outcomes are normally distributed
    ks = [kstest(df.loc[df[col]>-1, col], 'norm') for col in df]
    kpass = True if np.min(np.array(ks)[:,1])>0.05 else False
    test_name = "ttest" if kpass else "Mann-Whitney"
    posthoc = pd.DataFrame()
    
    if kpass:
        #Do ttest if vars normally distributed
        test = ttest(*[df.loc[-np.isnan(df[col]), col] for col in df])
        
    else:
        # Do KW test to test if medians different
        test = kruskal(*[df.loc[:, col] for col in df], nan_policy='omit')
        if test[1] < 0.05:
            # Do posthoc dunns if model is significant
            posthoc = ph.posthoc_dunn(ImgPost, group_col='variable', val_col='value', p_adjust='sidak')

    if verbose:
        print(f"Performed {test_name}: test statistic = {test[0]} p-value:{test[1]}")
    
    return {'test':test_name, 'result':test, 'posthoc':posthoc}