import os
import numpy as np
import pandas as pd
from dtw import dtw

import matplotlib as mpl
import matplotlib.pyplot as plt

import skfda
import skfda.misc.hat_matrix as hm
from skfda.inference.anova import oneway_anova
from skfda.representation import FDataGrid, FDataBasis
from skfda.ml.clustering import FuzzyCMeans, KMeans
from skfda.preprocessing.smoothing import KernelSmoother

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# create default cmap
def create_cmap(cmap_name, cmap_list):
    cmap_o = plt.get_cmap(cmap_name)
    for i in range(cmap_o.N):
        rgba = cmap_o(i)
        cmap_list.append(mpl.colors.rgb2hex(rgba))
    return cmap_list

cmap = []
for name in ['Paired', 'tab20b', 'tab20c']:
    cmap = create_cmap(name, cmap)

class FDA:
    '''
    An object of Functional Data for CellDrift derived perturbation responses.

    Parameters
    ----------
    df_zscore
        concatenated z scores file
    df_meta
        concatenated metadata file for all comparisons
    '''
    def __init__(self, df_zscore, df_meta):
        df_zscore.values[np.isnan(df_zscore.values)] = 0
        self.zscores = df_zscore
        self.meta = df_meta
        self.n_comps = self.meta.shape[0]
        self.n_times = len(self.meta['time'])
        self.n_types = len(self.meta['cell_type'])
        self.n_perts = len(self.meta['perturbation'])

    def create_fd_genes(
        self, 
        genes,
        cell_type,
        perturbation                    
    ):
        '''
        Create functional data for all genes in one cell type and perturbation combination.
        '''
        # initial
        df_meta = self.meta.copy()
        df_zscore = self.zscores.copy()

        df_meta_sub = df_meta.loc[(df_meta['cell_type'] == cell_type) & (df_meta['perturbation'] == perturbation), :].copy()
        df_meta_sub = df_meta_sub.sort_values(['time'], ascending = True)
        comps = list(df_meta_sub.index.values)
        df_zscore = df_zscore.loc[:, comps]

        timepoints = df_meta_sub['time']
        data_matrix = []
        for gene in genes:
            data_matrix.append(df_zscore.loc[gene, :].values)

        fd = skfda.FDataGrid(
            data_matrix = data_matrix,
            grid_points = timepoints
        )
        return fd, genes

    def create_fd_perts(
        self, 
        gene, 
        cell_types, 
        perturbations
    ):
        '''
        Create functional data object for cell types and perturbations of one gene
        '''
        # filtering
        df_meta = self.meta.copy()
        df_zscore = self.zscores.copy()

        df_meta = df_meta.loc[df_meta['cell_type'].isin(cell_type), :].copy()
        df_meta = df_meta.loc[df_meta['perturbation'].isin(perturbations), :].copy()

        df_zscore = df_zscore.loc[[gene], :].copy()
        
        # order
        data_matrix = []
        timepoints = []
        fd_groups = []
        samples = []

        for perturbation in perturbations:
            # preparing data matrix
            fd_groups.append(group_dict[perturbation])

            df_meta_sub = df_meta.loc[(df_meta['perturbation'] == perturbation), :].copy()
            df_meta_sub = df_meta_sub.sort_values(['time'], ascending = True)
            comps = list(df_meta_sub.index.values)
            df_zscore_sub = df_zscore.loc[:, comps].copy()
            timepoint = df_meta_sub['time'].to_numpy()

            # create raw data matrix
            zvalues = df_zscore_sub.loc[gene, :].values

            # ref time
            data_matrix.append(zvalues)
            timepoints.append(timepoint)

        # create representation of functional data
        fd = skfda.FDataGrid(
            data_matrix = data_matrix,
            grid_points = timepoints
        )
        return fd, fd_groups

    def align_create_fd_perts(self,
                            gene, 
                            cell_types, 
                            perturbations,
                            ref_perturbation, 
                            step_pattern = 'symmetric2', 
                            alternative_pattern = 'asymmetric'):
        '''
        Perform alignment using dynamic time warping (dtw) and create functional data
        '''

        df_meta = self.meta.copy()
        df_zscore = self.zscores.copy()

        df_meta = df_meta.loc[(df_meta['perturbation'].isin(perturbations)) & (df_meta['cell_type'].isin(cell_types)), :].copy()
        df_zscore = df_zscore.loc[[gene], :].copy()
        
        # order
        data_matrix = []
        timepoints = []
        fd_groups = []
        samples = []

        perturbations = list(np.setdiff1d(perturbations, ref_perturbation))
        perturbations = [ref_perturbation] + perturbations

        for perturbation in perturbations:
            for cell_type in cell_types:
                # preparing data matrix
                fd_groups.append(group_dict[perturbation])
                samples.append(perturbation)

                df_meta_sub = df_meta.loc[(df_meta['perturbation'] == perturbation), :].copy()
                df_meta_sub = df_meta_sub.sort_values(['time'], ascending = True)

                comps = list(df_meta_sub.index.values)
                df_zscore_sub = df_zscore.loc[:, comps].copy()
                timepoint = df_meta_sub['time'].to_numpy()

                # create raw data matrix
                zvalues = df_zscore_sub.loc[gene, :].values

                # ref time
                if perturbation == ref_perturbation:
                    ref_zvalues = zvalues
                    ref_timepoint = timepoint
                    # print(len(ref_timepoint))

                    data_matrix.append(ref_zvalues)
                    timepoints.append(ref_timepoint)

                # query time
                else:
                    try:
                        # print(step_pattern)
                        alignment = dtw(zvalues, ref_zvalues, keep_internals = False, step_pattern = step_pattern, window_type = 'none')
                    except:
                        # print(alternative_pattern)
                        alignment = dtw(zvalues, ref_zvalues, keep_internals = False, step_pattern = alternative_pattern, window_type = 'none')

                    aligned_zvalues = np.array(zvalues)[alignment.index1]
                    aligned_timepoint = np.array(ref_timepoint)[alignment.index2] # it contains multiple query -> one ref alignment

                    # calcuate avg values for duplicate times
                    updated_zvalues, updated_timepoint = make_unique_alignments(aligned_zvalues, aligned_timepoint)
                    data_matrix.append(updated_zvalues)
                    timepoints.append(updated_timepoint)

        # create a FDataGrid data
        fd = skfda.FDataGrid(
            data_matrix = data_matrix,
            grid_points = timepoints[0] # since all timepoints are identical, we only need one here.
        )
        return fd, fd_groups

def fda_cluster(fd, genes, n_clusters = 20, seed = 42, output_folder = 'fda_celldrift/'):
    
    # run clustering methods
    kmeans = KMeans(n_clusters = n_clusters, random_state = seed)
    kmeans.fit(fd)
    pred_labels = kmeans.predict(fd)
    
    fuzzy_kmeans = FuzzyCMeans(n_clusters = n_clusters, random_state = seed)
    fuzzy_kmeans.fit(fd)
    pred_labels_fuzzy = fuzzy_kmeans.predict(fd)

    df_clusters = pd.DataFrame()
    df_clusters['genes'] = genes
    df_clusters['clusters_kmeans'] = pred_labels
    df_clusters['clusters_fuzzy'] = pred_labels_fuzzy
    df_clusters.to_csv(output_folder + 'fda_clusters.txt', sep = '\t')

    return df_clusters   

def draw_smoothing(fd, df_clutser, genes, 
                    method = 'nw', n_neighbors = 2, bandwidth = 1,
                    cluster_key = 'cluster_fuzzy',
                    output_folder = 'fda_celldrift/figures/smoothing_percluster/'):
    
    genes_all = df_cluster.index.values
    sub_indices = list(np.where(np.in1d(genes_all, genes))[0])

    fd_sub = fd[sub_indices]

    if method == 'knn':
        smoother = KernelSmoother(kernel_estimator=hm.KNeighborsHatMatrix(n_neighbors = n_neighbors))
    elif method == 'llr':
        smoother = KernelSmoother(kernel_estimator=hm.LocalLinearRegressionHatMatrix(bandwidth = bandwidth))
    else:
        smoother = KernelSmoother(kernel_estimator=hm.NadarayaWatsonHatMatrix(bandwidth = bandwidth))
    
    smoother.fit(fd_sub)
    fd_smoothed = smoother.transform(fd_sub)
    
    fig, ax = plt.subplots()
    fd_smoothed.plot(axes = ax, group = [0] * fd_smoothed.size, group_colors = ['grey'], alpha = 0.5)
    return fig


def draw_smoothing_clusters(fd_whole, df_cluster, genes, 
                            n_neighbors = 2, bandwidth = 1, 
                            cluster_key = 'clusters_fuzzy', 
                            output_folder = 'fda_celldrift/figures/'):

    clusters = df_cluster[cluster_key].unique()
    fd_mean = fd_whole.mean()
    
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder) 

    fig, ax = plt.subplots()
    fd_mean.plot(axes = ax, group = [0], group_colors = ['black'], linewidth = 5)
    fig.savefig(output_folder + 'avg_curve.png', bbox_inches = 'tight')  

    for color, cluster in zip(cmap, clusters):
        # initial
        df_subcluster = df_cluster.loc[df_cluster[cluster_key] == cluster, :].copy()
        genes_sub = df_subcluster['genes']
        sub_indices = list(np.where(np.in1d(genes, genes_sub))[0])
        fd = fd_whole[sub_indices]

        # Basic visualization
        fig, axes = plt.subplots(1,2,figsize = (8,3))
        fd.plot(axes = axes[0])
        fd.scatter(s = 0.5, axes = axes[1])
        fig.savefig(output_folder + 'raw_' + str(cluster) + '.pdf', bbox_inches = 'tight')

        # do smoothing
        # K-nearest neighbours kernel smoothing.
        knn = KernelSmoother(kernel_estimator=hm.KNeighborsHatMatrix(n_neighbors = n_neighbors))
        knn.fit(fd)
        knn_fd = knn.transform(fd)

        # Local linear regression kernel smoothing.
        llr = KernelSmoother(kernel_estimator=hm.LocalLinearRegressionHatMatrix(bandwidth = bandwidth))
        llr.fit(fd)
        llr_fd = llr.transform(fd)

        # # Nadaraya-Watson kernel smoothing.
        nw = KernelSmoother(kernel_estimator=hm.NadarayaWatsonHatMatrix(bandwidth = bandwidth))
        nw.fit(fd)
        nw_fd = nw.transform(fd)

        # visualize fitted curves
        fig, ax = plt.subplots()
        knn_fd.plot(axes = ax, group = [0] * fd.size, group_colors = [color], alpha = 0.5)
        fd_mean.plot(axes = ax, group = [0], group_colors = ['black'])
        fig.savefig(output_folder + 'knn_smoothing_' + str(cluster) + '.pdf', bbox_inches = 'tight')    

        fig, ax = plt.subplots()
        llr_fd.plot(axes = ax, group = [0] * fd.size, group_colors = [color], alpha = 0.5)
        fd_mean.plot(axes = ax, group = [0], group_colors = ['black'])
        fig.savefig(output_folder + 'LR_smoothing_' + str(cluster) + '.pdf', bbox_inches = 'tight')    

        fig, ax = plt.subplots()
        nw_fd.plot(axes = ax, group = [0] * fd.size, group_colors = [color], alpha = 0.5)
        fd_mean.plot(axes = ax, group = [0], group_colors = ['black'])
        fig.savefig(output_folder + 'NW_smoothing_' + str(cluster) + '.pdf', bbox_inches = 'tight')    

def fpca(fd):
    return 1

def anova_test(fd, fd_groups):
    mild_idx = [i for i in range(len(fd_groups)) if fd_groups[i] == 'infection_mild']
    severe_idx = [i for i in range(len(fd_groups)) if fd_groups[i] == 'infection_severe']

    fd_mild = fd[mild_idx]
    fd_severe = fd[severe_idx]

    mild_mean = np.mean([np.mean(fd_mild.data_matrix[i]) for i in range(fd_mild.shape[0])])
    severe_mean = np.mean([np.mean(fd_severe.data_matrix[i]) for i in range(fd_severe.shape[0])])
    gap = severe_mean - mild_mean

    val, pval = oneway_anova(fd_mild, fd_severe)
    return val, pval, gap
