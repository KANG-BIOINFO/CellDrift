import os
import numpy as np
import pandas as pd
from dtw import dtw

import matplotlib as mpl
import matplotlib.pyplot as plt

import skfda
# import skfda.misc.hat_matrix as hm
from skfda.inference.anova import oneway_anova
from skfda.inference.hotelling import hotelling_t2
from skfda.representation import FDataGrid, FDataBasis
from skfda.ml.clustering import FuzzyCMeans, KMeans
# from skfda.preprocessing.smoothing import KernelSmoother
from skfda.preprocessing.smoothing.kernel_smoothers import KNeighborsSmoother, NadarayaWatsonSmoother, LocalLinearRegressionSmoother
from skfda.preprocessing.dim_reduction.feature_extraction import FPCA

from statsmodels.stats.multitest import fdrcorrection

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams.update({'font.size': 18})

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

#  remove duplicate time alignments
def make_unique_alignments(aligned_zvalues, aligned_timepoint):
    # get average z scores for multiple query -> one ref alignment
    df_align = pd.DataFrame({'aligned_timepoint': aligned_timepoint, 'aligned_zvalues': aligned_zvalues})
    df_avg = pd.DataFrame(df_align.groupby(by = ['aligned_timepoint']).mean()['aligned_zvalues'])

    updated_timepoint = list(df_avg.index.values)
    updated_zvalues = list(df_avg['aligned_zvalues'])
    return (updated_zvalues, updated_timepoint)

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
        self.n_times = len(self.meta['time'].unique())
        self.n_types = len(self.meta['cell_type'].unique())
        self.n_perts = len(self.meta['perturbation'].unique())
        self.n_reps = len(self.meta['rep'].unique())
        self.genes = df_zscore.index.values
        self.combs = df_meta.index.values

    def create_fd_genes(
        self, 
        genes,
        cell_type,
        perturbation,
        rep = None, 
    ):
        '''
        Create functional data for all genes in one cell type and perturbation combination.
        '''
        # initial
        df_meta = self.meta.copy()
        df_zscore = self.zscores.copy()

        if rep == None:
            rep = 1 if 1 in df_meta['rep'].unique() else 'rep1'

        df_meta_sub = df_meta.loc[(df_meta['cell_type'] == cell_type) & (df_meta['perturbation'] == perturbation), :].copy()
        df_meta_sub = df_meta_sub.loc[df_meta_sub['rep'] == rep, : ].copy()

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

    def align_create_fd_perts(self,
                            gene, 
                            cell_type, 
                            perturbations,
                            ref_perturbation, 
                            step_pattern = 'symmetric2', 
                            alternative_pattern = 'asymmetric',
                            reps = None):
        '''
        Perform alignment using dynamic time warping (dtw) and create functional data
        '''

        df_meta = self.meta.copy()
        df_zscore = self.zscores.copy()
        if reps == None:
            reps = self.n_reps

        df_meta = df_meta.loc[(df_meta['perturbation'].isin(perturbations)) & (df_meta['cell_type'] == cell_type), :].copy()
        df_zscore = df_zscore.loc[[gene], :].copy()
        
        # order
        data_matrix = []
        timepoints = []
        fd_groups = []
        samples = []

        if ref_perturbation == None:
            ref_perturbation = perturbation[0]

        perturbations = list(np.setdiff1d(perturbations, ref_perturbation))
        perturbations = [ref_perturbation] + perturbations
        
        for rep in range(1, reps + 1):
            for perturbation in perturbations:
                # preparing data matrix
                fd_groups.append(perturbation) # a better way of definition?
                samples.append(perturbation)

                df_meta_sub = df_meta.loc[(df_meta['perturbation'] == perturbation) & (df_meta['rep'] == rep), :].copy()
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

    def check_timeSeries_consistency(
        self,
        perturbations
    ):
        df_meta = self.meta.copy()
        df_zscore = self.zscores.copy()

        for idx, perturb in enumerate(perturbations):
            df_meta_sub = df_meta.loc[df_meta['perturbation'] == perturb, :].copy()
            times = df_meta_sub['time'].values
            if idx == 0:
                ref_times = times
            else:
                if (len(times) != len(ref_times)) or (not np.all(ref_times == times)):
                    return False
        return True

    def create_fd_types_perts(
        self, 
        gene, 
        cell_type, 
        perturbations,
        reps = None,
        ref_perturbation = None, 
        step_pattern = 'symmetric2', 
        alternative_pattern = 'asymmetric',
    ):
        '''
        Create functional data object for cell types and perturbations of one gene
        '''
        if self.check_timeSeries_consistency(perturbations):
            # filtering
            df_meta = self.meta.copy()
            df_zscore = self.zscores.copy()

            # need to add an alert to ensure timepoints are matching !!!
            perturbations = [perturbations] if type(perturbations) == 'str' else perturbations

            if reps == None:
                reps = self.meta['rep'].unique()
            df_meta = df_meta.loc[df_meta['cell_type'] == cell_type, :].copy()
            df_meta = df_meta.loc[df_meta['perturbation'].isin(perturbations), :].copy()
            df_meta = df_meta.loc[df_meta['rep'].isin(reps), :].copy()
            df_meta['comb'] = [i + '-' + j for (i, j) in zip(df_meta['cell_type'], df_meta['perturbation'])]

            df_zscore = df_zscore.loc[[gene], :].copy()
            
            # order
            data_matrix = []
            fd_groups = []
            samples = []
            timepoints_ref = None

            for comb in df_meta['comb'].unique():
                for rep in reps:
                    # preparing data matrix
                    fd_groups.append(comb)

                    df_meta_sub = df_meta.loc[(df_meta['comb'] == comb), :].copy()
                    df_meta_sub = df_meta_sub.loc[(df_meta_sub['rep'] == rep), :].copy()
                    
                    df_meta_sub = df_meta_sub.sort_values(['time'], ascending = True)
                    comps = list(df_meta_sub.index.values)
                    df_zscore_sub = df_zscore.loc[:, comps].copy()
                    timepoints = df_meta_sub['time'].to_numpy()
                    
                    if np.array_equal(timepoints_ref, None) and np.array_equal(timepoints, timepoints_ref):
                        raise Exception(
                            'Time series in different cell types / perturbations are not consistent.\
                            Please check function align_create_fd_perts to make sure time series are consistent.'
                        )
                    timepoints_ref = timepoints.copy()

                    # create raw data matrix
                    zvalues = df_zscore_sub.loc[gene, :].values
                    data_matrix.append(zvalues)

            fd = skfda.FDataGrid(
                data_matrix = data_matrix,
                grid_points = timepoints
            )
        
        else:
            fd, fd_groups = self.align_create_fd_perts(gene = gene, 
                                                        cell_type = cell_type, 
                                                        perturbations = perturbations,
                                                        ref_perturbation = ref_perturbation, 
                                                        step_pattern = step_pattern, 
                                                        alternative_pattern = alternative_pattern,
                                                        reps = reps)
        return fd, fd_groups
        

def fda_cluster(fd, genes, n_clusters = 20, seed = 42, suffix = '', output_folder = 'Temporal_CellDrift/'):
    
    '''
    Find functional data analysis clusters.

    :param fd: FDA object.
    :param genes: genes used to run clustering algorithms.
    :param n_clusters: number of clusters. Default is 20.
    :param seed: random seed. Default is 42.
    :param suffix: suffix in the name of output file.
    :param output_folder: location of output folder.
    '''

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
    df_clusters.to_csv(output_folder + 'FDA_clusters' + suffix + '.txt', sep = '\t')

    return df_clusters   

def draw_smoothing(fd, method = 'nw', n_neighbors = 2, bandwidth = 1,
                    output_folder = 'Figures_CellDrift/', fig_name_suffix = None):
    '''
    Draw a smoothed curve for functional data.

    :param fd: FDA object, such as the FDA object for all cell-type+perturbation combinations in a specific gene.
    :param method: method of smoothing, including knn, llr and nw. Default is nw
    :param n_neighbors: number of neighbors for Kernel Smoothing method KNeighborsHatMatrix. Default is 2.
    :param bandwidth: bandwidth of Kernel Smoothing method LocalLinearRegressionHatMatrix. Default is 1.
    :param output_folder: output folder for smoothing figures.
    :param fig_name_suffix: suffix in the output figure name.
    '''

    if method == 'knn':
        # smoother = KernelSmoother(kernel_estimator=hm.KNeighborsHatMatrix(n_neighbors = n_neighbors))
        smoother = KNeighborsSmoother(smoothing_parameter = 2)
    elif method == 'llr':
        # smoother = KernelSmoother(kernel_estimator=hm.LocalLinearRegressionHatMatrix(bandwidth = bandwidth))
        smoother = LocalLinearRegressionSmoother(smoothing_parameter = 2)
    else:
        # smoother = KernelSmoother(kernel_estimator=hm.NadarayaWatsonHatMatrix(bandwidth = bandwidth))
        smoother = NadarayaWatsonSmoother(smoothing_parameter = 2)
    
    smoother.fit(fd)
    fd_smoothed = smoother.transform(fd)
    
    fig, ax = plt.subplots()
    fd_smoothed.plot(axes = ax, group = [0] * fd_smoothed.size, group_colors = ['grey'], alpha = 0.5)

    if fig_name_suffix == None:
        fig_name_suffix = ''

    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    fig.savefig(output_folder + 'smoothing_' + method + '_' + fig_name_suffix + '.png', bbox_inches = 'tight')
    return fig

def draw_smoothing_clusters(fd_whole, df_cluster, 
                            n_neighbors = 2, bandwidth = 1, 
                            cluster_key = 'clusters_fuzzy', 
                            output_folder = 'Temporal_CellDrift/figures/',
                            ylim = None):
    '''
    Draw smoothed curves for all CellDrift derived clusters.

    :param fd_whole: FDA object, such as the FDA object for all genes in a specific cell-type+perturbation combination.
    :param df_cluster: dataframe of cluster output from function fda_cluster.
    :param n_neighbors: number of neighbors for Kernel Smoothing method KNeighborsHatMatrix. Default is 2.
    :param bandwidth: bandwidth of Kernel Smoothing method LocalLinearRegressionHatMatrix. Default is 1.
    :param cluster_key: the key of cluster in the data frame df_cluster
    :param output_folder: output folder for smoothing figures.
    '''

    clusters = df_cluster[cluster_key].unique()
    fd_mean = fd_whole.mean()
    genes = list(df_cluster['genes'])
    
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
        # knn = KernelSmoother(kernel_estimator=hm.KNeighborsHatMatrix(n_neighbors = n_neighbors))
        knn = KNeighborsSmoother(smoothing_parameter = 2)
        knn.fit(fd)
        knn_fd = knn.transform(fd)

        # Local linear regression kernel smoothing.
        # llr = KernelSmoother(kernel_estimator=hm.LocalLinearRegressionHatMatrix(bandwidth = bandwidth))
        llr = llr(smoothing_parameter = 2)
        llr.fit(fd)
        llr_fd = llr.transform(fd)

        # # Nadaraya-Watson kernel smoothing.
        # nw = KernelSmoother(kernel_estimator=hm.NadarayaWatsonHatMatrix(bandwidth = bandwidth))
        nw = NadarayaWatsonSmoother(smoothing_parameter = 2)
        nw.fit(fd)
        nw_fd = nw.transform(fd)

        # visualize fitted curves
        fig, ax = plt.subplots()
        knn_fd.plot(axes = ax, group = [0] * fd.size, group_colors = [color], alpha = 0.5)
        fd_mean.plot(axes = ax, group = [0], group_colors = ['black'])
        plt.xlabel('Timepoint'); plt.ylabel('Contrast Coefficient'); plt.title('Temporal pattern of cluster ' + str(cluster))
        if ylim != None: plt.ylim(ylim[0], ylim[1])
        fig.savefig(output_folder + 'knn_smoothing_' + str(cluster) + '.pdf', bbox_inches = 'tight')    

        fig, ax = plt.subplots()
        llr_fd.plot(axes = ax, group = [0] * fd.size, group_colors = [color], alpha = 0.5)
        fd_mean.plot(axes = ax, group = [0], group_colors = ['black'])
        plt.xlabel('Timepoint'); plt.ylabel('Contrast Coefficient'); plt.title('Temporal pattern of cluster ' + str(cluster))
        if ylim != None: plt.ylim(ylim[0], ylim[1])
        fig.savefig(output_folder + 'LR_smoothing_' + str(cluster) + '.pdf', bbox_inches = 'tight')    

        fig, ax = plt.subplots()
        nw_fd.plot(axes = ax, group = [0] * fd.size, group_colors = [color], alpha = 0.5)
        fd_mean.plot(axes = ax, group = [0], group_colors = ['black'])
        plt.xlabel('Timepoint'); plt.ylabel('Contrast Coefficient'); plt.title('Temporal pattern of cluster ' + str(cluster))
        if ylim != None: plt.ylim(ylim[0], ylim[1])
        fig.savefig(output_folder + 'NW_smoothing_' + str(cluster) + '.pdf', bbox_inches = 'tight')    

def run_fpca(fd, n_components = 3):
    # some issues in sklearn
    fpca_grid = FPCA(n_components = n_components)
    fpca_grid = fpca_grid.fit(fd)

    print(fpca_grid)
    return fpca_grid

def anova_test(fd, fd_groups, perturbations):
    fd_collections = []
    mean_vals = {}
    for perturbation in perturbations:
        idx = [i for i in range(len(fd_groups)) if perturbation in fd_groups[i]]
        fd_sub = fd[idx]
        mean_val = np.mean([np.mean(fd_sub.data_matrix[i]) for i in range(fd_sub.shape[0])])
        mean_vals[perturbation] = mean_val
        fd_collections.append(fd_sub)
        
    # run one-way anova
    val, pval = oneway_anova(*fd_collections)
    return val, pval, mean_vals

def run_anova_test(
    fda_obj, 
    genes, 
    cell_type, 
    perturbations, 
    reps = None, 
    ref_perturbation = None, 
    step_pattern = 'symmetric2', 
    alternative_pattern = 'asymmetric'
):
    df_anova = pd.DataFrame(columns = ['statistics', 'pval'] + ['mean_' + pert for pert in perturbations])
    for gene in genes:
        fd, fd_groups = fda_obj.create_fd_types_perts(gene, cell_type, perturbations, reps, ref_perturbation, step_pattern, alternative_pattern)
        val, pval, mean_vals = anova_test(fd, fd_groups, perturbations)
        mean_vals_flatten = [mean_vals[i] for i in perturbations]
        df_anova.loc[gene, :] = [val, pval] + mean_vals_flatten
    
    _, fdr = fdrcorrection(list(df_anova['pval']))
    df_anova['pval_adjusted'] = fdr 

    new_cols = ['statistics', 'pval', 'pval_adjusted'] + ['mean_' + pert for pert in perturbations]
    df_anova = df_anova.loc[:, new_cols]

    return df_anova