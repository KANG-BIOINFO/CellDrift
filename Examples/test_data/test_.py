import os
import pandas as pd 
import numpy as np
import scanpy as sc
import CellDrift as ct
# import skfda
# import skfda.misc.hat_matrix as hm
# from skfda.inference.anova import oneway_anova
# from skfda.representation import FDataGrid, FDataBasis
# from skfda.ml.clustering import FuzzyCMeans, KMeans
# from skfda.preprocessing.smoothing import KernelSmoother

adata = sc.read('simulation_n_times_10_rep0.h5ad')
adata.obs['size_factor'] = np.sum(adata.X, axis = 1)
adata.obs['batch'] = 0

adata = ct.setup_celldrift(adata, cell_type_key = 'cell_type', perturb_key = 'perturb', control_name = 'Control', time_key = 'time',
                        perturb_name = None, size_factor_key = 'size_factor', batch_key = 'batch', min_cells_perGene = 0)

adata = ct.model_timescale(adata, n_processes = 1, chunksize=60, pairwise_contrast_only = True, adjust_batch = False)

ct.organize_output(output_folder = 'output_celldrift/')

df_zscore = pd.read_csv('fda_celldrift/pairwise_zscores_combined_.txt', sep = '\t', header = 0, index_col = 0)
df_meta = pd.read_csv('fda_celldrift/pairwise_contrasts_metadata_.txt', sep = '\t', header = 0, index_col = 0)
fda = ct.FDA(df_zscore, df_meta)

fd, genes = fda.create_fd_genes(genes = df_zscore.index.values, cell_type = 'Type_0', perturbation = 'Perturb_0')
df_cluster = ct.fda_cluster(fd, genes, n_clusters = 3)

genes = df_zscore.index.values
ct.draw_smoothing_clusters(fd, df_cluster, genes = genes, 
                            n_neighbors = 2, bandwidth = 1, 
                            cluster_key = 'clusters_fuzzy', 
                            output_folder = 'fda_celldrift/figures/')
