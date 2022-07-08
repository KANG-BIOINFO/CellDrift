Tutorial
==================

.. _narrative_tutorial:

HIV Infection Study
------------------------------

This is an example of the CellDrift application in a longitudinal HIV-1 hyperacute infection study.
This scRNA-seq data covers peripheral blood mononuclear cells from untreated individuals with HIV infections before and after acute infection. Longitudinal samples were collected from these patients from week 0 to year 1 post-infection. Multiple cell types involved in immune responses were annotated and analyzed. The data was downloaded from single cell portal SCP256.

We first prepare the input data.
>>> import numpy as np
    import pandas as pd
    import scanpy as sc
    import CellDrift as ct

    adata = sc.read('../raw_v2.h5ad')
    kc = 'cell_type' # key of cell type
    kp = 'disease_v2' # key of perturbation
    kt = 'time_v2' # key of time points

    print(adata)
    print(adata.obs.head())

    >>> adata
    AnnData object with n_obs × n_vars = 59286 × 16980
        obs: 'cell_type', 'sample', 'nUMI', 'disease', 'time', 'time_v2', 'disease_v2'
        
    >>> adata.obs.head()
        cell_type      sample  nUMI                 disease     time  time_v2 disease_v2
    NAME                                                                                   
    S00001    B cell  P3_4 Weeks  1660  HIV infectious disease  4 Weeks       28        hiv
    S00002    B cell  P3_4 Weeks  1198  HIV infectious disease  4 Weeks       28        hiv
    S00003    B cell  P3_4 Weeks  1459  HIV infectious disease  4 Weeks       28        hiv
    S00004    B cell  P3_4 Weeks  1402  HIV infectious disease  4 Weeks       28        hiv
    S00005    B cell  P3_4 Weeks  1179  HIV infectious disease  4 Weeks       28        hiv

It is recommended to do the feature selection. The reason for feature selection is to investigate the most interesting genes and to reduce running time.
>>> # select highly variable genes
    def select_variable_genes(data, n_top_genes):
        # normalize raw data
        sc.pp.filter_genes(data, min_cells = 200)
        sc.pp.normalize_total(data, target_sum=1e4)
        sc.pp.log1p(data)
        # select top variable genes
        sc.pp.highly_variable_genes(data, n_top_genes = n_top_genes) 
        high_var_genes = data.var_names[data.var.highly_variable]

        return high_var_genes

    high_var_genes = select_variable_genes(adata.copy(), n_top_genes = 1200)
    adata = adata[:, high_var_genes].copy()

Then we set up the CellDrift object, which is a basic format with necessary information for downsteam analysis.
.. code-block:: 
    adata = ct.setup_celldrift(
        adata, 
        cell_type_key = kc, 
        perturb_key = kp, 
        time_key = kt, # the name of time covariate. Must be numeric
        control_name = 'ctrl', 
        perturb_name = None, 
        size_factor_key = 'nUMI', 
        batch_key = None, 
        n_reps = 3, 
        n_cells_perBlock = 100, 
        use_pseudotime = False, 
        min_cells_perGene = 200
    )

After we get the CellDrift object in a required format, we run the generalized linear model across input time points.
.. code-block:: 
    adata = ct.model_timescale(
        adata, 
        n_processes = 8, # number of processes for multiprocessing
        chunksize = 100, # number of genes in each chunk
        pairwise_contrast_only = True, 
        adjust_batch = False
    )
The output of generalized linear model is stored in the folder ``Coefficients_CellDrift``

We load the contrast coefficients (z scores) from the last step and then build up a ``FDA`` object.
.. code-block::
    # load data
    df_zscore = pd.read_csv('Temporal_CellDrift/Contrast_Coefficients_combined_zscores_.txt', sep = '\t', header = 0, index_col = 0)
    df_meta = pd.read_csv('Temporal_CellDrift/Contrast_Coefficients_combined_metadata_.txt', sep = '\t', header = 0, index_col = 0)

    # re-annotate discrete time points into continuous time points
    time_origin = [180, 365,  28,  21,   1,  14,   7]
    time_new = [6, 7, 5, 4, 1, 3, 2]
    time_dict = dict(zip(time_origin, time_new))
    df_meta['time'] = [time_dict[i] for i in df_meta['time']]

    # create object
    cell_type = 'monocyte'
    perturbations = ['hiv']
    perturbation = 'hiv'
    fda = ct.FDA(df_zscore, df_meta)

We then do the temporal clustering on the ``FDA`` object. 
.. code-block:: 
    # find clusters
    input_genes = df_zscore.index.values
    fd1, genes = fda.create_fd_genes(input_genes, cell_type = cell_type, perturbation = perturbation)
    df_cluster = ct.fda_cluster(fd1, genes, n_clusters = 8, seed = 42, output_folder = 'Temporal_CellDrift/')

    # visualize clusters
    ct.draw_smoothing_clusters(fd1, df_cluster, n_neighbors = 2, bandwidth = 1, 
                            cluster_key = 'clusters_fuzzy', output_folder = 'Temporal_CellDrift/cluster_fuzzy/')

    ct.draw_smoothing_clusters(fd1, df_cluster, n_neighbors = 2, bandwidth = 1, 
                            cluster_key = 'clusters_kmeans', output_folder = 'Temporal_CellDrift/cluster_kmeans/')
