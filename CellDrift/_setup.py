from typing import Optional

import os
import logging
import warnings
import random
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from ._trajectory import *

warnings.filterwarnings('ignore')
logger = logging.getLogger(__name__)

def prepare_item(input_item):
    input_item = input_item.replace('-', '_')
    input_item = input_item.replace('.', '_')
    input_item = input_item.replace(' ', '_')
    input_item = input_item.replace('(', '_')
    input_item = input_item.replace(')', '_')
    return input_item

def setup_celldrift(
    adata: AnnData, 
    cell_type_key: str, 
    perturb_key: str, 
    size_factor_key: str, 
    time_key: Optional[str],
    batch_key: Optional[str], 
    control_name: Optional[str], 
    perturb_name: Optional[list] = None, 
    use_pseudotime: bool = False,
    pseudotime_n_bins: int = 10,
    n_reps: int = 3,
    n_cells_perBlock: int = 50,
    min_cells_perGene: int = None
):
    '''
    Set up CellDrift with key words set.

    :param adata: AnnData as the input.
    :param cell_type_key: The key in adata.obs that indicates cell type assignments.
    :param perturb_key: The key in adata.obs that indicates perturbation assignments
    :param size_factor_key: The key of size factor for anndata. e.g. n_counts.
    :param time_key: The key of time point indicator. Values in this covariate must be numeric
    :param batch_key: The key of batch covariate.
    :param control_name: control value in the perturb_key.
    :param perturb_name: a list of perturbation labels of interest. Default is None. If None, all unique values apart from control_name from perturb_key will be used.
    :param use_pseudotime: whether time key is pseudotime or not. Default is False.
    :param pseudotime_n_bins: number of pseudotime bins if use_pseudotime is True. Default is 10
    :param n_reps: number of replicate runs of CellDrift. Default is 3.
    :param n_cells_perBlock: number of cells in each (cell type + perturbation + time point) combination. Default is 50
    :param min_cells_perGene: filtering criteria for genes. Default is None.

    Examples::

        import scanpy as sc
        import CellDrift as ct
        kc = 'Cell_Type'
        kp = 'Disease_State'
        kt = 'Time_Since_Onset'
        kp_p = 'COVID_SEVERE'
        kp_c = 'CTRL'
        adata = sc.read('test.h5ad')
        adata = ct.setup_celldrift(
            adata, 
            cell_type_key = kc, 
            perturb_key = kp, 
            time_key = kt,
            control_name = kp_c, 
            perturb_name = None, 
            size_factor_key = 'size_factor', 
            batch_key = None,
            min_cells_perGene = 200,
            use_pseudotime = False
        )
    '''
    # initialization
    logger.info('---Step1: Set up the celldrift object...')
    cell_type_key_n = prepare_item(cell_type_key)
    perturb_key_n = prepare_item(perturb_key)
    size_factor_key_n = prepare_item(size_factor_key)
    time_key_n = prepare_item(time_key)

    if not pd.api.types.is_numeric_dtype(adata.obs[time_key_n]):
        raise Exception('Values in ' + time_key_n + ' must be numeric (int or float).')
    
    if not pd.api.types.is_numeric_dtype(adata.obs[size_factor_key_n]):
        raise Exception('Values in ' + size_factor_key_n + ' must be numeric (int or float).')

    if batch_key != None:
        batch_key_n = prepare_item(batch_key)
        logger.info('---Step1: rename keys to ' + cell_type_key_n + ', ' + perturb_key_n + ', ' + size_factor_key_n + ', ' + batch_key_n)
    else:
        batch_key_n = batch_key
        logger.info('---Step1: rename keys to ' + cell_type_key_n + ', ' + perturb_key_n + ', ' + size_factor_key_n + ', ')

    perturb_name = list(np.setdiff1d(np.unique(adata.obs[perturb_key]), [control_name])) if perturb_name == None else perturb_name

    # set cell type / perturbation annotations as category
    adata.obs = adata.obs.astype({cell_type_key: 'category', perturb_key: 'category'})
    adata.obs[cell_type_key] = adata.obs[cell_type_key].cat.remove_unused_categories()
    adata.obs[perturb_key] = adata.obs[perturb_key].cat.remove_unused_categories()
    adata.obs[perturb_key].cat.reorder_categories([control_name] + perturb_name, inplace = True)
    
    logger.info('perturbation levels includes: ')
    logger.info(adata.obs[perturb_key].cat.categories)

    # rename columns
    adata.obs = adata.obs.rename({cell_type_key: cell_type_key_n,
                                perturb_key: perturb_key_n,
                                size_factor_key: size_factor_key_n,
                                batch_key: batch_key_n})                                
    
    # store attributes
    adata.uns['celldrift'] =  {'cell_type_key': cell_type_key_n,
                                'perturb_key': perturb_key_n,
                                'size_factor_key': size_factor_key_n,
                                'batch_key': batch_key_n,
                                'time_key': time_key_n,
                                'group_control': control_name,
                                'group_perturb': perturb_name,
                                'n_cell_types': len(adata.obs[cell_type_key_n].cat.categories),
                                'n_perturbations': len(adata.obs[perturb_key_n].cat.categories),
                                'n_coefficients': len(adata.obs[cell_type_key_n].cat.categories) * len(adata.obs[perturb_key_n].cat.categories),
                                'figure_dir': 'Figures_CellDrift/',
                                'output_dir': 'Coefficients_CellDrift/',
                                'temporal_dir': 'Temporal_CellDrift/'}
    # filter genes
    if min_cells_perGene != None:
        logger.info('---Step1: do the gene filtering')
        logger.info('---Step1: number of genes before filtering: ' + str(adata.shape[1]))

        sc.pp.filter_genes(adata, min_cells = 0)
        sc.pp.filter_genes(adata, min_cells = min_cells_perGene)

        logger.info('---Step1: number of genes after filtering: ' + str(adata.shape[1]))
    # reformat genes so that they won't confuse models
    logger.info('---Step1: reformat gene names to avoid confusing the glm model.')
    genes = adata.var.index.values
    adata.var['original_names'] = genes
    adata.var['new_index'] = [prepare_item(i) for i in genes]
    adata.var.set_index('new_index', inplace = True)
    # remove weird genes
    weird_genes = [g for g in adata.var.index.values if g[0].isdigit()] # genes start with digits
    filtered_genes = np.setdiff1d(adata.var.index.values, weird_genes)
    adata = adata[:, filtered_genes]

    # cretae pseudotime bins
    if use_pseudotime == True:
        adata = create_pseudotime_bins(
            adata = adata.copy(), 
            trajectory_key = cell_type_key_n, 
            pseudotime_key = time_key_n, 
            n_bins = pseudotime_n_bins, 
            min_cells_pert = min_cells_perGene
        )

    # select cells for replicates
    df_cells = adata.obs.copy()
    df_cells['combs'] = [(i + '-' + j + '-' + str(k)) for (i,j,k) in \
                        zip(df_cells[cell_type_key_n], df_cells[perturb_key_n], df_cells[adata.uns['celldrift']['time_key']])]
    combs_unique = df_cells['combs'].unique()
    
    for rep in range(1, n_reps + 1):
        cells_pooled = []    
        for comb in combs_unique:
            df_subset = df_cells.loc[df_cells['combs'] == comb, :].copy()
            max_cells = df_subset.shape[0]
            random.seed(rep)
            try:
                selected_cells = random.sample(list(df_subset.index.values), n_cells_perBlock)
            except:
                selected_cells = random.sample(list(df_subset.index.values), max_cells)
                logger.info('---Step1:' + comb + ' has less than ' + str(n_cells_perBlock) + 'cells.')

            cells_pooled += selected_cells
        rep_key = 'selection_rep_' + str(rep) 
        df_cells[rep_key] = False
        df_cells.loc[cells_pooled, rep_key] = True
    
    adata.obs = df_cells.copy()
    adata.uns['celldrift']['n_reps'] = n_reps
    adata.uns['celldrift']['rep_key'] = ['selection_rep_' + str(i) for i in range(1, n_reps + 1)]
        
    return adata

def update_celldrift(adata):
    cell_type_key = adata.uns['celldrift']['cell_type_key']
    perturb_key = adata.uns['celldrift']['perturb_key']
    control_group = adata.uns['celldrift']['group_control']

    perturb_groups = np.setdiff1d(adata.obs[perturb_key].unique(), [control_group])
    adata.obs[perturb_key].cat.categories = [control_group] + list(perturb_groups)
    adata.obs[cell_type_key].cat.categories = adata.obs[cell_type_key].unique()

    adata.uns['celldrift']['n_perturbations'] = len(adata.obs[perturb_key].unique())
    adata.uns['celldrift']['n_cell_types'] = len(adata.obs[perturb_key].unique())
    return adata

