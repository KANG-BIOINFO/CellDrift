import scanpy as sc
import numpy as np
import pandas as pd
import warnings
import os
import logging
from ._visualize import plot_gene_QC

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
    adata, 
    cell_type_key, 
    perturb_key, 
    size_factor_key, 
    time_key,
    batch_key, 
    control_name, 
    perturb_name = None, 
    min_cells_perGene = None
):
    '''
    Set up CellDrift with key words set.

    Parameters
    ----------
    adata
        AnnData as the input.
    cell_type_key
        The key in adata.obs that indicates cell type assignments.
    perturb_key
        The key in adata.obs that indicates perturbation assignments
    size_factor_key
        The key of size factor for anndata. e.g. n_counts
    time_key
        The key of time point indicator
    batch_key:
        The key of batch covariate
    control_name
        A pairwise list that indicates perturbations users are interested.
        control_name indicates perturbation group in control;
    min_cells_perGene
        filtering criteria for genes.
    '''
    # initialization
    logger.info('Set up the celldrift object...')
    cell_type_key_n = prepare_item(cell_type_key)
    perturb_key_n = prepare_item(perturb_key)
    size_factor_key_n = prepare_item(size_factor_key)
    batch_key_n = prepare_item(batch_key)

    logger.info('rename keys to ' + cell_type_key_n + ', ' + perturb_key_n + ', ' + size_factor_key_n + ', ' + batch_key_n)
    perturb_name = list(np.setdiff1d(np.unique(adata.obs[perturb_key]), [control_name])) if perturb_name == None else [perturb_name]

    # set cell type / perturbation annotations as category
    adata.obs = adata.obs.astype({cell_type_key: 'category', perturb_key: 'category'})
    adata.obs[cell_type_key] = adata.obs[cell_type_key].cat.remove_unused_categories()
    adata.obs[perturb_key] = adata.obs[perturb_key].cat.remove_unused_categories()
    adata.obs[perturb_key].cat.categories = [control_name] + perturb_name
    
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
                                'time_key': time_key,
                                'group_control': control_name,
                                'group_perturb': perturb_name,
                                'n_cell_types': len(adata.obs[cell_type_key_n].cat.categories),
                                'n_perturbations': len(adata.obs[perturb_key_n].cat.categories),
                                'n_coefficients': len(adata.obs[cell_type_key_n].cat.categories) * len(adata.obs[perturb_key_n].cat.categories),
                                'figure_dir': 'figures_celldrift/',
                                'output_dir': 'output_celldrift/'}

    # filter genes
    if min_cells_perGene != None:
        logger.info('do the gene filtering')
        logger.info('number of genes before filtering: ' + str(adata.shape[1]))

        sc.pp.filter_genes(adata, min_cells = 0)
        plot_gene_QC(adata, min_cells_perGene)
        sc.pp.filter_genes(adata, min_cells = min_cells_perGene)

        logger.info('number of genes after filtering: ' + str(adata.shape[1]))

    # reformat genes so that they won't confuse models
    logger.info('reformat gene names to avoid confusing the glm model.')
    genes = adata.var.index.values
    adata.var['original_names'] = genes
    adata.var['new_index'] = [prepare_item(i) for i in genes]
    adata.var.set_index('new_index', inplace = True)

    # remove weird genes
    weird_genes = [g for g in adata.var.index.values if g[0].isdigit()] # genes start with digits
    filtered_genes = np.setdiff1d(adata.var.index.values, weird_genes)
    adata = adata[:, filtered_genes]

    return adata
