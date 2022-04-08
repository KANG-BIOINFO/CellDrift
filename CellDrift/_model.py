import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import logging
import multiprocessing as mp
from tqdm import tqdm
from patsy import dmatrix
import scipy.stats as st
import pickle

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats.distributions import chi2

import matplotlib.pyplot as plt
from itertools import combinations

logger = logging.getLogger(__name__)

def extract_covar_name(idx_list, key_celltype, key_perturb):
    '''
    Get metadata for coefficients. Four levels are included now, including
    '''
    covar_name_describe = pd.DataFrame(index = idx_list, columns = ['level', 'celltype','perturb', 'x_label'])
    for i in idx_list:
        if i == 'Intercept':
            covar_name_describe.loc[i] = ['intercept',
                                          '',
                                          '',
                                          'intercept']
        elif (key_celltype in i) and (key_perturb not in i):
            covar_name_describe.loc[i] = ['celltype', 
                                          i.lstrip(key_celltype + '[T.').rstrip(']'), 
                                          '',
                                          i.lstrip(key_celltype + '[T.').rstrip(']')]
        elif (key_celltype not in i) and (key_perturb in i):
            covar_name_describe.loc[i] = ['perturb', 
                                          '', 
                                          i.split('T.')[1][:-1],
                                          i.split('T.')[1][:-1]]
        elif (key_celltype in i) and (key_perturb in i):
            covar_name_describe.loc[i] = ['celltype:perturb', 
                                          i.split(':')[0].lstrip(key_celltype + '[T.').rstrip(']'), 
                                          i.split(':')[1].split('T.')[1][:-1],
                                          i.split(':')[0].lstrip(key_celltype + '[T.').rstrip(']') + '-' + i.split(':')[1].split('T.')[1][:-1]]
    return covar_name_describe


def LRT(L1, L2):
    LR = abs(2*(L1-L2))
    p = chi2.sf(LR, 7)
    return p


def emmeans(df_input, betas, v, key_cell, key_perturb, ctrl):
    '''
    Applying emmeans package in GLM model for post-hoc analysis.
    Adapt from https://glennwilliams.me/blog/2021/09/07/estimating-marginal-means-and-pairwise-tests-by-hand-in-python/

    Parameters
    ----------
    df_input
        input of GLM model
    betas
        inferred coefficients from the model
    v
        variance-covariance matrix of coefficients
    '''
    # levels
    n_types = len(df_input[key_cell].cat.categories)
    n_perturbs = len(df_input[key_perturb].cat.categories)

    grid = np.array(np.meshgrid(
                    list(df_input[key_cell].cat.categories),
                    list(df_input[key_perturb].cat.categories)
    )).reshape(2, n_perturbs * n_types).T

    grid = pd.DataFrame(grid, columns = ['cell', 'stim']) # rename them into cell and stim purposely

    # matrix
    mat = dmatrix(
        'C(cell, Treatment(0))*C(stim, Treatment(0))', 
        grid, 
        return_type = 'matrix'
    )

    emmeans = grid
    # emmeans['means'] = mat @ betas # not necessary
    
    # variance-covariance matrix
    vcov = v
    vcov = vcov[~vcov.index.str.contains('Var|Cor')]
    vcov = vcov.loc[:,~vcov.columns.str.contains('Var|Cor')]

    # emmeans['SE'] = np.sqrt(np.diagonal(mat @ vcov) @ mat.T) # not necessary

    # making pairwise mean difference
    combo_names = emmeans.cell + '_' + emmeans.stim
    contrast_pairs = list(combinations(combo_names, 2))
    combo_names_dummy = emmeans.cell + ',' + emmeans.stim
    contrast_pairs_dummy = list(combinations(combo_names_dummy, 2)) # dummy names
    
    contrast_names = []
    cell_names = []
    perturb_names = []
    for names, names_dummy in zip(contrast_pairs, contrast_pairs_dummy):
        contrast_names.append('-'.join(names))
        cell_names.append((names_dummy[0].split(',')[0], names_dummy[1].split(',')[0]))
        perturb_names.append((names_dummy[0].split(',')[1], names_dummy[1].split(',')[1]))
    
    limits = list(combo_names.index)
    combo_cols = list(combinations(limits, 2))
    combo_cols = pd.DataFrame(combo_cols, columns = ['A', 'B'])

    model_pairwise_contrasts = mat[combo_cols.A, :] - mat[combo_cols.B, :]
    model_pairwise_contrasts = model_pairwise_contrasts.T

    # calculate contrast means
    pairwise = pd.DataFrame(np.array([contrast_names, cell_names, perturb_names]).T, columns = ['contrast', 'cell_type', 'perturbation'])
    pairwise['mean'] = np.transpose(betas @ model_pairwise_contrasts)

    # calculate contrast se 
    for i in pairwise.index:
        grad_g = model_pairwise_contrasts[:, i]
        pairwise.at[i, 'SE'] = np.sqrt(
        grad_g @ vcov @ grad_g.T
        )
    
    # calculate constrast z score and p value
    pairwise['z'] = pairwise['mean']/pairwise['SE']
    pairwise['p'] = 2*st.norm.sf(abs(pairwise['z']))

    confint = st.norm.interval(
        0.95, 
        loc=pairwise['mean'], 
        scale=pairwise['SE']
    )
    pairwise['lci'] = confint[0]
    pairwise['uci'] = confint[1]

    pairwise = pairwise[[
        'contrast', 
        'cell_type',
        'perturbation',
        'mean', 
        'lci', 
        'uci', 
        'SE', 
        'z', 
        'p'
    ]].round(4)

    # p value adjustment using FDR correction
    pairwise['p_fdr'] = fdrcorrection(pvals = pairwise['p'])[1]
    
    # subset
    selected_p = [i for i in np.unique(pairwise['perturbation']) if (i[0] == ctrl)]
    selected_c = [i for i in np.unique(pairwise['cell_type']) if (i[0] == i[1])]

    # reverse contrasts (make sure the ref is control)
    pairwise = pairwise.loc[(pairwise['perturbation'].isin(selected_p)) & (pairwise['cell_type'].isin(selected_c)), :]
    pairwise['mean'] = - pairwise['mean']
    pairwise['lci'] = - pairwise['lci']
    pairwise['uci'] = - pairwise['uci']
    pairwise['z'] = - pairwise['z']

    pairwise['perturbation'] = [x[::-1] for x in pairwise['perturbation']]
    pairwise['contrast'] = [(c[0] + '_' + p[0] + '-' + c[1] + '_' + p[1]) for (c, p) in zip(pairwise['cell_type'], pairwise['perturbation'])]
    pairwise.set_index('contrast', inplace = True)
    return pairwise

def glm_gene_chunk(genes, df_input, adjust_batch = True, add_dummy = True):
    '''
    Get GLM models for genes in a chunk during multiprocessing
    '''
    # add dummy values
    if add_dummy:
        1
    # xxx

    df_output_chunk = pd.DataFrame()
    df_pairwiseComp_chunk = pd.DataFrame()
    # build the glm model
    for gene in genes:
        if adjust_batch:
            formula = gene + ' ~ ' + key_celltype + ' * ' + key_perturb + ' + ' + key_batch
            formula_ref = gene + ' ~ ' + key_celltype + ' + ' + key_perturb + ' + ' + key_batch # reference glm model without interaction
        else:
            formula = gene + ' ~ ' + key_celltype + ' * ' + key_perturb
            formula_ref = gene + ' ~ ' + key_celltype + ' + ' + key_perturb # reference glm model without interaction
        
        try:
            # glm_result = smf.glm(formula = formula,
            #                     offset = np.log(df_input[key_size_factor]), 
            #                     data = df_input,
            #                     family = sm.families.NegativeBinomial()).fit(method='lbfgs')

            # glm_ref_result = smf.glm(formula = formula_ref,
            #                             offset = np.log(df_input[key_size_factor]), 
            #                             data = df_input,
            #                             family = sm.families.NegativeBinomial()).fit(method='lbfgs')

            glm_result = smf.glm(formula = formula,
                                offset = np.log(df_input[key_size_factor]), 
                                data = df_input,
                                family = sm.families.NegativeBinomial()).fit()

            glm_ref_result = smf.glm(formula = formula_ref,
                                        offset = np.log(df_input[key_size_factor]), 
                                        data = df_input,
                                        family = sm.families.NegativeBinomial()).fit()
            logger.info('Success for gene: ' + gene)
        except:
            logger.info('GLM can not be fitted for gene: ' + gene)
            continue
        
        # likelihood ratio test
        covariates_cp = [i for i in glm_result.params.index.values if (key_batch not in i)]

        L1 = glm_result.llf
        L2 = glm_ref_result.llf
        p = LRT(L1, L2)

        # add multiple comparison here
        betas = glm_result.params
        betas = betas[covariates_cp]
        try:
            v = glm_result.cov_params()
        except:
            logger.info('Value Error occurred for gene ' + gene + ' with error information: need covariance of parameters for computing (unnormalized) covariances.')
            df_input.to_csv('debug_input_table_' + gene + '.txt', sep = '\t')
            continue

        v = v.loc[covariates_cp,:]
        v = v[covariates_cp]

        df_pairwiseComp = emmeans(df_input, betas, v, key_cell = key_celltype, key_perturb = key_perturb, ctrl = group_control)
    
        # add two columns about pts information
        df_nonzero = df_input[gene].ne(0).groupby([df_input[key_celltype], df_input[key_perturb]]).sum().div(df_input.groupby([key_celltype, key_perturb]).size())
        contrast_tuples = [(c[0], p[0]) for (c, p) in zip(df_pairwiseComp['cell_type'], df_pairwiseComp['perturbation'])]
        reference_tuples = [(c[0], p[1]) for (c, p) in zip(df_pairwiseComp['cell_type'], df_pairwiseComp['perturbation'])]
        df_pairwiseComp['pts_contrast'] = list(df_nonzero.loc[contrast_tuples])
        df_pairwiseComp['pts_reference'] = list(df_nonzero.loc[reference_tuples])

        df_pairwiseComp['gene'] = gene
        df_pairwiseComp_chunk = pd.concat([df_pairwiseComp_chunk, df_pairwiseComp], axis = 0)

        # format output
        conf_int = glm_result.conf_int().loc[covariates_cp,:].copy()
        conf_int.rename({0: 'lower bound', 1: 'upper bound'}, axis = 1, inplace = True)
        
        pvals = pd.DataFrame(glm_result.pvalues[covariates_cp])
        pvals.rename({0: 'pvals'}, axis = 1, inplace = True)

        std_errs = pd.DataFrame(glm_result.bse[covariates_cp])
        std_errs.rename({0: 'standard err'}, axis = 1, inplace = True)

        mu = pd.DataFrame(glm_result.params[covariates_cp])
        mu.rename({0: 'mu'}, axis = 1, inplace = True)
        
        # add coefficient metadat onto the table
        df_meta_coeffs = extract_covar_name(conf_int.index.values, key_celltype, key_perturb)        
        
        df_output = pd.concat([df_meta_coeffs, conf_int, pvals, mu, std_errs], axis = 1)
        df_output['gene'] = gene
        df_output['LRT_pval'] = p
        df_output_chunk = pd.concat([df_output_chunk, df_output], axis = 0)

    return df_output_chunk, df_pairwiseComp_chunk


def model_genes(adata, gene_selection = 'all', n_processes = 16, chunksize = 100, adjust_batch = True, add_dummy = True, pairwise_contrast_only = False, output_suffix = None):
    '''
    Build a model for each gene and curate predictions from iterations.
    Multiprocessing is used here.

    Parameters
    ----------
    adata
        anndata
    gene_selection
        Whether to select all genes or not
    n_processes
        number of processes to be used
    '''

    # initialization

    logger.info('Start to run the model.')
    logger.info('Number of processes: ' + str(n_processes))
    logger.info('Chunksize: ' + str(chunksize))

    global key_celltype; global key_perturb; global key_size_factor; global key_batch; global group_control
    
    key_celltype = adata.uns['celldrift']['cell_type_key']
    key_perturb = adata.uns['celldrift']['perturb_key']
    key_size_factor = adata.uns['celldrift']['size_factor_key']    
    key_batch = adata.uns['celldrift']['batch_key']  
    group_control = adata.uns['celldrift']['group_control']
    n_coeffs = adata.uns['celldrift']['n_coefficients']  
    output_dir = adata.uns['celldrift']['output_dir']

    if gene_selection == 'all':
        genes = adata.var.index.values
        n_genes = len(genes)

    df_meta =  adata.obs[[key_celltype,key_perturb, key_size_factor, key_batch]].copy() # extract cell type and perturbation columns 

    # batch (batch preparation is to avoid large dense tables)
    batch_size = n_processes * chunksize # 800 or 1600
    final_batches = []
    for batch_idx in range(0, n_genes, batch_size):
        logger.info('Start to run batch: ' + str(batch_idx))
        
        # prepare huge expression table for all chunks in one batch
        genes_idx_batch = range(batch_idx, min(batch_idx + batch_size, n_genes))
        genes_batch = genes[genes_idx_batch]
        df_batch_expr = pd.DataFrame(data = adata[:, genes_idx_batch].X.toarray(),
                                    index = adata.obs.index.values,
                                    columns = genes_batch)
        df_batch_input = pd.concat([df_meta, df_batch_expr], axis = 1)
        logger.info('Size of input table: ' + str(sys.getsizeof(df_batch_input) / (1024*1024)) + ' Mb.')

        # separate batch table into chunks
        genes_chunks = []
        df_chunks = []

        for chunk_idx in range(0, batch_size, chunksize):
            genes_chunk = genes_batch[chunk_idx: min(chunk_idx + chunksize, batch_size)]
            genes_chunks.append(genes_chunk)
            cols = list(df_meta.columns) + list(genes_chunk)
            df_chunks.append(df_batch_input[cols].copy())
        
        # multiprocessing
        pool = mp.Pool(processes = n_processes)
        out = [pool.apply_async(glm_gene_chunk, args = (chunk, df_chunk, adjust_batch, add_dummy)) for (chunk, df_chunk) in zip(genes_chunks, df_chunks)]
        pool.close()
        pool.join()

        logger.info('Size of outcome of apply: ' + str(len(out)))
        final = []

        for collect in out:
            final.append(collect.get())
        
        '''
        ray.init(num_cpus = n_processes)
        futures = [glm_gene_chunk.remote(chunk, df_chunk) for (chunk, df_chunk) in zip(genes_chunks, df_chunks)]
        final = ray.get(futures)
        '''
        final_batches += final
        logger.info('Finish batch: ' + str(batch_idx))

    # concatenate results
    logger.info('Start to format modeling output')
    df_all = pd.DataFrame()
    df_pairwiseComp_all = pd.DataFrame()

    for (df_sub, df_pairwiseComp_sub) in final_batches:
        df_all = pd.concat([df_all, df_sub], axis = 0)
        df_pairwiseComp_all = pd.concat([df_pairwiseComp_all, df_pairwiseComp_sub], axis = 0)
    
    suffix = '' if output_suffix == None else output_suffix

    df_all.to_csv(output_dir + 'glm_predictions' + suffix + '.txt', sep = '\t')
    df_pairwiseComp_all.to_csv(output_dir + 'glm_predictions_pairwise_comparisons' + suffix + '.txt', sep = '\t')
    
    df_pred = df_all.copy()
    genes = [df_pred['gene'][i] for i in np.arange(df_pred.shape[0]) if i % n_coeffs == 0]
    coeffs = df_pred['x_label'].values[:n_coeffs]
    adata.uns['celldrift']['coefficient_names'] = coeffs
    adata = adata[:, genes].copy() # subset of genes which were fitted successfully using GLM.

    df_mu = pd.DataFrame(index = genes, columns = coeffs, data = np.array(df_pred['mu']).reshape([len(genes), len(coeffs)]))
    df_sigma = pd.DataFrame(index = genes, columns = coeffs, data = np.array(df_pred['standard err']).reshape([len(genes), len(coeffs)]))
    df_pvals = pd.DataFrame(index = genes, columns = coeffs, data = np.array(df_pred['pvals']).reshape([len(genes), len(coeffs)]))
    df_LRT_pval = pd.DataFrame(index = genes, columns = ['LRT_pval'], data = list(df_pred['LRT_pval'][list(range(0, df_pred.shape[0], len(coeffs)))]))

    if not pairwise_contrast_only:
        df_mu.to_csv(output_dir + 'glm_predictions_parameters' + suffix + '.txt', sep = '\t')
        df_sigma.to_csv(output_dir + 'glm_predictions_standard_errors' + suffix + '.txt', sep = '\t')
        df_pvals.to_csv(output_dir + 'glm_predictions_pvals' + suffix + '.txt', sep = '\t')
        df_LRT_pval.to_csv(output_dir + 'glm_predictions_LRT_pvals' + suffix + '.txt', sep = '\t')

    # select specific perturbation comparisons
    # terms = []
    # for celltype in np.unique(df_input[key_celltype]):
    #     term_name = celltype + '_' + group_perturb + '-' + celltype + '_' + group_control
    #     terms.append(term_name)
    n_pairwise_terms = adata.uns['celldrift']['n_cell_types'] * (adata.uns['celldrift']['n_perturbations'] - 1)
    terms = df_pairwiseComp_all.index.values[:n_pairwise_terms]
    adata.uns['celldrift']['multicomp_terms'] = terms

    df_multicomp_mu = pd.DataFrame(index = genes, columns = terms, data = np.array(df_pairwiseComp_all['mean']).reshape([len(genes), len(terms)]))
    df_multicomp_SE = pd.DataFrame(index = genes, columns = terms, data = np.array(df_pairwiseComp_all['SE']).reshape([len(genes), len(terms)]))
    df_multicomp_lci = pd.DataFrame(index = genes, columns = terms, data = np.array(df_pairwiseComp_all['lci']).reshape([len(genes), len(terms)]))
    df_multicomp_uci = pd.DataFrame(index = genes, columns = terms, data = np.array(df_pairwiseComp_all['uci']).reshape([len(genes), len(terms)]))
    df_multicomp_pvals = pd.DataFrame(index = genes, columns = terms, data = np.array(df_pairwiseComp_all['p_fdr']).reshape([len(genes), len(terms)]))

    # add information to anndata
    adata.varm['pred_mu'] = df_mu.loc[adata.var_names, :].values
    adata.varm['pred_sigma'] = df_sigma.loc[adata.var_names, :].values
    adata.varm['pred_pvals'] = df_pvals.loc[adata.var_names, :].values

    adata.varm['multicomp_mu'] = df_multicomp_mu.loc[adata.var_names, :].values
    adata.varm['multicomp_SE'] = df_multicomp_SE.loc[adata.var_names, :].values
    adata.varm['multicomp_lci'] = df_multicomp_lci.loc[adata.var_names, :].values
    adata.varm['multicomp_uci'] = df_multicomp_uci.loc[adata.var_names, :].values
    adata.varm['multicomp_pvals'] = df_multicomp_pvals.loc[adata.var_names, :].values
    
    adata.var['LRT_pvals'] = df_LRT_pval['LRT_pval']

    return adata


def model_selection(adata, pval_LRT = 0.05, pval_multicomp = 0.05):
    '''
    Select model for better downstream analysis
    '''
    # categorize genes using 
    df_genes = adata.var.copy()

    df_multicomp_pvals = pd.DataFrame(data = adata.varm['multicomp_pvals'], 
                                    index = adata.var_names, 
                                    columns = adata.uns['celldrift']['multicomp_terms'])

    genes_perturbed = df_multicomp_pvals.index.values[np.min(df_multicomp_pvals, axis = 1) < pval_multicomp]
    genes_perturbed_AcrossTypes = df_genes.loc[df_genes['LRT_pvals'] < pval_LRT, :].index.values
    logger.info('Number of perturbed genes: ' + str(len(genes_perturbed)))
    logger.info('Number of genes perturbed differently across cell types: ' + str(len(genes_perturbed_AcrossTypes)))

    adata.var['perturbed_genes'] = 0
    adata.var['differentially_perturbed_genes'] = 0
    adata.var.loc[genes_perturbed, 'perturbed_genes'] = 1
    adata.var.loc[genes_perturbed_AcrossTypes, 'differentially_perturbed_genes'] = 1

    return adata

def model_timescale(adata, **kwargs):
    '''
    run the GLM model across multiple time points. Also refer to model_genes to see model for one time point.

    Parameters
    ----------
    adata
        CellDrift anndata
    **kwargs
        model_genes proporties
    '''
    time_key = adata.uns['celldrift']['time_key']
    time_points = adata.obs[time_key].unique()
    output_dir = adata.uns['celldrift']['output_dir']

    for time in time_points:
        time_val = str(np.round(time,2))
        logger.info('Run GLM model for time point ' + time_val)
        adata_time = adata[adata.obs[time_key] == time, : ].copy()
        adata_time = model_genes(adata_time, output_suffix = '_time_' + time_val, **kwargs)
        adata_time = model_selection(adata_time)
        adata_time.write(output_dir + 'time_' + time_val + '.h5ad')
        

def organize_output(output_folder = 'output_celldrift/', suffix = ''):
    output_files = [i for i in os.listdir(output_folder) if (i.startswith('glm_predictions_pairwise_comparisons_'))]

    df_meta = pd.DataFrame(columns = ['time', 'perturbation', 'cell_type'])
    
    for idx, output_file in enumerate(tqdm(output_files)):    
        # get basic information
        df_pairwise = pd.read_csv(output_folder + output_file, sep = '\t', header = 0)
        time = output_file.split('glm_predictions_pairwise_comparisons_time_')[1].split('.txt')[0]
        
        n_types = len(np.unique(df_pairwise['cell_type']))
        n_perts = len(np.unique(df_pairwise['perturbation']))
        n_genes = len(np.unique(df_pairwise['gene']))

        # make graph per (time + pert)
        df_pairwise['zscore'] = df_pairwise['z']
        coeff_vals = df_pairwise['zscore'].values.reshape([n_genes, n_types * n_perts])

        genes = list(df_pairwise.sort_values(by = ['cell_type', 'perturbation', 'gene'])['gene'][ : n_genes])
        contrasts = [(time + '_' + i) for i in list(df_pairwise['contrast'][ : n_types * n_perts])]
        avail_types = [(i.split('\',')[0].split('(\'')[1]) for i in list(df_pairwise['cell_type'][:(n_types * n_perts)])]
        avail_perts = [(i.split('\',')[0].split('(\'')[1]) for i in list(df_pairwise['perturbation'][:(n_types * n_perts)])]
        
        df_graph = pd.DataFrame(data = coeff_vals, index = genes, columns = contrasts)

        if idx == 0:
            df_combined = df_graph.copy()
        else:
            df_combined = pd.concat([df_combined, df_graph], axis = 1, join = 'outer')
        
        for contrast, cell_type, pert in zip(contrasts, avail_types, avail_perts):
            df_meta.loc[contrast, : ] =  [time, pert, cell_type]

    df_combined.to_csv('fda_celldrift/pairwise_zscores_combined_' + suffix + '.txt', sep = '\t')
    df_meta.to_csv('fda_celldrift/pairwise_contrasts_metadata_' + suffix + '.txt', sep = '\t')

    return df_combined, df_meta
