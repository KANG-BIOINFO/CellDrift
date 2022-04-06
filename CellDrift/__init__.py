from ._model import model_genes, model_selection
from ._visualize import plot_coefficients, plot_cell_type_correlation, plot_patternScore, plot_coefficient_comparison, plot_specific_neighborhood, plot_correlation, plot_pairwise_coefficients, plot_top_perturbed_genes, plot_neigh_coeffs_umap, plot_neigh_coeffs_violin, extract_global_neighbor_coeffs_genes
from ._ensemble import ensemble_clustering
from ._local import make_nhoods, calculate_neigh_relative_likelihood, neighbor_GLM, neighbor_GLM_linear, outputDF_to_anndata, cluster_ndata
from ._setup import setup_celldrift
from ._pattern import prioritize_patterns, find_cooccuring_pattern
from ._decomposition import run_factor_analysis

import logging
import datetime
import os

if not os.path.isdir("logging_celldrift"):
    os.mkdir("logging_celldrift")

if not os.path.isdir("output_celldrift"):
    os.mkdir("output_celldrift")

if not os.path.isdir("figures_celldrift"):
    os.mkdir("figures_celldrift")
    
logging.basicConfig(level = logging.INFO, filename = "logging_celldrift/celldrift_" + datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + ".log",
                     filemode = 'w+', format='Date-Time : %(asctime)s : Line No. : %(lineno)d - %(message)s')


__all__ = ["model_genes", "plot_coefficients", "ensemble_clustering",
            "make_nhoods", "calculate_neigh_relative_likelihood", "setup_celldrift", "prioritize_patterns",
            "plot_cell_type_correlation", "plot_patternScore", "plot_coefficient_comparison",
            "plot_specific_neighborhood", "neighbor_GLM", "model_genes_linear", "neighbor_GLM_linear",
            "model_selection", "emmeans", "plot_correlation", "plot_pairwise_coefficients",
            "plot_top_perturbed_genes", "cluster_ndata", "plot_neigh_coeffs_umap",
            "plot_neigh_coeffs_violin", "find_cooccuring_pattern", "extract_global_neighbor_coeffs_genes",
            "run_factor_analysis"]
