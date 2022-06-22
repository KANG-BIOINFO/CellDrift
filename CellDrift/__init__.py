from ._model import model_genes, model_timescale
from ._setup import setup_celldrift
from ._temporal import FDA, draw_smoothing, draw_smoothing_clusters, fda_cluster, run_fpca, run_anova_test

import logging
import datetime
import os

if not os.path.isdir('Logging_CellDrift'):
    os.mkdir('Logging_CellDrift')

if not os.path.isdir('Coefficients_CellDrift'):
    os.mkdir('Coefficients_CellDrift')

if not os.path.isdir('Figures_CellDrift'):
    os.mkdir('Figures_CellDrift')

if not os.path.isdir('Temporal_CellDrift'):
    os.mkdir('Temporal_CellDrift')
    
logging.basicConfig(level = logging.INFO, filename = 'Logging_CellDrift/logging_' + datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '.log',
                     filemode = 'w+', format='Date-Time : %(asctime)s : Line No. : %(lineno)d - %(message)s')

__all__ = [
    'setup_celldrift',
    'model_genes', 
    'model_timescale',
    'FDA',
    'draw_smoothing',
    'draw_smoothing_clusters',
    'fda_cluster',
    'run_fpca',
    'run_anova_test',
]
