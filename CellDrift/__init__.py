from ._model import model_genes, model_selection, model_timescale, organize_output
from ._setup import setup_celldrift
from ._temporal import FDA, draw_smoothing, draw_smoothing_clusters, fda_cluster

import logging
import datetime
import os

if not os.path.isdir('logging_celldrift'):
    os.mkdir('logging_celldrift')

if not os.path.isdir('output_celldrift'):
    os.mkdir('output_celldrift')

if not os.path.isdir('figures_celldrift'):
    os.mkdir('figures_celldrift')

if not os.path.isdir('fda_celldrift'):
    os.mkdir('fda_celldrift')
    
logging.basicConfig(level = logging.INFO, filename = 'logging_celldrift/celldrift_' + datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '.log',
                     filemode = 'w+', format='Date-Time : %(asctime)s : Line No. : %(lineno)d - %(message)s')

__all__ = [
    'setup_celldrift',
    'model_genes', 
    'model_selection',
    'model_timescale',
    'organize_output',
    'FDA',
    'draw_smoothing',
    'draw_smoothing_clusters',
    'fda_cluster'
]
