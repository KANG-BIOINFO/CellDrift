from typing import Optional

import os
import logging
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from ._setup import *

warnings.filterwarnings('ignore')
logger = logging.getLogger(__name__)

# create pseudotime bins
def create_pseudotime_bins(adata, trajectory_key, pseudotime_key, n_bins = 10, min_cells_pert = 10):
    df_meta = adata.obs.copy()
    trajs = df_meta[trajectory_key].unique().tolist()
    num_trajs = len(trajs)

    df_meta_new = pd.DataFrame()
    for traj in trajs:
        df_sub = df_meta.loc[df_meta[trajectory_key] == traj, :].copy()

        pseudotime = np.array(df_sub[pseudotime_key])
        max_pseudo, min_pseudo = np.max(pseudotime), np.min(pseudotime)
        
        # assign bin ids
        bins = np.linspace(min_pseudo, max_pseudo, n_bins + 1)
        df_sub['pseudotime_bin_id'] = np.digitize(df_sub[pseudotime_key], bins[:-1])
        
        # assign median pseudotime
        median_pseudotime = {}
        for bin_id in range(1, n_bins + 1):
            med_value = np.round(np.median(df_sub.loc[df_sub['pseudotime_bin_id'] == bin_id, pseudotime_key]),3)
            median_pseudotime[bin_id] = med_value
        
        df_sub['pseudotime_bin_name'] = [median_pseudotime[i] for i in df_sub['pseudotime_bin_id']]
        df_meta_new = pd.concat([df_meta_new, df_sub], axis = 0)
    
    adata.obs = df_meta_new.copy()
    adata.uns['celldrift']['time_key'] = 'pseudotime_bin_name'
    return adata