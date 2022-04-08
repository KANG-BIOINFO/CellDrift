# CellDrift
CellDrift: temporal perturbation effects for single cell data

Perturbation effects on gene programs are commonly investigated in single-cell experiments. Existing models measure perturbation responses independently across time series, disregarding the temporal consistency of specific gene programs. We introduce CellDrift, a generalized linear model based functional data analysis approach to investigate temporal gene patterns in response to perturbations. 
![overview](Examples/overview.png)

### Reference
```
It has been presented as a poster on ProbGen22. The manuscript is in preparation and will come out soon.
```

### Installation
```python
git clone https://github.com/KANG-BIOINFO/CellDrift.git
cd CellDrift
pip install .
```

### Tutorial
- [Example on a simple simulated data]()
- [Example on COVID-19 Atlas (under construction)]()
- [Example on Gut Differentiation (under construction)]()

### Quick Start
```python
import numpy as np
import scanpy as sc
import CellDrift as ct
```

Load data and preparation
```python
adata = sc.read("example.h5ad")
adata.obs['size_factor'] = np.sum(adata.X, axis = 1)
```

Set up CellDrift object
```python
adata = ct.setup_celldrift(
    adata, 
    cell_type_key = 'cell_type',
    perturb_key = 'perturb', 
    time_key = 'time', # the name of time covariate. Must be numeric
    control_name = 'Control', 
    perturb_name = None, 
    size_factor_key = 'size_factor', 
    batch_key = 'batch', 
    min_cells_perGene = 0
)
```

Run GLM model 
```python
adata = ct.model_timescale(
    adata, 
    n_processes = 16, # number of processes for multiprocessing
    chunksize = 60, # number of genes in each chunk
    pairwise_contrast_only = True, 
    adjust_batch = False
)
```

Organize the output for functional data analysis (FDA)
```python
ct.organize_output(output_folder = 'output_celldrift/')
```

set up FDA object
```python
df_zscore = pd.read_csv('fda_celldrift/pairwise_zscores_combined_.txt', sep = '\t', header = 0, index_col = 0)
df_meta = pd.read_csv('fda_celldrift/pairwise_contrasts_metadata_.txt', sep = '\t', header = 0, index_col = 0)

fda = ct.FDA(df_zscore, df_meta)
```

temporal clustering
```python
fd, genes = fda.create_fd_genes(genes = df_zscore.index.values, cell_type = 'Type_0', perturbation = 'Perturb_0')
df_cluster = ct.fda_cluster(fd, genes, n_clusters = 3)
```

visualization for each temporal cluster
```python
genes = df_zscore.index.values
ct.draw_smoothing_clusters(
    fd, 
    df_cluster, 
    genes = genes, 
    n_neighbors = 2, 
    bandwidth = 1, 
    cluster_key = 'clusters_fuzzy', 
    output_folder = 'fda_celldrift/figures/'
)
```