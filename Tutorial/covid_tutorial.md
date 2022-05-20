# Example of COVID-19 Atlas

This is our favoriate example for the demostration of CellDrift pipeline. We used a large scale [blood atlas](https://doi.org/10.1016/j.cell.2022.01.012) of COVID-19 from Cell.
To focus on general cell types in the huge data, we selected a subset of cells from original single-cell dataset. The subset data contains 6 disease conditions, 6 major cell types and a series of time, which can be downloaded in this [link](https://zenodo.org/record/6466615#.Yl2GJdNKhTY), 

0. Load data and preparation
```python
import numpy as np
import scanpy as sc
import CellDrift as ct

adata = sc.read('COMBAT-CITESeq-DATA_selection_v2_time_raw.h5ad')
adata.obs['size_factor'] = np.sum(adata.X, axis = 1) # compuate size factor
adata

>>> adata
AnnData object with n_obs × n_vars = 116124 × 20807
    obs: 'Annotation_cluster_id', 'Annotation_cluster_name', 'Annotation_minor_subset', 'Annotation_major_subset', 'Annotation_cell_type', 'GEX_region', 'QC_ngenes', 'QC_total_UMI', 'QC_pct_mitochondrial', 'QC_scrub_doublet_scores', 'COMBAT_ID', 'scRNASeq_sample_ID', 'COMBAT_participant_timepoint_ID', 'Source', 'Age', 'Sex', 'Race', 'BMI', 'Hospitalstay', 'Death28', 'Institute', 'PreExistingHeartDisease', 'PreExistingLungDisease', 'PreExistingKidneyDisease', 'PreExistingDiabetes', 'PreExistingHypertension', 'PreExistingImmunocompromised', 'Smoking', 'Symptomatic', 'Requiredvasoactive', 'Respiratorysupport', 'SARSCoV2PCR', 'Outcome', 'TimeSinceOnset', 'Ethnicity', 'Tissue', 'DiseaseClassification', 'Pool_ID', 'Channel_ID'
    var: 'gene_ids', 'feature_types'
    uns: 'Institute', 'ObjectCreateDate', 'Source_colors', 'Technology', 'genome_annotation_version'
    obsm: 'X_umap', 'X_umap_source'
```

1. Preprocessing
```python
genes = [i for i in bdata.var.index.values if (not i.startswith('AB_'))]
adata = adata[:, genes].copy()
```

2. Set up CellDrift object
```python
adata = ct.setup_celldrift(
    adata, 
    cell_type_key = 'Annotation_major_subset',
    perturb_key = 'Source', 
    time_key = 'TimeSinceOnset', # the name of time covariate. Must be numeric
    control_name = 'Control', 
    perturb_name = None, 
    size_factor_key = 'size_factor', 
    batch_key = 'batch'
)
```

3. Run GLM model 
```python
adata = ct.model_timescale(
    adata, 
    n_processes = 16, # number of processes for multiprocessing
    chunksize = 60, # number of genes in each chunk
    pairwise_contrast_only = True, 
    adjust_batch = False
)
os.listdir('output_celldrift')

>>> ['time_1.0.h5ad', 'glm_predictions_time_0.33.txt', 'time_0.67.h5ad', 'glm_predictions_pairwise_comparisons_time_0.11.txt', 'glm_predictions_time_0.89.txt', 'glm_predictions_pairwise_comparisons_time_0.22.txt', 'time_0.33.h5ad', 'time_0.0.h5ad', 'glm_predictions_pairwise_comparisons_time_0.56.txt', 'glm_predictions_pairwise_comparisons_time_0.0.txt', 'glm_predictions_pairwise_comparisons_time_0.33.txt', 'glm_predictions_pairwise_comparisons_time_1.0.txt', 'glm_predictions_time_0.22.txt', 'time_0.11.h5ad', 'glm_predictions_time_0.44.txt', 'time_0.22.h5ad', 'time_0.89.h5ad', 'glm_predictions_pairwise_comparisons_time_0.44.txt', 'glm_predictions_time_1.0.txt', 'glm_predictions_time_0.11.txt', 'time_0.56.h5ad', 'glm_predictions_time_0.67.txt', 'glm_predictions_time_0.56.txt', 'glm_predictions_time_0.0.txt', 'time_0.44.h5ad', 'time_0.78.h5ad', 'glm_predictions_pairwise_comparisons_time_0.89.txt', 'glm_predictions_pairwise_comparisons_time_0.67.txt', 'glm_predictions_time_0.78.txt', 'glm_predictions_pairwise_comparisons_time_0.78.txt']
```

4. Organize the output for functional data analysis (FDA)
```python
ct.organize_output(output_folder = 'output_celldrift/')
os.listdir('fda_celldrift')

>>> ['pairwise_contrasts_metadata_.txt', 'pairwise_zscores_combined_.txt']
```


5. set up FDA object
```python
df_zscore = pd.read_csv('fda_celldrift/pairwise_zscores_combined_.txt', sep = '\t', header = 0, index_col = 0)
df_meta = pd.read_csv('fda_celldrift/pairwise_contrasts_metadata_.txt', sep = '\t', header = 0, index_col = 0)

fda = ct.FDA(df_zscore, df_meta)
```

6. temporal clustering
```python
fd, genes = fda.create_fd_genes(genes = df_zscore.index.values, cell_type = 'Type_0', perturbation = 'Perturb_0')
df_cluster = ct.fda_cluster(fd, genes, n_clusters = 3)
df_cluster.head()

>>>
     genes  clusters_kmeans  clusters_fuzzy
0   Gene_0                0               0
1   Gene_1                0               0
2  Gene_10                0               0
3  Gene_11                0               0
4  Gene_12                0               0
```

7. visualization for each temporal cluster
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

