# CellDrift
CellDrift: temporal perturbation effects for single cell data

Shifted gene programs are the key to understanding
perturbation responses in single-cell RNA sequencing
experiments. With the increasing complexity of
perturbational experiments, generative models, such as
scGen and CPA, have been used to interrogate
perturbation latent features utilizing the power of deep
neural networks. However, a lack of interpretability still
prevents biologists from straightforwardly understanding
perturbation responses. Here we present CellDrift, a
generalized linear model (GLM) that accounts for major
covariates, including perturbation groups, cell types, and
their interactions in perturbational single-cell data. We
applied Function Data Analysis (FDA) based on the
results of GLM for perturbational studies with time series
and identified temporal patterns of gene programs in
perturbation responses.

### Main functions
- Conduct differential expression analysis for super complicated biological metadata of single-cell analysis
- Visualization of results (gene modules) in hierarchical and comprehensive way.
- Large scale gene enrichment analysis supported by ToppGene.
- ToppCluster-based comparable gene enrichment analysis.

### Reference
```
Jin, Kang, et al. "An Interactive Single Cell Web Portal Identifies Gene and Cell Networks in COVID-19 Host Responses." Iscience (2021): 103115.
```

### Installation
```python
git clone https://github.com/KANG-BIOINFO/ToppCellPy.git
cd ToppCellPy
pip install .
```

### Tutorial
- [Example on COVID-19 data](https://nbviewer.jupyter.org/github/KANG-BIOINFO/ToppCellPy/blob/main/test/COVID-19%20example.ipynb)
- [Example on IFN-stimulated PBMC data](https://nbviewer.jupyter.org/github/KANG-BIOINFO/ToppCellPy/blob/main/test/IFN-stimulated%20PBMC%20example.ipynb)

### Quick Start
```python
import scanpy as sc
import numpy as np
import ToppCellPy as tp
```

Load example data
```python
adata = sc.read("batch2_all_normalized_filtered.h5ad")
```

Set up shred object run shred plan. shred_plan indicates how to run differential expression analysis and bin_group indicates how to make bins for heatmap visualization.
```python
# stim and cell are metadata of stimulation and cell type in the anndata
shred = tp.Shred(adata = adata,
                shred_plan = ["stim", "cell", "stim+cell|stim"],
                bin_group = ["stim", "cell"],
                order_bins = None,
                order_modules = None,
                method = "wilcoxon")
shred.do_shredplan()
```

Create heatmap data frame and draw heatmap
```python
shred.create_heatmap_matrix()
shred.draw_heatmap()
```

Create GCT file (GCT heatmap is ready to be loaded in Morpheus)
```python
shred.createGCT()
# load the output file of this code onto Morpheus (https://software.broadinstitute.org/morpheus/)
```

Do enrichment for all modules
```python
shred.enrich_modules(categories = ["GeneOntologyCellularComponent"])
```

Draw ToppCluster plot
```python
shred.toppcluster()
```

An alternative option for running the whole pipeline above is using the one-stop code below:
```python
shred = tp.Shred(adata = adata,
            shred_plan = ["stim", "cell", "stim+cell|stim"],
            bin_group = ["stim", "cell"],
            order_bins = None,
            order_modules = None,
            method = "wilcoxon")
shred.toppcell_batchRun(top_n_genes = 200,
                        enrich_categories = ["GeneOntologyCellularComponent"],
                        enrich_ranked = False,
                        toppcluster_run = True,
                        createGCT = True)
```