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

### Reference
```
It has been presented as a poster on ProbGen22. The manuscript is in preparation and will come out soon.
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