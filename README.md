# CellDrift
CellDrift: temporal perturbation effects for single cell data

Perturbation effects on gene programs are commonly investigated in single-cell experiments. Existing models measure perturbation responses independently across time series, disregarding the temporal consistency of specific gene programs. We introduce CellDrift, a generalized linear model based functional data analysis approach to investigate temporal gene patterns in response to perturbations. 
![overview](Examples/overview_CellDrift.png)

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