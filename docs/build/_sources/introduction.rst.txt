Introduction
==================

.. _narrative_background:

Why do we develop CellDrift?
-------------------------
Researchers have applied the single-cell RNA-sequencing technology in experiments with perturbation settings such as diseases, treatmentS, genetic mutations, and organ differentiation to explore transcriptional profiles across various biochemical states. 

The response to perturbation can vary over time, which is overlooked in many single cell studies. As a result, we develop CellDrift, a generalized linear model-based functional data analysis method capable of identifying covarying temporal patterns of various cell types in response to perturbations. It includes functions below:

1. Disentangle common and cell type specific perturbation effects across time;
2. Identify patterns of genes that have similar temporal perturbation responses;
3. Prioritize genes with distinct temporal perturbation responses between perturbations or cell types;
4. Infer differential genes of perturbational states in the pseudo-time trajectories.

![-](overview.png)

What is the appropriate input data?
-----------------------------------


What is the expected output?
----------------------------
