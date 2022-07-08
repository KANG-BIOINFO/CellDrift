Installation
============

Prerequisite
------------
It's recommended to create a new environment using conda (python 3.7 is recommended)::

    conda create -n celldrift_py python=3.7

Install prerequisite package scikit-fda (development version)::

    conda activate celldrift_py # activate celldrift environment
    pip install git+https://github.com/GAA-UAM/scikit-fda.git

Installation
------------
Install CellDrift package::

    git clone https://github.com/KANG-BIOINFO/CellDrift.git
    cd CellDrift
    pip install .
