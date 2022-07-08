import os
from setuptools import setup
from setuptools import find_packages

# build long description
base_dir = os.path.dirname(os.path.abspath(__file__))
long_description = '\n\n'.join([open(os.path.join(base_dir,'README.md'),'r').read()])

setup(
    name = 'CellDrift',
    version = '0.1.3',
    description = 'CellDrift: A Python Package to Infer Temporal Patterns of Peturbation Effects in Single Cell Data',
    author = 'Kang Jin',
    author_email = 'jinkg@mail.uc.edu',
    maintainer = 'Kang Jin',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    url = 'https://github.com/KANG-BIOINFO/CellDrift',
    project_urls={
        'Ducumentation':'https://celldrift.readthedocs.io/en/latest/index.html#',
    },
    packages = find_packages(),
    install_requires = [
        'numpy',
        'pandas',
        'scanpy>=1.6.0',
        'matplotlib',
        'seaborn',
        'scipy',
        'statsmodels',
        'dtw'
    ],
)