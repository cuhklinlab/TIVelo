## TIVelo
a tool to be published on Nature.

## Installation


## pytorch
conda create -n tivelo python=3.9
conda activate tivelo
conda install pytorch torchvision torchaudio pytorch-cuda=12.1 -c pytorch -c nvidia -y

## numba
conda install numba 
To enable CUDA GPU support for Numba, install the latest graphics drivers from NVIDIA for your platform. (Note that the open source Nouveau drivers shipped by default with many Linux distributions do not support CUDA.) Then install the CUDA Toolkit package.

For CUDA 12, cuda-nvcc and cuda-nvrtc are required:

$ conda install -c conda-forge cuda-nvcc cuda-nvrtc "cuda-version>=12.0" -y
For CUDA 11, cudatoolkit is required:

$ conda install -c conda-forge cudatoolkit "cuda-version>=11.2,<12.0"
You do not need to install the CUDA SDK from NVIDIA.

## scanpy
conda install -c conda-forge scanpy python-igraph leidenalg -y





## scVelo
pip install -U scvelo -y
Parts of scVelo (directed PAGA and Louvain modularity) require (optional):
pip install igraph louvain -y
Using fast neighbor search via hnswlib further requires (optional):
pip install pybind11 hnswlib -y







## Jupyterlab 
To run the tutorials in a notebook locally, please install:
conda install jupyterlab -y


