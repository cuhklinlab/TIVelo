# Installation

TIVelo requires Python 3.8 or later. We recommend using Miniconda for managing the environment.

### Step 1: Create and Activate the Conda Environment
First, create a new Conda environment with Python 3.9:
```bash
conda create -n tivelo python=3.9 -y
conda activate tivelo
```

### Step 2: Install Dependencies

We have published the TIVelo package on PyPI. To ensure a smooth and stable installation process, we recommend installing large dependencies separately before installing TIVelo in a Conda environment.

#### PyTorch
Install PyTorch along with torchvision, torchaudio, and CUDA support:
```bash
conda install pytorch torchvision torchaudio pytorch-cuda=12.1 -c pytorch -c nvidia -y
```

#### Numba
Install Numba:

To enable CUDA GPU support for Numba, install the latest NVIDIA graphics drivers for your platform (the open-source Nouveau drivers do not support CUDA). Then install the CUDA Toolkit package.

For CUDA 12, install the following:
```bash
conda install -c conda-forge cuda-nvcc cuda-nvrtc "cuda-version>=12.0" -y
```

For CUDA 11, install the following:
```bash
conda install -c conda-forge cudatoolkit "cuda-version>=11.2,<12.0" -y
```

Note: You do not need to install the CUDA SDK from NVIDIA.

Cpu version
```bash
conda install numba
```

#### Scanpy
Install Scanpy along with additional dependencies:
```bash
conda install -c conda-forge scanpy python-igraph leidenalg -y
```

#### scVelo
Install scVelo:
```bash
pip install -U scvelo
```

Optional dependencies for directed PAGA and Louvain modularity:
```bash
pip install igraph louvain
```

Optional dependencies for fast neighbor search via hnswlib:
```bash
pip install pybind11 hnswlib
```

### Step 3: Install TIVelo
Finally, install TIVelo:
```bash
pip install tivelo
```

## JupyterLab
To run the tutorials in a notebook locally, please install JupyterLab:
```bash
conda install jupyterlab -y
```

With these steps, TIVelo and its dependencies will be installed and ready for use.