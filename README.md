## TIVelo: RNA velocity estimation leveraging cluster-level trajectory inference

![workflow_compressed](https://github.com/aqlkzf/typoraimg/raw/main//imgmac/workflow_compressed.png)



RNA velocity inference is a valuable tool for understanding cell development, differentiation, and disease progression. However, existing RNA velocity inference methods typically rely on explicit assumptions of ordinary differential equations (ODE), which prohibits them to capture complex transcriptome expression patterns. In this study, we introduce TIVelo, a novel RNA velocity estimation approach that first determines the velocity direction at the cell cluster level based on trajectory inference, before estimating velocity for individual cells. TIVelo calculates an orientation score to infer the direction at the cluster level without an explicit ODE assumption, which effectively captures complex transcriptional patterns, avoiding potential inconsistencies in velocity estimation for genes that do not follow the simple ODE assumption. We validated the effectiveness of TIVelo by its application to 16 real datasets and the comparison with five benchmarking methods.

## Demo
The instructions for running TIVelo for different kinds of datasets can be found in [examples](https://github.com/cuhklinlab/TIVelo/tree/main/examples).

## Reproducibility 
To reproduce our results, please refer to the folder [notebooks](https://github.com/cuhklinlab/TIVelo/tree/main/docs/source/notebooks/notebooks) or our [tutorial website](https://tivelo.readthedocs.io/en/latest/).

## Datasets 
All datasets used in this research are openly accessible to the public. To achieve the datasets used in this research, please refer to our [Figshare](https://figshare.com/s/d95ebd3f89e991047c07) page. The datasets in Figshare are the version after scVelo standard preprocessing [pipeline](https://scvelo.readthedocs.io/en/stable/VelocityBasics.html).

## Installation

TIVelo requires Python 3.8 or later. We recommend using Miniconda for managing the environment. The typical time for installing our package is 3 minutes.

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
conda install pytorch==2.5.0 torchvision==0.20.0 torchaudio==2.5.0 pytorch-cuda=12.1 -c pytorch -c nvidia -y
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
pip install  scvelo==0.3.1
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

