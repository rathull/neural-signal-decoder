# Neural Signal Decoder
Information regarding project available in `report.pdf`.

This repo is still a work in progress. All results are avilable in `plots.ipynb`, but this will be updated to separate Python scripts soon.

# Getting Started

### Step 1: Creating the Environment
1. Create a new Conda environment with the necessary libraries:

```bash
conda create -n signals python=3.9 jupyter numpy pandas scipy scikit-learn matplotlib seaborn
conda activate signals
```

2. Install PyTorch with MPS support (on Mac) or CUDA support:
```bash
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/nightly/cpu
```
