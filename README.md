# Neural Signal Decoder
Work in progres...
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
