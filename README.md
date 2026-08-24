# PHOCI: Predictor of Higher-Order Chromatin Interactions

## Introduction

PHOCI (Predictor of Higher-Order Chromatin Interactions) is a deep learning-based predictive model designed to calculate the probability of chromatin multi-way interactions by integrating Hi-C data and epigenomic data. Using random walks to generate numerous candidate multi-way interaction samples and filtering them through probability predictions, PHOCI further explores chromatin multi-way interaction association rules related to specific genomic loci, such as enhancer-promoter multi-way association rules.

### Background

The 3D structure of chromatin and its multi-way interactions play a critical role in understanding gene regulation and genome functionality. However, traditional techniques like Hi-C have limitations in comprehensively capturing multi-way interactions. Emerging technologies such as Pore-C offer higher-resolution data but face challenges like data scarcity and the lack of systematic probabilistic models.

### Our Solution

PHOCI addresses these challenges through the following strategies:
1. **Graph Convolutional Network (GCN)-based Predictive Model**: Combines Hi-C and epigenomic data to directly predict chromatin multi-way interactions.
2. **Random Walk for Candidate Generation**: Performs random walks on chromatin contact topological graphs constructed from Hi-C data to generate candidate multi-way interactions.
3. **Probability Estimation and Validation**: Estimates probabilities for candidate samples to validate their plausibility.
4. **Association Rule Mining**: Utilizes the Apriori algorithm to discover significant multi-way association rules (e.g., enhancer-promoter relationships).

---

## File Structure

The following outlines the main files and directories in this project:
```
.
|-- annotation                    # Data annotation directory
|-- apriori_plots                 # Visualized results of Apriori rules
|-- apriori_rules                 # Extracted multi-way interaction association rules
|-- data                          # Folder for processed data
|-- models                        # Deep learning model files
|-- results                       # Data results for multi-way interaction association rules
|-- aggregator.py                 # Model aggregation module
|-- batch.py                      # Batch processing module
|-- bigwig_process.py             # BigWig data processing for epigenomics
|-- config.py                     # Configuration file
|-- data_process_for_prediction.py    # Data processing script
|-- encoder.py                    # Encoder module
|-- environment.yml               # Conda environment configuration file
|-- gene_position.py              # Gene position information utility
|-- hic_process.py                # Hi-C data processing script
|-- PHOCI_prediction.ipynb        # Main prediction Jupyter Notebook
|-- plot_results.py               # Result visualization tools
|-- random_walk.py                # Random walk generation script
|-- README.md                     # Project documentation
|-- run_data_process_for_prediction.py   # Main data processing script for prediction
|-- run_data_process_for_training.py   # Main data processing script for training
|-- training.py                   # Training script using GM12878 cell-line as example
|-- make_splits_val.py            # Negative sampling script
|-- sampler.py                    # Negative sampling script
`-- utils.py                      # General utility module
```

---

## Quick Start

### Prerequisites

Ensure you have [Conda](https://docs.conda.io/en/latest/miniconda.html) installed on your system.

### Environment Setup

To streamline the installation process, use the provided `environment.yml` file to quickly create and configure a virtual environment.

1. **Clone the repository**:
   ```bash
   git clone [https://github.com/yunyiwu/PHOCI.git](https://github.com/yunyiwu/PHOCI.git)
   cd PHOCI
   ```

2. **Create and activate the environment**:
   ```bash
   conda env create -f environment.yml
   conda activate phoci

   ```

> **Note**: If `conda env create fails due to some reasons, follow the **Alternative Setup** below.


3. **Run the project**:
Execute the required scripts or applications as per your project needs.

---

### Alternative Setup (Manual Installation)

If installation via `environment.yml` fails or hangs, you can manually build the environment step-by-step:

1. **Clone the repository**:
```bash
git clone https://github.com/yunyiwu/PHOCI.git
cd PHOCI

```


2. **Step-by-step Installation**:
* **Step 1: Create base Python 3.8 environment**
```bash
conda create -n phoci python=3.8 -y
conda activate phoci

```


* **Step 2: Upgrade basic build tools**
```bash
pip install --upgrade pip setuptools wheel

```


* **Step 3: Install PyTorch (CUDA 11.7)**
```bash
pip install torch==2.0.1 torchvision==0.15.2 torchaudio==2.0.2

```


* **Step 4: Install core data science packages**
```bash
pip install pandas==1.5.3 numpy==1.24.3 scipy==1.10.1 scikit-learn==1.3.2 matplotlib==3.7.2 seaborn==0.12.2 h5py==3.11.0 notebook==6.5.7

```


* **Step 5: Install PyTorch Geometric (PyG) and PyG C++ extensions**
```bash
pip install torch-geometric==2.5.2 \
  torch-scatter==2.1.2+pt20cu117 \
  torch-sparse==0.6.18+pt20cu117 \
  torch-cluster==1.6.3+pt20cu117 \
  torch-spline-conv==1.2.2+pt20cu117 \
  -f https://data.pyg.org/whl/torch-2.0.0+cu117.html

```


* **Step 6: Install 3D genomics / Hi-C analysis tools**
```bash
# Install general 3D genomics packages
pip install cooler==0.9.3 bioframe==0.3.3 pybigwig==0.3.22 pyfaidx==0.8.1.2 trackc==0.0.18

# Option A: Install hicstraw from PyPI
pip install hicstraw

# Option B: Install from a local pre-built wheel (if PyPI build fails, and for Linux only)
# pip install third_party_packages/hic_straw-1.3.1-cp38-cp38-linux_x86_64.whl
```

* **Step 7: Install utilities**
```bash
pip install torchmetrics==1.4.1 lightning-utilities==0.11.6 aiohttp==3.10.3 pyyaml==6.0.2 tqdm==4.66.5
pip install efficient-apriori==2.0.5
```
---

### Docker Containerization
To avoid potential system-level driver mismatches or C++ extension build errors, we provide a pre-configured `Dockerfile` based on NVIDIA CUDA 11.7.1.

1. **Build the Docker Image**:
```bash
docker build -t phoci:latest .

```

2. **Run Container with GPU Support**:
```bash
docker run --gpus all -it --rm -v $(pwd):/workspace/PHOCI phoci:latest
```

### Overview of Environment Dependencies

The `environment.yml` pins all essential packages for running PHOCI:

* **Python**: `3.8`
* **Deep Learning Stack**: PyTorch (2.0.1, CUDA 11.7) & PyTorch Geometric (PyG 2.5.2)
* **3D Genomics Analysis**: `cooler`, `bioframe`, `pyBigWig`, `pyfaidx`, `trackc`, `hicstraw`
* **Scientific Computing**: NumPy, Pandas, SciPy, Scikit-learn
* **Data Visualization**: Matplotlib, Seaborn

> **Note**: See `environment.yml` for the complete list of pinned dependencies.

---

## Data Preparation

1. **Configuration File Setup**:
   - Prepare paths for Hi-C files (with `.hic` extension) and epigenomic data files (with `.bigWig` extension).
   - Edit the `config.py` file with the following details:
     - Test cell line name: `test_cell_line`.
     - Hi-C file path: `config_test["hic_file"]`.
     - Folder path for BigWig files: `config_test["bigwig_dir"]`.
     - List of BigWig filenames: `config_test["bigwigs"]` (without file extensions). The list should follow this order:
       **H3K4me3, H3K27ac, H3K27me3, H3K4me1, H3K36me3, H3K9me3, H3K9ac, H3K4me2, H4K20me1, H2AFZ, H3K79me2, CTCF, POLR2A, RAD21, ATAC**.

   **Example Configuration**:
   ```python
   ###########################set by users########################
   test_cell_line = "H1"

   config_test["hic_file"] = "raw_data/H1_hg38/4DNFID162B9J.hic"
   config_test["bigwig_dir"] = "raw_data/H1_hg38/bigwigs/"
   config_test["bigwigs"] = [
       "ENCFF493QWY",
       "ENCFF314KQD",
       "ENCFF345VHG",
       "ENCFF088MXE",
       "ENCFF488THD",
       "ENCFF183MHJ",
       "ENCFF084JKQ",
       "ENCFF860NVB",
       "ENCFF156JZY",
       "ENCFF296IBP",
       "ENCFF401PZS",
       "ENCFF648BTZ",
       "ENCFF933YTR",
       "ENCFF002NBT",
       "4DNFICPNO4M5"
   ]
   ###########################set by users########################
   ```

2. **Run the Data Processing Script**:
   - Execute the `run_data_process_for_prediction.py` script to process raw data.
   - The processed files will be automatically saved in the `data` folder.

---

## Prediction and Association Rule Mining

1. **Set Parameters**:
   - Open `PHOCI_prediction.ipynb` in Jupyter Notebook.
   - Configure the following parameters:
     - **Cell Line**: `cell_line`.
     - **Model Directory**: `model_dir` (default: comprehensive model).
     - **Target Gene**: `gene_name`.
     - **Visualization File Paths**: `file_paths` for BigWig and BigBed files. Set `file_paths` to `False` if visualization is not needed.

   **Example Configuration**:
   ```python
   cell_line = 'K562'

   model_dir = 'models/comprehensive/'

   gene_name = 'MYB'

   file_paths = {
       'atac': "raw_data/K562_hg38/bigwigs/ENCFF754EAC.bigWig",
       'h3k27ac': "raw_data/K562_hg38/bigwigs/ENCFF381NDD.bigWig",
       'h3k4me1': "raw_data/K562_hg38/bigwigs/ENCFF761XBZ.bigWig",
       'ctcf': "raw_data/K562_hg38/bigwigs/ENCFF675GVW.bigWig",
       'rad21': "raw_data/K562_hg38/bigwigs/ENCFF652NKM.bigWig",
       'chrom_state': "annotation/ENCFF319VXX_K562.bigBed"
   }

   # If visualization is not needed:
   file_paths = False
   ```

2. **Run the Main Function**:
   - Execute the following command to run predictions and mine association rules:
     ```python
     run_phoci(cell_line, model_dir, file_paths)
     ```
   - The script will automatically:
     - Generate candidate multi-way interaction samples using random walks.
     - Score and filter the candidates using the trained model.
     - Mine association rules for multi-way interactions.

3. **View Results**:
   - **Multi-Way Interaction Data**: Stored in the `apriori_rules` folder.
   - **Visualization Results**: Found in the `apriori_plots` folder, displaying significant association rules based on the Apriori algorithm.

---

## Key Features

1. **Multi-Way Interaction Probability Prediction**:
   - Combines Hi-C and epigenomic data using a GCN-based approach.
   - Provides chromatin multi-way interaction predictions for various cell types.

2. **Sample Generation and Validation**:
   - Generates candidate multi-way interaction samples using random walks.
   - Validates predictions with Pore-C experimental data.

3. **Association Rule Mining**:
   - Analyzes functional chromatin interaction rules based on predicted data.
   - Offers detailed visualization and analysis of enhancer-promoter relationships and other multi-way interactions.

---

## Contributors

Special thanks to the team members for their valuable contributions to this project:

*   **Yunyi Wu**, **Zhi Yang**, **Yanqing Wang**,  **Lipeng Li**, and **Kai Huang** (Huang Lab, Institute of Systems and Physical Biology, Shenzhen Bay Laboratory)
*   **Xing Jiang** and **Chen Yu** (Yu Lab, Institute of Cancer Research, Shenzhen Bay Laboratory)
*   **Jinsheng Xu** and **Chunhui Hou** (Hou Lab, State Key Laboratory of Genetic Resources and Evolution, Kunming Institute of Zoology, Chinese Academy of Sciences)
---

## License

This project is open-source under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

---

Feel free to contact us with any questions or suggestions!
