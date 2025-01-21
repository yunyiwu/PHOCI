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
|-- run_data_process_for_prediction.py   # Main data processing script
`-- utils.py                      # General utility module
```

---

## Quick Start

### Environment Setup

To streamline the installation process, use the provided `environment.yml` file to quickly create and configure a virtual environment. Ensure [conda](https://docs.conda.io/en/latest/miniconda.html) is installed, then follow these steps:

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yunyiwu/PHOCI.git
   cd PHOCI
   ```

2. **Create and activate the environment**:
   ```bash
   conda env create -f environment.yml
   conda activate phoci
   ```

3. **Run the project**:
   Execute the required scripts or applications as per your project needs.

---

### Overview of Environment File

The `environment.yml` file includes the following key dependencies:
- **Python**: Version 3.8 or later.
- **Main Frameworks and Tools**:
  - PyTorch
  - NumPy, Pandas
  - Scikit-learn
  - Matplotlib, Seaborn
  - NetworkX
  - Other relevant libraries (see file details for more).

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

   config_test["hic_file"] = "/public/home/wuyy/H1_hg38/4DNFID162B9J.hic"
   config_test["bigwig_dir"] = "/public/home/wuyy/H1_hg38/bigwigs/"
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
       'atac': "/public/home/wuyy/K562_hg38/bigwigs/ENCFF754EAC.bigWig",
       'h3k27ac': "/public/home/wuyy/K562_hg38/bigwigs/ENCFF381NDD.bigWig",
       'h3k4me1': "/public/home/wuyy/K562_hg38/bigwigs/ENCFF761XBZ.bigWig",
       'ctcf': "/public/home/wuyy/K562_hg38/bigwigs/ENCFF675GVW.bigWig",
       'rad21': "/public/home/wuyy/K562_hg38/bigwigs/ENCFF652NKM.bigWig",
       'chrom_state': "/public/home/wuyy/K562_hg38/chrom_state/ENCFF319VXX.bigBed"
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

Special thanks to the following team members for their contributions:
- Yunyi Wu
- Xing Jiang
- Zhi Yang
- Lipeng Li
- Jinsheng Xu

---

## License

This project is open-source under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

---

Feel free to contact us with any questions or suggestions!