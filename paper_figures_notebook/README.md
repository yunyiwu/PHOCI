# Reproducibility Notebooks for PHOCI Figures

This directory contains Jupyter Notebooks for reproducing the figures and analysis presented in the **PHOCI** paper.

---

## 🛠️ Environment Requirements

Before running the notebooks, please ensure you have configured the appropriate environment:

1. **PHOCI Environment (Default)**:
   - Used for the majority of the analysis and visualization notebooks in this directory.
   - Please follow the installation instructions in the main repository root (`../README.md`) to set up the dependencies.

2. **AlphaGenome Environment**:
   - Required specifically for `Fig2(m).ipynb` and `Fig6(c,e,g)&FigS14.ipynb`.
   - Please set up the environment following the instructions from Google DeepMind's official repository:  
     👉 [alphagenome_research](https://github.com/google-deepmind/alphagenome_research)

---

## 📂 Notebook Summary & Requirements

| Notebook | Corresponding Figures | Required Environment | Dependencies / Prerequisites |
| :--- | :--- | :--- | :--- |
| `Fig1(c-f).ipynb` | Figure 1 (c–f) | `PHOCI` | Standard execution |
| `Fig2(a-c).ipynb` | Figure 2 (a–c) | `PHOCI` | Full pipeline for Fig 2(a); modify `model_eval` model parameter for Fig 2(b, c) |
| `Fig2(d-f).ipynb` | Figure 2 (d–f) | `PHOCI` | Standard execution |
| `Fig2(g-i).ipynb` | Figure 2 (g–i) | `PHOCI` | Standard execution |
| `Fig2(j).ipynb` | Figure 2 (j) | `PHOCI` | Standard execution |
| `Fig2(k).ipynb` | Figure 2 (k) | `PHOCI` | Standard execution |
| `Fig2(m).ipynb` | Figure 2 (m) | `AlphaGenome` | Requires AlphaGenome environment |
| `Fig3.ipynb` | Figure 3 | `PHOCI` | Requires TensorBoard export for UMAP; demonstrates a representative example |
| `Fig4.ipynb` | Figure 4 | `PHOCI` | Standard execution |
| `Fig5.ipynb` | Figure 5 | `PHOCI` | Plotting code only; requires pre-running PHOCI inference |
| `Fig6(a).ipynb` | Figure 6 (a) | `PHOCI` | Plotting code only; requires pre-running PHOCI inference |
| `Fig6(c,e,g)&FigS14.ipynb` | Figure 6 (c, e, g) & Fig S14 | `AlphaGenome` | Requires AlphaGenome environment |

---

## ⚠️ Important Usage Notes

### 1. Generating Figure 2(a–c) Panels (`Fig2(a-c).ipynb`)
- `Fig2(a-c).ipynb` contains the full sequential code for generating **Figure 2(a)**.
- To produce the results for **Figure 2(b)** and **Figure 2(c)**, you only need to modify the evaluation model parameter (`eval` model) in the notebook code.

### 2. Intermediate TensorBoard Export & Example Run (`Fig3.ipynb`)
- **TensorBoard Projection:** Running `Fig3.ipynb` requires exporting intermediate embeddings into **TensorBoard** (TensorBoard Projector) to perform the **UMAP** dimensionality reduction and spatial layout visualization.
- **Running Other Examples:** The provided notebook illustrates the complete pipeline using one representative example. You can reproduce results for other target cases or loci by modifying the corresponding configuration parameters within the notebook.

### 3. Pre-computed Results Required (`Fig5.ipynb` & `Fig6(a).ipynb`)
- These notebooks contain visualization and plotting scripts designed for specific cell lines and genomic locus prediction results.
- **Prerequisite:** You must first execute the main PHOCI prediction pipeline to generate the output prediction files for the corresponding cell lines and genes before executing these plotting scripts.

### 4. Data Availability & Contact
- Most required evaluation and plotting datasets are provided directly in this repository or associated Zenodo releases.
- Due to file size limitations, a small fraction of raw/intermediate output data may not be fully included in the repository.
- If you need additional raw data or intermediate files to reproduce specific figures, please feel free to reach out or open an issue.
