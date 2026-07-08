---
language: en
license: mit
tags:
- 3d-genomics
- graph-neural-networks
- epigenomics
- chromatin-interactions
---

# Model Card: PHOCI (Predictor of Higher-Order Chromatin Interactions)

PHOCI is a deep learning-based framework designed for the probabilistic modeling and prediction of high-dimensional, multi-way chromatin interactions[cite: 1]. By integrating widely available pairwise Hi-C data with cell-type-specific epigenomic features, PHOCI bypasses the high costs and spatial sparsity associated with direct multi-way sequencing techniques (such as Pore-C)[cite: 1].

## 1. Model Details

- **Developed by:** Yunyi Wu, Xing Jiang, Zhi Yang, Yanqing Wang, Lipeng Li, Jinsheng Xu, Pan Deng, Chunhui Hou, Chen Yu, Kai Huang (Shenzhen Bay Laboratory & SMART)[cite: 1]
- **Model Type:** Graph Neural Network (GNN) combined with Multi-layer Perceptron (MLP)[cite: 1]
- **Architecture Overview:** 
  - **Graph Encoder:** A 3-layer GraphSAGE architecture that extracts localized spatial features from a 5-Mb cis-topological graph built using routine 2D Hi-C data[cite: 1].
  - **Feature Integration:** Merges graph topology with 15 cell-type-specific epigenetic features[cite: 1].
  - **ESA Module:** An Embedding Sampling and Aggregation module that processes candidate multi-locus clusters derived via biased random walks using a dual-aggregation mechanism (max-min difference and mean pooling)[cite: 1].
  - **Classifier:** An MLP that outputs a final probability score reflecting the likelihood of simultaneous multi-way contact[cite: 1].

## 2. Intended Uses & Limitations

### Intended Uses
- Predicting cell-type-specific 3-way to 6-way high-order chromatin interactions in uncharacterized cell lines using only baseline genomic and epigenetic inputs[cite: 1].
- Identifying recurrent multi-locus enhancer-promoter hubs and deciphering complex transcriptional regulatory logic[cite: 1].
- Prioritizing functional genomic hubs for downstream validation (e.g., multiplexed CRISPR interference assays)[cite: 1].

### Limitations & Biases
- **Fragment-shifting:** Predicted multi-way interaction hubs may exhibit a subtle spatial offset of 1 to 2 bins (approximately 5–10 kb) relative to physical locations due to inherent cross-linking ambiguities in the training data[cite: 1].
- **Resolution and Span Constraints:** The current framework operates at a fixed 5-kb resolution within local 5-Mb cis-subgraph windows, meaning it may be less sensitive to sub-5kb fine structures or long-range trans-chromosomal hubs[cite: 1].
- **Static Modeling:** The architecture does not explicitly account for real-time transcription factor binding kinetics, eRNA activity, or the highly dynamic behavior of phase-separated transcriptional condensates[cite: 1].

## 3. Training & Epigenomic Features

The model utilizes 5-kb genomic bins as graph nodes[cite: 1]. Each node is encoded with a 15-dimensional feature vector capturing key epigenetic and chromatin accessibility states, including[cite: 1]:
- **Histone Modifications:** H3K27ac (active enhancers), H3K4me3 (active promoters), H3K9me3 (heterochromatin/repression), etc[cite: 1].
- **Chromatin Accessibility:** ATAC-seq data[cite: 1].
- **Structural Proteins:** CTCF binding profiles[cite: 1].

The underlying topological backbone is derived from standard 2D Hi-C contact maps, which guide the biased random walk sampling to compress the combinatorial search space[cite: 1].

## 4. Evaluation & Validation Results

### Robustness Against Adversarial Perturbations
To ensure PHOCI captures authentic high-order cooperative behavior rather than mere spatial proximity, the model was evaluated against three increasingly stringent negative sampling strategies[cite: 1]:
1. **Sized Negative Sampling (SNS):** Tests model baseline distance and genomic biases[cite: 1].
2. **Motif Negative Sampling (MNS):** Increases local cluster confusion[cite: 1].
3. **Clique Negative Sampling (CNS):** A stringent test where a true multi-way cluster node is replaced with a "decoy" node that preserves all pairwise contacts but lacks higher-order coherence[cite: 1]. PHOCI consistently achieved high AUC and AP scores under CNS, demonstrating superior discrimination[cite: 1].

### Biological and Experimental Validation
- **Emergent Macrostructures:** Unsupervised clustering (UMAP & K-means) of the learned latent embeddings spontaneously recapitulates macroscopic genomic architectures, distinguishing A/B compartments and hierarchical TAD/sub-TAD structures[cite: 1].
- **Wet-Lab Validation (MYB Locus in K562):** PHOCI correctly predicted multi-way regulatory modules linking the *MYB* promoter with distinct distal enhancers[cite: 1]. Multiplexed CRISPRi experiments confirmed that dual perturbation of these predicted modules leads to non-additive, synergistic transcriptional repression, providing functional proof of the model's predictions[cite: 1].

## 5. Citation / Reference

If you use this model or the PHOCI framework in your research, please cite our preprint[cite: 1]:

```bibtex
@article{wu2026phoci,
  title={Predicting Higher-Order Chromatin Interactions with PHOCI},
  author={Wu, Yunyi and Jiang, Xing and Yang, Zhi and Wang, Yanqing and Li, Lipeng and Xu, Jinsheng and Deng, Pan and Hou, Chunhui and Yu, Chen and Huang, Kai},
  journal={Preprint},
  year={2026}
}
