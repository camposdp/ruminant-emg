# Dictionary-Based Classification of Bovine Masseter sEMG for Ingestive Pattern Recognition

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI Sensors Journal](https://img.shields.io/badge/DOI-IEEE%20Sensors%20Journal-10.1109%2FJSEN.2020.2977768-blue)](https://ieeexplore.ieee.org/abstract/document/9020064)
[![DOI Sensors Letters](https://img.shields.io/badge/DOI-IEEE%20Sensors%20Letters-10.1109%2FLSENS.2024.3233549-blue)](https://ieeexplore.ieee.org/abstract/document/10591327)

## 🐄 Project Overview

This repository contains MATLAB code for the **classification of bovine masseter surface electromyography (sEMG) signals** to recognize ingestive behaviors, particularly **eating vs. rumination**.

The main contributions are:

- **Double Threshold Onset Segmentation (DTOS)** to detect individual chewing events.
- **Sparse Dictionary Learning (FDDL and others)** to classify signals without manual feature extraction.
- Full comparison with handcrafted feature-based pipelines (RMS, MAV, WL, ZC, etc.).

This work was part of my PhD at UTFPR and is documented in:

- [PhD Thesis (2019)](https://riut.utfpr.edu.br/jspui/handle/1/4637)
- [IEEE Sensors Journal (2020)](https://ieeexplore.ieee.org/abstract/document/9020064)
- [IEEE Sensors Letters (2024)](https://ieeexplore.ieee.org/abstract/document/10591327)

---

## 📖 Background and Motivation

Monitoring **feed intake and rumination** is essential in Precision Livestock Farming (PLF). sEMG of the jaw (masseter muscle) allows precise tracking of chewing behaviors.

Challenges addressed:

- Traditional pipelines rely on manual feature extraction, which can fail under noise or variable conditions.
- Dictionary Learning allows **end-to-end classification** directly from raw signal segments, improving robustness and reducing the need for handcrafted features.

We implemented and tested:

- **DTOS segmentation**: isolates chew events.
- **FDDL (Fisher Discriminant Dictionary Learning)**: builds class-specific dictionaries for classification.
- Other algorithms: LC-KSVD, DLSI, SRC for benchmarking.

---

## 🔬 Pipeline Summary

1. **Acquisition**: sEMG signals recorded at 2000 Hz from the cow's masseter muscle.
2. **Preprocessing**: Bandpass filter (20–450 Hz), rectification, smoothing (150 samples window).
3. **Segmentation (DTOS)**: Thresholding method to isolate chew onsets and offsets.
4. **Feature Extraction (optional)**: Standard EMG features (MAV, RMS, WL, ZC, etc.).
5. **Dictionary Learning**: Training with class-specific atoms using FDDL and others.
6. **Classification**: New segments are sparsely coded; classification is based on minimal reconstruction error.

This approach outperforms traditional LDA pipelines, especially under noise.

---

## 📁 Repository Structure

| Folder                                                  | Description                                                                                                                            |
|---------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------|
| `codigos_sparse/`                                       | Main scripts: preprocessing, segmentation, training/testing (FDDL, LDA, SRC).                                                          |
| `DICTOL/`                                               | [DICTOL Toolbox](https://github.com/tiepvupsu/DICTOL) – dictionary learning algorithms (FDDL, LC-KSVD, DLSI, SRC, etc.).                |
| `DFDL/`                                                 | [DFDL Toolbox](https://github.com/tiepvupsu/DFDL) – optional, for additional dictionary learning experiments.                           |
| `features/`                                             | Time-domain EMG feature extraction functions.                                                                                           |
| `ClassifierToolbox-master/`                             | [Classifier Toolbox](https://github.com/kasai-kiryu/ClassifierToolbox) – optional, used for extra classifiers like SVM/LDA.             |
| `sparseLDA_v2/`                                         | Sparse Linear Discriminant Analysis code (feature selection).                                                                           |
| [`metodo/`](https://github.com/camposdp/ruminant-emg/tree/master/metodo) | Figures from the thesis/articles (plots, confusion matrices, etc.).                                                                    |
| `.mat` files                                            | Data files (e.g., `EMG_base.mat`, `Segments_labels_Cpartition.mat`).                                                                    |

---

## 🌐 External Dependencies

- **DICTOL:** https://github.com/tiepvupsu/DICTOL
- **DFDL:** https://github.com/tiepvupsu/DFDL
- **Classifier Toolbox:** https://github.com/kasai-kiryu/ClassifierToolbox

These are bundled here for convenience but check their original licenses for reuse and citation requirements.

## ▶️ How to Run the Code

This section explains how to **reproduce the experiments** and results.

---

### 1️⃣ Preprocessing & Segmentation

**Script: `PROJETO_A_preprocess_bovinos_criarMatriz_classifica.m`**

- Loads raw EMG data (`EMG_base.mat`).
- Applies:
  - Bandpass filtering (20–450 Hz),
  - Rectification,
  - Smoothing (moving average),
  - **Double Threshold Onset Segmentation (DTOS)**.

- Extracts handcrafted features (if enabled).
- Saves:
  - `DATA_features_and_labels.mat`: feature matrix + labels.
  - `Segments_labels_Cpartition.mat`: segmented raw signals + labels + cross-validation partition.

---

### 2️⃣ Traditional Feature-Based Classification

**Script: `PROJETO_A_classificar_bovinos_variando_segment.m`**

- Loads feature matrix from `DATA_features_and_labels.mat`.
- Runs LDA classifier (default) or others if configured.
- Outputs:
  - Accuracy scores,
  - Confusion matrices,
  - Optionally plots performance.

This pipeline reproduces the **baseline** comparisons between DTOS vs. fixed windows vs. manual (visual) segmentation.

---

### 3️⃣ Dictionary Learning (FDDL)

**Script: `treinar_FDDL.m`**

- Loads segmented signal data (`Segments_labels_Cpartition.mat`).
- Defines:
  - `Ki`: dictionary size(s),
  - `lambda1/lambda2`: FDDL regularization,
  - `sparsitythres`: sparsity level.

- Runs **10-fold cross-validation**:
  - Trains dictionary,
  - Tests classification on each fold.

- Outputs:
  - Accuracy per fold,
  - Accuracy vs. dictionary size plot (if enabled).

---

### 4️⃣ Additional Experiments

**Available scripts:**

- `treinar_LC_kSVD.m`: Label-Consistent K-SVD.
- `treinar_DLSI.m`: Discriminative Local Sparse Coding.
- `treinar_SRC.m`: Sparse Representation-based Classifier.
- `treinar_FDDL_bestresults_with_noise.m`: Tests performance with added Gaussian noise.
- `treinar_FDDL_variando_dim_dic.m`: Tests multiple dictionary sizes.

---

### 🔎 Plotting & Visualization

- `plotar_confusion.m`: Confusion matrix plot.
- `plotar_resultados.m`: Accuracy vs. dictionary size or SNR.
- `plotar_exemplo_do_processo.m`: Pipeline illustration (signal -> segmentation -> classification).

Example figures are available in [`metodo/`](https://github.com/camposdp/ruminant-emg/tree/master/metodo).

---

## ⚙️ Important Parameters

| Parameter           | Description                                                  | Default                         |
|---------------------|--------------------------------------------------------------|---------------------------------|
| `fs`                | Sampling rate                                                | 2000 Hz                         |
| `Wsmooth`           | Smoothing window (samples)                                   | 150 samples                     |
| `threshold (k)`     | DTOS threshold multiplier                                    | 4.0 × std (baseline)            |
| `Wcrit`             | Min/max segment length                                       | [0.1×fs, 1.0×fs]                |
| `Ki`                | Dictionary sizes tested                                      | [100, 200, 300, ...]            |
| `lambda1 / lambda2` | FDDL regularization weights                                  | ~0.01 each                      |
| `sparsitythres`     | Sparsity level for sparse coding                             | Typically 0.05–0.1              |
| `folds`             | Number of folds in cross-validation                          | 10                              |

You can edit these at the top of each script as needed.

---

## 🔧 Installation & Setup

- **MATLAB:** R2018b / R2019a recommended (later versions likely work).
- **Toolboxes:** No special MATLAB toolboxes needed (external libraries bundled).
- **SPAMS / MEX:** Precompiled binaries included (for Windows/Linux). Recompile if using macOS.

**Add paths:**

```matlab
addpath(genpath('path_to_project/codigos_sparse'));
addpath(genpath('path_to_project/DICTOL'));
addpath(genpath('path_to_project/DFDL'));
addpath(genpath('path_to_project/features'));
addpath(genpath('path_to_project/ClassifierToolbox-master'));
addpath(genpath('path_to_project/sparseLDA_v2'));
savepath;

## 📜 Citation

If you use this repository or part of the code in your research, please cite the following works:

- **Campos, D.P., Abatti, P.J., Bertotti, F.L., Hill, J.A.G., Silveira, A.L.F.**  
  *Surface electromyography segmentation and feature extraction for ingestive behavior recognition in ruminants.*  
  Computers and Electronics in Agriculture, vol. 153, pp. 325–333, 2018.  
  [DOI: 10.1016/j.compag.2018.08.039](https://doi.org/10.1016/j.compag.2018.08.039)

- **Campos, D.P., Lazzaretti, A.E., Bertotti, F.L., Gomes, O.A., Hill, J.A.G., Silveira, A.L.F.**  
  *Single-Channel sEMG Dictionary Learning Classification of Ingestive Behavior in Cows.*  
  IEEE Sensors Journal, vol. 20, no. 13, pp. 7199–7207, 2020.  
  [DOI: 10.1109/JSEN.2020.2977768](https://ieeexplore.ieee.org/abstract/document/9020064)

- **Campos, D.P., Lazzaretti, A.E., Bertotti, F.L., Hill, J.A.G., Silveira, A.L.F.**  
  *Event-Driven sEMG Feature Extraction for Classifying Ingestive Behavior in Ruminants.*  
  IEEE Sensors Letters, vol. 8, no. 8, 2024.  
  [DOI: 10.1109/LSENS.2024.3233549](https://ieeexplore.ieee.org/abstract/document/10591327)

- **Campos, D.P.**  
  *PhD Thesis: Classificação de sinais eletromiográficos do músculo masseter de bovinos baseada em dicionários para reconhecimento de padrões ingestivos.*  
  Universidade Tecnológica Federal do Paraná (UTFPR), 2019.  
  [Link](https://riut.utfpr.edu.br/jspui/handle/1/4637)

👉 **For DICTOL and DFDL users:** please cite the original repositories:

- **DICTOL Toolbox:** https://github.com/tiepvupsu/DICTOL  
- **DFDL Toolbox:** https://github.com/tiepvupsu/DFDL

And acknowledge Tiep Vu and collaborators' original publications.

---

## 🙌 Acknowledgements

This research was supported by:

- **CAPES** and **CNPq** (Brazilian research funding agencies),
- **UTFPR** (Federal University of Technology – Paraná),
- Experimental and field support from **Instituto de Desenvolvimento Rural do Paraná (IDR-Paraná).**

Special thanks to:

- Prof. Paulo José Abatti (advisor),
- Prof. André Eugênio Lazzaretti (advisor),
- Collaborators: Fábio Luiz Bertotti, João Ari Gualberto Hill, André Luís Finkler da Silveira.

We also thank **Tiep Vu** and **Hiroyuki Kasai** for making their toolboxes available to the research community.

---

## 📬 Contact

For questions, suggestions, or collaboration opportunities, feel free to contact:

**Daniel Prado de Campos**  
✉️ danielcampos@utfpr.edu.br

Or open an issue in the repository.

---

## 🔒 License

This project’s original code is released under the **GPL v3 License** for academic and research use.  
For third-party dependencies (DICTOL, DFDL, Classifier Toolbox), check their respective licenses.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

---

