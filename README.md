# SRX-DTI

SRX-DTI is a machine learning framework for predicting Drug–Target Interactions (DTIs).  
This method integrates multiple protein sequence descriptors, molecular fingerprints, data balancing strategies, feature selection techniques, and gradient boosting to improve DTI prediction performance.

---

## Method Overview

The SRX-DTI framework consists of the following stages:

1. **Feature Extraction**
   - Protein sequence descriptors (e.g., AAC and others)
   - Drug encoding using FP2 molecular fingerprints

2. **Handling Imbalanced Data**
   - Proposed **One-SVM-US** technique for data balancing

3. **Feature Selection**
   - Proposed **FFS-RF** algorithm to remove irrelevant features
   - Selection of optimal feature subset

4. **Classification**
   - Final balanced dataset is trained using **XGBoost**
   - Performance evaluated using five-fold cross-validation

---

## Datasets

This study uses four benchmark (gold standard) datasets released by Yamanishi et al.:

- Enzymes (EN)
- G-protein-coupled receptors (GPCR)
- Ion Channels (IC)
- Nuclear Receptors (NR)

Datasets are publicly available at:  
http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/

---

## Framework Architecture

![SRX-DTI Framework](https://user-images.githubusercontent.com/72028345/204578716-30f41a3e-0f22-4881-82dc-f0af97e1eb52.png)

---

## Environment

- Python 3.9 or higher
- Required Python libraries (see project dependencies)

---

## How to Run

1️⃣ Feature extraction  
```
./feature extraction/00-AAC.py
```
(For other descriptors, modify the corresponding script.)

2️⃣ Data balancing and feature selection  
```
./NR-run/run.py
```

3️⃣ DTI prediction and evaluation  
```
./DTI prediction/DTI_predictor.py
```

---

## Published Paper

Full methodological details are available in our published paper:

Khojasteh, H., Pirgazi, J., & Ghanbari Sorkhi, A. (2023).  
**Improving prediction of drug-target interactions based on fusing multiple features with data balancing and feature selection techniques.**  
PLOS ONE, 18(8), e0288173.

🔗 https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0288173

---

## Contact

Hakimeh Khojasteh  
khojasteh.hb@gmail.com  

---

## Citation

If you use this code in your research, please cite:

```
Khojasteh, H., Pirgazi, J., & Ghanbari Sorkhi, A. (2023).
Improving prediction of drug-target interactions based on fusing multiple features with data balancing and feature selection techniques.
PLOS ONE, 18(8), e0288173.
```
