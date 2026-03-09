# RNA-seq Data Cleaning & Integration

## 📁 Project Structure

**Src Files**
- **`cleaning.py`**: Loads OHSU/TARGET expression and mapping tables, fixes orientation, maps Entrez→HUGO symbols, matches genes, applies variance filtering, log2-transforms and scales data, and saves `X_ohsu.npy`, `X_target.npy`, `ohsu_cleaned_expression.csv`, and `target_cleaned_expression.csv`. Run: `python src/cleaning.py`.

- **`fusions.py`**: Exploratory script that inspects the clinical TSV (`aml_ohsu_2022_clinical_data.tsv`) to list and count fusion-like labels (helpful for label curation). Run: `python src/fusions.py`.

- **`eln_classifier.py`**: Trains a ridge (L2) logistic regression to predict ELN favourable vs adverse using `cleaned/X_ohsu.npy` and `cleaned/ohsu_favorable_labels.npy`. Evaluates performance and saves `models/ridge_eln_model.joblib` and `models/ridge_scaler.joblib`. Run: `python src/eln_classifier.py`.

- **`predict_target.py`**: Loads `cleaned/X_target.npy`, the trained ELN model and scaler from `models/`, predicts ELN favourable probabilities for TARGET samples, and writes `results/TARGET_ELN_predictions.csv`. Run: `python src/predict_target.py`.

- **`predict_tcga.py`**: Loads `cleaned/X_tcga.npy`, the trained ELN model and scaler from `models/`, predicts ELN favourable probabilities for TCGA samples, and writes `results/TCGA_ELN_predictions.csv`. Run: `python src/predict_tcga.py`.

- **`target_fusion.py`**: Performs fusion type inference for TARGET samples using expression patterns and writes `results/fusion_inference/TARGET_fusion_inference.csv`. Run: `python src/target_fusion.py`.

- **`top_genes.py`**: Extracts and ranks top genes by model importance and expression patterns across datasets, generating `results/TARGET_top_genes.csv`, `results/TCGA_top_genes.csv`, and fusion-specific gene lists. Run: `python src/top_genes.py`.

- **`validation.py`**: Evaluates model performance on TCGA data using cytogenetics-derived ELN labels. Generates ROC curves, precision-recall curves, confusion matrices, and calibration plots in `results/tcga_validation/`. Run: `python src/validation.py`.


## ⚙️ Environment Setup (Conda)

Open a terminal in VS Code and run:

```bash
conda create -n rna_seq_env python=3.9 pandas numpy -y
conda activate rna_seq_env
```

Verify installation:

```bash
python -c "import pandas, numpy; print('Environment ready')"
```

---

## 📦 Install Dependencies (optional if using requirements.txt)

```bash
pip install -r requirements.txt
```

---

## ▶️ Run the Cleaning Script

From the project root or inside `src`:

```bash
cd src
python cleaning.py
```

---

## 🧠 What the Script Does

* Loads OHSU and TARGET RNA-seq datasets
* Matches common genes between datasets
* Transposes matrices (genes → columns)
* Log₂ transforms TARGET expression values
* Outputs aligned data ready for downstream analysis

---

## 📊 Input Data

The script expects the following files in the project root:

* `data_mrna_seq_rpkm.txt`
* `data_mrna_seq_tpm.txt`

Both must be tab-separated with gene names in the first column.

---

## ❗ Notes

* File paths in `cleaning.py` are set for **Windows**.
* If you move the project, update the paths accordingly.
* Large files (>100MB) are not recommended for GitHub without Git LFS.

---

## 🚀 Future Improvements

* Save cleaned datasets to `/output`
* Add exploratory data analysis
* Add model training pipeline
* Make paths OS-independent using `pathlib`

---

## 👤 Author

MSc Bioinformatics Project

---

