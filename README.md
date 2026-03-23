# RNA-seq Data Cleaning & Integration

## 📁 Project Structure

**Src Files** (execution order for a full pipeline):

1. **`cleaning.py`**: ingest OHSU/TARGET expression, map Entrez→HUGO, align genes, apply variance filter, log2 + scaling, save cleaned matrices.
   - output: `cleaned/ohsu_cleaned_expression.csv`, `cleaned/target_cleaned_expression.csv`, `cleaned/X_ohsu.npy`, `cleaned/X_target.npy`.
   - run: `python src/cleaning.py`

2. **`eln_classifier.py`**: train ridge logistic regression on OHSU with ELN labels, save model and scaler.
   - input: `cleaned/X_ohsu.npy`, `cleaned/ohsu_favorable_labels.npy`
   - output: `models/ridge_eln_model.joblib`, `models/ridge_scaler.joblib`
   - run: `python src/eln_classifier.py`

3. **`predict_oshu.py`**: apply trained model to OHSU 2018 cohort for validation.
   - input: `cleaned/X_validata.npy`, `models/*`
   - output: `results/OHSU_favourable_fusion_predictions.csv`
   - run: `python src/predict_oshu.py`

4. **`predict_target.py`**: apply trained model to TARGET cohort.
   - input: `cleaned/X_target.npy`, `models/*`
   - output: `results/TARGET_favourable_fusion_predictions.csv`
   - run: `python src/predict_target.py`

5. **`predict_tcga.py`**: apply trained model to TCGA cohort.
   - input: `cleaned/X_tcga.npy`, `models/*`
   - output: `results/TCGA_favourable_fusion_predictions.csv`
   - run: `python src/predict_tcga.py`

6. **`target_fusion.py`**: run fusion program inference pipeline for TARGET/TCGA.
   - output: `results/fusion_inference/TARGET_fusion_inference.csv`
   - run: `python src/target_fusion.py`

7. **`top_genes.py`**: analyze top gene expression differences and model coefficients for OHSU/TCGA/TARGET and VALIDATA (OHSU 2018), output CSVs and figures.
   - output:
     - `results/OHSU_2022_top_genes.csv`, `results/OHSU_2018_top_genes.csv`, `results/TCGA_top_genes.csv`, `results/TARGET_top_genes.csv`
     - gene overlap files + plots in `results/figures`
   - run: `python src/top_genes.py`

8. **`validation.py`**: benchmark the model on TCGA using cytogenetics-derived labels.
   - output: ROC/PR/confusion/calibration plots in `results/tcga_validation/`
   - run: `python src/validation.py`

---
or streamlit run app.py

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

Erin Wright

---

